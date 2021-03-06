#define _CRT_SECURE_NO_WARNINGS 1
#include <vector>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

#include <random>
#include <cmath>
#include <omp.h>
#include <string>
#include <iostream>
#include <stdio.h>
#include <algorithm>
#include <list>

#define _USE_MATH_DEFINES
#include <math.h>

#include "Vector.h"
#include <chrono>

// TODO : fix to automatically adjust for the number of CPU cores

std::vector<std::default_random_engine> rng; // Modifier par le nombre de coeurs de l'ordi

void init_rng(std::vector<std::default_random_engine> &rng)
{
    std::random_device r;
    for (int i = 0; i < omp_get_max_threads(); i++)
    {
        rng.push_back(std::default_random_engine(r()));
    };
};

std::uniform_real_distribution<double> uniform(0., 1.);

const double epsilon(0.01); // Pour corriger les bugs liés à la précision numérique
const int maxReflectionNumber(5);

struct Matter
{
    Vector albedo;
    bool mirror;
    bool transparent;
    double n;
};

struct Intersection
{
    bool exists;
    Vector point;
    double t;
    Vector normalVector;
    Vector albedo;
};

class Ray
{
public:
    explicit Ray(const Vector &origin, const Vector &direction) : origin(origin), direction(direction){};
    Vector origin, direction;
};

class Object
{
public:
    virtual Intersection intersect(const Ray &r) const = 0;

    Vector origin;
    Matter matter;
    double intensity;
};

class BBox
{
public:
    ~BBox(){};

    explicit BBox(Vector minPoint = Vector(0, 0, 0), Vector maxPoint = Vector(0, 0, 0))
    {
        this->minPoint = minPoint;
        this->maxPoint = maxPoint;
    };

    Intersection intersect(const Ray &r) const
    {
        double xMin = (minPoint[0] - r.origin[0]) / r.direction[0];
        double xMax = (maxPoint[0] - r.origin[0]) / r.direction[0];
        double yMin = (minPoint[1] - r.origin[1]) / r.direction[1];
        double yMax = (maxPoint[1] - r.origin[1]) / r.direction[1];
        double zMin = (minPoint[2] - r.origin[2]) / r.direction[2];
        double zMax = (maxPoint[2] - r.origin[2]) / r.direction[2];

        if (xMin > xMax)
            std::swap(xMin, xMax);
        if (yMin > yMax)
            std::swap(yMin, yMax);
        if (zMin > zMax)
            std::swap(zMin, zMax);

        double tMin = std::max(xMin, std::max(yMin, zMin));
        double tMax = std::min(xMax, std::min(yMax, zMax));

        bool hasInter = (tMax > 0 && tMin <= tMax);
        // On ne s'intéresse qu'au premier booléen, les autres champs sont remplis avec n'importe quoi
        return Intersection({hasInter, r.origin, 0, r.direction, Vector(0, 0, 0)});
    };

    Vector minPoint;
    Vector maxPoint;
};

class BVH
{
public:
    explicit BVH()
    {
    }

    BVH *lc, *rc;
    int lowerIndex, upperIndex;
    BBox boundingBox;
};

class Sphere : public Object
{
public:
    explicit Sphere(const Vector &origin, const double radius, const Matter &matter, const bool reverseNormal = false, const double intensity = 0)
    {
        this->origin = origin;
        this->radius = radius;
        this->matter = matter;
        this->reverseNormal = reverseNormal;
        this->intensity = intensity;
    }

    // Fonction qui calcule le vecteur normal étant donné un point
    Vector getNormalVector(const Vector &point) const
    {
        Vector normalVector = (point - origin);
        normalVector.normalize();
        return reverseNormal ? (-1) * normalVector : normalVector;
    }

    double radius;
    bool reverseNormal;
    bool isLight;

    virtual Intersection intersect(const Ray &r) const
    {
        double a = 1;
        double b = 2 * dot(r.direction, r.origin - origin);
        double c = dot(r.origin - origin, r.origin - origin) - radius * radius;

        double delta = b * b - 4 * a * c;

        Intersection intersection;
        Vector intersectionPoint;
        double t(0);

        // Cas sans racines réelles : pas d'intersection, on renvoie n'importe quoi
        if (delta < 0)
        {
            intersection = {false, intersectionPoint, t, Vector(0, 0, 0), Vector(0, 0, 0)};
            return intersection;
        }

        /*
        A partir de ici, delta >= 0. On calcule les deux valeurs de
        l'intersection. On s'intéresse à la plus petite intersection ayant t > 0
        */

        double t1 = (-b - sqrt(delta)) / (2 * a);
        double t2 = (-b + sqrt(delta)) / (2 * a);

        // Cas sans intersection
        if (t2 < 0)
        {
            intersection = {false, intersectionPoint, t, Vector(0, 0, 0), Vector(0, 0, 0)};
            return intersection;
        }

        t = t1 > 0 ? t1 : t2;                           // Soit t1 soit t2 en fonction de la positivité de t1
        intersectionPoint = r.origin + t * r.direction; // Le point associé
        intersection = {true, intersectionPoint, t, this->getNormalVector(intersectionPoint), matter.albedo};
        return intersection;
    }
};

class TriangleIndices
{
public:
    TriangleIndices(int vtxi = -1, int vtxj = -1, int vtxk = -1, int ni = -1, int nj = -1, int nk = -1, int uvi = -1, int uvj = -1, int uvk = -1, int group = -1, bool added = false) : vtxi(vtxi), vtxj(vtxj), vtxk(vtxk), uvi(uvi), uvj(uvj), uvk(uvk), ni(ni), nj(nj), nk(nk), group(group){};
    int vtxi, vtxj, vtxk; // indices within the vertex coordinates array
    int uvi, uvj, uvk;    // indices within the uv coordinates array
    int ni, nj, nk;       // indices within the normals array
    int group;            // face group
};

class TriangleMesh : public Object
{
public:
    ~TriangleMesh() {}

    explicit TriangleMesh(const Vector &origin, const Matter &matter)
    {
        this->origin = origin;
        this->matter = matter;
    }

    void readOBJ(const char *obj)
    {

        char matfile[255];
        char grp[255];

        FILE *f;
        f = fopen(obj, "r");
        int curGroup = -1;
        while (!feof(f))
        {
            char line[255];
            if (!fgets(line, 255, f))
                break;

            std::string linetrim(line);
            linetrim.erase(linetrim.find_last_not_of(" \r\t") + 1);
            strcpy(line, linetrim.c_str());

            if (line[0] == 'u' && line[1] == 's')
            {
                sscanf(line, "usemtl %[^\n]\n", grp);
                curGroup++;
            }

            if (line[0] == 'v' && line[1] == ' ')
            {
                Vector vec;

                Vector col;
                if (sscanf(line, "v %lf %lf %lf %lf %lf %lf\n", &vec[0], &vec[1], &vec[2], &col[0], &col[1], &col[2]) == 6)
                {
                    col[0] = std::min(1., std::max(0., col[0]));
                    col[1] = std::min(1., std::max(0., col[1]));
                    col[2] = std::min(1., std::max(0., col[2]));

                    vertices.push_back(vec);
                    vertexcolors.push_back(col);
                }
                else
                {
                    sscanf(line, "v %lf %lf %lf\n", &vec[0], &vec[1], &vec[2]);
                    vertices.push_back(vec);
                }
            }
            if (line[0] == 'v' && line[1] == 'n')
            {
                Vector vec;
                sscanf(line, "vn %lf %lf %lf\n", &vec[0], &vec[1], &vec[2]);
                normals.push_back(vec);
            }
            if (line[0] == 'v' && line[1] == 't')
            {
                Vector vec;
                sscanf(line, "vt %lf %lf\n", &vec[0], &vec[1]);
                uvs.push_back(vec);
            }
            if (line[0] == 'f')
            {
                TriangleIndices t;
                int i0, i1, i2, i3;
                int j0, j1, j2, j3;
                int k0, k1, k2, k3;
                int nn;
                t.group = curGroup;

                char *consumedline = line + 1;
                int offset;

                nn = sscanf(consumedline, "%u/%u/%u %u/%u/%u %u/%u/%u%n", &i0, &j0, &k0, &i1, &j1, &k1, &i2, &j2, &k2, &offset);
                if (nn == 9)
                {
                    if (i0 < 0)
                        t.vtxi = vertices.size() + i0;
                    else
                        t.vtxi = i0 - 1;
                    if (i1 < 0)
                        t.vtxj = vertices.size() + i1;
                    else
                        t.vtxj = i1 - 1;
                    if (i2 < 0)
                        t.vtxk = vertices.size() + i2;
                    else
                        t.vtxk = i2 - 1;
                    if (j0 < 0)
                        t.uvi = uvs.size() + j0;
                    else
                        t.uvi = j0 - 1;
                    if (j1 < 0)
                        t.uvj = uvs.size() + j1;
                    else
                        t.uvj = j1 - 1;
                    if (j2 < 0)
                        t.uvk = uvs.size() + j2;
                    else
                        t.uvk = j2 - 1;
                    if (k0 < 0)
                        t.ni = normals.size() + k0;
                    else
                        t.ni = k0 - 1;
                    if (k1 < 0)
                        t.nj = normals.size() + k1;
                    else
                        t.nj = k1 - 1;
                    if (k2 < 0)
                        t.nk = normals.size() + k2;
                    else
                        t.nk = k2 - 1;
                    indices.push_back(t);
                }
                else
                {
                    nn = sscanf(consumedline, "%u/%u %u/%u %u/%u%n", &i0, &j0, &i1, &j1, &i2, &j2, &offset);
                    if (nn == 6)
                    {
                        if (i0 < 0)
                            t.vtxi = vertices.size() + i0;
                        else
                            t.vtxi = i0 - 1;
                        if (i1 < 0)
                            t.vtxj = vertices.size() + i1;
                        else
                            t.vtxj = i1 - 1;
                        if (i2 < 0)
                            t.vtxk = vertices.size() + i2;
                        else
                            t.vtxk = i2 - 1;
                        if (j0 < 0)
                            t.uvi = uvs.size() + j0;
                        else
                            t.uvi = j0 - 1;
                        if (j1 < 0)
                            t.uvj = uvs.size() + j1;
                        else
                            t.uvj = j1 - 1;
                        if (j2 < 0)
                            t.uvk = uvs.size() + j2;
                        else
                            t.uvk = j2 - 1;
                        indices.push_back(t);
                    }
                    else
                    {
                        nn = sscanf(consumedline, "%u %u %u%n", &i0, &i1, &i2, &offset);
                        if (nn == 3)
                        {
                            if (i0 < 0)
                                t.vtxi = vertices.size() + i0;
                            else
                                t.vtxi = i0 - 1;
                            if (i1 < 0)
                                t.vtxj = vertices.size() + i1;
                            else
                                t.vtxj = i1 - 1;
                            if (i2 < 0)
                                t.vtxk = vertices.size() + i2;
                            else
                                t.vtxk = i2 - 1;
                            indices.push_back(t);
                        }
                        else
                        {
                            nn = sscanf(consumedline, "%u//%u %u//%u %u//%u%n", &i0, &k0, &i1, &k1, &i2, &k2, &offset);
                            if (i0 < 0)
                                t.vtxi = vertices.size() + i0;
                            else
                                t.vtxi = i0 - 1;
                            if (i1 < 0)
                                t.vtxj = vertices.size() + i1;
                            else
                                t.vtxj = i1 - 1;
                            if (i2 < 0)
                                t.vtxk = vertices.size() + i2;
                            else
                                t.vtxk = i2 - 1;
                            if (k0 < 0)
                                t.ni = normals.size() + k0;
                            else
                                t.ni = k0 - 1;
                            if (k1 < 0)
                                t.nj = normals.size() + k1;
                            else
                                t.nj = k1 - 1;
                            if (k2 < 0)
                                t.nk = normals.size() + k2;
                            else
                                t.nk = k2 - 1;
                            indices.push_back(t);
                        }
                    }
                }

                consumedline = consumedline + offset;

                while (true)
                {
                    if (consumedline[0] == '\n')
                        break;
                    if (consumedline[0] == '\0')
                        break;
                    nn = sscanf(consumedline, "%u/%u/%u%n", &i3, &j3, &k3, &offset);
                    TriangleIndices t2;
                    t2.group = curGroup;
                    if (nn == 3)
                    {
                        if (i0 < 0)
                            t2.vtxi = vertices.size() + i0;
                        else
                            t2.vtxi = i0 - 1;
                        if (i2 < 0)
                            t2.vtxj = vertices.size() + i2;
                        else
                            t2.vtxj = i2 - 1;
                        if (i3 < 0)
                            t2.vtxk = vertices.size() + i3;
                        else
                            t2.vtxk = i3 - 1;
                        if (j0 < 0)
                            t2.uvi = uvs.size() + j0;
                        else
                            t2.uvi = j0 - 1;
                        if (j2 < 0)
                            t2.uvj = uvs.size() + j2;
                        else
                            t2.uvj = j2 - 1;
                        if (j3 < 0)
                            t2.uvk = uvs.size() + j3;
                        else
                            t2.uvk = j3 - 1;
                        if (k0 < 0)
                            t2.ni = normals.size() + k0;
                        else
                            t2.ni = k0 - 1;
                        if (k2 < 0)
                            t2.nj = normals.size() + k2;
                        else
                            t2.nj = k2 - 1;
                        if (k3 < 0)
                            t2.nk = normals.size() + k3;
                        else
                            t2.nk = k3 - 1;
                        indices.push_back(t2);
                        consumedline = consumedline + offset;
                        i2 = i3;
                        j2 = j3;
                        k2 = k3;
                    }
                    else
                    {
                        nn = sscanf(consumedline, "%u/%u%n", &i3, &j3, &offset);
                        if (nn == 2)
                        {
                            if (i0 < 0)
                                t2.vtxi = vertices.size() + i0;
                            else
                                t2.vtxi = i0 - 1;
                            if (i2 < 0)
                                t2.vtxj = vertices.size() + i2;
                            else
                                t2.vtxj = i2 - 1;
                            if (i3 < 0)
                                t2.vtxk = vertices.size() + i3;
                            else
                                t2.vtxk = i3 - 1;
                            if (j0 < 0)
                                t2.uvi = uvs.size() + j0;
                            else
                                t2.uvi = j0 - 1;
                            if (j2 < 0)
                                t2.uvj = uvs.size() + j2;
                            else
                                t2.uvj = j2 - 1;
                            if (j3 < 0)
                                t2.uvk = uvs.size() + j3;
                            else
                                t2.uvk = j3 - 1;
                            consumedline = consumedline + offset;
                            i2 = i3;
                            j2 = j3;
                            indices.push_back(t2);
                        }
                        else
                        {
                            nn = sscanf(consumedline, "%u//%u%n", &i3, &k3, &offset);
                            if (nn == 2)
                            {
                                if (i0 < 0)
                                    t2.vtxi = vertices.size() + i0;
                                else
                                    t2.vtxi = i0 - 1;
                                if (i2 < 0)
                                    t2.vtxj = vertices.size() + i2;
                                else
                                    t2.vtxj = i2 - 1;
                                if (i3 < 0)
                                    t2.vtxk = vertices.size() + i3;
                                else
                                    t2.vtxk = i3 - 1;
                                if (k0 < 0)
                                    t2.ni = normals.size() + k0;
                                else
                                    t2.ni = k0 - 1;
                                if (k2 < 0)
                                    t2.nj = normals.size() + k2;
                                else
                                    t2.nj = k2 - 1;
                                if (k3 < 0)
                                    t2.nk = normals.size() + k3;
                                else
                                    t2.nk = k3 - 1;
                                consumedline = consumedline + offset;
                                i2 = i3;
                                k2 = k3;
                                indices.push_back(t2);
                            }
                            else
                            {
                                nn = sscanf(consumedline, "%u%n", &i3, &offset);
                                if (nn == 1)
                                {
                                    if (i0 < 0)
                                        t2.vtxi = vertices.size() + i0;
                                    else
                                        t2.vtxi = i0 - 1;
                                    if (i2 < 0)
                                        t2.vtxj = vertices.size() + i2;
                                    else
                                        t2.vtxj = i2 - 1;
                                    if (i3 < 0)
                                        t2.vtxk = vertices.size() + i3;
                                    else
                                        t2.vtxk = i3 - 1;
                                    consumedline = consumedline + offset;
                                    i2 = i3;
                                    indices.push_back(t2);
                                }
                                else
                                {
                                    consumedline = consumedline + 1;
                                }
                            }
                        }
                    }
                }
            }
        }
        fclose(f);
    }

    void loadTexture(const char *textureFile)
    {
        useTextures = true;
        int width;
        int heigth;
        int channels;
        unsigned char *texture = stbi_load(textureFile, &width, &heigth, &channels, 3);
        widthList.push_back(width);
        heigthList.push_back(heigth);
        textures.push_back(texture);
    }

    BBox buildBB(int lowerIndex, int upperIndex)
    {
        TriangleIndices tri = indices[lowerIndex];
        double xMin = vertices[tri.vtxi][0];
        double xMax = vertices[tri.vtxi][0];
        double yMin = vertices[tri.vtxi][1];
        double yMax = vertices[tri.vtxi][1];
        double zMin = vertices[tri.vtxi][2];
        double zMax = vertices[tri.vtxi][2];

        for (int k = lowerIndex; k < upperIndex; k++)
        {
            TriangleIndices tri = indices[k];

            xMin = std::min(std::min(xMin, vertices[tri.vtxi][0]), std::min(vertices[tri.vtxj][0], vertices[tri.vtxk][0]));
            xMax = std::max(std::max(xMax, vertices[tri.vtxi][0]), std::max(vertices[tri.vtxj][0], vertices[tri.vtxk][0]));

            yMin = std::min(std::min(yMin, vertices[tri.vtxi][1]), std::min(vertices[tri.vtxj][1], vertices[tri.vtxk][1]));
            yMax = std::max(std::max(yMax, vertices[tri.vtxi][1]), std::max(vertices[tri.vtxj][1], vertices[tri.vtxk][1]));

            zMin = std::min(std::min(zMin, vertices[tri.vtxi][2]), std::min(vertices[tri.vtxj][2], vertices[tri.vtxk][2]));
            zMax = std::max(std::max(zMax, vertices[tri.vtxi][2]), std::max(vertices[tri.vtxj][2], vertices[tri.vtxk][2]));
        }

        Vector minPoint(xMin, yMin, zMin);
        Vector maxPoint(xMax, yMax, zMax);

        return BBox(minPoint, maxPoint);
    }

    void buildBVH(BVH *node, int lowerIndex, int upperIndex)
    {
        node->boundingBox = buildBB(lowerIndex, upperIndex);
        node->lowerIndex = lowerIndex;
        node->upperIndex = upperIndex;
        node->lc = NULL;
        node->rc = NULL;

        Vector diag = node->boundingBox.maxPoint - node->boundingBox.minPoint;

        int dim(-1);

        // On cherche l'axe ayant la plus grande diagonale
        if (diag[0] >= std::max(diag[1], diag[2]))
        {
            dim = 0;
        }

        else if (diag[1] >= std::max(diag[0], diag[2]))
        {
            dim = 1;
        }

        else
        {
            dim = 2;
        }

        double splitValue = 0.5 * (node->boundingBox.minPoint[dim] + node->boundingBox.maxPoint[dim]);

        int pivot = lowerIndex;

        for (int k = lowerIndex; k < upperIndex; k++)
        {
            TriangleIndices tri = indices[k];
            double midPointValue = (vertices[tri.vtxi][dim] + vertices[tri.vtxj][dim] + vertices[tri.vtxk][dim]) / 3;
            if (midPointValue < splitValue)
            {
                std::swap(indices[k], indices[pivot]);
                pivot++;
            }
        }

        if ((pivot <= lowerIndex) || (pivot >= upperIndex - 1) || (upperIndex - lowerIndex < 5))
            return;

        node->lc = new BVH;
        node->rc = new BVH;
        buildBVH(node->lc, lowerIndex, pivot + 1);
        buildBVH(node->rc, pivot + 1, upperIndex);
    }

    void initBVH()
    {
        buildBVH(&bvh, 0, indices.size());
    }

    Intersection intersectWithOneFace(const Ray &r, const TriangleIndices &triangle) const
    {
        // Points du triangle
        Vector verticeI = vertices[triangle.vtxi];
        Vector verticeJ = vertices[triangle.vtxj];
        Vector verticeK = vertices[triangle.vtxk];

        Vector e1 = verticeJ - verticeI;
        Vector e2 = verticeK - verticeI;

        Vector normalVector = cross(e1, e2);

        Vector toOrigin = (r.origin - verticeI);
        Vector tangentialVector = cross(toOrigin, r.direction);
        double denominator = dot(r.direction, normalVector);

        double beta = -dot(e2, tangentialVector) / denominator;
        double gamma = +dot(e1, tangentialVector) / denominator;
        double alpha = 1. - beta - gamma;

        double t = -dot(toOrigin, normalVector) / denominator;
        bool valid = (0 <= alpha && alpha <= 1 && 0 <= beta && beta <= 1 && 0 <= gamma && gamma <= 1 && t >= 0);

        if (!valid)
            return Intersection({false, verticeI, 0, normalVector, matter.albedo});

        Vector intersectionPoint = verticeI + e1 * beta + e2 * gamma;

        Vector ni = normals[triangle.ni];
        Vector nj = normals[triangle.nj];
        Vector nk = normals[triangle.nk];
        normalVector = alpha * ni + beta * nj + gamma * nk;
        normalVector.normalize();

        // Texture
        Vector albedo;

        if (useTextures)
        {
            // Case with one texture file or multiple files
            int group = (textures.size() == 1) ? 0 : triangle.group;
            int height = heigthList[group];
            int width = widthList[group];
            Vector uv = alpha * uvs[triangle.uvi] + beta * uvs[triangle.uvj] + gamma * uvs[triangle.uvk];
            int uvx = fabs(fmod(uv[0], 1.) * width);
            int uvy = height - fabs(fmod(uv[1], 1.) * height);
            albedo = Vector(
                std::pow(textures[group][(uvy * width + uvx) * 3 + 0] / 255., 2.2),
                std::pow(textures[group][(uvy * width + uvx) * 3 + 1] / 255., 2.2),
                std::pow(textures[group][(uvy * width + uvx) * 3 + 2] / 255., 2.2));
        }
        else
        {
            albedo = matter.albedo;
        }

        return Intersection{true, intersectionPoint, t, normalVector, albedo};
    }

    virtual Intersection intersect(const Ray &r) const
    {
        Intersection finalIntersection{false, Vector(0, 0, 0), 0, Vector(0, 0, 0), matter.albedo};
        double minimalTSoFar(std::numeric_limits<double>::max());

        if (!bvh.boundingBox.intersect(r).exists)
            return finalIntersection;

        std::list<const BVH *> nodes;
        nodes.push_back(&bvh);

        while (!nodes.empty())
        {
            const BVH *currentBVH = nodes.front();
            nodes.pop_front();

            if (currentBVH->lc)
            {
                if (currentBVH->lc->boundingBox.intersect(r).exists)
                {
                    nodes.push_front(currentBVH->lc);
                }

                if (currentBVH->rc->boundingBox.intersect(r).exists)
                {
                    nodes.push_front(currentBVH->rc);
                }
            }

            else
            {
                Intersection currentIntersection{false, Vector(0, 0, 0), 0, Vector(0, 0, 0), Vector(0, 0, 0)};

                for (int k = currentBVH->lowerIndex; k < currentBVH->upperIndex; k++)
                {
                    currentIntersection = intersectWithOneFace(r, indices[k]);

                    if (currentIntersection.exists && currentIntersection.t < minimalTSoFar)
                    {
                        finalIntersection = currentIntersection;
                        minimalTSoFar = currentIntersection.t;
                    }
                }
            }
        }
        return finalIntersection;
    }

    void move(double scale, Vector offset)
    {
        for (int k = 0; k < vertices.size(); k++)
        {
            vertices[k] = scale * vertices[k] + offset;
        }
    }

    void swapAxis(int axis1, int axis2)
    {
        for (int k = 0; k < normals.size(); k++)
        {
            std::swap(normals[k][axis1], normals[k][axis2]);
        }

        for (int k = 0; k < vertices.size(); k++)
        {
            std::swap(vertices[k][axis1], vertices[k][axis2]);
        }
    }

    void rotate(int axis, double angle)
    {
        int axis1 = (axis + 1) % 3;
        int axis2 = (axis + 2) % 3;

        angle = angle * M_PI / 180;

        double cosAngle = cos(angle);
        double sinAngle = sin(angle);

        double temp1;
        double temp2;

        for (int k = 0; k < normals.size(); k++)
        {
            temp1 = normals[k][axis1];
            temp2 = normals[k][axis2];

            normals[k][axis1] = temp1 * cosAngle + temp2 * sinAngle;
            normals[k][axis2] = temp1 * (-sinAngle) + temp2 * cosAngle;
        }

        for (int k = 0; k < vertices.size(); k++)
        {
            temp1 = vertices[k][axis1];
            temp2 = vertices[k][axis2];

            vertices[k][axis1] = temp1 * cosAngle + temp2 * sinAngle;
            vertices[k][axis2] = temp1 * (-sinAngle) + temp2 * cosAngle;
        }
    }

    std::vector<TriangleIndices> indices;
    std::vector<Vector> vertices;
    std::vector<Vector> normals;
    BVH bvh;

    bool useTextures;
    std::vector<unsigned char *> textures;
    std::vector<int> widthList;
    std::vector<int> heigthList;
    std::vector<Vector> uvs;
    std::vector<Vector> vertexcolors;
};

struct IntersectionWithScene
{
    bool exists;
    Vector point;
    double t;
    int objectNumber;
    Vector normalVector;
    Vector albedo;
};

Vector randomCos(const Vector &normalVector)
{
    int threadNumber = omp_get_thread_num();

    Vector firstTangentialVector;

    if (abs((normalVector[2]) < abs(normalVector[0])) && (abs(normalVector[2]) < abs(normalVector[1])))
    {
        firstTangentialVector = Vector(-normalVector[1], normalVector[0], 0.);
    }
    else if (abs((normalVector[1]) < abs(normalVector[0])) && (abs(normalVector[1]) < abs(normalVector[2])))
    {
        firstTangentialVector = Vector(normalVector[2], 0., -normalVector[0]);
    }
    else
    {
        firstTangentialVector = Vector(0., normalVector[2], -normalVector[1]);
    }
    firstTangentialVector.normalize();
    Vector secondTangentialVector = cross(normalVector, firstTangentialVector);

    // Génération des composantes aléatoires
    double u = uniform(rng[threadNumber]);
    double v = uniform(rng[threadNumber]);

    double sqrtV = sqrt(1 - v);
    double twoPiU = 2 * M_PI * u;

    double x = cos(twoPiU) * sqrtV;
    double y = sin(twoPiU) * sqrtV;
    double z = sqrt(v);

    // Génération de la direction du vecteur aléatoire (déjà normalisée par construction)
    return x * firstTangentialVector + y * secondTangentialVector + z * normalVector;
}

class Scene
{
public:
    explicit Scene() : n(1.){};

    void addObject(Object *object)
    {
        this->objects.push_back(object);
    }

    IntersectionWithScene intersect(const Ray &r) const
    {
        double minimalTSoFar(std::numeric_limits<double>::max());
        bool hasIntersection(false);
        Vector intersectionPoint;
        int objectNumber;
        Intersection currentIntersection;
        Vector normalVector;
        Vector albedo;

        for (int k = 0; k < objects.size(); k++)
        {
            currentIntersection = objects[k]->intersect(r);

            if (currentIntersection.exists && (currentIntersection.t < minimalTSoFar))
            {
                hasIntersection = true;
                minimalTSoFar = currentIntersection.t;
                intersectionPoint = currentIntersection.point;
                objectNumber = k;
                normalVector = currentIntersection.normalVector;
                albedo = currentIntersection.albedo;
            }
        }
        return IntersectionWithScene{hasIntersection, intersectionPoint, minimalTSoFar, objectNumber, normalVector, albedo};
    }

    Vector getColor(const Ray &r, int reflectionNumber, bool indirectLighting, bool softShadows, bool fresnel, bool showLight = false) const
    {
        Vector color(0, 0, 0);

        int threadNumber = omp_get_thread_num();

        if (reflectionNumber == maxReflectionNumber)
            return color;

        IntersectionWithScene currentIntersection;
        currentIntersection = this->intersect(r);

        if (currentIntersection.exists)
        {
            Vector intersectionPoint = currentIntersection.point;
            Object *intersectingObject = objects[currentIntersection.objectNumber];
            Vector normalVector = currentIntersection.normalVector;
            Vector intersectingAlbedo = currentIntersection.albedo;

            if (intersectingObject->intensity > epsilon && showLight)
            {
                return intersectingObject->intensity * light->matter.albedo;
            }

            if (intersectingObject->matter.mirror)
            {
                Vector reflectedDirection = r.direction + (-2) * dot(r.direction, normalVector) * normalVector;
                Ray reflectedRay(intersectionPoint + epsilon * normalVector, reflectedDirection);
                return getColor(reflectedRay, reflectionNumber + 1, indirectLighting, softShadows, fresnel, showLight);
            }

            else if (intersectingObject->matter.transparent)
            {
                double cosAngle = dot(r.direction, normalVector);
                double n1;
                double n2;

                if (cosAngle < 0) // On rentre dans la sphère
                {
                    n1 = this->n;
                    n2 = intersectingObject->matter.n;
                }

                else // On sort de la sphère
                {
                    n1 = intersectingObject->matter.n;
                    n2 = this->n;
                    normalVector = (-1) * normalVector;
                    cosAngle = -cosAngle;
                }

                double insideSquareRoot = 1 - pow(n1 / n2, 2) * (1 - pow(cosAngle, 2));

                double k0 = pow((n1 - n2) / (n1 + n2), 2);
                double R = k0 + (1 - k0) * pow(1 - abs(cosAngle), 5);
                bool fresnelCase = fresnel && (uniform(rng[threadNumber]) < R);

                if ((insideSquareRoot < 0) || fresnelCase)
                {
                    Vector reflectedDirection = r.direction + (-2) * dot(r.direction, normalVector) * normalVector;
                    Ray reflectedRay(intersectionPoint + epsilon * normalVector, reflectedDirection);
                    return getColor(reflectedRay, reflectionNumber + 1, indirectLighting, softShadows, fresnel);
                }

                else
                {
                    Vector normalComponent;
                    Vector tangentialComponent;
                    tangentialComponent = n1 / n2 * (r.direction - cosAngle * normalVector);
                    normalComponent = -sqrt(insideSquareRoot) * normalVector;

                    Vector refractedDirection = normalComponent + tangentialComponent;
                    Ray refractedRay(intersectionPoint - epsilon * normalVector, refractedDirection);
                    return getColor(refractedRay, reflectionNumber + 1, indirectLighting, softShadows, fresnel);
                }
            }

            else
            {
                if (softShadows)
                {
                    // Vecteur reliant l'intersection et la source
                    Vector goingFromLightVector = intersectionPoint - light->origin;
                    goingFromLightVector.normalize(); // Normalisé

                    Vector randomVector = randomCos(goingFromLightVector);
                    randomVector.normalize();

                    Vector randomPoint = light->origin + light->radius * randomVector; // Point aléatoire à la surface de la source de lumière
                    Vector randomDirection = randomPoint - intersectionPoint;
                    double dSquared = randomDirection.norm2();
                    randomDirection.normalize();

                    // Calcul de la fonction de visibilité

                    // Rayon partant de l'intersection vers la lumière, légèrement décollé pour pallier les problèmes de précision numérique
                    // En direction du point aléatoire à la surface de la sphère
                    Ray rayGoingToLight(intersectionPoint + epsilon * normalVector, randomDirection);

                    // Intersection avec la scène
                    IntersectionWithScene intersectionGoingToLight;
                    intersectionGoingToLight = this->intersect(rayGoingToLight);

                    /* On teste la présence d'une ombre : s'il y a une intersection,
                    est-elle plus proche que la source de lumière ? */
                    bool hasShadow;
                    // L'évaluation du && est paresseuse
                    hasShadow = intersectionGoingToLight.exists && (std::pow(intersectionGoingToLight.t, 2) < 0.99 * dSquared);
                    ;
                    double V(hasShadow ? 0. : 1.); // 0 si ombre, 1 sinon

                    double probability = std::max(0., dot(goingFromLightVector, randomVector)) / M_PI / (light->radius * light->radius);
                    double J = std::max(0., dot(randomVector, (-1) * randomDirection)) / dSquared;

                    Vector BRDF = currentIntersection.albedo / M_PI;
                    color = V * BRDF * light->intensity / (4 * M_PI * M_PI * light->radius * light->radius) * std::max(0., dot(randomDirection, normalVector)) * J / probability;
                }

                else
                {
                    // Vecteur reliant l'intersection et la source
                    Vector goingToLightVector = light->origin - intersectionPoint;
                    double dSquared = goingToLightVector.norm2();
                    goingToLightVector.normalize(); // Normalisé

                    // Calcul de la fonction de visibilité

                    // Rayon partant de l'intersection vers la lumière, légèrement décollé pour pallier les problèmes de précision numérique
                    Ray rayGoingToLight(intersectionPoint + epsilon * normalVector, goingToLightVector);

                    // Intersection avec la scène
                    IntersectionWithScene intersectionGoingToLight;
                    intersectionGoingToLight = this->intersect(rayGoingToLight);
                    /* On teste la présence d'une ombre : s'il y a une intersection,
                    est-elle plus proche que la source de lumière ? */
                    bool hasShadow;
                    // L'évaluation du && est paresseuse
                    hasShadow = intersectionGoingToLight.exists && (intersectionGoingToLight.t < 0.99 * (std::sqrt(dSquared) - light->radius));
                    double V(hasShadow ? 0. : 1.); // 0 si ombre, 1 sinon

                    Vector BRDF = currentIntersection.albedo / M_PI;
                    color = V * BRDF * this->light->intensity / (4 * M_PI * dSquared) * std::max(0., dot(goingToLightVector, normalVector));
                }

                // Eclairage indirect
                if (indirectLighting)
                {
                    // Génération de la direction du vecteur aléatoire (déjà normalisée par construction)
                    Vector randomDirection = randomCos(normalVector);
                    Ray randomRay(intersectionPoint + epsilon * normalVector, randomDirection);

                    color = color + currentIntersection.albedo * getColor(randomRay, reflectionNumber + 1, indirectLighting, softShadows, fresnel);
                }
                return color;
            }
        }

        // Cas sans intersection, color est par défaut à (0, 0, 0)
        return color;
    }

    std::vector<Object *> objects;
    Sphere *light;
    double lightIntensity;
    double n;
};

int main()
{
    init_rng(rng);

    // Configuration de la scène
    bool ANTI_ALIASING(false);
    bool INDIRECT_LIGHTING(false);
    bool SOFT_SHADOWS(false);
    bool FRESNEL(false);
    int MAX_RAYS_MONTE_CARLO = (ANTI_ALIASING || INDIRECT_LIGHTING || SOFT_SHADOWS || FRESNEL) ? 64 : 1; // Nombre de rayons par pixel

    double fov = 60 * M_PI / 180;
    double tanfov2 = tan(fov / 2);
    double stdAntiAliasing = 0.3;

    const int width = 512;
    const int height = 512;

    Vector originPoint(0, 0, 0);  // Position de l'image
    Vector cameraPoint(0, 0, 55); // Position de la caméra

    Vector originLight(-10, 20, 40); // Position de la source de lumière
    double radiusLight(8);
    double lightIntensity(5E9); // Intensité de la source de lumière
    Vector rhoLight(1, 1, 1);
    Matter matterLight({rhoLight, false, false, 1});
    Sphere sLight(originLight, radiusLight, matterLight, false, lightIntensity);

    // Sphère principale
    Vector rhoMain(0.1, 0.1, 1); // Albédo de la sphère cible
    Matter matterMain({rhoMain, false, false, 1.5});
    Vector originMain(0, -5, 15);
    double radiusMain(5); // Rayon de la sphère cible
    Sphere sMain(originMain, radiusMain, matterMain);

    // Sphère devant la caméra (magenta)
    Vector rhoFront(1, 0, 1); // Albédo de la boule magenta
    Matter matterFront({rhoFront, false, false, 1.3});
    Vector originFront(0, 0, -1000);
    double radiusFront(940); // Rayon de la boule magenta
    Sphere sFront(originFront, radiusFront, matterFront);

    // Sphère derrière la caméra (cyan)
    Vector rhoBack(0, 1, 1); // Albédo de la boule cyan
    Matter matterBack({rhoBack, false, false, 1.3});
    Vector originBack(0, 0, 1000);
    double radiusBack(940); // Rayon de la boule cyan
    Sphere sBack(originBack, radiusBack, matterBack);

    // Sphère en haut (blanche)
    Vector rhoUp(1, 1, 1); // Albédo de la boule blanche
    Matter matterUp({rhoUp, false, false, 1.3});
    Vector originUp(0, 1000, 0);
    double radiusUp(940); // Rayon de la boule blanche
    Sphere sUp(originUp, radiusUp, matterUp);

    // Sphère en bas (blanche)
    Vector rhoDown(1, 1, 1); // Albédo de la boule blanche
    Matter matterDown({rhoDown, false, false, 1.3});
    Vector originDown(0, -1000, 0);
    double radiusDown(990); // Rayon de la boule blanche
    Sphere sDown(originDown, radiusDown, matterDown);

    // Sphère à gauche (vert)
    Vector rhoLeft(0, 1, 0); // Albédo de la boule verte
    Matter matterLeft({rhoLeft, false, false, 1.3});
    Vector originLeft(-1000, 0, 0);
    double radiusLeft(940); // Rayon de la boule verte
    Sphere sLeft(originLeft, radiusLeft, matterLeft);

    // Sphère à droite (rouge)
    Vector rhoRight(1, 0, 0); // Albédo de la boule rouge
    Matter matterRight({rhoRight, false, false, 1.3});
    Vector originRight(1000, 0, 0);
    double radiusRight(940); // Rayon de la boule rouge
    Sphere sRight(originRight, radiusRight, matterRight);

    // Mesh
    Vector rhoMesh(1, 1, 1); // Albédo de la boule rouge
    Matter matterMesh({rhoMesh, false, false, 1.3});
    Vector originMesh(0, 0, 0);
    TriangleMesh mesh(originMesh, matterMesh);

    // Goodboi
    mesh.readOBJ("./mesh/dog/dog.obj");
    mesh.loadTexture("./mesh/dog/dog.jpg");
    mesh.rotate(0, 90);
    mesh.rotate(1, -45);
    mesh.move(1, Vector(0, -10, 0));

    std::cout << std::endl
              << "Building BVH..." << std::endl;
    mesh.initBVH();
    std::cout << "Done!" << std::endl;

    // Création de la scène
    Scene scene;
    scene = Scene();
    scene.addObject(&sLight);
    scene.light = &sLight;

    // scene.addObject(&sMain);
    scene.addObject(&sFront);
    scene.addObject(&sBack);
    scene.addObject(&sUp);
    scene.addObject(&sDown);
    scene.addObject(&sLeft);
    scene.addObject(&sRight);
    scene.addObject(&mesh);

    std::vector<unsigned char> image(width * height * 3, 0);

    auto start = std::chrono::high_resolution_clock::now();

#pragma omp parallel for schedule(dynamic, 1)
    for (int i = 0; i < height; i++)
    {

        for (int j = 0; j < width; j++)
        {
            Vector color(0, 0, 0);

            for (int rayNumber = 0; rayNumber < MAX_RAYS_MONTE_CARLO; rayNumber++)
            {
                double x = 0;
                double y = 0;

                // Anti-aliasing
                if (ANTI_ALIASING)
                {
                    int threadNumber = omp_get_thread_num();

                    double r1 = uniform(rng[threadNumber]);
                    double r2 = uniform(rng[threadNumber]);

                    double twoPiR = 2 * M_PI * r1;
                    double sqrtLogR = sqrt((-2) * log(r2)) * stdAntiAliasing;

                    x = cos(twoPiR) * sqrtLogR;
                    y = sin(twoPiR) * sqrtLogR;
                }

                // Direction du rayon lancé depuis la caméra vers le pixel en (i, j)
                Vector u(j - width / 2 + 0.5 + x, height / 2 - i + 0.5 + y, -width / (2 * tanfov2)); // Facteur 0.5 pour être au centre d'un pixel
                u.normalize();                                                                       // On normalise la direction
                Ray r(cameraPoint, u);                                                               // On crée le rayon

                color = color + scene.getColor(r, 0, INDIRECT_LIGHTING, SOFT_SHADOWS, FRESNEL, true);
            }

            color = color / MAX_RAYS_MONTE_CARLO;

            for (int k = 0; k < 3; k++)
            {
                image[(i * width + j) * 3 + k] = std::min(255., std::pow(color[k], 1. / 2.2));
            }
        }
    }

    stbi_write_png("output.png", width, height, 3, &image[0], 0);

    auto end = std::chrono::high_resolution_clock::now();
    auto diff = end - start;
    auto diff_sec = std::chrono::duration_cast<std::chrono::milliseconds>(diff);
    std::cout << std::endl
              << "Run time : " << diff_sec.count() << "ms (" << diff_sec.count() / 1000. << "s)" << std::endl;

    return 0;
}