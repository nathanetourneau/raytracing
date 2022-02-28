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

#define _USE_MATH_DEFINES
#include <math.h>

#include "Vector.h"
#include "progressbar.h"

// TODO : fix to automatically adjust for the number of CPU cores
std::default_random_engine rng[8]; // Modifier par le nombre de coeurs de l'ordi
std::uniform_real_distribution<double> uniform(0., 1.);

const double epsilon(0.01); // Pour corriger les bugs liés à la précision numérique
const int maxReflectionNumber(5);

struct LightSource
{
    double intensity;
    Vector origin;
};

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

        bool hasInter = (tMin < tMax);
        // On ne s'intéresse qu'au premier booléen, les autres champs sont remplis avec n'importe quoi
        return Intersection({hasInter, r.origin, 0, r.direction});
    };

    Vector minPoint;
    Vector maxPoint;
};

class Sphere : public Object
{
public:
    explicit Sphere(const Vector &origin, const double radius, const Matter &matter, const bool reverseNormal = false)
    {
        this->origin = origin;
        this->radius = radius;
        this->matter = matter;
        this->reverseNormal = reverseNormal;
    };

    // Fonction qui calcule le vecteur normal étant donné un point
    Vector getNormalVector(const Vector &point) const
    {
        Vector normalVector = (point - origin);
        normalVector.normalize();
        return reverseNormal ? (-1) * normalVector : normalVector;
    }

    double radius;
    bool reverseNormal;

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
            intersection = {false, intersectionPoint, t, Vector(0, 0, 0)};
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
            intersection = {false, intersectionPoint, t, Vector(0, 0, 0)};
            return intersection;
        }

        t = t1 > 0 ? t1 : t2;                           // Soit t1 soit t2 en fonction de la positivité de t1
        intersectionPoint = r.origin + t * r.direction; // Le point associé
        intersection = {true, intersectionPoint, t, this->getNormalVector(intersectionPoint)};
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
        this->boundingBox = BBox();
    };

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

    void buildBB()
    {
        double xMin, xMax, yMin, yMax, zMin, zMax;

        for (int k = 0; k < vertices.size(); k++)
        {
            xMin = std::min(xMin, vertices[k][0]);
            xMax = std::max(xMax, vertices[k][0]);
            yMin = std::min(yMin, vertices[k][1]);
            yMax = std::max(yMax, vertices[k][1]);
            zMin = std::min(zMin, vertices[k][2]);
            zMax = std::max(zMax, vertices[k][2]);
        }

        Vector minPoint(xMin, yMin, zMin);
        Vector maxPoint(xMax, yMax, zMax);

        this->boundingBox = BBox(minPoint, maxPoint);
    };

    Intersection intersectWithOneFace(const Ray &r, const Vector &A, const Vector &B, const Vector &C, const Vector &ni, const Vector &nj, const Vector &nk) const
    {
        Vector e1 = B - A;
        Vector e2 = C - A;
        Vector AO = r.origin - A;
        Vector normalVector = cross(e1, e2);
        double denominator = dot(r.direction, normalVector);
        Vector tangentialVector = cross(AO, r.direction);

        double beta = -dot(e2, tangentialVector) / denominator;
        double gamma = dot(e1, tangentialVector) / denominator;
        double alpha = 1 - beta - gamma;
        double t = (-1) * dot(AO, normalVector) / denominator;

        bool hasInter = (alpha >= 0 && alpha <= 1 && beta >= 0 && beta <= 1 && gamma >= 0 && gamma <= 1);

        if ((t < 0) || (!hasInter))
        {
            return Intersection({false, A, t, normalVector});
        }

        Vector intersectionPoint = r.origin + t * r.direction;
        normalVector = alpha * ni + beta * nj + gamma * nk;
        normalVector.normalize();
        return Intersection{true, intersectionPoint, t, normalVector};
    }

    virtual Intersection intersect(const Ray &r) const
    {
        if (!this->boundingBox.intersect(r).exists)
        {
            return Intersection({false, r.origin, 0, r.direction});
        }

        double minimalTSoFar(std::numeric_limits<float>::max());

        Intersection currentIntersection{false, Vector(0, 0, 0), 0, Vector(0, 0, 0)};
        Intersection finalIntersection{false, Vector(0, 0, 0), 0, Vector(0, 0, 0)};

        for (int k = 0; k < indices.size(); k++)
        {
            TriangleIndices triangle = indices[k];
            currentIntersection = intersectWithOneFace(
                r,
                vertices[triangle.vtxi],
                vertices[triangle.vtxj],
                vertices[triangle.vtxk],
                normals[triangle.vtxi],
                normals[triangle.vtxj],
                normals[triangle.vtxk]);

            if ((currentIntersection.exists) && (currentIntersection.t < minimalTSoFar))
            {
                finalIntersection = currentIntersection;
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

    void invertNormals()
    {
        for (int k = 0; k < normals.size(); k++)
        {
            normals[k] = (-1) * normals[k];
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

    std::vector<TriangleIndices> indices;
    std::vector<Vector> vertices;
    std::vector<Vector> normals;
    std::vector<Vector> uvs;
    std::vector<Vector> vertexcolors;
    BBox boundingBox;
};

struct IntersectionWithScene
{
    bool exists;
    Vector point;
    double t;
    int objectNumber;
    Vector normalVector;
};

Vector randomCos(const Vector &normalVector)
{
    std::default_random_engine &engine = rng[omp_get_thread_num()];

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
    double u = uniform(engine);
    double v = uniform(engine);

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

    void addLight(const LightSource &light)
    {
        this->light = light;
    }

    IntersectionWithScene intersect(const Ray &r) const
    {
        double minimalTSoFar(std::numeric_limits<float>::max());
        bool hasIntersection(false);
        Vector intersectionPoint;
        int objectNumber;
        Intersection currentIntersection;
        Vector normalVector;

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
            }
        }
        return IntersectionWithScene{hasIntersection, intersectionPoint, minimalTSoFar, objectNumber, normalVector};
    }

    Vector getColor(const Ray &r, int reflectionNumber, bool indirectLighting, bool fresnel) const
    {
        Vector color(0, 0, 0);

        if (reflectionNumber == maxReflectionNumber)
            return color;

        IntersectionWithScene currentIntersection;
        currentIntersection = this->intersect(r);

        if (currentIntersection.exists)
        {
            Vector intersectionPoint = currentIntersection.point;
            Object *intersectingObject = objects[currentIntersection.objectNumber];
            Vector normalVector = currentIntersection.normalVector;

            if (intersectingObject->matter.mirror)
            {
                Vector reflectedDirection = r.direction + (-2) * dot(r.direction, normalVector) * normalVector;
                Ray reflectedRay(intersectionPoint + epsilon * normalVector, reflectedDirection);
                return getColor(reflectedRay, reflectionNumber + 1, indirectLighting, fresnel);
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
                std::default_random_engine &engine = rng[omp_get_thread_num()];
                bool fresnelCase = fresnel && (uniform(engine) < R);

                if ((insideSquareRoot < 0) || fresnelCase)
                {
                    Vector reflectedDirection = r.direction + (-2) * dot(r.direction, normalVector) * normalVector;
                    Ray reflectedRay(intersectionPoint + epsilon * normalVector, reflectedDirection);
                    return getColor(reflectedRay, reflectionNumber + 1, indirectLighting, fresnel);
                }

                else
                {
                    Vector normalComponent;
                    Vector tangentialComponent;
                    tangentialComponent = n1 / n2 * (r.direction - cosAngle * normalVector);
                    normalComponent = -sqrt(insideSquareRoot) * normalVector;

                    Vector refractedDirection = normalComponent + tangentialComponent;
                    Ray refractedRay(intersectionPoint - epsilon * normalVector, refractedDirection);
                    return getColor(refractedRay, reflectionNumber + 1, indirectLighting, fresnel);
                }
            }

            else
            {
                // Vecteur reliant l'intersection et la source
                Vector goingToLightVector = light.origin - intersectionPoint;
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
                hasShadow = intersectionGoingToLight.exists && (std::pow(intersectionGoingToLight.t, 2) < dSquared);

                double visibility(hasShadow ? 0. : 1.); // 0 si ombre, 1 sinon

                color = visibility * intersectingObject->matter.albedo * this->light.intensity * std::max(0., dot(goingToLightVector, normalVector)) / dSquared / M_PI;

                // Eclairage indirect
                if (indirectLighting)
                {
                    // Génération de la direction du vecteur aléatoire (déjà normalisée par construction)
                    Vector randomDirection = randomCos(normalVector);
                    Ray randomRay(intersectionPoint + epsilon * normalVector, randomDirection);

                    color = color + intersectingObject->matter.albedo * getColor(randomRay, reflectionNumber + 1, indirectLighting, fresnel);
                }
                return color;
            }
        }

        // Cas sans intersection, color est par défaut à (0, 0, 0)
        return color;
    }

    std::vector<Object *> objects;
    LightSource light;
    double n;
};

int main()
{
    // Configuration de la scène
    bool ANTI_ALIASING(false);
    bool INDIRECT_LIGHTING(false);
    bool FRESNEL(false);
    int MAX_RAYS_MONTE_CARLO = (ANTI_ALIASING || INDIRECT_LIGHTING || FRESNEL) ? 16 : 1; // Nombre de rayons par pixel

    double fov = 60 * M_PI / 180;
    double tanfov2 = tan(fov / 2);
    double stdAntiAliasing = 0.3;

    const int width = 1024;
    const int height = 1024;

    Vector originPoint(0, 0, 0);  // Position de l'image
    Vector cameraPoint(0, 0, 55); // Position de la caméra

    Vector lightPoint(-10, 20, 40);    // Position de la source de lumière
    double lightIntensity(3000000000); // Intensité de la source de lumière
    LightSource lightSource({lightIntensity, lightPoint});

    // Sphère principale
    Vector rhoMain(1, 1, 1); // Albédo de la sphère cible
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

    // Triangle (rouge)
    Vector rhoTriangle(1, 0, 0); // Albédo de la boule rouge
    Matter matterTriangle({rhoTriangle, false, false, 1.3});
    Vector originTriangle(0, 0, 0);
    TriangleMesh triangle(originTriangle, matterTriangle);
    triangle.readOBJ("beautifulgirl.obj");
    triangle.swapAxis(1, 2);
    triangle.move(20, Vector(0, -10, 0));
    triangle.invertNormals();
    triangle.buildBB();

    // Création de la scène
    Scene scene;
    scene = Scene();
    scene.addLight(lightSource);

    // scene.addObject(&sMain);
    scene.addObject(&sFront);
    scene.addObject(&sBack);
    scene.addObject(&sUp);
    scene.addObject(&sDown);
    scene.addObject(&sLeft);
    scene.addObject(&sRight);
    scene.addObject(&triangle);

    std::vector<unsigned char> image(width * height * 3, 0);

    ProgressBar pg;
    pg.start(height);

    auto start = std::chrono::high_resolution_clock::now();

#pragma omp parallel for
    for (int i = 0; i < height; i++)
    {
        for (int j = 0; j < width; j++)
        {

            std::default_random_engine &engine = rng[omp_get_thread_num()];

            Vector color(0, 0, 0);

            for (int rayNumber = 0; rayNumber < MAX_RAYS_MONTE_CARLO; rayNumber++)
            {
                double x = 0;
                double y = 0;

                // Anti-aliasing
                if (ANTI_ALIASING)
                {
                    double r1 = uniform(engine);
                    double r2 = uniform(engine);

                    double twoPiR = 2 * M_PI * r1;
                    double sqrtLogR = sqrt((-2) * log(r2)) * stdAntiAliasing;

                    x = cos(twoPiR) * sqrtLogR;
                    y = sin(twoPiR) * sqrtLogR;
                }

                // Direction du rayon lancé depuis la caméra vers le pixel en (i, j)
                Vector u(j - width / 2 + 0.5 + x, height / 2 - i + 0.5 + y, -width / (2 * tanfov2)); // Facteur 0.5 pour être au centre d'un pixel
                u.normalize();                                                                       // On normalise la direction
                Ray r(cameraPoint, u);                                                               // On crée le rayon

                color = color + scene.getColor(r, 0, INDIRECT_LIGHTING, FRESNEL);
            }

            color = color / MAX_RAYS_MONTE_CARLO;

            for (int k = 0; k < 3; k++)
            {
                image[(i * width + j) * 3 + k] = std::min(255., std::pow(color[k], 1. / 2.2));
            }
        }
        pg.update(i);
    }

    stbi_write_png("output.png", width, height, 3, &image[0], 0);

    auto end = std::chrono::high_resolution_clock::now();
    auto diff = end - start;
    auto diff_sec = std::chrono::duration_cast<std::chrono::milliseconds>(diff);
    std::cout << std::endl
              << "Run time : " << diff_sec.count() << "ms (" << diff_sec.count() / 1000. << "s)" << std::endl;

    return 0;
}