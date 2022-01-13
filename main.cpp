#define _CRT_SECURE_NO_WARNINGS 1
#include <vector>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

class Vector
{

public:
    explicit Vector(double x = 0, double y = 0, double z = 0)
    {
        data[0] = x;
        data[1] = y;
        data[2] = z;
    }
    double norm2() const
    {
        return data[0] * data[0] + data[1] * data[1] + data[2] * data[2];
    }
    double norm() const
    {
        return sqrt(norm2());
    }
    void normalize()
    {
        double n = norm();
        data[0] /= n;
        data[1] /= n;
        data[2] /= n;
    }
    double operator[](int i) const { return data[i]; };
    double &operator[](int i) { return data[i]; };
    double data[3];
};

Vector operator+(const Vector &a, const Vector &b)
{
    return Vector(a[0] + b[0], a[1] + b[1], a[2] + b[2]);
}
Vector operator-(const Vector &a, const Vector &b)
{
    return Vector(a[0] - b[0], a[1] - b[1], a[2] - b[2]);
}
Vector operator*(const double a, const Vector &b)
{
    return Vector(a * b[0], a * b[1], a * b[2]);
}
Vector operator*(const Vector &a, const double b)
{
    return Vector(a[0] * b, a[1] * b, a[2] * b);
}
Vector operator/(const Vector &a, const double b)
{
    return Vector(a[0] / b, a[1] / b, a[2] / b);
}
double dot(const Vector &a, const Vector &b)
{
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}
Vector cross(const Vector &a, const Vector &b)
{
    return Vector(a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]);
}

class Ray
{
public:
    explicit Ray(const Vector &origin, const Vector &direction) : C(origin), u(direction){};
    Vector C, u;
};

struct Intersection
{
    bool exists;
    Vector intersectionPoint;
    double t;
};

class Sphere
{
public:
    explicit Sphere(const Vector &origin, double radius, Vector &albedo) : O(origin), R(radius), rho(albedo){};

    Intersection intersect(const Ray &r) const
    {
        double a = 1;
        double b = 2 * dot(r.u, r.C - O);
        double c = dot(r.C - O, r.C - O) - R * R;

        double delta = b * b - 4 * a * c;

        Intersection intersection;
        Vector intersectionPoint;
        double t;

        // Cas sans racines réelles : pas d'intersection, on renvoie n'importe quoi
        if (delta < 0)
        {
            intersection = {false, intersectionPoint, t};
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
            intersection = {false, intersectionPoint, t};
            return intersection;
        }

        t = t1 > 0 ? t1 : t2;              // Soit t1 soit t2 en fonction de la positivité de t1
        intersectionPoint = r.C + t * r.u; // Le point associé
        intersection = {true, intersectionPoint, t};
        return intersection;
    }

    // Fonction qui calcule le vecteur normal étant donné un point
    Vector getNormalVector(Vector &point) const
    {
        Vector normalVector = (point - O);
        normalVector.normalize();
        return normalVector;
    }

    Vector O;
    double R;
    Vector rho;
};

struct IntersectionWithScene
{
    bool exists;
    Vector intersectionPoint;
    double t;
    int sphereNumber;
};

class Scene
{
public:
    explicit Scene(){};

    void addSphere(Sphere sphere)
    {
        spheres.push_back(sphere);
    }

    IntersectionWithScene intersect(const Ray &r)
    {
        double minimalTSoFar(std::numeric_limits<int>::max());
        bool hasIntersection(false);
        Vector intersectionPoint;
        int sphereNumber;
        Intersection currentIntersection;

        for (int k = 0; k < spheres.size(); k++)
        {
            currentIntersection = spheres[k].intersect(r);
            hasIntersection = hasIntersection || currentIntersection.exists;

            if (currentIntersection.exists && (currentIntersection.t < minimalTSoFar))
            {
                minimalTSoFar = currentIntersection.t;
                intersectionPoint = currentIntersection.intersectionPoint;
                sphereNumber = k;
            }
        }
        return IntersectionWithScene{hasIntersection, intersectionPoint, minimalTSoFar, sphereNumber};
    }

    std::vector<Sphere> spheres;
};

int main()
{
    int W = 1024;
    int H = 1024;

    double fov = 60 * M_PI / 180;
    double tanfov2 = tan(fov / 2);

    Vector O(0, 0, 0);              // Position de l'image
    Vector C(0, 0, 55);             // Position de la caméra
    Vector lightPoint(-10, 20, 40); // Position de la source de lumière
    double I(1000000);              // Intensité de la source de lumière

    // Création de la scène
    Scene scenery;
    scenery = Scene();

    // Sphère principale
    Vector rhoMain(0.8, 0.8, 0.8); // Albédo de la sphère cible
    double RMain(10);              // Rayon de la sphère cible
    Sphere sMain(O, RMain, rhoMain);

    // Sphère devant la caméra (verte)
    Vector rhoFront(0.5, 1., 0.5); // Albédo de la boule verte
    Vector originFront(0, 0, -1000);
    double RFront(940); // Rayon de la boule verte
    Sphere sFront(originFront, RFront, rhoFront);

    // Sphère derrière la caméra (magenta)
    Vector rhoBack(1., 0.5, 1.); // Albédo de la boule magenta
    Vector originBack(0, 0, 1000);
    double RBack(940); // Rayon de la boule magenta
    Sphere sBack(originBack, RBack, rhoBack);

    // Sphère en haut (rouge)
    Vector rhoUp(1., 0.5, 0.5); // Albédo de la boule rouge
    Vector originUp(0, 1000, 0);
    double RUp(940); // Rayon de la boule rouge
    Sphere sUp(originUp, RUp, rhoUp);

    // Sphère en bas (bleue)
    Vector rhoDown(0.5, 0.5, 1.); // Albédo de la boule bleue
    Vector originDown(0, -1000, 0);
    double RDown(990); // Rayon de la boule verte
    Sphere sDown(originDown, RDown, rhoDown);

    scenery.addSphere(sMain);
    scenery.addSphere(sFront);
    scenery.addSphere(sBack);
    scenery.addSphere(sUp);
    scenery.addSphere(sDown);

    Sphere intersectingSphere(sMain); // TODO: fix the bug to remove (sMain) from the constructor

    std::vector<unsigned char> image(W * H * 3, 0);
    for (int i = 0; i < H; i++)
    {
        for (int j = 0; j < W; j++)
        {
            Vector u(j - W / 2 + 0.5, H - i - H / 2 + 0.5, -W / (2 * tanfov2));
            u.normalize();
            Ray r(C, u);

            image[(i * W + j) * 3 + 0] = 0;
            image[(i * W + j) * 3 + 1] = 0;
            image[(i * W + j) * 3 + 2] = 0;

            IntersectionWithScene intersection;
            intersection = scenery.intersect(r);

            if (intersection.exists)
            {
                // On calcule le vecteur normal au point d'intersection
                intersectingSphere = scenery.spheres[intersection.sphereNumber];

                Vector intersectionPoint = intersection.intersectionPoint;
                Vector normalVector = intersectingSphere.getNormalVector(intersectionPoint); // Normalisé

                // Vecteur reliant l'intersection et la source
                Vector goingToLightVector = lightPoint - intersectionPoint;
                double dSquared = goingToLightVector.norm2();
                goingToLightVector.normalize(); // Normalisé

                // Calcul de la fonction de visibilité

                // Rayon partant de l'intersection vers la lumière
                Ray rayGoingToLight(intersectionPoint, goingToLightVector);

                // Intersection avec la scène
                IntersectionWithScene intersectionGoingToLight;
                intersectionGoingToLight = scenery.intersect(rayGoingToLight);

                /* On teste la présence d'une ombre : s'il y a une intersection, 
                est-elle plus proche que la source de lumière ? */
                bool hasShadow;
                // L'évaluation du && est paresseuse
                hasShadow = intersectionGoingToLight.exists && (pow(intersectionGoingToLight.t, 2) < dSquared);

                double visibility(hasShadow ? 0. : 1.); // 0 si ombre, 1 sinon

                Vector scale;
                scale = intersectingSphere.rho * I * std::max(0., dot(goingToLightVector, normalVector)) / dSquared / M_PI;
                scale = visibility * scale;

                for (int k = 0; k < 3; k++)
                {
                    image[(i * W + j) * 3 + k] = std::min(255., scale[k]);
                }
            }
        }
    }

    stbi_write_png("image.png", W, H, 3, &image[0], 0);
    return 0;
}