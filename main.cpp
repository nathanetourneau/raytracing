#define _CRT_SECURE_NO_WARNINGS 1
#include <vector>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

#include "Vector.h"

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
};

class Sphere
{
public:
    explicit Sphere(const Vector &origin, const double radius, const Matter &matter) : origin(origin), radius(radius), matter(matter){};

    Intersection intersect(const Ray &r) const
    {
        double a = 1;
        double b = 2 * dot(r.u, r.C - origin);
        double c = dot(r.C - origin, r.C - origin) - radius * radius;

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
        Vector normalVector = (point - origin);
        normalVector.normalize();
        return normalVector;
    }

    Vector origin;
    double radius;
    Matter matter;
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

    void addSphere(const Sphere &sphere)
    {
        spheres.push_back(sphere);
    }

    IntersectionWithScene intersect(const Ray &r) const
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

    LightSource light;
};

int main()
{
    int W = 1024;
    int H = 1024;

    double fov = 60 * M_PI / 180;
    double tanfov2 = tan(fov / 2);

    Vector originPoint(0, 0, 0);  // Position de l'image
    Vector cameraPoint(0, 0, 55); // Position de la caméra

    Vector lightPoint(-10, 20, 40);   // Position de la source de lumière
    double lightIntensity(100000000); // Intensité de la source de lumière
    LightSource lightSource({lightIntensity, lightPoint});

    const double epsilon(0.001); // Pour corriger les bugs liés à la précision numérique

    // Création de la scène
    Scene scene;
    scene = Scene();
    scene.light = lightSource;

    // Sphère principale
    Vector rhoMain(0.8, 0.8, 0.8); // Albédo de la sphère cible
    Matter matterMain({rhoMain, false, false});
    double radiusMain(10); // Rayon de la sphère cible
    Sphere sMain(originPoint, radiusMain, matterMain);

    // Sphère devant la caméra (verte)
    Vector rhoFront(0.5, 1., 0.5); // Albédo de la boule verte
    Matter matterFront({rhoFront, false, false});
    Vector originFront(0, 0, -1000);
    double radiusFront(940); // Rayon de la boule verte
    Sphere sFront(originFront, radiusFront, matterFront);

    // Sphère derrière la caméra (magenta)
    Vector rhoBack(1., 0.5, 1.); // Albédo de la boule magenta
    Matter matterBack({rhoBack, false, false});
    Vector originBack(0, 0, 1000);
    double radiusBack(940); // Rayon de la boule magenta
    Sphere sBack(originBack, radiusBack, matterBack);

    // Sphère en haut (rouge)
    Vector rhoUp(1., 0.5, 0.5); // Albédo de la boule rouge
    Matter matterUp({rhoUp, false, false});
    Vector originUp(0, 1000, 0);
    double radiusUp(940); // Rayon de la boule rouge
    Sphere sUp(originUp, radiusUp, matterUp);

    // Sphère en bas (bleue)
    Vector rhoDown(0.5, 0.5, 1.); // Albédo de la boule bleue
    Matter matterDown({rhoDown, false, false});
    Vector originDown(0, -1000, 0);
    double radiusDown(990); // Rayon de la boule verte
    Sphere sDown(originDown, radiusDown, matterDown);

    scene.addSphere(sMain);
    scene.addSphere(sFront);
    scene.addSphere(sBack);
    scene.addSphere(sUp);
    scene.addSphere(sDown);

    Sphere intersectingSphere(sMain); // TODO: fix the bug to remove (sMain) from the constructor

    std::vector<unsigned char> image(W * H * 3, 0);
    for (int i = 0; i < H; i++)
    {
        for (int j = 0; j < W; j++)
        {
            Vector u(j - W / 2 + 0.5, H - i - H / 2 + 0.5, -W / (2 * tanfov2));
            u.normalize();
            Ray r(cameraPoint, u);

            image[(i * W + j) * 3 + 0] = 0;
            image[(i * W + j) * 3 + 1] = 0;
            image[(i * W + j) * 3 + 2] = 0;

            IntersectionWithScene intersection;
            intersection = scene.intersect(r);

            if (intersection.exists)
            {
                // On calcule le vecteur normal au point d'intersection
                intersectingSphere = scene.spheres[intersection.sphereNumber];

                Vector intersectionPoint = intersection.intersectionPoint;
                Vector normalVector = intersectingSphere.getNormalVector(intersectionPoint); // Normalisé

                // Vecteur reliant l'intersection et la source
                Vector goingToLightVector = lightPoint - intersectionPoint;
                double dSquared = goingToLightVector.norm2();
                goingToLightVector.normalize(); // Normalisé

                // Calcul de la fonction de visibilité

                // Rayon partant de l'intersection vers la lumière, légèrement décollé pour pallier les problèmes de précision numérique
                Ray rayGoingToLight(intersectionPoint + epsilon * goingToLightVector, goingToLightVector);

                // Intersection avec la scène
                IntersectionWithScene intersectionGoingToLight;
                intersectionGoingToLight = scene.intersect(rayGoingToLight);

                /* On teste la présence d'une ombre : s'il y a une intersection, 
                est-elle plus proche que la source de lumière ? */
                bool hasShadow;
                // L'évaluation du && est paresseuse
                hasShadow = intersectionGoingToLight.exists && (std::pow(intersectionGoingToLight.t, 2) < dSquared);

                double visibility(hasShadow ? 0. : 1.); // 0 si ombre, 1 sinon

                Vector scale;
                scale = visibility * intersectingSphere.matter.albedo * scene.light.intensity * std::max(0., dot(goingToLightVector, normalVector)) / dSquared / M_PI;

                for (int k = 0; k < 3; k++)
                {
                    image[(i * W + j) * 3 + k] = std::min(255., std::pow(scale[k], 1. / 2.2));
                }
            }
        }
    }

    stbi_write_png("output.png", W, H, 3, &image[0], 0);
    return 0;
}