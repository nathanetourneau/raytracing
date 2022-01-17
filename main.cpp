#define _CRT_SECURE_NO_WARNINGS 1
#include <vector>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

#include "Vector.h"

const double epsilon(0.001); // Pour corriger les bugs liés à la précision numérique

const int maxReflectionNumber(5);

const int W = 1024;
const int H = 1024;

class Ray
{
public:
    explicit Ray(const Vector &origin, const Vector &direction) : origin(origin), direction(direction){};
    Vector origin, direction;
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

struct Intersection
{
    bool exists;
    Vector point;
    double t;
};

class Sphere
{
public:
    explicit Sphere(const Vector &origin, const double radius, const Matter &matter) : origin(origin), radius(radius), matter(matter){};

    Intersection intersect(const Ray &r) const
    {
        double a = 1;
        double b = 2 * dot(r.direction, r.origin - origin);
        double c = dot(r.origin - origin, r.origin - origin) - radius * radius;

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

        t = t1 > 0 ? t1 : t2;                           // Soit t1 soit t2 en fonction de la positivité de t1
        intersectionPoint = r.origin + t * r.direction; // Le point associé
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
    Vector point;
    double t;
    int sphereNumber;
};

class Scene
{
public:
    explicit Scene(){};

    void addSphere(const Sphere &sphere)
    {
        this->spheres.push_back(sphere);
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
        int sphereNumber;
        Intersection currentIntersection;

        for (int k = 0; k < spheres.size(); k++)
        {
            currentIntersection = spheres[k].intersect(r);

            if (currentIntersection.exists && (currentIntersection.t < minimalTSoFar))
            {
                hasIntersection = true;
                minimalTSoFar = currentIntersection.t;
                intersectionPoint = currentIntersection.point;
                sphereNumber = k;
            }
        }
        return IntersectionWithScene{hasIntersection, intersectionPoint, minimalTSoFar, sphereNumber};
    }

    Vector getColor(Ray r, int reflectionNumber) const
    {
        if (reflectionNumber == 0)
            return Vector(0, 0, 0);
        Vector color(0, 0, 0);

        IntersectionWithScene currentIntersection;
        currentIntersection = this->intersect(r);
        if (currentIntersection.exists)
        {
            Sphere intersectingSphere(spheres[currentIntersection.sphereNumber]);

            if (intersectingSphere.matter.mirror)
            {
                Vector normalVector = intersectingSphere.getNormalVector(currentIntersection.point);
                Vector reflectedDirection = r.direction + (-2) * dot(r.direction, normalVector) * normalVector;
                Ray reflectedRay(r.origin + epsilon * normalVector, reflectedDirection);
                color = getColor(reflectedRay, reflectionNumber - 1);
            }

            else
            {
                Vector intersectionPoint = currentIntersection.point;
                Vector normalVector = intersectingSphere.getNormalVector(intersectionPoint); // Normalisé

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

                Vector scale;
                color = visibility * intersectingSphere.matter.albedo * this->light.intensity * std::max(0., dot(goingToLightVector, normalVector)) / dSquared / M_PI;
            }
        }
        return color;
    }

    std::vector<Sphere> spheres;

    LightSource light;
};

int main()
{
    Vector originPoint(0, 0, 0);  // Position de l'image
    Vector cameraPoint(0, 0, 55); // Position de la caméra

    Vector lightPoint(-10, 20, 40);   // Position de la source de lumière
    double lightIntensity(600000000); // Intensité de la source de lumière
    LightSource lightSource({lightIntensity, lightPoint});

    // Sphère principale
    Vector rhoMain(1, 1, 1); // Albédo de la sphère cible
    Matter matterMain({rhoMain, true, false});
    Vector originMain(15, 5, 0);
    double radiusMain(10); // Rayon de la sphère cible
    Sphere sMain(originMain, radiusMain, matterMain);

    // Sphère secondaire
    Vector rhoBis(1, 1, 1); // Albédo de la sphère cible
    Matter matterBis({rhoMain, true, false});
    Vector originBis(-15, 5, 0);
    double radiusBis(10); // Rayon de la sphère cible
    Sphere sBis(originBis, radiusBis, matterBis);

    // Sphère devant la caméra (verte)
    Vector rhoFront(0, 1, 0); // Albédo de la boule verte
    Matter matterFront({rhoFront, false, false});
    Vector originFront(0, 0, -1000);
    double radiusFront(940); // Rayon de la boule verte
    Sphere sFront(originFront, radiusFront, matterFront);

    // Sphère derrière la caméra (magenta)
    Vector rhoBack(1, 0, 1); // Albédo de la boule magenta
    Matter matterBack({rhoBack, false, false});
    Vector originBack(0, 0, 1000);
    double radiusBack(940); // Rayon de la boule magenta
    Sphere sBack(originBack, radiusBack, matterBack);

    // Sphère en haut (rouge)
    Vector rhoUp(1, 0, 0); // Albédo de la boule rouge
    Matter matterUp({rhoUp, false, false});
    Vector originUp(0, 1000, 0);
    double radiusUp(940); // Rayon de la boule rouge
    Sphere sUp(originUp, radiusUp, matterUp);

    // Sphère en bas (bleue)
    Vector rhoDown(0, 0, 1); // Albédo de la boule bleue
    Matter matterDown({rhoDown, false, false});
    Vector originDown(0, -1000, 0);
    double radiusDown(990); // Rayon de la boule verte
    Sphere sDown(originDown, radiusDown, matterDown);

    // Création de la scène
    Scene scene;
    scene = Scene();
    scene.addLight(lightSource);

    scene.addSphere(sMain);
    scene.addSphere(sBis);
    scene.addSphere(sFront);
    //scene.addSphere(sBack);
    scene.addSphere(sUp);
    scene.addSphere(sDown);

    Sphere intersectingSphere(sMain); // TODO: fix the bug to remove (sMain) from the constructor

    double fov = 60 * M_PI / 180;
    double tanfov2 = tan(fov / 2);

    std::vector<unsigned char> image(W * H * 3, 0);
    for (int i = 0; i < H; i++)
    {
        for (int j = 0; j < W; j++)
        {
            // Direction du rayon lancé depuis la caméra vers le pixel en (i, j)
            Vector u(j - W / 2 + 0.5, H - i - H / 2 + 0.5, -W / (2 * tanfov2)); // Facteur 0.5 pour être au centre d'un pixel
            u.normalize();                                                      // On normalise la direction
            Ray r(cameraPoint, u);                                              // On crée le rayon

            // On cherche ensuite si le rayon intersecte la scène
            IntersectionWithScene intersection;

            Vector color;
            color = scene.getColor(r, maxReflectionNumber);

            for (int k = 0; k < 3; k++)
            {
                image[(i * W + j) * 3 + k] = std::min(255., std::pow(color[k], 1. / 2.2));
            }
        }
    }

    stbi_write_png("output.png", W, H, 3, &image[0], 0);
    return 0;
}