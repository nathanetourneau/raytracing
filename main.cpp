#define _CRT_SECURE_NO_WARNINGS 1
#include <vector>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

#include <random>
#include <cmath>
#include <omp.h>

#define _USE_MATH_DEFINES
#include <math.h>

#include "Vector.h"

std::default_random_engine rng[8]; // Modifier par le nombre de coeurs de l'ordi
std::uniform_real_distribution<double> uniform(0., 1.);

const double epsilon(0.01); // Pour corriger les bugs liés à la précision numérique
const int maxReflectionNumber(5);

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
    double n;
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
        double t(0);

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

Vector random_cos(Vector &normalVector)
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
        Vector color(0, 0, 0);

        if (reflectionNumber == maxReflectionNumber)
            return color;

        IntersectionWithScene currentIntersection;
        currentIntersection = this->intersect(r);

        if (currentIntersection.exists)
        {
            Vector intersectionPoint = currentIntersection.point;
            Sphere intersectingSphere(spheres[currentIntersection.sphereNumber]);
            Vector normalVector = intersectingSphere.getNormalVector(intersectionPoint); // Normalisé

            if (intersectingSphere.matter.mirror)
            {
                Vector reflectedDirection = r.direction + (-2) * dot(r.direction, normalVector) * normalVector;
                Ray reflectedRay(intersectionPoint + epsilon * normalVector, reflectedDirection);
                return getColor(reflectedRay, reflectionNumber + 1);
            }

            else if (intersectingSphere.matter.transparent)
            {
                double cosAngle = dot(r.direction, normalVector);
                double n1;
                double n2;

                if (cosAngle < 0) // On rentre dans la sphère
                {
                    n1 = this->n;
                    n2 = intersectingSphere.matter.n;
                }

                else // On sort de la sphère
                {
                    n1 = intersectingSphere.matter.n;
                    n2 = this->n;
                    normalVector = (-1) * normalVector;
                    cosAngle = -cosAngle;
                }

                double insideSquareRoot = 1 - pow(n1 / n2, 2) * (1 - pow(cosAngle, 2));

                double k0 = pow((n1 - n2) / (n1 + n2), 2);
                double R = k0 + (1 - k0) * pow(1 - abs(cosAngle), 5);
                std::default_random_engine &engine = rng[omp_get_thread_num()];

                if ((insideSquareRoot < 0) || (uniform(engine) < R))
                {
                    Vector reflectedDirection = r.direction + (-2) * dot(r.direction, normalVector) * normalVector;
                    Ray reflectedRay(intersectionPoint + epsilon * normalVector, reflectedDirection);
                    return getColor(reflectedRay, reflectionNumber + 1);
                }

                else
                {
                    Vector normalComponent;
                    Vector tangentialComponent;
                    tangentialComponent = n1 / n2 * (r.direction - cosAngle * normalVector);
                    normalComponent = -sqrt(insideSquareRoot) * normalVector;

                    Vector refractedDirection = normalComponent + tangentialComponent;
                    Ray refractedRay(intersectionPoint - epsilon * normalVector, refractedDirection);
                    return getColor(refractedRay, reflectionNumber + 1);
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

                color = visibility * intersectingSphere.matter.albedo * this->light.intensity * std::max(0., dot(goingToLightVector, normalVector)) / dSquared / M_PI;

                // Eclairage indirect

                // Génération de la direction du vecteur aléatoire (déjà normalisée par construction)
                Vector randomDirection = random_cos(normalVector);
                Ray randomRay(intersectionPoint + epsilon * normalVector, randomDirection);

                color = color + intersectingSphere.matter.albedo * getColor(randomRay, reflectionNumber + 1);

                return color;
            }
        }

        // Cas sans intersection, color est par défaut à (0, 0, 0)
        return color;
    }

    std::vector<Sphere> spheres;
    LightSource light;
    double n;
};

int main()
{
    // Configuration de la scène

    int maxRaysForMonteCarlo(128); // Nombre de rayons par pixel

    double fov = 60 * M_PI / 180;
    double tanfov2 = tan(fov / 2);

    const int width = 1024;
    const int height = 1024;

    Vector originPoint(0, 0, 0);  // Position de l'image
    Vector cameraPoint(0, 0, 55); // Position de la caméra

    Vector lightPoint(-10, 20, 40);    // Position de la source de lumière
    double lightIntensity(1000000000); // Intensité de la source de lumière
    LightSource lightSource({lightIntensity, lightPoint});

    // Sphère principale
    Vector rhoMain(1, 1, 1); // Albédo de la sphère cible
    Matter matterMain({rhoMain, false, true, 1.5});
    Vector originMain(10, 0, 15);
    double radiusMain(10); // Rayon de la sphère cible
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

    // Création de la scène
    Scene scene;
    scene = Scene();
    scene.addLight(lightSource);

    scene.addSphere(sMain);
    scene.addSphere(sFront);
    scene.addSphere(sBack);
    scene.addSphere(sUp);
    scene.addSphere(sDown);
    scene.addSphere(sLeft);
    scene.addSphere(sRight);

    Sphere intersectingSphere(sMain); // TODO: fix the bug to remove (sMain) from the constructor

    std::vector<unsigned char> image(width * height * 3, 0);

#pragma omp parallel for
    for (int i = 0; i < height; i++)
    {
        for (int j = 0; j < width; j++)
        {

            std::default_random_engine &engine = rng[omp_get_thread_num()];

            Vector color(0, 0, 0);

            for (int rayNumber = 0; rayNumber < maxRaysForMonteCarlo; rayNumber++)
            {
                // Anti-aliasing
                double r1 = uniform(engine);
                double r2 = uniform(engine);

                double twoPiR = 2 * M_PI * r1;
                double sqrtLogR = sqrt((-2) * log(r2));

                double x = cos(twoPiR) * sqrtLogR;
                double y = sin(twoPiR) * sqrtLogR;

                // Direction du rayon lancé depuis la caméra vers le pixel en (i, j)
                Vector u(j - width / 2 + 0.5 + x, height / 2 - i + 0.5 + y, -width / (2 * tanfov2)); // Facteur 0.5 pour être au centre d'un pixel
                u.normalize();                                                                       // On normalise la direction
                Ray r(cameraPoint, u);                                                               // On crée le rayon

                color = color + scene.getColor(r, 0);
            }

            color = color / maxRaysForMonteCarlo;

            for (int k = 0; k < 3; k++)
            {
                image[(i * width + j) * 3 + k] = std::min(255., std::pow(color[k], 1. / 2.2));
            }
        }
    }

    stbi_write_png("output.png", width, height, 3, &image[0], 0);

    return 0;
}