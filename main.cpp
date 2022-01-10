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
    Ray(const Vector &origin, const Vector &direction) : C(origin), u(direction){};
    Vector C, u;
};

class Sphere
{
public:
    Sphere(const Vector &origin, double radius) : O(origin), R(radius){};

    bool intersect(const Ray &r) const
    {
        double a = 1;
        double b = 2 * dot(r.u, r.C - O);
        double c = dot(r.C - O, r.C - O) - R * R;
        double delta = b * b - 4 * a * c;

        bool intersection;

        // Pas de racines réelles : pas d'intersection
        if (delta < 0)
            return false;

        /* 
		A partir de ici, delta >= 0

		On calcule les deux valeurs de l'intersection. Si t est négatif, 
		l'intersection obtenue est incorrecte : on s'intérese à la demi-droite
		t >= 0. On va donc effectuer un premier test : si t2 < 0, il n'y a
		aucune intersection.

		Ensuite, si t2 >= 0 mais t1 < 0, l'intersection valable est celle
		fournie par t2.

		Si t1 >= 0, et t2 >= 0, l'intersection valable est celle qui est la
		plus proche de la caméra, c'est à dire celle fournie par t1.

		*/

        double t1 = (-b - sqrt(delta)) / (2 * a);
        double t2 = (-b + sqrt(delta)) / (2 * a);

        // Cas t1 < t2 < 0 : les deux intersections sont derrière la caméra
        if (t2 < 0)
            return false;

        else
            return true;
    }

    // Fonction qui calcule le vecteur normal étant donné un point
    Vector getNormalVector(const Vector &point) const
    {
        Vector normalVector = (point - O);
        normalVector.normalize();
        return normalVector;
    }

    Vector O;
    double R;
};

int main()
{
    int W = 1024;
    int H = 1024;

    Vector O(0, 0, 0), C(0, 0, 55);
    double R(10);

    double fov = 60 * M_PI / 180;
    double tanfov2 = tan(fov / 2);

    Sphere s(O, R);

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

            if (s.intersect(r))
            {
                for (int k = 0; k < 3; k++)
                {
                    image[(i * W + j) * 3 + k] = 255;
                }
            }
        }
    }
    stbi_write_png("image.png", W, H, 3, &image[0], 0);

    return 0;
}