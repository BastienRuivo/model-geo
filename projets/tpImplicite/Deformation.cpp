#include "Deformation.h"
#include <cmath>

#include "mat.h"


Point GlobalDeformation::TorsionHelicoidal::GetPoint(const Point & p, int axis, const float T) {
    Transform r = axis == 0 ? 
        RotationX(360 * p(axis) / T) : 
        axis == 1 ? RotationY(360 * p(axis) / T) : RotationZ(360 * p(axis) / T);
    return r(p);
}

Mesh GlobalDeformation::TorsionHelicoidal::Warp(const Mesh& mesh, int axis, float T) {
    Mesh res = Mesh(GL_TRIANGLES);
    assert(axis >= 0 && axis <= 2);
    assert(T > 0);
    auto positions = mesh.positions();

    for (size_t i = 0; i < positions.size(); i++)
    {
        res.vertex(GetPoint(positions[i], axis, T));
    }

    auto index = mesh.indices();

    for (size_t i = 0; i < index.size(); i += 3)
    {
        res.triangle(index[i], index[i + 1], index[i + 2]);
    }

    return res;
}

float LocalDeformation::Sphere::attenuation(const Point & p, const Point & center, float radius) {
    Vector PC = p - center;
    float d = length(PC);
    if(d < radius) {
        return 1 - d / radius;
    } else {
        return 0;
    }    
}

Point LocalDeformation::Sphere::GetPoint(const Point & p, const Vector & dir, const Point & center, float radius) {
    return p + dir * attenuation(p, center, radius);
}

Mesh LocalDeformation::Sphere::Warp(const Mesh& mesh, const Vector & dir, const Point & center, float radius) {
    Mesh res = Mesh(GL_TRIANGLES);
    auto positions = mesh.positions();

    for (size_t i = 0; i < positions.size(); i++)
    {
        res.vertex(GetPoint(positions[i], dir, center, radius));
    }

    auto index = mesh.indices();

    for (size_t i = 0; i < index.size(); i += 3)
    {
        res.triangle(index[i], index[i + 1], index[i + 2]);
    }

    return res;
}