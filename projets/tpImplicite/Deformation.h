#pragma once

#include "mesh.h"

class GlobalDeformation {
    public:
    class TorsionHelicoidal {
        public:
        static Point GetPoint(const Point & p, int axis, float T);
        static Mesh Warp(const Mesh& mesh, int axis, float T);
    };
};

class LocalDeformation
{
public:
    class Sphere
    {
    private:
    static float attenuation(const Point & p, const Point & center, float radius);
    public:
    static Point GetPoint(const Point & p, const Vector & dir, const Point & center, float radius);
    static Mesh Warp(const Mesh& mesh, const Vector & dir, const Point & center, float radius);
    };
    
};