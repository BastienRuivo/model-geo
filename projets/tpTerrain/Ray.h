#pragma once

#include <vec.h>


struct Ray {
    Point origin;
    Vector direction;
    Ray() {}
    Ray(const Point & origin, const Vector & direction) : origin(origin), direction(direction) {}
    Point operator()(float t) const {
        return origin + t * direction;
    }
};