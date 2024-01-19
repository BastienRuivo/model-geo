#pragma once

#include "vec.h"
#include <assert.h>

class Box2
{
protected:
    Point pMin;
    Point pMax;
public:
    Box2() : pMin(Point(0.0, 0.0, 0.0)), pMax(Point(0.0, 0.0, 0.0)) {}
    Box2(const Point & pMin, const Point & pMax) : pMin(pMin), pMax(pMax) {}
    Box2(const Point & pMin, float width, float height) : pMin(pMin), pMax(pMin + Point(width, height)) {}

    Point Min() const { return pMin; }
    Point Max() const { return pMax; }

    void setMin(const Point & pMin) { this->pMin = pMin; }
    void setMax(const Point & pMax) { this->pMax = pMax; }

    bool inside(const Point & p) const;
    bool intersect(const Box2 & box) const;

    Box2 operator+ (const Box2 & box) const;

    friend std::ostream & operator << (std::ostream & os, const Box2 & box);

    ~Box2() {}
};

