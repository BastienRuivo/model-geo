#include "Box2.h"

bool Box2::inside(const Point & p) const
{
    //printf("p =(%f, %f, %f) = %d\n", p.x, p.y, p.z, p.x >= pMin.x && p.x <= pMax.x && p.z >= pMin.y && p.z <= pMax.z);
    return (p.x >= pMin.x && p.x < pMax.x && p.z >= pMin.y && p.z < pMax.z && p.y >= pMin.y && p.y < pMax.y);
}

bool Box2::intersect(const Box2 & box) const
{
    return (pMin.x <= box.pMax.x && pMax.x >= box.pMin.x && pMin.y <= box.pMax.y && pMax.y >= box.pMin.y);
}

Box2 Box2::operator+ (const Box2 & box) const
{
    return Box2(min(pMin, box.pMin), max(pMax, box.pMax));
}

std::ostream & operator << (std::ostream & os, const Box2 & box)
{
    os << "Box2 {min: " << box.Min() << " max: " << box.Max() << "}";
    return os;
}