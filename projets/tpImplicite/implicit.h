#pragma once

#include "vec.h"
#include "mesh.h"
#include "box.h"

#include <algorithm>
#include <omp.h>
#include <string>
#include <fstream>
#include <iostream>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/export.hpp>
#include <boost/serialization/version.hpp>


#include "FastNoiseLite.h"

const float epsilon = 0.0001f;

// TEMPLATES
class Node {
public:
    friend class boost::serialization::access;
    Node() {}
    virtual float distance(const Vector& p) const = 0;
    Vector grad(const Vector& p) const;
    bool isInside(const Vector& p) const;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version) {
        if(version == 0) {
            Box bbox;
            ar & bbox;
        }
        //delete bbox;
    }
    virtual ~Node() {}
};

class BinaryOperator : public Node {
public:
    Node* left;
    Node* right;

    BinaryOperator(Node* left, Node* right) : left(left), right(right) {}   
    BinaryOperator() : left(nullptr), right(nullptr) {}
    void merge(Node* l, Node * r);

    template<class Archive>
    void serialize(Archive & ar, const unsigned int version) {
        ar & boost::serialization::base_object<Node>(*this);
        ar & left;
        ar & right;
    }

    ~BinaryOperator() {
        delete left;
        delete right;
    }
};

class Tree : public Node {
public:
    Node* root;

    Box SDFBox;

    Tree(Node* root) : root(root) {
        SDFBox.a = Vector(-2.5, -2.5, -2.5);
        SDFBox.b = Vector(2.5, 2.5, 2.5);
    }

    Tree(const std::string & filename) {
        SDFBox.a = Vector(-2.5, -2.5, -2.5);
        SDFBox.b = Vector(2.5, 2.5, 2.5);
        load(filename);
    }
    Tree() : root(nullptr) {}

    float distance(const Vector& p) const override;
    void Polygonize(int n, Mesh & m, const Box &b) const;
    Vector Dichotomy(Vector a, Vector b, double va, double vb, double length) const;

    bool intersectRayMarching(const Point& p, const Vector& d, Vector & hitPoint, float e = 0.0001) const;
    bool intersectRayMarchingBBOX(const Point& p, const Vector& d, Vector & hitPoint, float e = 0.0001) const;

    bool intersectSphereTracing(const Point& p, const Vector& d, Vector & hitPoint, float e = 0.00001) const;
    bool intersectSphereTracingBBOX(const Point& p, const Vector& d, Vector & hitPoint, float e = 0.00001) const;

    int countPrimitives() {
        int res = 0;
        std::vector<Node*> nodes;
        nodes.push_back(root);
        while(!nodes.empty()) {
            Node* n = nodes.back();
            nodes.pop_back();
            if(auto u = dynamic_cast<BinaryOperator*>(n)) {
                nodes.push_back(u->left);
                nodes.push_back(u->right);
            }
            else {
                res++;
            }
        }
        return res;
    }

    template<class Archive>
    void serialize(Archive & ar, const unsigned int version) {
        ar & boost::serialization::base_object<Node>(*this); // serialize base class
        ar & root; // serialize the root pointer
        if(version >= 1) {
            ar & SDFBox;
        }
    }

    void save(const std::string & filename) const {
        std::ofstream file(filename);
        if (file.is_open()) {
            boost::archive::text_oarchive oa(file);
            oa << *this;
        }
        else {
            std::cout << "Unable to open file " << filename << std::endl;
        }
    }

    void load(const std::string & filename) {
        std::cout<<"Loading "<<filename<<std::endl;
        std::ifstream file(filename);
        delete root; // delete the old root
        if (file.is_open()) {
            boost::archive::text_iarchive ia(file);
            ia >> *this;

        }
        else {
            std::cout << "Unable to open file " << filename << std::endl;
        }
    }

    ~Tree() {
        delete root;
    }
};

// OPERATORS

class Union : public BinaryOperator {
public:
    Union(Node* left, Node* right) : BinaryOperator(left, right) {}
    Union() : Union(nullptr, nullptr) {}

    float distance(const Vector& p) const override;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version) {
        ar & boost::serialization::base_object<BinaryOperator>(*this); // serialize base class
    }
};

class UnionSmooth : public BinaryOperator {
    float smooth;
public:
    UnionSmooth(Node* left, Node* right, float smooth = 1.f) : BinaryOperator(left, right), smooth(smooth) {}
    UnionSmooth() : UnionSmooth(nullptr, nullptr, 1.f) {}

    // float distance(const Vector& p) const override {
    //     return distance(p, false);
    // }

    float distance(const Vector& p) const override;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version) {
        ar & boost::serialization::base_object<BinaryOperator>(*this); // serialize base class
        ar & smooth;
    }
};

class Difference : public BinaryOperator {
public:
    Difference(Node* left, Node* right) : BinaryOperator(left, right) {}
    Difference() : Difference(nullptr, nullptr) {}

    float distance(const Vector& p) const override;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version) {
        ar & boost::serialization::base_object<BinaryOperator>(*this); // serialize base class
    }
};

class DifferenceSmooth : public BinaryOperator {
    float smooth;
public:
    DifferenceSmooth(Node* left, Node* right, float smooth = 1.f) : BinaryOperator(left, right), smooth(smooth) {
    }
    
    DifferenceSmooth() : DifferenceSmooth(nullptr, nullptr, 1.f) {}

    float distance(const Vector& p) const override;template<class Archive>
    void serialize(Archive & ar, const unsigned int version) {
        ar & boost::serialization::base_object<BinaryOperator>(*this); // serialize base class
        ar & smooth;
    }
};

class Intersection : public BinaryOperator {
public:
    Intersection(Node* left, Node* right) : BinaryOperator(left, right) {}
    Intersection() : Intersection(nullptr, nullptr) {}

    float distance(const Vector& p) const override;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version) {
        ar & boost::serialization::base_object<BinaryOperator>(*this); // serialize base class
    }
};

class IntersectionSmooth : public BinaryOperator {
    float smooth;
public:
    IntersectionSmooth(Node* left, Node* right, float smooth = 1.f) : BinaryOperator(left, right), smooth(smooth) {}
    IntersectionSmooth() : IntersectionSmooth(nullptr, nullptr, 1.f) {}

    float distance(const Vector& p) const override;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version) {
        ar & boost::serialization::base_object<BinaryOperator>(*this); // serialize base class
        ar & smooth;
    }
};

// PRIMITIVES

class Sphere : public Node {
public:
    float radius;
    Vector center;

    Sphere(const Vector& center, float radius) : center(center), radius(radius) {}
    Sphere() : center(Vector(0, 0, 0)), radius(1) {}

    float distance(const Vector& p) const override;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version) {
        ar & boost::serialization::base_object<Node>(*this); // serialize base class
        ar & center;
        ar & radius;
    }
};

class Cube : public Node {
public:
    float size;
    Vector center;

    Cube(const Vector& center, float size) : center(center), size(size) {}
    Cube() : Cube(Vector(0, 0, 0), 1) {}

    float distance(const Vector& p) const override;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version) {
        ar & boost::serialization::base_object<Node>(*this); // serialize base class
        ar & center;
        ar & size;
    }
};

class Tor : public Node {
public:
    Vector center;
    float r;
    float R;

    Tor(const Vector& center, float r, float R) : center(center), r(r), R(R) {}
    Tor() : Tor(Vector(0, 0, 0), 0.5, 1) {}
    float distance(const Vector& p) const override;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version) {
        ar & boost::serialization::base_object<Node>(*this); // serialize base class
        ar & center;
        ar & r;
        ar & R;
    }
};

class Cylinder : public Node {
public:
    Vector center;
    float radius;
    float height;

    Cylinder(const Vector& center, float radius) : center(center), radius(radius) {}
    Cylinder() : Cylinder(Vector(0, 0, 0), 1) {}
    float distance(const Vector& p) const override;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version) {
        ar & boost::serialization::base_object<Node>(*this); // serialize base class
        ar & center;
        ar & radius;
        ar & height;
    }
};

class Capsule : public Node {
public:
    Vector center;
    float radius;
    float height;

    Capsule(const Vector& center, float radius, float height) : center(center), radius(radius), height(height) {}
    Capsule() : Capsule(Vector(0, 0, 0), 1, 1) {}
    float distance(const Vector& p) const override;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version) {
        ar & boost::serialization::base_object<Node>(*this); // serialize base class
        ar & center;
        ar & radius;
        ar & height;
    }
};
