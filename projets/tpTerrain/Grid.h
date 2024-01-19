#pragma once

#include "Box2.h"
#include <vector>

class Grid : public Box2
{
protected:
    std::vector<double> values;
    int width;
    int height;
public:
    Grid() : Box2(), width(0), height(0) {}
    Grid(const Point & _min, const Point & _max, int gridWidth, int gridHeight) : Box2(_min, _max), width(gridWidth), height(gridHeight)
    {
        values.resize(width * height);
    }
    Grid(const Point & _min, float _width, float _height, int gridWidth, int gridHeight) : Box2(_min, _width, _height), width(gridWidth), height(gridHeight)
    {
        values.resize(gridWidth * gridHeight);
    }

    Grid(int width, int height) : Box2(), width(width), height(height)
    {
        values.resize(width * height);
    }

    Grid(const Grid & grid) : Box2(grid.pMin, grid.pMax), width(grid.width), height(grid.height) 
    {
        values.resize(width * height);
        for (int i = 0; i < width * height; ++i)
            values[i] = grid.values[i];
    }
    
    int getWidth() const { return width; }
    int getHeight() const { return height; }

    int index(int i, int j) const { 
        assert(i >= 0 && i < width && j >= 0 && j < height);
        return j * width + i; 
    }
    double operator[] (int i) const { 
        assert(i >= 0 && i < width * height);
        return values[i]; 
    }
    double & operator[] (int i) { 
        assert(i >= 0 && i < width * height);
        return values[i]; 
    }

    double operator() (int i, int j) const { 
        //std::cout << " i = " << i << " / " << width << " j = " << j << " / " << height << std::endl;
        if(i < 0) {
            std::cout<<"i = "<<i<<std::endl;
            assert(i >= 0);
        }
        if(i >= width) {
            std::cout<<"i = "<<i<<std::endl;
            assert(i < width);
        }
        if(j < 0) {
            std::cout<<"j = "<<j<<std::endl;
            assert(j >= 0);
        }
        if(j >= height) {
            std::cout<<"j = "<<j<<std::endl;
            assert(j < height);
        }
        return values[index(i, j)]; 
    }
    double & operator() (int i, int j) { 
        if(i < 0) {
            std::cout<<"i = "<<i<<std::endl;
            assert(i >= 0);
        }
        if(i >= width) {
            std::cout<<"i = "<<i<<std::endl;
            assert(i < width);
        }
        if(j < 0) {
            std::cout<<"j = "<<j<<std::endl;
            assert(j >= 0);
        }
        if(j >= height) {
            std::cout<<"j = "<<j<<std::endl;
            assert(j < height);
        }
        return values[index(i, j)]; 
    }

    double minValue() const;
    double maxValue() const;

    int minValueIndex() const;
    int maxValueIndex() const;

    void set(int i, int j, double value) { values[index(i, j)] = value; }

    void Convolute(const Grid & grid);
    void Normalize();
    void Clamp(double min, double max);

    void operator+= (const Grid & grid);
    void operator*= (const Grid & grid);
    void operator/= (const Grid & grid);
    void operator= (const Grid & grid);
    void operator+= (double value);
    void operator*= (double value);
    void operator/= (double value);
    void operator= (double value);

    // operator+
    Grid operator+ (const Grid & grid) const;
    Grid operator+ (double value) const;

    // operator-
    Grid operator- (const Grid & grid) const;
    Grid operator- (double value) const;

    // operator*
    Grid operator* (const Grid & grid) const;
    Grid operator* (double value) const;

    // operator/
    Grid operator/ (const Grid & grid) const;
    Grid operator/ (double value) const;

    // operator==
    bool operator== (const Grid & grid) const;
    bool operator== (double value) const;

    // operator!=
    bool operator!= (const Grid & grid) const;
    bool operator!= (double value) const;

    ~Grid();
};