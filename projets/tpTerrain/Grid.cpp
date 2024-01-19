#include "Grid.h"

void Grid::Convolute(const Grid & grid)
{
    int kernelWidth = grid.width;
    int kernelHeight = grid.height;
    for (int i = 0; i < width * height; ++i) {
        int x = i % width;
        int y = i / width;
        double sum = 0;
        for (int j = 0; j < kernelWidth * kernelHeight; ++j) {
            int kx = j % kernelWidth;
            int ky = j / kernelWidth;
            int x2 = x + kx - kernelWidth / 2;
            int y2 = y + ky - kernelHeight / 2;
            if (x2 >= 0 && x2 < width && y2 >= 0 && y2 < height)
                sum += values[y2 * width + x2] * grid.values[j];
        }
        values[i] = sum;
    }
}

double Grid::minValue() const {
    double min = values[0];
    for (size_t i = 1; i < width * height; i++)
    {
        if (values[i] < min)
            min = values[i];
    }
    return min;
}

double Grid::maxValue() const {
    double max = values[0];
    for (size_t i = 1; i < width * height; i++)
    {
        if (values[i] > max)
            max = values[i];
    }
    return max;
}

int Grid::minValueIndex() const {
    double min = values[0];
    int index = 0;
    for (size_t i = 1; i < width * height; i++)
    {
        if (values[i] < min) {
            min = values[i];
            index = i;
        }
    }
    return index;
}

int Grid::maxValueIndex() const {
    double max = values[0];
    int index = 0;
    for (size_t i = 1; i < width * height; i++)
    {
        if (values[i] > max) {
            max = values[i];
            index = i;
        }
    }
    return index;
}


void Grid::Normalize() {
    *this /= maxValue();
}

void Grid::Clamp(double min, double max) {
    for (size_t i = 0; i < width * height; i++)
    {
        if (values[i] < min)
            values[i] = min;
        else if (values[i] > max)
            values[i] = max;
    }
}

void Grid::operator+= (const Grid & grid)
{
    assert(width == grid.width && height == grid.height);
    for (int i = 0; i < width * height; ++i)
        values[i] += grid.values[i];
}

void Grid::operator*= (const Grid & grid)
{
    assert(width == grid.width && height == grid.height);
    for (int i = 0; i < width * height; ++i)
        values[i] *= grid.values[i];
}

void Grid::operator/= (const Grid & grid)
{
    assert(width == grid.width && height == grid.height);
    for (int i = 0; i < width * height; ++i)
        values[i] /= grid.values[i];
}

void Grid::operator= (const Grid & grid)
{
    height = grid.height;
    width = grid.width;
    values.resize(width * height);
    setMin(grid.Min());
    setMax(grid.Max());
    for (int i = 0; i < width * height; ++i)
        values[i] = grid.values[i];
}

void Grid::operator+= (double value)
{
    for (int i = 0; i < width * height; ++i)
        values[i] += value;
}

void Grid::operator*= (double value)
{
    for (int i = 0; i < width * height; ++i)
        values[i] *= value;
}

void Grid::operator/= (double value)
{
    for (int i = 0; i < width * height; ++i)
        values[i] /= value;
}

void Grid::operator= (double value)
{
    for (int i = 0; i < width * height; ++i)
        values[i] = value;
}

// operator+

Grid Grid::operator+ (const Grid & grid) const
{
    assert(width == grid.width && height == grid.height);
    Grid result(*this);
    result += grid;
    return result;
}

Grid Grid::operator+ (double value) const
{
    Grid result(*this);
    result += value;
    return result;
}

// operator-

Grid Grid::operator- (const Grid & grid) const
{
    assert(width == grid.width && height == grid.height);
    Grid result(*this);
    for (int i = 0; i < width * height; ++i)
        result.values[i] -= grid.values[i];
    return result;
}

Grid Grid::operator- (double value) const
{
    Grid result(*this);
    for (int i = 0; i < width * height; ++i)
        result.values[i] -= value;
    return result;
}

// operator*

Grid Grid::operator* (const Grid & grid) const
{   
    assert(width == grid.width && height == grid.height);
    Grid result(*this);
    for (int i = 0; i < width * height; ++i)
        result.values[i] *= grid.values[i];
    return result;
}

Grid Grid::operator* (double value) const
{
    Grid result(*this);
    for (int i = 0; i < width * height; ++i)
        result.values[i] *= value;
    return result;
}

// operator/

Grid Grid::operator/ (const Grid & grid) const
{
    assert(width == grid.width && height == grid.height);
    Grid result(*this);
    for (int i = 0; i < width * height; ++i)
        result.values[i] /= grid.values[i];
    return result;
}

Grid Grid::operator/ (double value) const
{
    Grid result(*this);
    for (int i = 0; i < width * height; ++i)
        result.values[i] /= value;
    return result;
}

// operator==

bool Grid::operator== (const Grid & grid) const
{
    assert(width == grid.width && height == grid.height);
    for (int i = 0; i < width * height; ++i)
        if (values[i] != grid.values[i])
            return false;
    return true;
}

bool Grid::operator== (double value) const
{
    for (int i = 0; i < width * height; ++i)
        if (values[i] != value)
            return false;
    return true;
}

// operator!=

bool Grid::operator!= (const Grid & grid) const
{
    assert(width == grid.width && height == grid.height);
    for (int i = 0; i < width * height; ++i)
        if (values[i] != grid.values[i])
            return true;
    return false;
}

bool Grid::operator!= (double value) const
{
    for (int i = 0; i < width * height; ++i)
        if (values[i] != value)
            return true;
    return false;
}

Grid::~Grid()
{
}
