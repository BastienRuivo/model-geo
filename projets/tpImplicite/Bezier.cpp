#include "Bezier.h"
#include <cmath>
#include <assert.h>

int Bezier::binomialCoefficient(int n, int k) const
{
    if (k == 0 || k == n) {
        return 1;
    }
    return binomialCoefficient(n - 1, k - 1) + binomialCoefficient(n - 1, k);
}

Point Bezier::getPoint(float u, float v) const
{
    Point p(0, 0, 0);
    int n = controlPoints.size() - 1;
    int m = controlPoints[0].size() - 1;

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            double coeff = bCoeffs[n][i] * bCoeffs[m][j] * pow(u, i) * pow(1 - u, n - i) * pow(v, j) * pow(1 - v, m - j);
            p = p + controlPoints[i][j] * coeff;
        }
    }
    return p;
}

Point Bezier::getPointCurve(float t) const
{
    assert(controlPoints.size() == 1);

    std::vector<Point> points(controlPoints[0].begin(), controlPoints[0].end());

    while (points.size() > 1)
    {
        std::vector<Point> newPoints(points.size() - 1);
        for (size_t i = 0; i < points.size() - 1; i++) {
            newPoints[i] = (1 - t) * points[i] + t * points[i + 1];
        }
        points = newPoints;
    }
    
    return points[0];
}

void Bezier::getLinesMesh(std::vector<Mesh> & lines) const {
    for (size_t i = 0; i < controlPoints.size(); i++)
    {
        for (size_t j = 0; j < controlPoints[i].size() - 1; j++)
        {
            assert(j < controlPoints[i].size() - 1);
            Mesh m = Mesh(GL_LINES);

            Color c = Color::random();
            m.color(c);
            m.vertex(controlPoints[i][j]);
            m.color(c);
            m.vertex(controlPoints[i][j+1]);
            lines.push_back(m);
        }
    }
}


Mesh Bezier::polygonize(int resolution) const
{
    Mesh m = Mesh(GL_TRIANGLES);
    std::cout<<"resolution = " << resolution <<std::endl;
    float step = 1.0f / (resolution - 1);
    for (int i = 0; i < resolution; i++)
    {
        float u = i * step;
        for (int j = 0; j < resolution; j++)
        {
            float v = j * step;
            Point p = getPoint(u, v);
            m.vertex(p);


            if(j > 0 && i > 0) {
                // triangle i-1, j-1, to i, j
                int prevL = (i - 1) * resolution;
                int curL = i * resolution;
                std::cout<<"prevL = " << prevL << " curL = " << curL << " i = " << i << " j = " << j << " vertex_count = " << m.vertex_count() << std::endl;
                m.triangle(prevL + j - 1, prevL + j, curL + j - 1);
                m.triangle(curL + j - 1, prevL + j, curL + j);
            }
        }
        std::cout<<m.vertex_count()<<std::endl;

    }
    return m;
}

Mesh Revolution::polygonize(int resolution) const
{
    Mesh m = Mesh(GL_TRIANGLES);

    float step = 1.0f / (resolution);

    for (size_t i = 0; i < resolution; i++)
    {
        Point po = bezier.getPointCurve(i * step);
        Point poi = bezier.getPointCurve((i + 1) * step);
        for (size_t j = 0; j < resolution; j++)
        {
            // rotate around axis

            float a0 = j * step * 360.0f;

            Transform t0 = Rotation(axis, a0);

            Point p_i_j = t0(po);
            Point p_i_j1 = t0(poi);

            m.vertex(p_i_j);
            m.vertex(p_i_j1);

            if(j > 0) {
                int s = m.vertex_count() - 4;
                //m.triangle(s, s + 1, s + 2);
                //m.triangle(s + 2, s + 1, s + 3);
                m.triangle(s, s + 2, s + 1);
                m.triangle(s + 2, s + 3, s + 1);
            }

            if(j == resolution - 1) {
                int s = m.vertex_count() - 2;
                int s0 = i * resolution * 2;
                //m.triangle(s, s + 1, s0);
                //m.triangle(s0, s + 1, s0 + 1); 
                m.triangle(s, s0, s + 1);
                m.triangle(s0, s0 + 1, s + 1); 
            }
        }
    }

    return m;
}