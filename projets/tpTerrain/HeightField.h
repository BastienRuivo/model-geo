#pragma once

#include "ScalarField.h"
#include "mesh.h"
#include "Ray.h"
#include <chrono>

const double epsilon = 0.01;
static std::vector<Point> pointsOnHemisphere = {Point(0.991311, 0.131538, 4.87444e-05), Point(0.0312899, 0.45865, -0.888066), Point(-0.955127, 0.218959, -0.199471), Point(0.702419, 0.678865, 0.213894), Point(-0.152766, 0.934693, -0.320954), Point(-0.635639, 0.519416, 0.571113), Point(0.486769, 0.0345721, -0.872846), Point(0.800781, 0.5297, 0.279585), Point(-0.475399, 0.00769819, -0.879737), Point(-0.741828, 0.0668422, 0.667251), Point(-0.631352, 0.686773, 0.360191), Point(-0.310663, 0.930436, -0.194362), Point(0.482861, 0.526929, -0.699422), Point(0.633735, 0.653919, 0.413243), Point(-0.615953, 0.701191, 0.359073), Point(0.547269, 0.762198, -0.345762), Point(-0.0780766, 0.0474645, 0.995817), Point(-0.0824997, 0.328234, -0.940987), Point(-0.439794, 0.75641, -0.484174), Point(0.929399, 0.365339, -0.0523934), Point(0.00346031, 0.98255, 0.185965), Point(-0.11241, 0.753356, -0.647934), Point(-0.578505, 0.0726859, -0.812434), Point(-0.315593, 0.884707, -0.343066), Point(-0.127951, 0.436411, 0.890603), Point(0.0908851, 0.477732, -0.873792), Point(0.0737833, 0.274907, 0.958636), Point(-0.62501, 0.166507, 0.762652), Point(-0.439116, 0.897656, 0.0372884), Point(0.840107, 0.0605643, -0.53903), Point(0.71304, 0.504523, -0.486857), Point(-0.942782, 0.319033, -0.096847), Point(0.866415, 0.493977, -0.0728893), Point(-0.100847, 0.0907329, 0.990756), Point(0.944044, 0.0737491, -0.321468), Point(-0.923265, 0.384142, -0.00410189), Point(-0.0687732, 0.913817, 0.40026), Point(-0.870177, 0.464446, -0.164564), Point(0.930856, 0.050084, -0.361938), Point(0.046102, 0.770205, -0.636129), Point(0.465985, 0.125365, -0.875866), Point(0.721677, 0.688455, 0.0721903), Point(0.525599, 0.629543, -0.572207), Point(-0.0595021, 0.725412, -0.685738), Point(0.458734, 0.888572, -0.00156261), Point(0.100327, 0.306322, 0.946626), Point(-0.508871, 0.513274, 0.691087), Point(-0.448196, 0.845982, -0.288852), Point(-0.459886, 0.841511, 0.283487), Point(-0.110136, 0.415395, 0.902949), Point(-0.859607, 0.467917, -0.205254), Point(-0.227974, 0.178328, 0.957198), Point(0.466629, 0.571655, 0.674884), Point(0.32318, 0.0330538, -0.94576), Point(-0.846672, 0.49848, -0.186183), Point(0.637447, 0.748293, -0.183629), Point(-0.428049, 0.890737, -0.152844), Point(-0.381785, 0.84204, -0.381063), Point(0.524764, 0.212752, 0.824232), Point(-0.218043, 0.130427, -0.967185), Point(0.808652, 0.274588, 0.520272), Point(0.909982, 0.414293, 0.0171525), Point(0.694364, 0.70982, 0.118384), Point(0.897822, 0.239911, -0.369267)};
class HeightField : public ScalarField
{
    protected:

    public:
    Mesh mesh;
    std::vector<Mesh> lines;
    HeightField() : ScalarField() {}
    HeightField(const Point & _min, const Point & _max, int gridWidth, int gridHeight) : ScalarField(_min, _max, gridWidth, gridHeight) {
    }
    HeightField(const Point & _min, float width, float height, int gridWidth, int gridHeight) : ScalarField(_min, width, height, gridWidth, gridHeight) {
    }  
    float laplacian(int i, int j) const;
    double Height(const Point & p) const;
    double Height(double x, double z) const;
    ImageData Slope(int sampleWidth, int sampleHeight) const;
    ImageData Laplacian(int sampleWidth, int sampleHeight) const;
    ImageData GradientMagnitude(int sampleWidth, int sampleHeight) const;
    ImageData NormalColor(int sampleWidth, int sampleHeight) const;
    double AverageSlope(int x, int z) const;
    float skyViewFixed(const Point & p);
    float accesibility(const Point & p);

    bool spheretrace(const Ray & ray, float & t) const;

    // TODO : IMPLEMENT INTERSECTION RAY / HEIGHTFIELD WITH SPHERE TRACING

    /// TODO : IMPLEMENT ACCESIBILITY 
    /// sampling n points on an hemisphere at a given position with raycast
    ImageData Accessible(int sampleWidth, int sampleHeight);
    /// TODO : IMPLEMENT LAPLACIAN

    Point Vertex(int i, int j, int w, int h) const;
    Vector Normal(int i, int j, int w, int h) const;

    void polygonize(int resolutionX, int resolutionY) {
        mesh = Mesh(GL_TRIANGLES);
        for (int i = 0; i < resolutionX; i++)
        {
            for (int j = 0; j < resolutionY; j++)
            {
                mesh.normal(0, 0, 1);
                mesh.color(1, 1, 1);
                mesh.texcoord((float)i / (float)resolutionX, (float)j / (float)resolutionY);
                mesh.vertex(Vertex(i, j, resolutionX, resolutionY));
                mesh.normal(mesh.vertex_count() - 1, Normal(i, j, resolutionX, resolutionY));
                //mesh.normal(mesh.vertex_count() - 1, Vector(0, 1, 0));
                //auto pos = mesh.positions()[mesh.vertex_count() - 1];
                //mesh.color(mesh.vertex_count() - 1, Color(pos.z, pos.z, pos.z));

                if(i > 0 && j > 0 && j < resolutionY - 1) {
                    mesh.triangle(i * resolutionY + j - 1, (i - 1) * resolutionY + j - 1, (i - 1) * resolutionY + j);
                    mesh.triangle((i - 1) * resolutionY + j, i * resolutionY + j, i * resolutionY + j - 1);
                }
            }
        }
    }

};

