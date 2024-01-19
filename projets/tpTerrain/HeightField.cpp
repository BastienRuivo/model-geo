#include "HeightField.h"
#include <algorithm>

double HeightField::Height(const Point & p) const {
    // p to local coordinates
    Vector p2;
    p2.x = (p.x - Min().x) / (Max().x - Min().x);
    p2.z = (p.z - Min().z) / (Max().z - Min().z);

    float cp = (p.x - Min().x) / (Max().x - Min().x);

    // cell coordinates
    int x = std::floor(p2.x * width);
    int z = std::floor(p2.z * height);

    x = std::max(0, std::min(x, width - 1));
    z = std::max(0, std::min(z, height - 1));

    // local cell coordinates
    p2.x = p2.x*width - x;
    p2.z = p2.z*height - z;

    int x1 = std::min(x + 1, width - 1);
    int z1 = std::min(z + 1, height - 1);

    if(x < 0 || z < 0 || x >= width || z >= height) {
        printf("x = %d z = %d\n", x, z);
    }


    // bilinear interpolation
    if(p2.x + p2.z < 1.0) {
        return ((1.0 - p2.x - p2.z) * (*this)(x, z) + p2.x * (*this)(x1, z) + p2.z * (*this)(x, z1)) * pMax.y;
    } else {
        return ((p2.x + p2.z - 1) * (*this)(x1, z1) + (1 - p2.z) * (*this)(x1, z) + (1 - p2.x) * (*this)(x, z1)) * pMax.y;
    }
    
}


double HeightField::Height(double x, double z) const {
    return Height(Point(x, 0.0, z));
}

ImageData HeightField::NormalColor(int sampleWidth, int sampleHeight) const {
    ImageData img(width, height, 3);
    for(int i = 0; i < width; i++) {
        for(int j = 0; j < height; j++) {
            Vector n = Normal(i, j, width, height);
            img.pixels[(i + j * width) * 3 + 0] = (n.x + 1) * 0.5 * 255;
            img.pixels[(i + j * width) * 3 + 1] = (n.y + 1) * 0.5 * 255;
            img.pixels[(i + j * width) * 3 + 2] = (n.z + 1) * 0.5 * 255;
        }
    }

    return img;
}

double HeightField::AverageSlope(int x, int z) const {
    double sum = 0;
    int nb = 0;
    double h1 = (*this)(x, z);
    if(x > 0) {
        sum += std::abs(h1 - (*this)(x - 1, z));
        nb++;
    }
    if(x < width - 1) {
        sum += std::abs(h1 - (*this)(x + 1, z));
        nb++;
    }
    if(z > 0) {
        sum += std::abs(h1 - (*this)(x, z - 1));
        nb++;
    }
    if(z < height - 1) {
        sum += std::abs(h1 - (*this)(x, z + 1));
        nb++;
    }
    return sum / nb;
}

ImageData HeightField::Slope(int sampleWidth, int sampleHeight) const {
    std::vector<float> slopes(width * height);
    float maxSlope = std::numeric_limits<float>::min();
    float minSlope = std::numeric_limits<float>::max();
    for (int i = 0; i < width; i++)
    {
        for (int j = 0; j < height; j++) {
            slopes[i + j * width] = AverageSlope(i, j);
            minSlope = std::min(minSlope, slopes[i + j * width]);
            maxSlope = std::max(maxSlope, slopes[i + j * width]);
        }
    }
    // normalize
    for (int i = 0; i < width * height; i++)
    {
        slopes[i] = (slopes[i] - minSlope) / (maxSlope - minSlope);
    }

    ImageData img(width, height, 3);
    for (int i = 0; i < width * height; i++)
    {
        img.pixels[i * 3 + 0] = 255 - slopes[i] * 255;
        img.pixels[i * 3 + 1] = 255 - slopes[i] * 255;
        img.pixels[i * 3 + 2] = 255;// - slopes[i] * 255;
    }

    std::cout<<"DONE SLOPE"<<std::endl;
    return img;
}

float HeightField::laplacian(int i, int j) const {
    float h = gradient(i, j).y;
    float dx, dz;
    dx = dz = 0;
    int sumx, sumz;
    sumx = sumz = 0;

    if(i > 0) {
        dx += gradient(i - 1, j).y - h;
        sumx++;
    }

    if(i < width - 1) {
        dx += gradient(i + 1, j).y - h;
        sumx++;
    }

    if(j > 0) {
        dz += gradient(i, j - 1).y - h;
        sumz++;
    }

    if(j < height - 1) {
        dz += gradient(i, j + 1).y - h;
        sumz++;
    }

    return (dx / sumx + dz / sumz) / 2;
}


ImageData HeightField::Laplacian(int sampleWidth, int sampleHeight) const {
    std::vector<float> res(sampleWidth * sampleHeight);
    float maxLaplacian = std::numeric_limits<float>::min();
    float minLaplacian = std::numeric_limits<float>::max();
    for (int i = 0; i < sampleWidth; i++)
    {
        for (int j = 0; j < sampleHeight; j++) {
            res[i + j * sampleWidth] = laplacian(i, j);
            minLaplacian = std::min(minLaplacian, res[i + j * sampleWidth]);
            maxLaplacian = std::max(maxLaplacian, res[i + j * sampleWidth]);
        }
    }

    ImageData img(sampleWidth, sampleHeight, 3);
    for (int i = 0; i < sampleWidth * sampleHeight; i++)
    {
        // calmp res to [0, 1]
        img.pixels[i * 3 + 0] = 255 - (res[i] - minLaplacian) / (maxLaplacian - minLaplacian) * 255;
        img.pixels[i * 3 + 1] = 255 - (res[i] - minLaplacian) / (maxLaplacian - minLaplacian) * 255;
        img.pixels[i * 3 + 2] = 255 - (res[i] - minLaplacian) / (maxLaplacian - minLaplacian) * 255;
    }

    std::cout<<"DONE res"<<std::endl;
    return img;
}

ImageData HeightField::GradientMagnitude(int sampleWidth, int sampleHeight) const {
    ImageData img(width, height, 3);
    Vector min = Vector(std::numeric_limits<float>::max(), std::numeric_limits<float>::max(), std::numeric_limits<float>::max());
    Vector max = Vector(std::numeric_limits<float>::min(), std::numeric_limits<float>::min(), std::numeric_limits<float>::min());
    for(int i = 0; i < width; i++) {
        for(int j = 0; j < height; j++) {
            Vector n = gradient(i, j);
            if(n.x < min.x) min.x = n.x;
            else if (n.x > max.x) max.x = n.x;
            if(n.z < min.z) min.z = n.z;
            else if (n.z > max.z) max.z = n.z;
        }
    }
    for(int i = 0; i < width; i++) {
        for(int j = 0; j < height; j++) {
            Vector n = gradient(i, j);
            float xc = (n.x - min.x) / (max.x - min.x);
            float zc = (n.z - min.z) / (max.z - min.z);

            xc = std::max(0.f, std::min(xc, 1.f));
            zc = std::max(0.f, std::min(zc, 1.f));
            float Magnitude = sqrt(xc * xc + zc * zc);
            img.pixels[(i + j * width) * 3 + 0] = xc * 255;
            img.pixels[(i + j * width) * 3 + 1] = 255 - Magnitude * 255;
            img.pixels[(i + j * width) * 3 + 2] = zc * 255;
        }
    }

    return img;
}

float HeightField::skyViewFixed(const Point & p) {
    float nbSky = 0;
    for(size_t i = 0; i < pointsOnHemisphere.size(); i++) {
        Ray ray(p, normalize(Vector(pointsOnHemisphere[i])));

        float t;
        if(!spheretrace(ray, t)) {
            nbSky = nbSky + 1.f;
        }
    }
    float r = nbSky / pointsOnHemisphere.size();
    return r;
}

float HeightField::accesibility(const Point & p) {
    float nbSky = 0;
    //lines.clear();
    Vector N = Normal(p.x, p.z, width, height);
    Vector T, B;
    T = normalize(Vector(-N.z, 0, N.x));
    B = normalize(cross(N, T));

    bool printLine = false;
    for(size_t i = 0; i < pointsOnHemisphere.size(); i++) {
        Point coordOnHemi = pointsOnHemisphere[i];

        Vector dirOnHemi = Vector(coordOnHemi.x, coordOnHemi.y, coordOnHemi.z);
        //Vector dirOnHemi = coordOnHemi.x * T + coordOnHemi.y * N + coordOnHemi.z * B;
        Ray ray(p, normalize(dirOnHemi));
        float t;
        if(printLine) 
        {
            Mesh line = Mesh(GL_LINES);
            if(!spheretrace(ray, t)) {
                nbSky = nbSky + 1.f;
                line.color(0, 1, 0);
            } 
            else {
                line.color(1, 0, 0);
            }
            line.vertex(p);
            line.vertex(ray(t));
            lines.push_back(line);
        } else  {
            if(!spheretrace(ray, t)) {
                nbSky = nbSky + 1.f;
            }
        }
        
    }
    return nbSky / pointsOnHemisphere.size();
}

bool HeightField::spheretrace(const Ray & ray, float & t) const {
    // maximal slope
    double lambda = pMax.y;
    t = epsilon;

    double maxDist = ((pMax.x - pMin.x) / width) * 64.0;
    int i = 0;
    while(t < maxDist && i < 200) {

        Point p = ray(t);
        if(!inside(p)) {
            return false;
        }
        double h = Height(p.x, p.z);
        if(p.y <= h) {
            return true;
        }
        t += std::max(abs(h - p.y) / lambda, epsilon);
        i++;
    }
    return false;
}

ImageData HeightField::Accessible(int sampleWidth, int sampleHeight) {
    ImageData img(sampleWidth, sampleHeight, 3);
    auto now = std::chrono::system_clock::now();
    lines.clear();
    #pragma omp parallel for
    for(int i = 0; i < sampleWidth; i++) {
        //printf("%d / %d\n", (i + 1), sampleWidth);
        for(int j = 0; j < sampleHeight; j++) {
            Point p = Vertex(i, j, sampleWidth, sampleHeight);
            float sky = accesibility(p);
            img.pixels[(i + j * sampleWidth) * 3 + 0] = sky * 255;
            img.pixels[(i + j * sampleWidth) * 3 + 1] = sky * 255;
            img.pixels[(i + j * sampleWidth) * 3 + 2] = sky * 255;
        }
    }
    auto end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - now;
    std::cout<<"Elapsed time : "<<elapsed_seconds.count()<<std::endl;
    return img;
}

Point HeightField::Vertex(int i, int j, int w, int h) const{
    
    float x = ((float)i / (float)w) * pMax.x;
    float z = ((float)j / (float)h) * pMax.z;

    //std::cout<< "Vertex :: i = " << i << " j = " << j << " w = " << w << " h = " << h << " x = " << x << " z = " << z << std::endl;

    float height = Height(x, z);
    

    //std::cout<<height<<std::endl;

    return Point(x, height, z);
}

Vector HeightField::Normal(int i, int j, int w, int h) const{
    Point t, b, l, r;

    t = Vertex(i, j - 1, w, h);
    b = Vertex(i, j + 1, w, h);
    l = Vertex(i - 1, j, w, h);
    r = Vertex(i + 1, j, w, h);

    Vector n = normalize(cross(t - b, l - r));

    return n;
}