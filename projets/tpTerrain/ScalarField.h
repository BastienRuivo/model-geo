#pragma once

#include "Grid.h"
#include "../../src/gKit/texture.h"
#include <string>

class ScalarField : public Grid
{
private:
public:
    ScalarField() : Grid() {}
    ScalarField(const Point & min, const Point & max, int gridWidth, int gridHeight) : Grid(min, max, gridWidth, gridHeight) {
    }
    ScalarField(const Point & min, float width, float height, int gridWidth, int gridHeight) : Grid(min, width, height, gridWidth, gridHeight) {
    }
    

    Vector gradient(int x, int z) const {
        float height = (*this)(x, z);
        float dx, dz;
        dx = dz = 0;
        int sumx, sumz;
        sumx = sumz = 0;

        if(x > 0) {
            dx += (*this)(x - 1, z) - height;
            sumx++;
        }

        if(x < width - 1) {
            dx += (*this)(x + 1, z) - height;
            sumx++;
        }

        if(z > 0) {
            dz += (*this)(x, z - 1) - height;
            sumz++;
        }

        if(z < height - 1) {
            dz += (*this)(x, z + 1) - height;
            sumz++;
        }
        

        auto res = Vector(dx / sumx, 0, dz / sumz);
        float hei = length(res);
        //normalize(res);
        res.y = hei;

        return res;
        
    }
    
    void blur();
    void smooth();

    void UpdateFromImage(const ImageData & img) {
        width = img.width;
        height = img.height;
        values.resize(width * height);
        for (int i = 0; i < width * height; i++)
        {
            auto value = img.pixels[i * 3 + 0];
            values[i] = (255 - value) / 255.0;
        }
    }

    void loadFromFile(const std::string & path, const Point & pmin = Point(0, 0, 0), const Point & pmax = Point(1, 1, 0), float maxHeight = 1.0) {
        auto img = read_image_data(path.c_str());
        if(!img.data()) {
            std::cerr << "Error loading image " << path << std::endl;
            return;
        }

        std::cout<< "Loading image " << path << " (" << img.width << "x" << img.height << " channels=" << img.channels << " total=" << img.pixels.size() << " / " << img.width * img.height * img.channels << ")" << std::endl;
        width = img.width;
        height = img.height;
        values.resize(width * height);
        for (int i = 0; i < width * height * img.channels; i+=img.channels)
        {
            auto value = img.pixels[i];
            int x = (i / img.channels) % width;
            int z = (i / img.channels) / width;
            values[i / img.channels] = value / 255.0;

        }

        pMin = pmin;
        pMax = pmax;

        std::cout << "Loaded image " << path << " (" << width << "x" << height << ")" << std::endl;
    }

    ImageData toImage() const {
        ImageData img(width, height, 3);
        for (int i = 0; i < width * height; i++)
        {
            img.pixels[i * 3 + 0] = values[i] * 255;
            img.pixels[i * 3 + 1] = values[i] * 255;
            img.pixels[i * 3 + 2] = values[i] * 255;
        }
        return img;
    }



    ~ScalarField() {
    }
};