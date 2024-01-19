#include "ScalarField.h"

void ScalarField::blur()
{
    Grid kernel(3, 3);
    kernel = 1.0 / 9.0;
    Convolute(kernel);
}

void ScalarField::smooth()
{
    Grid kernel(3, 3);
    kernel(0, 0) = 1.0;
    kernel(1, 0) = 2.0;
    kernel(2, 0) = 1.0;

    kernel(0, 1) = 2.0;
    kernel(1, 1) = 4.0;
    kernel(2, 1) = 2.0;

    kernel(0, 2) = 1.0;
    kernel(1, 2) = 2.0;
    kernel(2, 2) = 1.0;

    kernel = 1.0 / 16.0;

    Convolute(kernel);
}