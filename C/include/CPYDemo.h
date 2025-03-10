#ifndef CPYDEMO_H
#define CPYDEMO_H

#include "openSimplexNoise/OpenSimplexNoise/OpenSimplexNoise.h"
#include <string>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

namespace py = pybind11;

/**
 * @brief CPYDemo class for demonstrating basic C++ functionality
 */
class CPYDemo 
{
    public:
        CPYDemo() = default;
        ~CPYDemo() = default;
        void printWorld();
        long long count(long long n);
        void printHello();

        std::vector<std::vector<double>>  generate_noise_map(int width, int height, double scale, double offsetX = 0.0, double offsetY = 0.0);
        double testFunction(double n);
    private:
    protected:
};

#endif // CPYDEMO_H