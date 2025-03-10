#include "CPYDemo.h"
#include <iostream>
#include <vector>

namespace py = pybind11;

void CPYDemo::printWorld()
{
    std::cout << "Hello World" << std::endl; 
}

long long CPYDemo::count(long long n) {
    long long sum = 0;
    for (long long i = 0; i < n; i++) {
        sum += i;
    }
    return sum;
}

void CPYDemo::printHello()
{
    std::cout << "goodbye" << std::endl;
}

std::vector<std::vector<double>> CPYDemo::generate_noise_map(int width, int height, double scale, double offsetX, double offsetY)
{
    // Create a 2D vector to hold the results
    std::vector<std::vector<double>> result(height, std::vector<double>(width));
    
    // Create noise generator
    OpenSimplexNoise::Noise noise(321314);
    
    // Fill the vector with noise values
    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            double sX = (x + offsetX) * scale;
            double sY = (y + offsetY) * scale;
            
            // Calculate noise and store in vector
            result[y][x] = noise.eval(sX, sY);
        }
    }
    return result;
}

double CPYDemo::testFunction(double n)
{
    return n;
}