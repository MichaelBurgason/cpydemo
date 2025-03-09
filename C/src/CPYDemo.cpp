#include "CPYDemo.h"
#include <iostream>

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