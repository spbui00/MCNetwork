#pragma once

#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>
#include <iterator>
#include <ostream>
#include <random>
#include <string>
#include <sys/stat.h>
#include <type_traits>

namespace enhance {
extern unsigned int seed;
extern std::mt19937_64 rand_engine;

double random_triangle(double min, double peak, double max);

float fastExp(float x);
double mediumFastExp(double x);
double sigmoid(std::vector<double> &coeff, double x);
double polynom(std::vector<double> &coeff, double x);

unsigned long long fastExp2(int x);

std::string multipliplyString(std::string in, int n);

/*! random double from [a,b) */
inline double random_double(double a, double b) {
    std::uniform_real_distribution<double> distribution(a, b);
    return distribution(rand_engine);
}

/*! random int from [a,b] */
inline int random_int(int a, int b) {
    std::uniform_int_distribution<int> intdistribution(a, b);
    return intdistribution(rand_engine);
}

/*! random int from [a,b] */
inline unsigned int random_int(unsigned int a, unsigned int b) {
    std::uniform_int_distribution<unsigned int> intdistribution(a, b);
    return intdistribution(rand_engine);
}
} // namespace enhance
