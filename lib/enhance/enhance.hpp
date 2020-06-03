#pragma once

#include <ostream>
#include <iostream>
#include <cstdint>
#include <algorithm>
#include <random>
#include <iterator>
#include <type_traits>
#include <string>
#include <array>
#include <sys/stat.h>

namespace enhance
{

    extern unsigned int     seed;
    extern std::mt19937_64  rand_engine;

    inline double random_double(double, double);
    inline unsigned int    random_uns_int(int, int);
    inline int    random_int(int, int);

    double random_triangle(double min, double peak, double max);

    float fastExp(float x);
    double mediumFastExp(double x);
    double sigmoid(std::vector<double>& coeff, double x);
    double polynom(std::vector<double>& coeff, double x);


    unsigned long long fastExp2(int x);



    std::string multipliplyString(std::string in, int n);

    



    /*! random double from [a,b) */
    double random_double(double a, double b)
    {
        std::uniform_real_distribution<double> distribution(a,b);
        return distribution(rand_engine);
    }

    /*! random int from [a,b] */
    unsigned int random_uns_int(int a, int b)
    {
        std::uniform_int_distribution<int> intdistribution(a,b);
        return intdistribution(rand_engine);
    }

    
    /*! random int from [a,b] */
    int random_int(int a, int b)
    {
        std::uniform_int_distribution<int> intdistribution(a,b);
        return intdistribution(rand_engine);
    }



    // float fastExp(float x)
    // // !!! WARNING yielding completly wrong results if not x << n (here n=256) !!!
    // {
    //     // std::cout<<"x= "<<x;
    //     x = 1.0 + x / 256.0;
    //     x *= x; x *= x; x *= x; x *= x;
    //     x *= x; x *= x; x *= x; x *= x;
    //     // std::cout<<" exp(x) = "<<x<<std::endl;
    //     return x;        
    // }

    // float mediumFastExp(float x)
    // // !!! WARNING yielding completly wrong results if not x << n (here n=4096) !!!
    // {
    //     // std::cout<<"x= "<<x;
    //     x = 1.0 + x / 4096.0;
    //     x *= x; x *= x; x *= x; x *= x;
    //     x *= x; x *= x; x *= x; x *= x;
    //     x *= x; x *= x; x *= x; x *= x;
    //     // std::cout<<" exp(x) = "<<x<<std::endl;
    //     return x;        
    // }



}
namespace enh = enhance;
