
#include "util.hpp"



namespace util 
{
    void sphToCart(double Long, double cos, double sin, double &x, double &y, double &z)
    {
        x = cos * std::cos(Long);
        y = cos * std::sin(Long);
        z = sin;
    }
}
