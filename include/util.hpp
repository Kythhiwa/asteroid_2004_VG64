#pragma once

#include <cmath>
#include <string>

namespace util 
{
    void sphToCart(double Long, double cos, double sin, double &x, double &y, double &z);
    
    void convertToRaDec(const std::string &ra, const std::string &dec, double &Ra, double &Dec);

    void dateToJd(int year, int month, double day, double &jd);
}

