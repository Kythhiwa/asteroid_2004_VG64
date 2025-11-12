#pragma once

#include <string>
#include <vector>
#include <map>

struct dataObs {
    double Long;
    double cos;
    double sin;
    dataObs(double l, double c, double s)
        :   Long(l), cos(c), sin(s) {}
        
};


std::map<std::string, dataObs> read_obs(std::string filename);

