#pragma once

#include <string>
#include <vector>
#include <map>
#include <set>
#include <fstream>
#include <iostream>
#include <stdexcept>

#include "../external/sofa/20231011/c/src/sofa.h"
#include "../external/sofa/20231011/c/src/sofam.h"

#include "geometry.hpp"
#include "util.hpp"

/*
 * /note MJD UTC 
 */

struct EOPEntry
{
    double mjd;
    double ut1_utc;
    double xp;
    double yp;
    
    EOPEntry(double mjd, 
             double ut1_utc, 
             double xp,
             double yp)
        :   mjd(mjd), ut1_utc(ut1_utc), xp(xp), yp(yp) {}

    bool operator<(const EOPEntry &other) const
    {
         return mjd < other.mjd;
    }
};


std::map<std::string, Vector3D> read_obs(std::string filename);

std::set<EOPEntry> read_eop(std::string filename);
