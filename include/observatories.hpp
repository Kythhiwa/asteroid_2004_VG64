#pragma once

#include <algorithm>

#include "data_reader.hpp"
#include "ephemeris.hpp"
#include "util.hpp"
#include "../external/sofa/20231011/c/src/sofa.h"


class Observatories
{
public:
    Observatories(Ephemeris &eph);
    
    void loadFileObs(std::string filename);
    
    void loadFileEop(std::string filename);
    
    std::map<std::string, Vector3D> getObs() const;
    
    std::vector<std::string> getObsCodes() const;
    
    // JD (TDB)
    Vector3D getConvertObsForJD(double JD, std::string code);
    

    std::set<EOPEntry> getEop() const;

private:
    Ephemeris &eph;
    std::map<std::string, Vector3D>  obs;
    std::set<EOPEntry> eop;
    
};
