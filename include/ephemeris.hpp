#pragma once

#include <string>
#include <vector> 

#include "ephaccess.h"
#include "geometry.hpp"

class Ephemeris
{
public:
    enum class CelestialBody {
        SSB = 0,              
        Mercury = 1,
        Venus = 2,
        EMB = 3,              
        Mars = 4,             
        Jupiter = 5,          
        Saturn = 6,           
        Uranus = 7,           
        Neptune = 8,          
        Pluto = 9,            
        Sun = 10,
        Moon = 301,
        Earth = 399
    };
    

    Ephemeris();
    
    void loadFile(std::string filename);

    StateVector getStateVector(int body, int reference, double jd);
    
    void setDistance(int dist_type);
    void setTime(int time_type);

    ~Ephemeris();
private:
    EphAccess *eph;

};
