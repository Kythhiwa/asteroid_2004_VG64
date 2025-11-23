#include <stdexcept>

#include "ephemeris.hpp"



Ephemeris::Ephemeris()
{
    eph = ephCreate();
    if (!eph)
    {
        throw std::runtime_error("initialization error\n");
    }
    // 1 -  AU, 2 - KM
    ephSetDistanceUnits(eph, 2);
    // 3 - SEC, 4 - DAY
    ephSetTimeUnits(eph, 3);
}

void Ephemeris::loadFile(std::string filename)
{
    int result = ephLoadFile(eph, filename.c_str());
    if (result != 0)
    {
        ephDestroy(eph);
        eph = nullptr;
        throw std::runtime_error("load file error\n");
    }

}

StateVector Ephemeris::getStateVector(int body, int reference, double jd)
{
    double date0 = std::floor(jd);
    double date1 = jd - date0;

    std::vector<double> x = {0,0,0};
    std::vector<double> vx = {0,0,0};

    int result = ephCalculateRectangular(eph, 
                                         body, 
                                         reference,
                                         date0,
                                         date1,
                                         x.data(),
                                         vx.data());
    if (result != 0)
    {
        throw std::runtime_error("Failed to calculate ephemeris\n");
    }

    return {
        {x[0], x[1], x[2]},
        {vx[0], vx[1], vx[2]}
    };
}

 
void Ephemeris::setDistance(int dist_type)
{
    ephSetDistanceUnits(eph, dist_type);
}

void Ephemeris::setTime(int time_type)
{
    ephSetTimeUnits(eph, time_type);
}


Ephemeris::~Ephemeris()
{
    if (eph)
    {
        ephDestroy(eph);
        eph = nullptr;
    }
}
