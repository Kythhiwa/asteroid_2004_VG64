#pragma once

#include "ephemeris.hpp"


class Rk4Integrator
{
public:
    Rk4Integrator(Ephemeris &eph) : eph(eph) {};
    
    BodyVector rk4step(const BodyVector &state, double jd, double dt);

private: 
    Ephemeris &eph;
    
    BodyVector computeDer(const BodyVector &state, double jd);

    std::vector<std::pair<int, double>> planetsGM = {
        {(int)Ephemeris::CelestialBody::Sun, 132712440043.17},
        {(int)Ephemeris::CelestialBody::Mercury, 22031.78000},
        {(int)Ephemeris::CelestialBody::Venus,324858.59200},
        {(int)Ephemeris::CelestialBody::Earth,398600.43629 },
        {(int)Ephemeris::CelestialBody::Mars, 42828.37521},
        {(int)Ephemeris::CelestialBody::Jupiter, 126712764.13345},
        {(int)Ephemeris::CelestialBody::Saturn,37940585.20000},
        {(int)Ephemeris::CelestialBody::Uranus, 5794556.46575},
        {(int)Ephemeris::CelestialBody::Neptune, 6836527.10058},
        {(int)Ephemeris::CelestialBody::Pluto, 975.50118}
    };
};
