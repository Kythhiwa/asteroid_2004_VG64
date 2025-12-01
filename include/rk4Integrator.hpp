#pragma once

#include <algorithm>

#include <iomanip>
#include "ephemeris.hpp"
#include "../external/sofa/20231011/c/src/sofa.h"

class Rk4Integrator
{
public:
    struct State 
    {
        double jd;
        BodyVector x;
        State(double jd, BodyVector &x)
            : jd(jd), x(x) {}
    };

    Rk4Integrator(Ephemeris &eph) : eph(eph) {};
    
    BodyVector rk4step(const BodyVector &state, 
                       double jd, 
                       double dt);

    void cartToRaDec(const Vector3D& pos, 
                     double& ra, 
                     double& dec, 
                     double& dist);
    
    void computeObservedRaDec(double jd_tdb, 
                              const Vector3D& obs_pos, 
                              double& ra, 
                              double& dec, 
                              double& dist);

    Vector3D interpolatePosition(double jd) const;

    // dt (sec)
    void integrateOrbit(
            const BodyVector &state,
            double start_jd, 
            double duration_days,
            double dt = 3600);
   
  
    std::vector<State> getOrbit() const;

    
    Vector3D applyAstrometricCorrections(
            double jd, 
            const Vector3D &obs,
            const StateVector &earth,
            const StateVector &Sun);

private: 
     // JD (TDB)
    Vector3D lightTimeCorrection(double jd, const Vector3D &obs) const;
    
    Vector3D lightAbberation(const Vector3D &obs_ast, 
                             const Vector3D &obs,
                             const StateVector &Earth, 
                             const StateVector &Sun);

    Vector3D lightDeflection(const Vector3D &obs_ast,
                             const Vector3D &obs,
                             const StateVector &Earth,
                             const StateVector &Sun);

    BodyVector computeDer(const BodyVector &state, double jd);
    
    Ephemeris &eph;
    std::vector<State> states;

    double start_jd_;
    double step_jd_;

    std::vector<std::pair<int, double>> planetsGM = {
        {(int)Ephemeris::CelestialBody::Sun, 132712440043.17},
        {(int)Ephemeris::CelestialBody::Mercury, 22031.78000},
        {(int)Ephemeris::CelestialBody::Venus,324858.59200},
        {(int)Ephemeris::CelestialBody::Earth,398600.43629},
        {(int)Ephemeris::CelestialBody::Mars, 42828.37521},
        {(int)Ephemeris::CelestialBody::Jupiter, 126712764.13345},
        {(int)Ephemeris::CelestialBody::Saturn,37940585.20000},
        {(int)Ephemeris::CelestialBody::Uranus, 5794556.46575},
        {(int)Ephemeris::CelestialBody::Neptune, 6836527.10058},
        {(int)Ephemeris::CelestialBody::Pluto, 975.50118}
    };
};
