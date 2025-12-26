#pragma once

#include <array>

#include "rk4Integrator.hpp"
#include "data_reader.hpp"
#include "observatories.hpp"

class GaussNewton
{
public:
    GaussNewton(std::vector<RaDec> &a,
                Ephemeris &eph,
                Observatories &obs,
                Rk4Integrator &r) 
    : a(a), eph(eph), obs(obs), r(r) {}

    int start(std::size_t max_iters = 50, double eps = 1e-8);

    void setX0(BodyVector x);
    void setJD(double jd);

    BodyVector getXbest() const;
    

private:

    static constexpr std::array<double, 6> SCALE = { 
                0x1.0p27,  
                0x1.0p27,  
                0x1.0p27,  
                0x1.0p6,  
                0x1.0p6,  
                0x1.0p6
    };

    std::vector<RaDec> &a;   
    Ephemeris &eph;
    Observatories &obs;
    Rk4Integrator &r;

    BodyVector x0;
    BodyVector x_best;

    double jd_start;
    
    BodyVector scaledToPhysical(const BodyVector &x);
    BodyVector physicalToScaled(const BodyVector &x);
    
    
    // x - physical
    std::vector<std::vector<double>> computeJ(const BodyVector &x);

    // x - scaled
    std::vector<double> computeResiduals(const BodyVector &x);

    bool solveSystem(const std::vector<std::vector<double>> &A, 
                     const std::vector<double> &b, 
                     std::vector<double> &x);
};

