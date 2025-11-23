#include "rk4Integrator.hpp"

BodyVector Rk4Integrator::rk4step(const BodyVector &state, double jd, double dt)
{
    // k1 = f(t, y)
    BodyVector k1 = computeDer(state, jd);
    
    // k2 = f(t + dt/2, y + dt/2 * k1)
    BodyVector state_k2;
    state_k2.x = state.x + k1.v * (dt * 0.5);
    state_k2.v = state.v + k1.a * (dt * 0.5);
    BodyVector k2 = computeDer(state_k2, jd + dt * 0.5);
    
    // k3 = f(t + dt/2, y + dt/2 * k2)  
    BodyVector state_k3;
    state_k3.x = state.x + k2.v * (dt * 0.5);
    state_k3.v = state.v + k2.a * (dt * 0.5);
    BodyVector k3 = computeDer(state_k3, jd + dt * 0.5);
    
    // k4 = f(t + dt, y + dt * k3)
    BodyVector state_k4;
    state_k4.x = state.x + k3.v * dt;
    state_k4.v = state.v + k3.a * dt;
    BodyVector k4 = computeDer(state_k4, jd + dt);
    
    // y_new = y + dt/6 * (k1 + 2*k2 + 2*k3 + k4)
    BodyVector new_state;
    new_state.x = state.x + (k1.v + k2.v * 2.0 + k3.v * 2.0 + k4.v) * (dt / 6.0);
    new_state.v = state.v + (k1.a + k2.a * 2.0 + k3.a * 2.0 + k4.a) * (dt / 6.0);
    
    return new_state;
}

BodyVector Rk4Integrator::computeDer(const BodyVector &state, double jd)
{
    BodyVector deriv;
    deriv.v = state.v;  // dx/dt = v
    deriv.a = {0, 0, 0};  // dv/dt = a 
    
    for (auto &[b, gm] : planetsGM)
    {
        StateVector planet = eph.getStateVector(b, (int)Ephemeris::CelestialBody::SSB, jd);
        
        Vector3D r_vec = planet.x - state.x;
        double r = r_vec.len();
        
        if (r < 1e-10) continue;
        
        // a = GM * (r_vec) / r^3
        double factor = gm / (r * r * r);
        deriv.a.x += factor * r_vec.x;
        deriv.a.y += factor * r_vec.y;
        deriv.a.z += factor * r_vec.z;
    }
    
    return deriv;
}
