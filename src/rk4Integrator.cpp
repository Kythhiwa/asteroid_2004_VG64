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



void Rk4Integrator::integrateOrbit(
            const BodyVector &state,
            double start_jd, 
            double duration_days,
            double dt)
{
    start_jd_ = start_jd;

    BodyVector cur = state;
    double step_days = dt / 86400;
    double step = duration_days / step_days;
    
    step_jd_ = step_days;

    states.clear();
    states.reserve(step);

    for (std::size_t i = 0; i < step; ++i)
    {
        double cur_jd = start_jd + i * step_days;
        cur = rk4step(cur, cur_jd, dt);
        states.emplace_back(cur_jd, cur);
    }
}

std::vector<Rk4Integrator::State> Rk4Integrator::getOrbit() const 
{
    return states;
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

Vector3D Rk4Integrator::interpolatePosition(double jd) const
{
    if (jd < start_jd_)
    {
        throw std::runtime_error("jd < start_jd");
    }
    std::size_t idx = (jd - start_jd_) / step_jd_;

    if (idx == 0)
    {
        idx = 1;
    }
    if (idx >= states.size())
    {
        return states.back().x.x;
    }
    
    auto cur =  states[idx];
    auto prev = states[idx - 1];
    double jd1 = prev.jd;
    double jd2 = cur.jd;
    Vector3D pos1 = prev.x.x;
    Vector3D pos2 = cur.x.x;

    double t = (jd - jd1) / (jd2 - jd1);
    return pos1 + (pos2 - pos1) * t;
}

Vector3D Rk4Integrator::lightTimeCorrection(double jd, const Vector3D &obs) const
{
    double delta = 0.0;
    const double c = 299792.458; 
    const std::size_t max_iters = 20;
    const double eps = 1e-18;

    for (std::size_t i = 0; i < max_iters; ++i) 
    {
        double cur_jd = jd - delta / 86400;
        Vector3D obj_pos = interpolatePosition(cur_jd);
        double dist = (obj_pos - obs).len();

        double cur_delta = dist / c;
        
        if (std::abs(cur_delta - delta) < eps)
        {
            break;
        }
        
        delta = cur_delta;
    }

    return interpolatePosition(jd - delta / 86400);
}


void Rk4Integrator::cartToRaDec(const Vector3D& pos, double& ra, double& dec, double& dist)
{

    double c[3];
    c[0] = pos.x;
    c[1] = pos.y; 
    c[2] = pos.z;
    
    double theta, phi;
    iauC2s(c, &phi, &theta); 
    
    ra = iauAnp(phi) * 180.0 / M_PI;  
    dec = theta * 180.0 / M_PI;      
    dist = std::sqrt(c[0] * c[0] + c[1] * c[1] + c[2] * c[2]);
}

Vector3D Rk4Integrator::lightDeflection(const Vector3D &obs_ast,
                                        const Vector3D &obs,
                                        const StateVector &Earth,
                                        const StateVector &Sun)
{
    const double au_to_km = 149597870.7;
    
    Vector3D obs_pos_au = obs * (1.0 / au_to_km);
    Vector3D sun_pos_au = Sun.x * (1.0 / au_to_km);
    Vector3D ast_pos_au = (obs + obs_ast) * (1.0 / au_to_km); 
    Vector3D p_vec = obs_ast * (1.0 / au_to_km); 

    double p_len = p_vec.len();
    Vector3D p_dir = p_vec.normalized();
    double p[3] = {p_dir.x, p_dir.y, p_dir.z};

    Vector3D body_to_obs = obs_pos_au - sun_pos_au; 
    double dist_body_obs = body_to_obs.len();           
    Vector3D e_dir = body_to_obs.normalized();      
    double e[3] = {e_dir.x, e_dir.y, e_dir.z};

    double em = dist_body_obs; 

    Vector3D body_to_source = ast_pos_au - sun_pos_au; 
    Vector3D q_dir = body_to_source.normalized();
    
        
    double q[3] = {q_dir.x, q_dir.y, q_dir.z};

    double bm = 1.0;
    double dlim = 1e-6;
    double p1[3];       
    
    iauLd(bm, p, q, e, em, dlim, p1);
    
    Vector3D result_dir(p1[0], p1[1], p1[2]);
    return result_dir * p_len * au_to_km;
}


Vector3D Rk4Integrator::lightAbberation(
        const Vector3D &obs_ast, 
        const Vector3D &obs, 
        const StateVector &Earth, 
        const StateVector &Sun)
{
    const double au_to_km = 149597870.7;
    const double c_km_s = 299792.458; 
    
    Vector3D natural_dir = obs_ast.normalized(); 
    double pnat[3] = {natural_dir.x, natural_dir.y, natural_dir.z};
    
    Vector3D velocity_c = Earth.v * (1.0 / c_km_s);
    double v[3] = {velocity_c.x, velocity_c.y, velocity_c.z};
    
    Vector3D obs_pos_au = obs * (1.0 / au_to_km);
    Vector3D sun_pos_au = Sun.x * (1.0 / au_to_km);
    Vector3D sun_to_obs = obs_pos_au - sun_pos_au;
    double s = sun_to_obs.len();     

    double v_squared = velocity_c.x * velocity_c.x + 
                      velocity_c.y * velocity_c.y + 
                      velocity_c.z * velocity_c.z;
    double bm1 = sqrt(1.0 - v_squared);
    
    double ppr[3];
    iauAb(pnat, v, s, bm1, ppr);
    
    Vector3D aberrated_dir(ppr[0], ppr[1], ppr[2]);
    double distance = obs_ast.len();
    
    return aberrated_dir * distance;
}

Vector3D Rk4Integrator::applyAstrometricCorrections(
            double jd, 
            const Vector3D &obs,
            const StateVector &earth,
            const StateVector &Sun)
{
    Vector3D ast_ltc = lightTimeCorrection(jd, obs);

    Vector3D obs_ast = ast_ltc - obs;
    obs_ast = lightAbberation(obs_ast, obs, earth, Sun);
    
    obs_ast = lightDeflection(obs_ast, obs, earth, Sun);
    
    return obs_ast;
}
