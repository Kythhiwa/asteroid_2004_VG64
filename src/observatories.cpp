#include "observatories.hpp" 



Observatories::Observatories(Ephemeris &eph)
    : eph(eph)
{
}


void Observatories::loadFileObs(std::string filename)
{
    obs = read_obs(filename);
}

void Observatories::loadFileEop(std::string filename)
{
    eop = read_eop(filename);
}

std::set<EOPEntry> Observatories::getEop() const
{
    return eop;
}

std::map<std::string, Vector3D> Observatories::getObs() const
{
    return obs;
}


std::vector<std::string> Observatories::getObsCodes() const 
{
    std::vector<std::string> codes;
    for (auto c : obs)
    {
        codes.emplace_back(c.first);
    }
    return codes;
}


Vector3D Observatories::getConvertObsForJD(double JD_tdb, std::string code)
{
    double itrs[3] = {obs[code].x, obs[code].y, obs[code].z};
    double elong = atan2(itrs[1], itrs[0]);
    double u = sqrt(itrs[0]*itrs[0] + itrs[1]*itrs[1]);
    double v = itrs[2];

    // TDB -> TT -> UTC
    double tdb1 = floor(JD_tdb);
    double tdb2 = JD_tdb - tdb1;
    
    double dtr_approx = iauDtdb(tdb1, tdb2, 0.0, elong, u, v);
    double tt1_approx, tt2_approx;
    iauTdbtt(tdb1, tdb2, dtr_approx, &tt1_approx, &tt2_approx);
    
    double tai1, tai2;
    iauTttai(tt1_approx, tt2_approx, &tai1, &tai2);
    
    double utc1, utc2;
    int status = iauTaiutc(tai1, tai2, &utc1, &utc2);
    if (status != 0) {
        throw std::runtime_error("[ERROR] TAI -> UTC conversion failed");
    }
    
    double jd_utc = utc1 + utc2;
    double mjd_utc = jd_utc - 2400000.5;

    auto it = eop.lower_bound(EOPEntry{mjd_utc, 0, 0, 0});
    if (it == eop.end()) 
    {
        it = --eop.end(); 
    }
    else if (it != eop.begin()) 
    {
        auto prev_it = std::prev(it);
        if (fabs(mjd_utc - prev_it->mjd) < fabs(mjd_utc - it->mjd)) 
        {
            it = prev_it;
        }
    }

    // UTC -> UT1
    double ut1_jd = jd_utc + (it->ut1_utc / 86400.0);
    double ut11 = floor(ut1_jd);
    double ut12 = ut1_jd - ut11;

    // TDB -> TT
    double dtr_exact = iauDtdb(tdb1, tdb2, ut12, elong, u, v);
    double tt1_exact, tt2_exact;
    iauTdbtt(tdb1, tdb2, dtr_exact, &tt1_exact, &tt2_exact);

    //ITRS -> GCRS
    double rc2t[3][3];
    iauC2t06a(tt1_exact, tt2_exact, ut11, ut12, it->xp, it->yp, rc2t);

    // ITRS -> GCRS
    double gcrs[3];
    iauTrxp(rc2t, itrs, gcrs);
    
    StateVector earth_state = eph.getStateVector(
        (int)Ephemeris::CelestialBody::Earth, 
        (int)Ephemeris::CelestialBody::SSB, 
        JD_tdb
    );

    return Vector3D{gcrs[0], gcrs[1], gcrs[2]} + earth_state.x;
}


double Observatories::convertUtcToTdb(double jd_utc, std::string code)
{
    double itrs[3] = {obs[code].x, obs[code].y, obs[code].z};
    double elong = atan2(itrs[1], itrs[0]);
    double u = sqrt(itrs[0] * itrs[0] + itrs[1] * itrs[1 ]);
    double v = itrs[2];

    double mjd_utc =  jd_utc - 2400000.5;
    auto it = eop.lower_bound(EOPEntry{mjd_utc, 0, 0, 0});
    if (it == eop.end()) 
    {
        it = --eop.end(); 
    }
    else if (it != eop.begin()) 
    {
        auto prev_it = std::prev(it);
        if (fabs(mjd_utc - prev_it->mjd) < fabs(mjd_utc - it->mjd)) 
        {
            it = prev_it;
        }
    }
    double ut1_jd = jd_utc + (it->ut1_utc / 86400.0);
    double ut11 = floor(ut1_jd);
    double ut12 = ut1_jd - ut11;
    // UTC -> TAI -> TT
    double utc1 = floor(jd_utc);
    double utc2 = jd_utc - utc1;
    
    double tai1, tai2;
    int status = iauUtctai(utc1, utc2, &tai1, &tai2);
    if (status != 0)
    {
        throw std::runtime_error("[ERROR] UTC -> TAI conversion failed");
    }
    
    double tt1_target, tt2_target;
    iauTaitt(tai1, tai2, &tt1_target, &tt2_target);
    double jd_tt_target = tt1_target + tt2_target;

    // TDB iter 
    double jd_tdb = jd_tt_target; 
    
    for (std::size_t iter = 0; iter < 5; ++iter)
    {
        double tdb1 = floor(jd_tdb);
        double tdb2 = jd_tdb - tdb1;
        
        double dtr = iauDtdb(tdb1, tdb2, ut12, elong, u, v);
        
        // TDB -> TT
        double tt1_calc, tt2_calc;
        iauTdbtt(tdb1, tdb2, dtr, &tt1_calc, &tt2_calc);
        double JD_tt_calc = tt1_calc + tt2_calc;
        
        double delta = jd_tt_target - JD_tt_calc;
        jd_tdb += delta;
        
        if (fabs(delta) < 1e-18) 
        {
            break;
        }
    }
    
    return jd_tdb;
}
