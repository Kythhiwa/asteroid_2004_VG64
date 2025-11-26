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



Vector3D Observatories::getConvertObsForJD(double JD, std::string code)
{

    
    double itrs[3] = {obs[code].x, obs[code].y, obs[code].z};
    double elong = atan2(itrs[1], itrs[0]);
    double u = sqrt(itrs[0]*itrs[0] + itrs[1] * itrs[1]);
    double v = itrs[2];

    // TDB -> TT 
    double tdb1 = floor(JD);
    double tdb2 = JD - tdb1;
    double dtr = iauDtdb(tdb1, tdb2, 0.0, 0.0, 0.0, 0.0);
    double tt1, tt2;
    iauTdbtt(tdb1, tdb2, dtr, &tt1, &tt2);

    // TT -> UTC
    double tai1, tai2;
    iauTttai(tt1, tt2, &tai1, &tai2); // TT -> TAI
    double utc1, utc2;
    int status = iauTaiutc(tai1, tai2, &utc1, &utc2);
    if (status != 0) {
        throw std::runtime_error("[ERROR]TAI -> UTC");
    }

    double jd_utc_full = utc1 + utc2;
    double mjd_utc_full = jd_utc_full - 2400000.5;

    auto it = eop.lower_bound(EOPEntry{mjd_utc_full, 0, 0, 0});
    if (it == eop.end()) {
        --it; 
    }

    // UTC -> UT1
    double ut1_jd = jd_utc_full + (it->ut1_utc / 86400.0);
    double ut_frac = fmod(ut1_jd, 1.0);

    // TDB -> TT 
    double dtr_exact = iauDtdb(tdb1, tdb2, ut_frac, elong, u, v);
    double tt1_exact, tt2_exact; 
    iauTdbtt(tdb1, tdb2, dtr_exact, &tt1_exact, &tt2_exact);

    double ut11 = floor(ut1_jd);
    double ut12 = ut1_jd - ut11;
    
    //Matrix GCRS ->ITRS
    double rc2t[3][3];
    iauC2t06a(tt1_exact, tt2_exact, ut11, ut12, it->xp, it->yp, rc2t);

    // ITRS -> GCRS
    double gcrs[3];
    iauTrxp(rc2t, itrs, gcrs);
    
    StateVector planet = eph.getStateVector((int)Ephemeris::CelestialBody::Earth, (int)Ephemeris::CelestialBody::SSB, JD);

    return Vector3D{gcrs[0], gcrs[1], gcrs[2]} + planet.x;
}
