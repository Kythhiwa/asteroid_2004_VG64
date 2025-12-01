#include <iostream>
#include "include/data_reader.hpp"
#include "ephemeris.hpp"
#include "rk4Integrator.hpp"
#include "observatories.hpp"

#include <fstream>
#include <iomanip>


void solve(Ephemeris &eph, Rk4Integrator &r, Observatories &obs, std::vector<RaDec> &a)
{
    double start_jd = 2453320;
    BodyVector asteroid = {
        {1.367372443402559E+08, 6.602611574848668E+07, 3.033666479027805E+07},

        {8.299501903651528E+00, 2.704928307360539E+01,-1.899669293113739E+00}
    };
    
    auto Sun1 = eph.getStateVector((int)Ephemeris::CelestialBody::Sun, (int) Ephemeris::CelestialBody::SSB, start_jd);
    r.integrateOrbit(asteroid, start_jd, 12 * 365, 3600);
    {
        std::size_t i = 1;
        std::ofstream file("result.txt");

        file << std::fixed << std::setprecision(10);
        for (auto [cur_jd_utc, Ra, Dec, code] : a)
        {
            file << "observations " << i << " JD " << cur_jd_utc  <<  ": \n";
            double cur_jd_tdb = obs.convertUtcToTdb(cur_jd_utc, code);
            auto comp_pos_obs  = obs.getConvertObsForJD(cur_jd_tdb, code);
        
            auto Earth = eph.getStateVector((int)Ephemeris::CelestialBody::Earth, (int) Ephemeris::CelestialBody::SSB, cur_jd_tdb);

            auto Sun = eph.getStateVector((int)Ephemeris::CelestialBody::Sun, (int) Ephemeris::CelestialBody::SSB, cur_jd_tdb);

            auto ltc = r.applyAstrometricCorrections(cur_jd_tdb, comp_pos_obs, Earth, Sun);
            double Ra1, Dec1, dist1;
            r.cartToRaDec(ltc, Ra1, Dec1, dist1);

            file << "Exact   Ra = " << Ra << " Dec = " << Dec << "\n";
            file << "Compute Ra = " << Ra1 << " Dec = " << Dec1 << "\n";
            file << "Diff    Ra = " << Ra - Ra1 << " Dec = "  << Dec - Dec1 << "\n\n"; 

            ++i;
        }
        file.close();
    }
}

int main()
{
    Ephemeris eph;
    eph.loadFile("data/de440s.bsp");

    Observatories obs(eph);
    obs.loadFileObs("data/ObsCodes.html");
    obs.loadFileEop("data/finals.data.iau2000.txt");
    
    Rk4Integrator r(eph);

    auto a = read_observations("data/observations.txt");

    solve(eph, r, obs, a);
}
