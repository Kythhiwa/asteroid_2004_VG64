#include <iostream>
#include "include/data_reader.hpp"
#include "ephemeris.hpp"
#include "rk4Integrator.hpp"
#include "observatories.hpp"

#include <fstream>
#include <iomanip>

#include "gauss_newton.hpp"


void write(std::vector<RaDec> &a,
          Ephemeris &eph,
          Observatories &obs,
          Rk4Integrator &r,
          double jd,
          BodyVector &x,
          std::string filename)
{

    r.integrateOrbit(x, jd, 2 * 365, 3600);
    std::size_t i = 1;
    std::ofstream file(filename);

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

void test(std::vector<RaDec> &a,
          Ephemeris &eph,
          Observatories &obs,
          Rk4Integrator &r)
{
    GaussNewton G(a, eph, obs, r);

    double start_jd = 2453320;       
    G.setJD(start_jd);
    BodyVector asteroid = {
        {1.367372443402559E+08, 6.602611574848668E+07, 3.033666479027805E+07},

        {8.299501903651528E+00, 2.704928307360539E+01,-1.899669293113739E+00}
    };
    G.setX0(asteroid);
    auto asteroid_best = asteroid;
    if (G.start())
    {
        asteroid_best = G.getXbest();
    }
    std::cout << "#########################################################\n";
    std::cout << std::fixed << std::setprecision(8);
    std::cout << "x0:\n";
    std::cout << "Position: (" << asteroid.x.x << ", "
                  << asteroid.x.y << ", " << asteroid.x.z << ")\n";
    std::cout << "Velocity: (" << asteroid.v.x << ", "
              << asteroid.v.y << ", " << asteroid.v.z << ")\n";


    std::cout << "Best x0:\n";
    std::cout << "Position: (" << asteroid_best.x.x << ", "
                  << asteroid_best.x.y << ", " << asteroid_best.x.z << ")\n";
    std::cout << "Velocity: (" << asteroid_best.v.x << ", "
              << asteroid_best.v.y << ", " << asteroid_best.v.z << ")\n";
    write(a, eph, obs, r, start_jd, asteroid, "result_x0.txt");
    write(a, eph, obs, r, start_jd, asteroid_best, "result_xbest_GN.txt");
    std::cout << "Результаты интегрирования сохранены в файлы:\n";
    std::cout << "result_x0.txt - с начальными условиями\n";
    std::cout << "result_xbest_GN.txt - после применение GaussNewton\n";
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
    test(a, eph, obs, r);
    // solve(eph, r, obs, a);
}
