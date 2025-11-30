#include <iostream>
#include "include/data_reader.hpp"
#include "ephemeris.hpp"
#include "rk4Integrator.hpp"
#include "observatories.hpp"

#include <fstream>
#include <iomanip>

void test1()
{
    std::cout << "TEST_1 ##################################\n";

    Ephemeris eph;
    eph.loadFile("data/de440s.bsp");

    Observatories obs(eph);
    obs.loadFileObs("data/ObsCodes.html");
    obs.loadFileEop("data/finals.data.iau2000.txt");
    
    Rk4Integrator r(eph);

    double start_jd = 2459172.5;

    BodyVector asteroid = {
        {-4.985939512102024E+07, -2.426095591496699E+07, -1.027706231518984E+07},

        {3.914031781284947E+01, -4.393811128598539E+01, 1.898563558065064E+01}
    };


    r.integrateOrbit(asteroid, start_jd, 120, 100);
    
    double cur_jd = 2459272.5;

    auto comp_pos_ast_1 = r.interpolatePosition(cur_jd);
    auto comp_pos_obs_1 = obs.getConvertObsForJD(cur_jd, "145"); // 145 - 's-Gravenwezel (observatory)

    auto v = comp_pos_ast_1 - comp_pos_obs_1;
    double Ra, Dec, dist;
    r.cartToRaDec(v, Ra, Dec, dist);
    auto ltc = r.applyAstrometricCorrections(cur_jd, comp_pos_obs_1);
    double Ra1, Dec1, dist1;
    r.cartToRaDec(ltc, Ra1, Dec1, dist1);
    std::cout << std::fixed << std::setprecision(8);

    std::cout << "Start :" << start_jd << "\n";
    std::cout << "Cur : " << cur_jd << "\n";
    std::cout << "Exact:                        X = 2.811419190134476E+08 Y = 8.917364246930613E+07 Z =-1.221687550216738E+06\n";
    std::cout << "Compute:                      X = " << v.x << " Y = " << v.y << " Z = " << v.z << "\n";
    std::cout << "Compute(corrections)          X = " << ltc.x << " Y = " << ltc.y << " Z = " << ltc.z << "\n\n";

    std::cout << "Exact                        RA = 17.62° Dec = -0.2373° Distance = 294,950,000\n";
    std::cout << "Compute                      RA = " << " " << Ra << " Dec = " << Dec << " Distance = " << dist << "\n";  
    std::cout << "Compute(corrections)         RA = " << " " << Ra1 << " Dec = " << Dec1 << " Distance = " << dist1 << "\n";  


}

void test2()
{
    std::cout << "TEST_2 ##################################\n";
    Ephemeris eph;
    eph.loadFile("data/de440s.bsp");

    Observatories obs(eph);
    obs.loadFileObs("data/ObsCodes.html");
    obs.loadFileEop("data/finals.data.iau2000.txt");
    
    Rk4Integrator r(eph);

    double start_jd = 2459172.5;

    BodyVector asteroid = {
        {-4.985939512102024E+07, -2.426095591496699E+07, -1.027706231518984E+07},

        {3.914031781284947E+01, -4.393811128598539E+01, 1.898563558065064E+01}
    };


    r.integrateOrbit(asteroid, start_jd, 120, 100);
    
    double cur_jd = 2459180.5;

    auto comp_pos_ast_1 = r.interpolatePosition(cur_jd);

    auto comp_pos_obs_1 = obs.getConvertObsForJD(cur_jd, "145"); // 145 - 's-Gravenwezel (observatory)
    std::cout << std::fixed << std::setprecision(10);
    std::cout << 2453322.61435  << " " << obs.convertUtcToTdb(2453322.61435, "991") << "\n"; 
    auto v = comp_pos_ast_1 - comp_pos_obs_1;
    double Ra, Dec, dist;
    r.cartToRaDec(v, Ra, Dec, dist);
    auto ltc = r.applyAstrometricCorrections(cur_jd, comp_pos_obs_1);
    double Ra1, Dec1, dist1;
    r.cartToRaDec(ltc, Ra1, Dec1, dist1);
    std::cout << std::fixed << std::setprecision(8);
    std::cout << "Start :" << start_jd << "\n";
    std::cout << "Cur : " << cur_jd << "\n";
    std::cout << "Exact:                        X = -7.627297554771891E+07 Y = -1.705950656977384E+08 Z =-4.948204867554787E+07\n";
    std::cout << "Compute:                      X = " << v.x << " Y = " << v.y << " Z = " << v.z << "\n";
    std::cout << "Compute(corrections)          X = " << ltc.x << " Y = " << ltc.y << " Z = " << ltc.z << "\n\n";

    std::cout << "Exact                        RA = 245.20° Dec =-15.04° Distance = 190,710,000\n";
    std::cout << "Compute                      RA = " << " " << Ra << " Dec = " << Dec << " Distance = " << dist << "\n";  
    std::cout << "Compute(corrections)         RA = " << " " << Ra1 << " Dec = " << Dec1 << " Distance = " << dist1 << "\n";
}

void solve(Rk4Integrator &r, Observatories &obs, std::vector<RaDec> &a)
{
    double start_jd = 2453320;
    BodyVector asteroid = {
        {1.367372443402559E+08, 6.602611574848668E+07,3.033666479027805E+07},

        {8.299501903651528E+00, 2.704928307360539E+01,-1.899669293113739E+00}
    };

    r.integrateOrbit(asteroid, start_jd, 12 * 365, 3600);
    {
        std::size_t i = 1;
        std::ofstream file("result.txt");
        for (auto [cur_jd_utc, Ra, Dec, code] : a)
        {
            file << "observations " << i << " JD " << cur_jd_utc  <<  ": \n";
            double cur_jd_tdb = obs.convertUtcToTdb(cur_jd_utc, code);
            auto comp_pos_obs  = obs.getConvertObsForJD(cur_jd_tdb, code);
            auto ltc = r.applyAstrometricCorrections(cur_jd_tdb, comp_pos_obs);
            double Ra1, Dec1, dist1;
            r.cartToRaDec(ltc, Ra1, Dec1, dist1);

            file << std::fixed << std::setprecision(8);
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
    solve(r, obs, a);
}
