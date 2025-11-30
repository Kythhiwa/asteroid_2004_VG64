

#include "data_reader.hpp"



std::map<std::string, Vector3D> read_obs(std::string filename)
{
    std::ifstream file(filename);
    if (!file.is_open())
    {
        throw std::runtime_error("[ERROR DATA OBS] open file\n");
    }
    std::string trash;
    if (!std::getline(file, trash))
    {
        throw std::runtime_error("[ERROR DATA OBS] getline\n");
    }
    std::string names;
    if (!std::getline(file, names))
    {
        throw std::runtime_error("[ERROR DATA OBS] getline\n");
    }

    std::map<std::string, Vector3D> data;

    std::string cur;
    int i = 0;
    while(std::getline(file, cur))
    {
        std::string code;
        double Long, cos, sin;
        try 
        {
            code = cur.substr(0, 3);
            Long = std::stod(cur.substr(3, 10)) * M_PI / 180;
            cos = std::stod(cur.substr(13, 8));
            sin = std::stod(cur.substr(21, 9));
            double x, y, z;
            util::sphToCart(Long, cos, sin, x, y, z);
            data.emplace(code, Vector3D{x * 6371.0, y * 6371.0, z * 6371.0});
        }
        catch (...)
        {
            // std::cerr << "[SKIP LINE DATA OBS]\n";
        }
    }
    file.close();
    return data;
}


std::set<EOPEntry> read_eop(std::string filename)
{
    std::ifstream file(filename);
    if (!file.is_open())
    {
       throw std::runtime_error("[ERROR DATA EOP] open file\n");
    }

    std::set<EOPEntry> data;
    
    std::string cur;
    int i = 0;
    while(std::getline(file, cur))
    {
        double mjd;
        double ut1_utc;
        double xp;
        double yp;
        try
        {
            mjd = std::stod(cur.substr(7, 8));
            ut1_utc = std::stod(cur.substr(58, 10));
            xp = std::stod(cur.substr(18, 9));
            yp = std::stod(cur.substr(37, 9));
        }
        catch (...)
        {
            // std::cerr << "[SKIP LINE DATA EOP]\n";
            continue;
        }
        data.insert(EOPEntry{mjd, ut1_utc, xp * DAS2R, yp * DAS2R});
    }
    file.close();
    return data;
}


std::vector<RaDec> read_observations(std::string filename)
{
    std::ifstream file(filename);
    if (!file.is_open())
    {
       throw std::runtime_error("[ERROR DATA OBSERVATIONS] open file\n");
    }

    std::vector<RaDec> data;

    std::string cur;
    int i = 0;
    while(std::getline(file, cur))
    {
        double jd_utc;
        double Ra;
        double Dec;
        std::string code;
        try 
        {
            std::string date = cur.substr(15, 16);
            int year = std::stoi(date.substr(0, 4));
            if (year > 2007)
            {
                continue;
            }
            int month = std::stoi(date.substr(5, 2));
            double day = std::stod(date.substr(8, 8));
            util::dateToJd(year, month, day, jd_utc);
            code = cur.substr(77, 3);
            util::convertToRaDec(cur.substr(32, 11), cur.substr(44,11), Ra, Dec);
        }
        catch (...)
        {
            continue;
        }

        data.emplace_back(jd_utc, Ra, Dec, code);
    }

    return data;
}
