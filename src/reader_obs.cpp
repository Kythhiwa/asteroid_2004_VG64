#include <fstream>

#include <iostream>

#include "reader_obs.h"



std::map<std::string, dataObs> read_obs(std::string filename)
{
    std::ifstream file(filename);
    if (!file.is_open())
    {
        std::cerr << "[ERROR] open file\n";
    }
    std::string trash;
    if (!std::getline(file, trash))
    {
        std::cerr << "[ERROR] getline\n";
    }
    std::string names;
    if (!std::getline(file, names))
    {
        std::cerr << "[ERROR] getline\n";
    }

    std::map<std::string, dataObs> data;

    std::string cur;
    while(std::getline(file, cur))
    {
        std::string code;
        double Long, cos, sin;
        try 
        {
            code = cur.substr(0, 3);
            Long = std::stod(cur.substr(3, 10));
            cos = std::stod(cur.substr(13, 8));
            sin = std::stod(cur.substr(21, 9));
        }
        catch (...)
        {
            std::cerr << "[SKIP LINE]\n";
            continue;
        }
        data.emplace(code, dataObs{Long, cos, sin});
        
    }
    file.close();
    return data;
}
