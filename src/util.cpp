
#include "util.hpp"



namespace util 
{
    void sphToCart(double Long, double cos, double sin, double &x, double &y, double &z)
    {
        x = cos * std::cos(Long);
        y = cos * std::sin(Long);
        z = sin;
    }

    void convertToRaDec(const std::string &ra,const std::string &dec, double &Ra, double &Dec)
    {
        double ra_h = std::stod(ra.substr(0, 2));
        double ra_m = std::stod(ra.substr(3, 2));  
        double ra_s = std::stod(ra.substr(6, 6));            
             
        double dec_d = std::stod(dec.substr(0, 3)); 
        double dec_m = std::stod(dec.substr(4, 2)); 
        double dec_s = std::stod(dec.substr(7, 5));          
        
    Ra = (ra_h + ra_m/60.0 + ra_s/3600.0) * 15.0;
    
    if (dec_d < 0)
    {
        Dec = dec_d - dec_m/60.0 - dec_s/3600.0;
    }
    else
    {
        Dec = dec_d + dec_m/60.0 + dec_s/3600.0;
    }
    
    if (Ra >= 360.0) Ra -= 360.0;
    if (Ra < 0.0) Ra += 360.0;    
    }

    void dateToJd(int year, int month, double day, double &jd)
    {
        if (month <= 2)
        {
            year -= 1;
            month += 12;
        }
        
        int A = year / 100;
        int B = 2 - A + A/4;
        
        int day_int = (int)day;
        double day_frac = day - day_int;
        
        jd = (int)(365.25 * (year + 4716)) + 
                    (int)(30.6001 * (month + 1)) + 
                    day_int + B - 1524.5;
        
        jd += day_frac;
    }
}



