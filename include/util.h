#pragma once


namespace util 
{
    /**
     *  \brief Converting spherical coordinates to Cartesian coordinates
     *
     *  \param[in] Long - longitude (in degrees east of Greenwich)
     *  \param[in] cos - rho cos phi`
     *  \param[in] sin - rho sin phi'
     *
     *  where phi' is the geocentric latitude and rho is the geocentric distance in earth radii
     *  
     *  Cartesian coordinates
     *  \param[out] x 
     *  \param[out] y
     *  \param[out] z
     */

    void sphToCart(double Long, double cos, double sin, double &x, double &y, double &z);







}

