// @HEADER
// *****************************************************************************
//     Compadre: COMpatible PArticle Discretization and REmap Toolkit
//
// Copyright 2018 NTESS and the Compadre contributors.
// SPDX-License-Identifier: BSD-2-Clause
// *****************************************************************************
// @HEADER
#ifndef _GMLS_MANIFOLD_HPP_
#define _GMLS_MANIFOLD_HPP_

#include <Kokkos_Core.hpp>
#include <cmath>

#define PI 3.14159265358979323846

KOKKOS_INLINE_FUNCTION
double device_max(double d1, double d2) {
    return (d1 > d2) ? d1 : d2;
}

KOKKOS_INLINE_FUNCTION
double atan4(const double y, const double x) { 
    double result = 0.0;
    if (x == 0.0)
    {
        if (y > 0.0) 
            result = 0.5 * PI;
        else if ( y < 0.0 )
            result = 1.5 * PI;
        else if ( y == 0.0 )
            result = 0.0;
    }
    else if (y == 0)
    {
        if (x > 0.0)
            result = 0.0;
        else if ( x < 0.0 )
            result = PI;
    }
    else
    {
        double theta = std::atan2( std::abs(y), std::abs(x) );
        if (x > 0.0 && y > 0.0)
            result = theta;
        else if ( x < 0.0 && y > 0.0 ) 
            result = PI - theta;
        else if ( x < 0.0 && y < 0.0 ) 
            result = PI + theta;
        else if ( x > 0.0 && y < 0.0 )
            result = 2.0 * PI - theta;
    }
    return result;
}

KOKKOS_INLINE_FUNCTION
double latitude(double x, double y, double z) {
    return std::atan2(z, std::sqrt( x*x + y*y));
}

KOKKOS_INLINE_FUNCTION
double longitude(double x, double y, double z) {
    return atan4(y, x);
}

KOKKOS_INLINE_FUNCTION
double legendre54(double z) {
    return z * (  z * z - 1.0 ) * ( z * z - 1.0 );
}

KOKKOS_INLINE_FUNCTION
double sphere_harmonic54(double x, double y, double z) {
    const double lon = longitude(x, y, z);
    return std::cos(4.0 * lon) * legendre54(z);
}

KOKKOS_INLINE_FUNCTION
void curl_sphere_harmonic54(double *curl, double x, double y, double z) {
    const scalar_type lon = longitude(x, y, z); // theta
    const scalar_type lat = acos(z); // phi
    const scalar_type sigma_lon_comp = std::pow(sin(lat), 2) * (5.0* std::pow(cos(lat), 2) - 1.0) * cos(4.0 *lon);
    const scalar_type sigma_lat_comp = 4*cos(lat) * std::pow(sin(lat), 3) * sin(4.0 * lon);
    // solution oriented for inward normal, so we flip sign for outward
    curl[0] = -(-sin(lat)*sin(lon)*sigma_lon_comp + cos(lat)*cos(lon)*sigma_lat_comp);
    curl[1] = -(sin(lat)*cos(lon)*sigma_lon_comp + cos(lat)*sin(lon)*sigma_lat_comp);
    curl[2] = -(-sin(lat)*sigma_lat_comp);
}


KOKKOS_INLINE_FUNCTION
double laplace_beltrami_sphere_harmonic54(double x, double y, double z) {
    const double lon = longitude(x, y, z);
    return -30 * std::cos(4.0 * lon) * legendre54(z);
}

KOKKOS_INLINE_FUNCTION
void gradient_sphereHarmonic54_local(double *gradient, double x, double y, double z) {
    const double lat = latitude(x, y, z); // phi
    const double lon = longitude(x, y, z); // lambda

    const double A = -4.0 * std::pow(std::cos(lat),3) * std::sin(4.0 * lon) * std::sin(lat);
    const double B = 0.5* std::cos(4.0 * lon) * std::pow(std::cos(lat),3) * ( 5 * std::cos(2.0 * lat) - 3.0 );

    gradient[0] = A;
    gradient[1] = B;
}

KOKKOS_INLINE_FUNCTION
void gradient_sphereHarmonic54_ambient(double *gradient, double x, double y, double z) {
    const double lat = latitude(x, y, z); // phi
    const double lon = longitude(x, y, z); // lambda

    const double A = -4.0 * std::pow(std::cos(lat),3) * std::sin(4.0 * lon) * std::sin(lat);
    const double B = 0.5* std::cos(4.0 * lon) * std::pow(std::cos(lat),3) * ( 5 * std::cos(2.0 * lat) - 3.0 );

    gradient[0] = -A * std::sin(lon) - B * std::sin(lat) * std::cos(lon);
    gradient[1] = A * std::cos(lon) - B * std::sin(lat) * std::sin(lon);
    gradient[2] = B * std::cos(lat);
}

KOKKOS_INLINE_FUNCTION
void velocity_sphereHarmonic54_ambient(double *velocity, double x, double y, double z) {
    const double lat = latitude(x, y, z); // phi
    const double lon = longitude(x, y, z); // lambda

    const double U = 0.5* std::cos(4.0 * lon) * std::pow(std::cos(lat),3) * ( 5 * std::cos(2.0 * lat) - 3.0 );
    const double V = 4.0 * std::pow(std::cos(lat),3) * std::sin(4.0 * lon) * std::sin(lat);

    velocity[0] = -U * std::sin(lon) - V * std::sin(lat) * std::cos(lon);
    velocity[1] = U * std::cos(lon) - V * std::sin(lat) * std::sin(lon);
    velocity[2] = V * std::cos(lat);
}



/** Manifold GMLS Example 
 *
 *  Exercises GMLS operator evaluation with data over various orders and numbers of targets for targets including point evaluation, Laplace-Beltrami, gradient and gradient on a manifold.
 */
int main (int argc, char* args[]);

/**
 * \example "Manifold GMLS Tutorial" based on GMLS_Manifold.cpp
 * \section ex GMLS Example with Device Views
 *
 * This tutorial sets up a batch of GMLS problems, solves the minimization problems, and applies the coefficients produced to data.
 * 
 * \section ex1a Parse Command Line Arguments
 * \snippet GMLS_Manifold.cpp Parse Command Line Arguments
 *
 * \section ex1b Setting Up The Point Cloud
 * \snippet GMLS_Manifold.cpp Setting Up The Point Cloud
 *
 * \section ex1c Performing Neighbor Search
 * \snippet GMLS_Manifold.cpp Performing Neighbor Search
 *
 * \section ex2 Creating The Data
 * \snippet GMLS_Manifold.cpp Creating The Data
 *
 * \section ex3 Setting Up The GMLS Object
 * \snippet GMLS_Manifold.cpp Setting Up The GMLS Object
 * 
 * \section ex4 Apply GMLS Alphas To Data
 * \snippet GMLS_Manifold.cpp Apply GMLS Alphas To Data
 *
 * \section ex5 Check That Solutions Are Correct
 * \snippet GMLS_Manifold.cpp Check That Solutions Are Correct
 *
 * \section ex6 Finalize Program
 * \snippet GMLS_Manifold.cpp Finalize Program
 */ 

#endif
