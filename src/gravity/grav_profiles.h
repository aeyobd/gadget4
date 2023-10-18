#ifndef GRAV_PROFILES_H
#define GRAV_PROFILES_H

#include "gadgetconfig.h"
#include "math.h"
// note: all of these are normalized to GM = 1
// accelerations are normalized and positive so need multiplied by 
// unnormalized position vector



inline double square(double x) {
    return x*x;
}
inline double cube(double x) {
    return x*x*x;
}

inline double norm(double x, double y) {
    return sqrt(square(x) + square(y));
}

inline double norm(double x, double y, double z) {
    return sqrt(square(x) + square(y) + square(z));
}

double HernquistAcceleration(double r, double a) {
    return 1 / r / square(r + a);
}

double HernquistPotential(double r, double a) {
    return 1.0 / (r + a);
}

double NFW_g(double c) {
    return 1 / (log(1+c) - c/(1+c));
}
// these all assume GM = 1 for simplicity
//
double NFWAcceleration(double r, double R, double c) {
    return  1/square(r) * NFW_g(c) * ( 1/r * log(r/R + 1) - 1/(r + R) );
}


double NFWPotential(double r, double R, double c) {
    return  1/r * NFW_g(c) * log( 1 + r/R);
}

double DiskScale_r(double R, double z, double a, double b) {
    return norm(R, a + norm(z, b));
}

double DiskPotential(double R, double z, double a, double b) {
    return 1 / DiskScale_r(R, z, a, b);
}

double DiskAcceleration_R(double R, double z, double a, double b) {
    return 1 / cube(DiskScale_r(R, z, a, b));
}

double DiskAcceleration_Z(double R, double z, double a, double b) {
    return  ( 1 + a / norm(z, b) ) / cube(DiskScale_r(R, z, a, b));
}


#endif
