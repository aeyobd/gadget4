#ifndef GRAV_MILKYWAY_H
#define GRAV_MILKYWAY_H

#include "gadgetconfig.h"

#ifdef EXTERNALGRAVITY_MW

#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../data/allvars.h"
#include "../data/simparticles.h"




vector<double> MilkyWayAcceleration(vector<double> pos);
double MilkyWayPotential(vector<double> pos);
double HernquistAcceleration(double r, double a);
double HernquistPotential(double r, double a);
double NFW_g(double c);// these all assume GM = 1 for simplicity
//
double NFWAcceleration(double r, double R, double c);

double NFWPotential(double r, double R, double c);
double DiskScale_r(double R, double z, double a, double b);
double DiskPotential(double R, double z, double a, double b);
double DiskAcceleration_R(double R, double z, double a, double b);
double DiskAcceleration_Z(double R, double z, double a, double b);


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

#endif // GRAVEXTERNAL_MW

#endif /* GRAV_MILKYWAY_H */
