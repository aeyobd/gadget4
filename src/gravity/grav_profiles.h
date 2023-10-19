#ifndef GRAV_PROFILES_H
#define GRAV_PROFILES_H

#include "gadgetconfig.h"
#include <math.h>
#include <mpi.h>
#include <stdlib.h>
#include <string.h>



#include "../data/allvars.h"
#include "../data/dtypes.h"

#include "../data/intposconvert.h"
#include "../data/mymalloc.h"
#include "../domain/domain.h"
#include "../fmm/fmm.h"
#include "../gravtree/gravtree.h"
#include "../logs/logs.h"
#include "../logs/timer.h"
#include "../main/simulation.h"
#include "../mpi_utils/mpi_utils.h"
#include "../pm/pm.h"
#include "../system/system.h"
#include "../time_integration/timestep.h"



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

inline double r_sphere(vector<double> pos) {
    return sqrt(pos.r2());
}

inline double r_cyl(vector<double> pos) {
    double x = pos[0];
    double y = pos[1];
    return norm(x, y);
}

inline double z_cyl(vector<double> pos) {
    return pos[2];
}



double HernquistScalarAcceleration(double r, double a) {
    if (r > 0)
        return 1. / r / square(r + a);
    else
        return 0;
}

double HernquistPotential(double r, double a) {
    return -1.0 / (r + a);
}

double HernquistPotential(vector<double> pos, double a) {
    return HernquistPotential(r_sphere(pos), a);
}

double NFW_g(double c) {
    return 1. / (log(1.+c) - c/(1.+c));
}
// these all assume GM = 1 for simplicity
//
double NFWScalarAcceleration(double r, double R, double c) {
    if (r > 0)
        return  1. / square(r) * NFW_g(c) * ( 1/r * log(r/R + 1) - 1/(r + R) );
    else
        return 0;
}



double NFWPotential(double r, double R, double c) {
    return  -1/r * NFW_g(c) * log( 1 + r/R);
}

double NFWPotential(vector<double> pos, double R, double c) {
    return NFWPotential(r_sphere(pos), R, c);
}

double DiskScale_r(double R, double z, double a, double b) {
    return norm(R, a + norm(z, b));
}

double DiskPotential(double R, double z, double a, double b) {
    return -1 / DiskScale_r(R, z, a, b);
}

double DiskPotential(vector<double> pos, double a, double b) {
    return DiskPotential(r_cyl(pos), z_cyl(pos), a, b);
}

double DiskScalarAcceleration_R(double R, double z, double a, double b) {
    if (DiskScale_r(R, z, a, b) > 0)
        return 1 / cube(DiskScale_r(R, z, a, b));
    else
        return 0;
}

double DiskScalarAcceleration_z(double R, double z, double a, double b) {
    if (norm(z, b) > 0)
        return  ( 1 + a / norm(z, b) ) / cube(DiskScale_r(R, z, a, b));
    else
        return 0;
}



vector<double> NFWAcceleration(vector<double> pos, double R, double c) {
    double r = r_sphere(pos);
    double a = -NFWScalarAcceleration(r, R, c);
    return  a * pos;
}


vector<double> HernquistAcceleration(vector<double> pos, double a) {
    double r = sqrt(pos.r2());
    double acc = -HernquistScalarAcceleration(r, a);
    return  acc * pos;
}


vector<double> DiskAcceleration(vector<double> pos, double a, double b) {
    double x = pos[0];
    double y = pos[1];
    double z = pos[2];
    double R = norm(x, y);

    vector<double> z_vec = {0,0,z};
    vector<double> R_vec = {x,y,0};
    double a_R = -DiskScalarAcceleration_R(R, z, a, b);
    double a_z = -DiskScalarAcceleration_z(R, z, a, b);
    return  a_R*R_vec + a_z*z_vec;
}

#endif
