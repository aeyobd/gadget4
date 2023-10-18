/*******************************************************************************
 * Milkyway potential
 * Defined to have minimum dependencies
 *******************************************************************************/

#include "gadgetconfig.h"

#ifdef EXTERNALGRAVITY
#ifdef EXTERNALGRAVITY_MW

#include "grav_milkyway.h"




#include <math.h>
#include <mpi.h>
#include <stdlib.h>

#include "../data/allvars.h"
#include "../data/dtypes.h"
#include "../data/intposconvert.h"
#include "../data/mymalloc.h"
#include "../fmm/fmm.h"
#include "../logs/logs.h"
#include "../logs/timer.h"
#include "../main/simulation.h"
#include "../mpi_utils/mpi_utils.h"
#include "../pm/pm.h"
#include "../system/system.h"
#include "../time_integration/timestep.h"



vector<double> MilkyWayAcceleration(vector<double> pos)
{
    double x = pos[0];
    double y = pos[1];
    double z = pos[2];
    double r = sqrt(pos.r2());
    double R = sqrt(x*x + y*y);
    vector<double> z_vec{0, 0, z};
    vector<double> R_vec{x, y, 0};
    vector<double> acceleration{0., 0., 0.};


    if(r > 0) 
    {
      double a_bulge = - All.G * All.MWBulgeMass * HernquistAcceleration(r, All.MWBulge_A);
      double a_halo = - All.G * All.MWHaloMass * NFWAcceleration(r, All.MWHalo_R, All.MWHalo_C);
      double a_thin_R = - All.G * All.MWThinMass * DiskAcceleration_R(R, z, All.MWThin_A, All.MWThin_B);
      double a_thin_z = - All.G * All.MWThinMass * DiskAcceleration_Z(R,  z,All.MWThin_A, All.MWThin_B);
      double a_thick_R = - All.G * All.MWThickMass * DiskAcceleration_R(R, z, All.MWThick_A, All.MWThick_B);
      double a_thick_z = - All.G * All.MWThickMass * DiskAcceleration_Z(R, z, All.MWThick_A, All.MWThick_B);

      acceleration =  a_bulge*pos + a_halo *pos 
          + a_thin_R * R_vec + a_thin_z * z_vec
          + a_thick_R * R_vec + a_thick_z * z_vec;
    }

    return acceleration;

}


double MilkyWayPotential(vector<double> pos)
{
    double x = pos[0];
    double y = pos[1];
    double z = pos[2];
    double r = sqrt(pos.r2());
    double R = sqrt(x*x + y*y);

      double BulgePotential = - All.G * All.MWBulgeMass * HernquistPotential(r, All.MWBulge_A);
      double HaloPotential = - All.G * All.MWHaloMass * NFWPotential(r, All.MWHalo_R, All.MWHalo_C);
      double ThinPotential = - All.G * All.MWThinMass * DiskPotential(R, z, All.MWThin_A, All.MWThin_B);
      double ThickPotential = - All.G * All.MWThickMass * DiskPotential(R, z, All.MWThick_A, All.MWThick_B);
    return BulgePotential + HaloPotential + ThinPotential + ThickPotential;
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
#endif
