/*******************************************************************************
 * Milkyway potential
 * Defined to have minimum dependencies
 *******************************************************************************/

#include "gadgetconfig.h"

#ifdef EXTERNALGRAVITY
#ifdef EXTERNALGRAVITY_MW

#include "grav_milkyway.h"
#include "grav_profiles.h"




#include <math.h>
// #include <mpi.h>
#include <stdlib.h>

#include "../data/allvars.h"
#include "../data/dtypes.h"
// #include "../data/intposconvert.h"
// #include "../data/mymalloc.h"
// #include "../fmm/fmm.h"
// #include "../logs/logs.h"
// #include "../logs/timer.h"
// #include "../main/simulation.h"
// #include "../mpi_utils/mpi_utils.h"
// #include "../pm/pm.h"
// #include "../system/system.h"
// #include "../time_integration/timestep.h"



vector<double> MilkyWayAcceleration(vector<double> pos)
{
    vector<double> acceleration{0., 0., 0.};


    vector<double> a_bulge = All.MWBulgeMass * HernquistAcceleration(pos, All.MWBulge_A);
    vector<double> a_halo = All.MWHaloMass * NFWAcceleration(pos, All.MWHalo_R, All.MWHalo_C);
    vector<double> a_thin = All.MWThinMass * DiskAcceleration(pos, All.MWThin_A, All.MWThin_B);
    vector<double> a_thick = All.MWThickMass * DiskAcceleration(pos, All.MWThick_A, All.MWThick_B);

    acceleration =  All.G * (a_bulge + a_halo + a_thin + a_thick);

    return acceleration;

}


double MilkyWayPotential(vector<double> pos)
{
    double BulgePotential = All.MWBulgeMass * HernquistPotential(pos, All.MWBulge_A);
    double HaloPotential = All.MWHaloMass * NFWPotential(pos, All.MWHalo_R, All.MWHalo_C);
    double ThinPotential = All.MWThinMass * DiskPotential(pos, All.MWThin_A, All.MWThin_B);
    double ThickPotential =All.MWThickMass * DiskPotential(pos, All.MWThick_A, All.MWThick_B);

    return All.G * ( 
            BulgePotential + HaloPotential + ThinPotential + ThickPotential
            );
}





#endif
#endif
