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

#endif // GRAVEXTERNAL_MW

#endif /* GRAV_MILKYWAY_H */
