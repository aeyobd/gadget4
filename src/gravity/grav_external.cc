/*******************************************************************************
 * \copyright   This file is part of the GADGET4 N-body/SPH code developed
 * \copyright   by Volker Springel. Copyright (C) 2014-2020 by Volker Springel
 * \copyright   (vspringel@mpa-garching.mpg.de) and all contributing authors.
 *******************************************************************************/

/*! \file gravity_external.cc
 *
 * \brief can add an optional external gravity field
 */

#include "gadgetconfig.h"

#ifdef EXTERNALGRAVITY

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

void sim::gravity_external(void)
{
#ifdef PERIODIC
  // pick middle of (stretched) box, but could also choose other point
  vector<double> pos_center{0.5 * All.BoxSize / LONG_X, 0.5 * All.BoxSize / LONG_Y, 0.5 * All.BoxSize / LONG_Z};
#else
  // here pick origin
  vector<double> pos_center{0, 0, 0};
#endif
  MyIntPosType intpos_center[3];
  Sp.pos_to_intpos(pos_center.da, intpos_center);

  for(int i = 0; i < Sp.TimeBinsGravity.NActiveParticles; i++)
    {
      int target = Sp.TimeBinsGravity.ActiveParticleList[i];

#if defined(EVALPOTENTIAL) || defined(OUTPUT_POTENTIAL)
      Sp.P[target].ExtPotential = 0;
#endif

#ifdef EXTERNALGRAVITY_MW
      {
        vector<double> pos;
        Sp.nearest_image_intpos_to_pos(Sp.P[target].IntPos, intpos_center, pos.da);

        double r = sqrt(pos.r2());
        double x = pos[0];
        double y = pos[1];
        double R = sqrt(x*x + y*y);
        double z = pos[2];
        vector<double> z_hat (0,0,1);
        vector<double> R_hat (x/R, y/R, 0);
        vector<double> r_hat = 1/r * pos;
        double thin_S = sqrt(All.MWThin_B*All.MWThin_B * z*z);
        double thin_D = sqrt(R*R + pow(All.MWThin_A + thin_S, 2));
        double thick_S = sqrt(All.MWThick_B*All.MWThick_B * z*z);
        double thick_D = sqrt(R*R + pow(All.MWThick_A + thick_S, 2));


        if(r > 0) 
        {
          double a_bulge = - All.G * All.MWBulgeMass  
              / pow((r + All.MWBulge_A), 2);

          double a_halo = -1/r * All.G * All.MWHaloMass 
              / (log(2) - 0.5) 
              * ( 1/r * log(r/All.MWHalo_R + 1) - 1/(r+R) );

          double a_thin_r = - All.G * All.MWThinMass * R 
              / (thin_D*thin_D*thin_D);

          double a_thin_z = - All.G * All.MWThinMass * z 
              / (thin_D*thin_D*thin_D) * (All.MWThin_A/thin_S + 1);


          double a_thick_r = - All.G * All.MWThickMass * R 
              / (thick_D*thick_D*thick_D);
          double a_thick_z = -All.G * All.MWThickMass * z 
              / (thick_D*thick_D*thick_D) * (All.MWThick_A/thick_S + 1);

          Sp.P[target].GravAccel += a_bulge*r_hat + a_halo *r_hat 
              + a_thin_r * R_hat + a_thin_z * z_hat
              + a_thick_r * R_hat + a_thick_z * z_hat;
        }

#if defined(EVALPOTENTIAL) || defined(OUTPUT_POTENTIAL)
        double BulgePotential = - All.G * All.MWBulgeMass / (r + All.MWBulge_A);
        double HaloPotential = - All.G * All.MWHaloMass / (r * (log(2) - 0.5)) * log(1 + r/All.MWHalo_R);
        double ThinPotential = - All.G * All.MWThinMass / thin_D;
        double ThickPotential = - All.G * All.MWThickMass / thick_D;
        Sp.P[target].ExtPotential += (-All.G * All.MWBulgeMass / (r + All.MWBulge_A));
#endif
      }

#endif

#ifdef EXTERNALGRAVITY_STATICHQ
      {
        vector<double> pos;
        Sp.nearest_image_intpos_to_pos(Sp.P[target].IntPos, intpos_center, pos.da);

        double r = sqrt(pos.r2());

        double m = All.Mass_StaticHQHalo * pow(r / (r + All.A_StaticHQHalo), 2);

        if(r > 0)
          Sp.P[target].GravAccel += (-All.G * m / (r * r * r)) * pos;

#if defined(EVALPOTENTIAL) || defined(OUTPUT_POTENTIAL)
        Sp.P[target].ExtPotential += (-All.G * All.Mass_StaticHQHalo / (r + All.A_StaticHQHalo));
#endif
      }
#endif
    }
}

#endif
