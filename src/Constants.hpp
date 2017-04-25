/******************************************************************************
 *  File: Constants.hpp
 *  A bunch of constants used in the radiation solver
 *
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 ******************************************************************************/
#ifndef FLIP_SIM_CONSTANTS_H
#define FLIP_SIM_CONSTANTS_H

#include <sys/time.h>
#include "PrecisionDef.h"

// Physical constants
const _T GRAVITY                 = 9.8;
const _T WATER_DENSITY           = 1000;
const _T AIR_DENSITY             = 1.184;
const _T ATM                     = 101325;                      /* 1 atmosphere presure 101325 (pa) */
const _T GAS_CONSTANT            = 8.314472;
const _T HEAT_RATIO              = 1.4;                         /* specific heat ratio (gamma) */
const _T SURFACE_TENSION_COEFF   = 0.0726;                      /* N/m water at 293K */
const _T SURFACE_TENSION_COEFF_2 = SURFACE_TENSION_COEFF * 2.0;
const _T AIR_THERMAL_DIFFUSIVITY = 2.122169236803534e-05;       /* m^2/s at 293K */
const _T SOUND_SPEED_WATER       = 1496;                        /* m/s speed in fresh water */
const _T SOUND_SPEED_AIR         = 343;
const _T WATER_VISCOSITY         = 8.9e-04;                     /* water viscosity coefficient in 25'C */
const _T AIR_MOLAR_MASS          = 0.029;                       /* molar mass of air 29.0 g/molar */
const _T DRAG_COEFFI             = 1.;

// the max number of sources for multipole approximation
#define D_ORDER_IN   2
#define D_ORDER_OUT  2

const int MAX_SRC_IN             = 40;          /* for interior solver */
const int MAX_SRC_OUT            = 145;         /* for exterior solver */

const struct timeval RPC_CONFIG_TIMEOUT = {180, 0};
const struct timeval RPC_SOLVE_TIMEOUT  = {7200, 0}; /* two hours */

// ****************************************************************************
// Parameters need to be changed for different simulation configuration
const unsigned int LATTICE_NX = 70;
const unsigned int LATTICE_NY = 90;
const unsigned int LATTICE_NZ = 70;
const _T            CELL_SIZE = 0.002;

/*
 * This value is only used for adaptive sampling.
 * See sec. 5.7 in the paper and the implementation in 
 * file Bubble.hpp
 */
const _T BUB_AVG_TIMESTEP = 0.0005;

#endif

