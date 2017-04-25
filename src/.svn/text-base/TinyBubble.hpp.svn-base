/******************************************************************************
 *  File: TinyBubble.hpp
 *  Implement the vibration model of the acoustic bubble
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
#ifndef TINY_BUBBLE_H
#   define TINY_BUBBLE_H

#include <math.h>
#include "utils/math.hpp"
#include "Constants.hpp"
#include "linearalgebra/Vector3.hpp"

/*!
 * Tiny bubble vibration model
 */
namespace TinyBubble
{
    /*! 
     * get the Minnaert resonance frequency for a given 
     * bubble radius
     * v = sqrt(3*k*p0/rho)/r
     */
    inline double minnaert_resonance_freq(double r, 
            double k=HEAT_RATIO, 
            double p0=ATM)
    {
        return sqrt(3.*k*p0/WATER_DENSITY)/r;
    }

    /*!
     * get the resonance frequencey taking into account the effects of thermal
     * conduction and surface tension
     * NOTE: The returned value is actually the angular velocity, w.
     * frequency: f = w / (2*PI)
     */
    inline double resonance_angular_velocity(double r, 
            double k=HEAT_RATIO, 
            double p0=ATM)
    {
        return sqrt((3*k*p0*
                (1.+(SURFACE_TENSION_COEFF_2 / (p0*r))) / 
                WATER_DENSITY) - 
                (SURFACE_TENSION_COEFF_2/(WATER_DENSITY*r))) / r;
    }

    /*!
     * Get the thermal dimensionless damping constant according to Equ.(3.210)
     * in Acoustic Bubble book
     */
    inline double thermal_dimensionless_damping_constant(double r, 
            double k=HEAT_RATIO, double p0=ATM)
    {
        double freq= resonance_angular_velocity(r) * 0.5 * M_1_PI;
        double Gth = 0.75*M_1_PI*HEAT_RATIO*p0 / 
            (WATER_DENSITY*AIR_THERMAL_DIFFUSIVITY);
        double g   = 1.+(1. - 1./(3.*k))*SURFACE_TENSION_COEFF_2/(p0*r);
        double t   = (16.*Gth*g) / (9.*M_SQR(HEAT_RATIO-1.)*freq);
        return 2.*(sqrt(t-3.) - 
                (((3.*HEAT_RATIO)-1.)/(3.*(HEAT_RATIO-1.)))) / (t - 4.);
    }
    
    /*!
     * Get the radiation dimensionless damping constant according to 
     * Equ.(3.218) in the Acoustic Bubble book
     */
    inline double radiation_dimensionless_damping_constant(double r
            /*double k=HEAT_RATIO, double p0=ATM*/)
    {
        return resonance_angular_velocity(r)*r / SOUND_SPEED_WATER;
    }

    /*!
     * get the ratio of thermal conductivity to volumetric heat capacity.
     * This formula is got from:
     *     http://users.wpi.edu/~ierardi/FireTools/air_prop.html
     */
    inline double air_thermal_diffusivity(double T)
    {
        return 9.1018e-11*M_SQR(T)+8.8197e-8*T-1.0654e-5;
    }

    /*!
     * NOTE: the returned "av" is acutally the angular velocity
     * Return the resonance angular velocity and dimensionless damping constant
     */
    inline void get_pulsation_parameters(double& av, double& damping,
            double r, double k=HEAT_RATIO, double p0=ATM)
    {
        av = resonance_angular_velocity(r);
        double f = av * 0.5 * M_1_PI;

        damping = av*r / SOUND_SPEED_WATER;

        double Gth = 0.75*M_1_PI*HEAT_RATIO*p0 / 
                (WATER_DENSITY*AIR_THERMAL_DIFFUSIVITY);
        double g   = 1.+(1. - 1./(3.*k))*SURFACE_TENSION_COEFF_2/(p0*r);
        double t   = (16.*Gth*g) / (9.*M_SQR(HEAT_RATIO-1.)*f);
        damping += (2.*(sqrt(t-3.) - 
                (((3.*HEAT_RATIO)-1.)/(3.*(HEAT_RATIO-1.)))) / (t - 4.));
    }

    /*!
     * Get the parameter for damping vibration given the dimensionless 
     * damping contant (delta) and resonance angular velocity (omega0)
     */
    inline void get_damping_parameters(double& omega_d, double& beta, 
                double delta, double omega0)
    {
        double deltaSqr = M_SQR(delta);
        double omegaSqr = M_SQR(omega0);
        double betaSqr = omegaSqr * deltaSqr / (deltaSqr+4.);
        beta = sqrt(betaSqr);
        omega_d = sqrt(omegaSqr - betaSqr);
    }

    template<typename T>
    inline T reynolds_number_3d(const Vector3<T>& v_bub, 
                                const Vector3<T>& v_liq,
                                T r_bub)
    {
        return (T)2*r_bub*WATER_DENSITY*(v_bub-v_liq).length() /
               WATER_VISCOSITY;
    }

    template<typename T>
    inline void drag_force_3d(const Vector3<T>& v_bub, const Vector3<T>& v_liq, 
                              T r_bub, Vector3<T>& f)
    {
        T rn = reynolds_number_3d(v_bub, v_liq, r_bub);
        f = v_liq - v_bub;
        if ( rn <= (T)1 )
            f *= (T)6 * M_PI * r_bub * WATER_VISCOSITY;
        else
            f *= (T)0.5 * DRAG_COEFFI * WATER_DENSITY * M_PI * M_SQR(r_bub) * f.length();
    }
}

#endif
