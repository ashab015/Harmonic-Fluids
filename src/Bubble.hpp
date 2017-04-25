/******************************************************************************
 *  File: Bubble.hpp
 *  Data structure for acoustic bubble
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
#ifndef BUBBLESIM_BUBBLE_HPP
#   define BUBBLESIM_BUBBLE_HPP

#include <complex>
#include "TinyBubble.hpp"
#include "linearalgebra/Vector3.hpp"
#include "geometric/Point3.hpp"

/*!
 * Bubble information
 */
template<typename T = double>
struct Bubble
{
    uint32_t    id;
    T           rad;                    // radius
    T           rho;                    // density
    T           invRho;                 // 1 / density
    T           invMass;                // 1 / mass
    Point3<T>   pos;                    // current position
    Vector3<T>  vel;                    // current velocity

    // bubble pulsation parameters
    T           beta;                   // damping coefficient (exp(-beta*t))
    T           omega0; 
    T           omegad;
    T           waveNumAir;             // wave number \omega / c in air
    T           waveNumFluid;

    T           iso;
    T           sampleDt;
    T           lastSampleTs;
    T           sampleDensityCoeffi;
    T           creationTs;             // bubble creation timestep
    bool        sampled;
    bool        killed;

    std::complex<T>*   coeffi;           // coefficients from least-square solver
    int                numSrc[2];        // number of sources associated with this bubble
    const Point3<T>**  regularSrc[2];
    
    /*!
     * Construct a new bubble
     * 1. set the bubble related parameters
     * 2. compute the initial frequency
     * 3. allocate memory to store radiation estimation data.
     */
    Bubble(uint32_t i, T r, const Point3<T>& pos,
           const Vector3<T>& vel):id(i), rad(r), 
            rho((ATM + SURFACE_TENSION_COEFF_2/r)*AIR_MOLAR_MASS /
                (GAS_CONSTANT * 298.15)),
            invRho((T)1 / rho), 
            invMass(3. / (4. * M_PI * M_TRI(r) * rho)),
            pos(pos), vel(vel), 
            sampleDt(0), lastSampleTs(0),
            killed(false), 
            coeffi(NULL)
    {
        double delta;
        // compute omega0 and delta
        TinyBubble::get_pulsation_parameters(omega0, delta, (double)rad);
        // compute omegad and beta
        TinyBubble::get_damping_parameters(omegad, beta, delta, (double)omega0);

        waveNumAir   = omegad / SOUND_SPEED_AIR;
        waveNumFluid = omegad / SOUND_SPEED_WATER;
        sampleDensityCoeffi = 0.55 * 2 * BUB_AVG_TIMESTEP;

        numSrc[0] = int((20 - (1000. - omegad * 0.5 * M_1_PI) / 75) + 0.5);
        numSrc[0] = int(std::max(20, std::min(40, numSrc[0]))*1.5);
        numSrc[1] = std::min(122, numSrc[0]*2 + 40);

        coeffi = new std::complex<T>[std::max(numSrc[0], numSrc[1]) * 9];

        regularSrc[0] = new const Point3<T>*[numSrc[0]];
        regularSrc[1] = new const Point3<T>*[numSrc[1]];
    }

    ~Bubble()
    {
        if ( coeffi ) 
        {
            delete []coeffi;
            delete []regularSrc[0];
            delete []regularSrc[1];
        }
    }

    // decide if we need to sample this bubble, if yes, 
    // update the corresponding state variables
    bool should_sample(T ts)
    {
        if ( killed ) return (sampled = false);

        T gd = exp(-beta*(ts-creationTs));
        if ( gd < 0.0000001 )               // vibration amplitude is small enough
        {
            killed = true;                  // indicate it would not sample the bubble any more
            return (sampled = false);
        }
        return (sampled = true);

        if ( lastSampleTs + sampleDt > ts )
            return (sampled = false);

        lastSampleTs = ts;
        sampleDt = fmin(sampleDensityCoeffi / gd, 0.0008);
        return (sampled = true);
    }
};

#endif
