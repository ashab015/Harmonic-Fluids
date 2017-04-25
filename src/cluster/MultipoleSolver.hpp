/******************************************************************************
 *  File: MultipoleSolver.hpp
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
#ifndef CLUSTER_MULTIPOLE_SOLVER_HPP
#   define CLUSTER_MULTIPOLE_SOLVER_HPP

#include <map>

#include "HelmholtzSolutions.hpp"
#include "harmonic_fluids.hpp"

/*!
 * service of multipole solver running on clusters
 */
class MultipoleSolver
{
    public:
        MultipoleSolver():logfd(stderr)
        {
            sample_boundary();
        }

        /*
         * parallel solve using OpenMP
         */
        const hf_results* solve(
                const hf_samples* samples, 
                const hf_bubbles* bubs);

    private:
        void sample_boundary();

        void compute_free_boundary_pressure(
                std::complex<double> Uw[], std::complex<double> b[],
                const hf_bubble* bub, const int nrow, const int width);
        void solid_1pass(const hf_solid_sample* sample,
                const Vector3d& normal, const hf_bubble* bub,
                const int nrow, const int width, const double wn,
                std::complex<double> Uw[], std::complex<double> b[]);
        void compute_solid_normal_derivative_1pass(
                std::complex<double> Uw[], std::complex<double> b[],
                const hf_bubble* bub, const int nrow, const int width);
        void compute_fluid_normal_derivative(
                std::complex<double> Uw[], std::complex<double> b[],
                const std::complex<double> coeff[], const hf_bubble* bub,
                const int nrow, const int width);
        void solid_2pass(const hf_solid_sample* sample,
                const Vector3d& normal, const hf_bubble* bub,
                const int nrow, const int width, const double wn,
                std::complex<double> Uw[]);
        void compute_solid_normal_derivative_2pass(
                std::complex<double> Uw[], 
                const hf_bubble* bub, 
                const int nrow, 
                const int width);

    public:
        void monopole_derivative(const hf_bubble* bub, const double wn, 
                const Point3d& pt, const Vector3d& normal, 
                std::complex<double>& out) const;

    private:
        const hf_bubbles*           m_bubs;
        const hf_samples*           m_samples;

        hf_results                  m_ret;
        std::vector<hf_result>      m_results;      // result for each bubble
        std::vector< std::vector< std::complex<double> > >  m_data;
        std::vector< std::vector< std::complex<double> > >  m_d0;   // $SPLASH-6

        //// $FITTING-ERROR
        //std::map< int, std::complex<double>* >              m_coeffMap; // map from bubble id to the result coefficients
        std::map< int, hf_result* >              m_coeffMap; // map from bubble id to the result coefficients

    public:
        std::vector<Point3d>        m_bottom;

        SphericalBesselHarmonics<double, D_ORDER_IN>    m_sbh;
        SphericalHankelHarmonics<double, D_ORDER_OUT>   m_shh;
        FILE* logfd;

    public:
        static const int ORDER_IN;
        static const int ORDER_OUT;
        static const int M_STRIDE_IN;
        static const int M_STRIDE_OUT;
        static const int NUM_FIXED_SRC;

        static const struct timeval CONFIG_TIMEOUT;
};
#endif
