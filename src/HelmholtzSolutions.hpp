/******************************************************************************
 *  File: HelmholtzSolution.hpp
 *  Fundamental solution of Helmholtz equation in free space
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
#ifndef HELMHOLTZ_SOLUTION_HPP
#   define HELMHOLTZ_SOLUTION_HPP

#include "SpecialFunc.hpp"
#include "linearalgebra/SphericalFunc.hpp"
#ifdef HAVE_BOOST
#   include <boost/static_assert.hpp>
#endif

/*
 * Regular solution
 * compute j_n * sph(m, n), the basic solution of helmholtz equation
 */
template <typename T, int Order>
class SphericalBesselHarmonics
{
    public:
        /*
         * \param wn    the wave number
         * \param src   the source position
         * \param pos   the receiver position
         * \param output store the output data
         * \param incl  increment of the output pointer
         */
        void eval(T wn, const Point3<T>& src, const Point3<T>& pos, 
                  std::complex<T>* output, int incl)
        {
            Point3<T> scoord;

            cartesian_to_spherical(src, pos, scoord);
            const T kr = wn * scoord.r;

            // compute j0(kr) * sph(0, 0)
            output[0] = m_spb.j0(kr) * m_sph.Y0();
            int idx = incl;

            if ( Order > 0 )
            {
                const T j1 = m_spb.j1(kr);
                for(int i = -1;i <= 1;++ i, idx += incl)
                    output[idx] = j1 * m_sph.Y1(i, scoord.theta, scoord.phi);
            }

            if ( Order > 1 )
            {
                const T j2 = m_spb.j2(kr);
                for(int i = -2;i <= 2;++ i, idx += incl)
                    output[idx] = j2 * m_sph.Y2(i, scoord.theta, scoord.phi);
            }
        }

        /* 
         * compute directional derivative 
         */
        void dir_deriv(T wn, const Point3<T>& src, const Point3<T>& pos,
                const Vector3<T>& dir, std::complex<T>* output, int incl)
        {
            Point3<T> scoord;
            const Vector3<T> rv = pos - src;

            const T xysqr = cartesian_to_spherical(rv, scoord);
            const T xyd   = sqrt(xysqr);
            const T invxysqr = 1. / xysqr;
            const T invxyd   = 1. / xyd;
            const T kr = wn * scoord.r;

            const T invr = 1. / scoord.r;
            const T invr2 = M_SQR(invr);

            const T prpx = rv.x * invr, prpy = rv.y * invr, prpz = rv.z * invr;
            const T pthetapx = rv.x*rv.z*invr2*invxyd;
            const T pthetapy = rv.y*rv.z*invr2*invxyd;
            const T pthetapz = -xyd*invr2;
            const T pphipx   = -rv.y * invxysqr;
            const T pphipy   =  rv.x * invxysqr;
            const T pphipz   = 0;

            const T sfr     = prpx*dir.x + prpy*dir.y + prpz*dir.z;
            const T sftheta = pthetapx*dir.x + pthetapy*dir.y + pthetapz*dir.z; 
            const T sfphi   = pphipx*dir.x + pphipy*dir.y + pphipz*dir.z;

            // order 0:
            output[0] = sfr*wn*m_spb.j0_deriv(kr)*m_sph.Y0();
            int idx = incl;

            if ( Order > 0 )
            {
                const T j1      = m_spb.j1(kr);
                const T j1deriv = wn * m_spb.j1_deriv(kr);

                for(int i = -1;i <= 1;++ i, idx += incl)
                    output[idx] =  sfr * j1deriv * m_sph.Y1(i, scoord.theta, scoord.phi) +
                                   sftheta * j1 * m_sph.Y1_deriv_theta(i, scoord.theta, scoord.phi) +
                                   sfphi * j1 * m_sph.Y1_deriv_phi(i, scoord.theta, scoord.phi);
            }

            if ( Order > 1 )
            {
                const T j2      = m_spb.j2(kr);
                const T j2deriv = wn * m_spb.j2_deriv(kr);
                for(int i = -2;i <= 2;++ i, idx += incl)
                    output[idx] = sfr * j2deriv * m_sph.Y2(i, scoord.theta, scoord.phi) + 
                                  sftheta * j2 * m_sph.Y2_deriv_theta(i, scoord.theta, scoord.phi) + 
                                  sfphi * j2 * m_sph.Y2_deriv_phi(i, scoord.theta, scoord.phi);
            }
        }

    private:
#ifdef HAVE_BOOST
        //// validate ORDER at compile time
        BOOST_STATIC_ASSERT(Order <= 2 && Order >= 0);
#endif

        SphericalBessel<T>     m_spb;
        SphericalHarmonics<T>  m_sph;
};

/*
 * Singular solution of Helmholtz equation
 */
template <typename T, int Order>
class SphericalHankelHarmonics
{
    public:
        /*
         * \param wn    the wave number
         * \param src   the source position
         * \param pos   the receiver position
         * \param output store the output data
         * \param incl  increment of the output pointer
         */
        void eval_2nd(T wn, const Point3<T>& src, const Point3<T>& pos, 
                     std::complex<T>* output, int incl)
        {
            Point3<T> scoord;
            cartesian_to_spherical(src, pos, scoord);
            const T kr = wn * scoord.r;

            // compute h0(kr)*sph(0,0)
            output[0] = m_spk.h2nd0(kr) * m_sph.Y0();
            int idx = incl;

            if ( Order > 0 )
            {
                const std::complex<T> h1 = m_spk.h2nd1(kr);
                for(int i = -1;i <=1;++ i, idx += incl)
                    output[idx] = h1 * m_sph.Y1(i, scoord.theta, scoord.phi);
            }

            if ( Order > 1 )
            {
                const std::complex<T> h2 = m_spk.h2nd2(kr);
                for(int i = -2;i <= 2;++ i, idx += incl)
                    output[idx] = h2 * m_sph.Y2(i, scoord.theta, scoord.phi);
            }
        }

        void dir_deriv_2nd(T wn, const Point3<T>& src, const Point3<T>& pos,
                const Vector3<T>& dir, std::complex<T>* output, int incl)
        {
            Point3<T> scoord;
            const Vector3<T> rv = pos - src;

            const T xysqr = cartesian_to_spherical(rv, scoord);
            const T xyd   = sqrt(xysqr);
            const T invxysqr = 1. / xysqr;
            const T invxyd   = 1. / xyd;
            const T kr = wn * scoord.r;

            const T invr = 1. / scoord.r;
            const T invr2 = M_SQR(invr);

            const T prpx = rv.x * invr, prpy = rv.y * invr, prpz = rv.z * invr;
            const T pthetapx = rv.x*rv.z*invr2*invxyd;
            const T pthetapy = rv.y*rv.z*invr2*invxyd;
            const T pthetapz = -xyd*invr2;
            const T pphipx   = -rv.y * invxysqr;
            const T pphipy   =  rv.x * invxysqr;
            const T pphipz   = 0;

            const T sfr     = prpx*dir.x + prpy*dir.y + prpz*dir.z;
            const T sftheta = pthetapx*dir.x + pthetapy*dir.y + pthetapz*dir.z; 
            const T sfphi   = pphipx*dir.x + pphipy*dir.y + pphipz*dir.z;

            // order 0:
            output[0] = sfr*wn*m_spk.h2nd0_deriv(kr) * m_sph.Y0();
            int idx = incl;

            if ( Order > 0 )
            {
                const std::complex<T> h1      = m_spk.h2nd1(kr);
                const std::complex<T> h1deriv = wn * m_spk.h2nd1_deriv(kr);

                for(int i = -1;i <= 1;++ i, idx += incl)
                    output[idx] = sfr * h1deriv * m_sph.Y1(i, scoord.theta, scoord.phi) +
                                  sftheta * h1 * m_sph.Y1_deriv_theta(i, scoord.theta, scoord.phi) +
                                  sfphi * h1 * m_sph.Y1_deriv_phi(i, scoord.theta, scoord.phi);
            }

            if ( Order > 1 )
            {
                const std::complex<T> h2      = m_spk.h2nd2(kr);
                const std::complex<T> h2deriv = wn * m_spk.h2nd2_deriv(kr);

                for(int i = -2;i <= 2;++ i, idx += incl)
                    output[idx] = sfr * h2deriv * m_sph.Y2(i, scoord.theta, scoord.phi) + 
                                  sftheta * h2 * m_sph.Y2_deriv_theta(i, scoord.theta, scoord.phi) +
                                  sfphi * h2 * m_sph.Y2_deriv_phi(i, scoord.theta, scoord.phi);
            }
        }

    private:
#ifdef HAVE_BOOST
        //// validate ORDER at compile time
        BOOST_STATIC_ASSERT(Order <= 2 && Order >= 0);
#endif

        SphericalHankel<T>      m_spk;
        SphericalHarmonics<T>   m_sph;
};
#endif
