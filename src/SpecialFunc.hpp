/******************************************************************************
 *  File: SpecialFunc.hpp
 *  Computing the spherical Bessel/Hankel function and spherical harmonics with
 *  the order less than 3.
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
#ifndef SPECIAL_FUNC_HPP
#   define SPECIAL_FUNC_HPP

#include <math.h>
#include <complex>
#include <stdlib.h>
#include <stdio.h>
#include "utils/math.hpp"

/*
 * Computing spherical Bessel function with the order less than 3
 */
template <typename T>
class SphericalBessel
{
    public:
        T j(const int order, T x) const
        {
            switch (order)
            {
                case 0:
                    return j0(x);
                case 1:
                    return j1(x);
                case 2:
                    return j2(x);
                default:
                    fprintf(stderr, "SphericalBessel-j: Order %d hasn't been supported!\n",
                            order);
                    return 0;
            }
        }

        T j0(T x) const
        {
            return sin(x) / x;
        }

        T j1(T x) const
        {
            T invx = 1. / x;
            return sin(x)*invx*invx - cos(x)*invx;
        }

        T j2(T x) const
        {
            T invx = 1. / x;
            return (3*invx*invx - 1)*sin(x)*invx - 3*cos(x)*invx*invx;
        }

        /* analytic derivative of j0 */
        T j0_deriv(T x) const
        {
            T invx = 1. / x;
            return cos(x)*invx - sin(x)*invx*invx;
        }

        /* analytic derivative of j1 */
        T j1_deriv(T x) const
        {
            T invx = 1. / x;
            T invx2 = invx*invx;
            T invx3 = invx*invx2;
            T sinx = sin(x);
            return 2*cos(x)*invx2 - 2*sinx*invx3 + sinx*invx;
        }

        /* analytic derivative of j2 */
        T j2_deriv(T x) const
        {
            T t1 = x*x;      
            T t2 = t1*t1;      
            T t4 = sin(x);      
            T t7 = 1/t1;      
            T t9 = 3.0*t7-1.0;      
            T t10 = cos(x);      
            return -6.0/t2*t4+t9*t10/x-t9*t4*t7+3.0*t4*t7+6.0*t10/t1/x;
        }

        // ===================================================== //
        T y(const int order, T x) const
        {
            switch (order)
            {
                case 0:
                    return y0(x);
                case 1:
                    return y1(x);
                case 2:
                    return y2(x);
                default:
                    fprintf(stderr, "SphericalBessel-y: Order %d hasn't been supported!\n",
                            order);
                    return 0;
            }
        }

        T y0(T x) const
        {
            return -cos(x) / x;
        }

        T y1(T x) const
        {
            T invx = 1. / x;
            return -(cos(x)*invx*invx + sin(x)*invx);
        }

        T y2(T x) const
        {
            T invx = 1. / x;
            return (-3*invx*invx + 1)*cos(x)*invx - 3*sin(x)*invx*invx;
        }

        /* analytic derivative of y0 */
        T y0_deriv(T x) const
        {
            T invx = 1. / x;
            return sin(x)*invx + cos(x)*invx*invx;
        }

        T y1_deriv(T x) const
        {
            T invx = 1. / x;
            T invx2 = invx*invx;
            T invx3 = invx2*invx;
            T cosx = cos(x);
            return -cosx*invx + 2*sin(x)*invx2 + 2*cosx*invx3;
        }

        T y2_deriv(T x) const
        {
            T t1 = x*x;      
            T t2 = t1*t1;      
            T t4 = cos(x);      
            T t7 = 1/t1;      
            T t9 = -3.0*t7+1.0;      
            T t10 = sin(x);      
            return 6.0/t2*t4-t9*t10/x-t9*t4*t7-3.0*t4*t7+6.0*t10/t1/x;
        }
};

template <typename T>
class SphericalHankel
{
    public:
        /* zero order hankel func of first kind */
        std::complex<T> h1st0(T x) const
        {
            return std::complex<T>(m_spb.j0(x), m_spb.y0(x));
        }

        std::complex<T> h1st1(T x) const
        {
            return std::complex<T>(m_spb.j1(x), m_spb.y1(x));
        }

        std::complex<T> h1st2(T x) const
        {
            return std::complex<T>(m_spb.j2(x), m_spb.y2(x));
        }

        /* derivative of zero order hankel func of first kind */
        std::complex<T> h1st0_deriv(T x) const
        {
            return std::complex<T>(m_spb.j0_deriv(x), m_spb.y0_deriv(x));
        }

        /* derivative of first order hankel func of first kind */
        std::complex<T> h1st1_deriv(T x) const
        {
            return std::complex<T>(m_spb.j1_deriv(x), m_spb.y1_deriv(x));
        }

        std::complex<T> h1st2_deriv(T x) const
        {
            return std::complex<T>(m_spb.j2_deriv(x), m_spb.y2_deriv(x));
        }

        // =========================================================== //
        /* zero order hankel func of 2nd kind */
        std::complex<T> h2nd0(T x) const
        {
            return std::complex<T>(m_spb.j0(x), -m_spb.y0(x));
        }

        std::complex<T> h2nd1(T x) const
        {
            return std::complex<T>(m_spb.j1(x), -m_spb.y1(x));
        }

        std::complex<T> h2nd2(T x) const
        {
            return std::complex<T>(m_spb.j2(x), -m_spb.y2(x));
        }

        /* derivative of zero order hankel func of 2nd kind */
        std::complex<T> h2nd0_deriv(T x) const
        {
            return std::complex<T>(m_spb.j0_deriv(x), -m_spb.y0_deriv(x));
        }

        std::complex<T> h2nd1_deriv(T x) const
        {
            return std::complex<T>(m_spb.j1_deriv(x), -m_spb.y1_deriv(x));
        }

        std::complex<T> h2nd2_deriv(T x) const
        {
            return std::complex<T>(m_spb.j2_deriv(x), -m_spb.y2_deriv(x));
        }

    private:
        SphericalBessel<T>  m_spb;
};

template <typename T>
class SphericalHarmonics
{
    public:
        std::complex<T> Y0() const
        {  return 0.5*sqrt(M_1_PI); }

        std::complex<T> Y0_deriv() const
        {  return 0; }

        std::complex<T> Y1(const int m, T theta, T phi) const
        {
            double a;
            switch (m)
            {
                case -1:
                    a = 0.5*sqrt(1.5*M_1_PI)*sin(theta);
                    return std::complex<T>(a*cos(-phi), a*sin(-phi));
                case 0:
                    return std::complex<T>(0.5*sqrt(3*M_1_PI)*cos(theta), 0);
                case 1:
                    a = -0.5*sqrt(1.5*M_1_PI)*sin(theta);
                    return std::complex<T>(a*cos(phi), a*sin(phi));
                default:
                    fprintf(stderr, "ERROR: %d for Y1\n", m);
                    exit(1);
            }
        }

        std::complex<T> Y1_deriv_theta(const int m, T theta, T phi) const
        {
            double a;
            switch (m)
            {
                case -1:
                    a = 0.5*sqrt(1.5*M_1_PI)*cos(theta);
                    return std::complex<T>(a*cos(-phi), a*sin(-phi));
                case 0:
                    return std::complex<T>(-0.5*sqrt(3*M_1_PI)*sin(theta), 0);
                case 1:
                    a = -0.5*sqrt(1.5*M_1_PI)*cos(theta);
                    return std::complex<T>(a*cos(phi), a*sin(phi));
                default:
                    fprintf(stderr, "ERROR: %d for Y1_deriv_theta\n", m);
                    exit(1);
            }
        }

        std::complex<T> Y1_deriv_phi(const int m, T theta, T phi) const
        {
            double a;
            switch (m)
            {
                case -1:
                    a = 0.5*sqrt(1.5*M_1_PI)*sin(theta);
                    break;
                case 0:
                    return std::complex<T>(0, 0);
                case 1:
                    a = -0.5*sqrt(1.5*M_1_PI)*sin(theta);
                    break;
                default:
                    fprintf(stderr, "ERROR: %d for Y1_deriv_phi\n", m);
                    exit(1);
            }
            return a*std::complex<T>(-m*sin(m*phi), m*cos(m*phi));
        }

        std::complex<T> Y2(const int m, T theta, T phi) const
        {
            double a;
            double costheta;
            double sintheta;

            switch (M_ABS(m))
            {
                case 0:
                    costheta = cos(theta);
                    return std::complex<T>(
                            0.25*sqrt(5*M_1_PI)*(3*M_SQR(costheta) - 1),
                            0);
                case 1:
                    a = 0.5*sqrt(7.5*M_1_PI)*sin(theta)*cos(theta);
                    return -m*a*std::complex<T>(cos(m*phi), sin(m*phi));
                case 2:
                    sintheta = sin(theta);
                    a = 0.25*sqrt(7.5*M_1_PI)*M_SQR(sintheta);
                    return a*std::complex<T>(cos(m*phi), sin(m*phi));
                default:
                    fprintf(stderr, "ERROR: %d for Y2\n", m);
                    exit(1);
            }
        }

        std::complex<T> Y2_deriv_theta(const int m, T theta, T phi) const
        {
            double costheta;
            double sintheta;
            double a;

            switch (M_ABS(m))
            {
                case 0:
                    return std::complex<T>(
                            0.25*sqrt(5*M_1_PI)*(-6)*cos(theta)*sin(theta),
                            0);
                case 1:
                    sintheta = sin(theta);
                    costheta = cos(theta);
                    a = 0.5*sqrt(7.5*M_1_PI)*(M_SQR(costheta)-M_SQR(sintheta));
                    return -m*a*std::complex<T>(cos(m*phi), sin(m*phi));
                case 2:
                    a = 0.5*sqrt(7.5*M_1_PI)*sin(theta)*cos(theta);
                    return a*std::complex<T>(cos(m*phi), sin(m*phi));
                default:
                    fprintf(stderr, "ERROR: %d for Y2_deriv_theta\n", m);
                    exit(1);
            }
        }

        std::complex<T> Y2_deriv_phi(const int m, T theta, T phi) const
        {
            double a;
            double sintheta;

            switch (M_ABS(m))
            {
                case 0:
                    return std::complex<T>(0, 0);
                case 1:
                    a = -m*0.5*sqrt(7.5*M_1_PI)*sin(theta)*cos(theta);
                    break;
                case 2:
                    sintheta = sin(theta);
                    a = 0.25*sqrt(7.5*M_1_PI)*M_SQR(sintheta);
                    break;
                default:
                    fprintf(stderr, "ERROR: %d for Y2\n", m);
                    exit(1);
            }
            return a*std::complex<T>(-m*sin(m*phi), m*cos(m*phi));
        }
};

#endif
