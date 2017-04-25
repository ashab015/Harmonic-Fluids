/******************************************************************************
 *  File: math.h
 *  Utilities for some math functions
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
#ifndef MATH_UTILS_H
#   define MATH_UTILS_H

#ifndef M_PI
#   define M_PI           3.14159265358979323846  /* pi */
#endif

#define M_ABS(a)      ((a)>0?(a):-(a))
#define M_NEG(a)      (-(a))
#define M_SQR(a)      ((a)*(a))
#define M_TRI(a)      ((a)*(a)*(a))
#define M_MAX(a, b)   ((a)>(b)?(a):(b))
#define M_MIN(a, b)   ((a)<(b)?(a):(b))
#define M_DEG2RAD(x)  (((x)*M_PI) / 180.0)
#define M_RAD2DEG(x)  (((x)*180.) / M_PI)

inline int gcd(int a, int b)
{
    for(int c;b;c=a,a=b,b=c%b) ;
    return a;
}

inline int lcm(int a, int b)
{
    return a/gcd(a,b)*b;
}

template<typename T>
inline T clamp(T a, T minv, T maxv)
{
    return a <= minv ? minv : (a > maxv ? maxv : a);
}

#endif
