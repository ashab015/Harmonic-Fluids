/******************************************************************************
 *  File: SphericalFunc.hpp
 *  Computing in spherical coordinates
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
#ifndef SPHERICAL_HANKEL_HPP
#   define SPHERICAL_HANKEL_HPP

#include <complex>
#include "Constants.hpp"
#include "geometric/Point3.hpp"
#include "Vector3.hpp"

/*!
 * out.x ==> r
 * out.y ==> theta
 * out.z ==> phi
 */
template<typename T>                /* $TESTED$ */
inline void cartesian_to_spherical(const Point3<T>& from,
        const Point3<T>& to, Point3<T>& out)
{
    Vector3<T> v = to - from;
    out.r = v.length();

    if ( out.r < 0.8*CELL_SIZE )
    {
        fprintf(stderr, "ERROR: <cartesian_to_spherical> r is less than CELL_SIZE! %lf [%lf,%lf,%lf]==>[%lf,%lf,%lf]\n", 
                out.r, from.x, from.y, from.z, to.x, to.y, to.z);
        if ( out.r < EPS )
        {
            out.zero();
            fprintf(stderr, "ERROR: <cartesian_to_spherical> r is too small! %lf\n", out.r);
            return;
        }
    }

    out.theta = acos(v.z / out.r); 
    out.phi   = atan2(v.y, v.x);
}

template<typename T>                /* $TESTED$ */
inline T cartesian_to_spherical(const Vector3<T>& v, Point3<T>& out)
{
    out.r = v.length();

    if ( out.r < 0.8*CELL_SIZE )
    {
        fprintf(stderr, "ERROR: <cartesian_to_spherical> r is less than CELL_SIZE! [%lf,%lf,%lf] LEN=%lf\n", 
                v.x, v.y, v.z, out.r);
        if ( out.r < EPS )
        {
            out.zero();
            fprintf(stderr, "ERROR: <cartesian_to_spherical> r is too small! %lf\n", out.r);
            return 0;
        }
    }

    T ret;
    if ( (ret = v.x*v.x + v.y*v.y) < EPS )
    {
        fprintf(stderr, "ERROR: <cartesian_to_spherical> v.x^2+v.y^2 is too small! %lf\n", ret);
    }

    out.theta = acos(v.z / out.r); 
    out.phi   = atan2(v.y, v.x);

    return ret;
}

#endif
