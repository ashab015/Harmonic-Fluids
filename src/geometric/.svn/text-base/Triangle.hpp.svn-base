/******************************************************************************
 *  File: Triangle.hpp
 *  simple triangle class in 3D
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
#ifndef GEOMETRIC_TRIANGLE_HPP
#   define GEOMETRIC_TRIANGLE_HPP

#include "Point3.hpp"

#ifdef USE_NAMESPACE
namespace carbine
{
#endif

template <typename T>
class Triangle
{
    public:
        static T area(const Point3<T>& v0, const Point3<T>& v1, const Point3<T>& v2)
        {
            return 0.5 * (v2 - v0).crossProduct(v1 - v0).length();
        }
};

#ifdef USE_NAMESPACE
}
#endif
#endif
