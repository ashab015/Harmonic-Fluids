/******************************************************************************
 *  File: trivial_type.hpp
 *  Extract the trivial type from a more complicated type
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

#ifndef CARBINE_TYPE_TRAIT_HPP
#   define CARBINE_TYPE_TRAIT_HPP

#include <complex>
#include "null_type.hpp"

namespace carbine
{
    template<typename T>
    struct TrivialType
    { typedef NullType    type; };

    template<>
    struct TrivialType<double>
    { typedef double      type; };

    template<>
    struct TrivialType<float>
    { typedef float       type; };

    template<>
    struct TrivialType< std::complex<double> >
    { typedef double      type; };

    template<>
    struct TrivialType< std::complex<float> >
    { typedef float       type; };
}

#endif
