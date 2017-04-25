/******************************************************************************
 *  File: PrecisionDef.hpp
 *  Switch computation precision
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
#ifndef PRECISION_DEF
#   define PRECISION_DEF

#ifdef USE_FLOAT

#   define _T  float
const _T EPS  = 1e-8;
const _T EPS2 = 1e-16;

#else

#   define _T  double
const _T EPS  = 1e-12;
const _T EPS2 = 1e-24;

#endif

#endif
