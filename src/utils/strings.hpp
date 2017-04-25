/******************************************************************************
 *  File: string.hpp
 *  Manipulate stl::string
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
#ifndef CARBINE_STRING_UTILS_H
#   define CARBINE_STRING_UTILS_H

#include <vector>
#include <string>
#include <sstream>

//! string to int function
#define S2D(s,d)  {std::istringstream(s)>>d;}

namespace carbine
{

//! Trim the string
std::string trim(const std::string text);
//! Return whether the string a starts with string b
bool start_with(const char* pre, const char* str);

bool end_with(const char* suf, const char* str);

template<typename T> 
inline int Int(const T &t)
{
    int r;
    std::stringstream s;
    s << t;
    s >> r;
    return r;
}

template<typename T> 
inline double Double(const T &t)
{
    double r;
    std::stringstream s;
    s << t;
    s >> r;
    return r;
}

template<typename T>
inline float Float(const T& t)
{
    float r;
    std::stringstream s;
    s << t;
    s >> r;
    return r;
}

}

#endif
