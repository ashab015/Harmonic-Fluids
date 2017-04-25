/******************************************************************************
 *  File: IdGenerator.hpp
 *  Generate the sequential IDs
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
#ifndef UTILS_ID_GENERATOR_H
#define     UTILS_ID_GENERATOR_H

#include <assert.h>

namespace carbine
{

class CircularIdGenerator
{
    public:
        CircularIdGenerator(int maxId):m_maxId(maxId), m_nextId(0)
        {
            assert(maxId > 0);
        }

        inline int next_id()
        {
            m_nextId = (m_nextId + 1) % m_maxId;
            return m_nextId;
        }

    private:
        const int   m_maxId;
        int         m_nextId;
};

}

#endif
