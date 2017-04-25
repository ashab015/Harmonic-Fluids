/******************************************************************************
 *  File: macros.h
 *  Define some utility macros
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

/*
 * The returned value should be zero, if an operation is succeed.
 * otherwise, non-zero is returned.
 */
#ifndef MACRO_DEF_H
#   define MACRO_DEF_H

#ifndef ERROR_RETURN
#   define ERROR_RETURN        -1
#endif

#ifndef SUCC_RETURN
#   define SUCC_RETURN          0
#endif

#ifndef _FAILED
#   define _FAILED(p)      ((p) < 0)
#endif

#ifndef _SUCCEEDED
#   define _SUCCEEDED(p)  ((p) == 0)
#endif

#include <stdio.h>

#if defined(DEBUG) | defined(_DEBUG)
#   ifndef V
#       define V(x)                                                     \
        {                                                               \
            int hr = x;                                                 \
            if( _FAILED(hr) ) {                                         \
                fprintf(stderr, "DEBUG Failed: line %d, file \"%s\"\n", \
                        __LINE__, __FILE__);                            \
                fflush(stderr);                                         \
            }                                                           \
        }
#   endif
#   ifndef V_RETURN
#       define V_RETURN(x)                                              \
        {                                                               \
            int hr = x;                                                 \
            if( _FAILED(hr) ) {                                         \
                fprintf(stderr, "DEBUG Failed: line %d, file \"%s\"\n", \
                        __LINE__, __FILE__);                            \
                fflush(stderr);                                         \
                return hr;                                              \
            }                                                           \
        }
#   endif
#else
#   ifndef V
#       define V(x)         { x; }
#   endif
#   ifndef V_RETURN
#       define V_RETURN(x)  { int hr = x; if( _FAILED(hr) ) { return hr; } }
#   endif
#endif

/*
 * Safely delete/release a reference
 */
#ifdef __cplusplus

#ifndef SAFE_DELETE
#   define SAFE_DELETE(p)       { if(p) { delete (p);     (p)=NULL; } }
#endif    
#ifndef SAFE_DELETE_ARRAY
#   define SAFE_DELETE_ARRAY(p) { if(p) { delete[] (p);   (p)=NULL; } }
#endif    
#ifndef SAFE_RELEASE
#   define SAFE_RELEASE(p)      { if(p) { (p)->Release(); (p)=NULL; } }
#endif

#endif

#endif

