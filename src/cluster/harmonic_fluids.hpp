/******************************************************************************
 *  File: harmonic_fluids.hpp
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
#ifndef _HARMONIC_FLUIDS_H_RPCGEN
#define _HARMONIC_FLUIDS_H_RPCGEN

#include <rpc/rpc.h>
#include <vector>
#include <complex>
#include "Bubble.hpp"
#include "geometric/TriangleMesh.hpp"

struct cluster_config_arg;

struct hf_solid_sample
{
    int         len;
    Point3d*    vtx;
    double*     w;
};

struct hf_surf_sample
{
    int         len;
    Point3d*    vtx;
    Vector3d*   nml;
    double*     w;
};

struct hf_samples
{
    double          ts;
    hf_surf_sample  surf_samples[3];
    hf_solid_sample solid_samples[4];
};

//=========================================================
struct hf_bubble 
{
    int         id;
    double      freq;
    double      rad;
    Point3d     pos;

    struct 
    {
        u_int    len;
        Point3d* val;
    } regular_src[2];
};

struct hf_bubbles
{
    double      ts;
    u_int       len;
    hf_bubble*  val;
};

//=========================================================
struct hf_result
{
    int     id;
    int     len;
    std::complex<double>*   c0;
    std::complex<double>*   c;
    double  fittingError[2];    // $FITTING-ERROR
};

struct hf_results
{
    int         len;
    hf_result*  val;

    hf_results():len(0), val(NULL) {}
    hf_results(int l):len(l), val(NULL) {}
};

#define HARMONIC_FLUIDS_PROGRAM 0x20000586
#define HARMONIC_FLUIDS_VERSION_1 1

#define HF_CONFIG 1
#define HF_SOLVE  2

int harmonic_fluids_program_1_freeresult(SVCXPRT *, xdrproc_t, caddr_t);

/* the xdr functions */
bool_t xdr_solid_sample(XDR *xdrs, hf_solid_sample* objp);
bool_t xdr_surf_sample(XDR* xdrs, hf_surf_sample* objp);
bool_t xdr_samples(XDR *xdrs, hf_samples* objp);

bool_t xdr_tuple3d(XDR *xdrs, Tuple3d *objp);
bool_t xdr_bubble(XDR *xdrs, hf_bubble* bub);
bool_t xdr_bubbles(XDR *xdrs, hf_bubbles* objp);

bool_t xdr_complex(XDR *xdrs, std::complex<double>* c);
bool_t xdr_result(XDR *xdrs, hf_result* ret);
bool_t xdr_results(XDR *xdrs, const hf_results* objp);

bool_t xdr_solid_sample_enc(XDR *xdrs, const TriangleMesh<double>*);
bool_t xdr_samples_enc(XDR *xdrs, cluster_config_arg* objp);
bool_t xdr_samples(XDR *xdrs, hf_samples* objp);
bool_t xdr_solid_sample(XDR *xdrs, hf_solid_sample* objp);

#endif /* !_HARMONIC_FLUIDS_H_RPCGEN */
