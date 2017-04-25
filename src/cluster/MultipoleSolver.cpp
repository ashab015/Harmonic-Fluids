#include "MultipoleSolver.hpp"
#include <fstream>
#include <unistd.h>

#include "parallel/TaskSchedule.hpp"
#include "Constants.hpp"
#include "linearalgebra/Tuple3.hpp"
#include "linearalgebra/LeastSquareSolver.hpp"
#include "utils/strings.hpp"

const int MultipoleSolver::ORDER_IN     = D_ORDER_IN;
const int MultipoleSolver::ORDER_OUT    = D_ORDER_OUT;
const int MultipoleSolver::M_STRIDE_IN  = (D_ORDER_IN+1)*(D_ORDER_IN+1);
const int MultipoleSolver::M_STRIDE_OUT = (D_ORDER_OUT+1)*(D_ORDER_OUT+1);
const int MultipoleSolver::NUM_FIXED_SRC= 9;

const hf_results* MultipoleSolver::solve(
        const hf_samples* samples, 
        const hf_bubbles* bubs)
{
    using namespace std;

    m_samples = samples;
    m_bubs    = bubs;

    fprintf(logfd, "MSG: solve for %u bubbles ...\n", m_bubs->len);

    //// allocate memory
    m_results.resize(m_bubs->len);
    m_data.resize(m_bubs->len);
    m_d0.resize(m_bubs->len);
    m_coeffMap.clear();
    for(size_t i = 0;i < m_bubs->len;++ i)
    {
        m_results[i].len = m_bubs->val[i].regular_src[1].len * M_STRIDE_OUT;

        m_data[i].resize(m_results[i].len);
        m_d0[i].resize(int(m_bubs->val[i].regular_src[0].len * M_STRIDE_IN));

        m_results[i].id  = m_bubs->val[i].id;
        m_results[i].c   = &m_data[i][0];
        m_results[i].c0  = &m_d0[i][0];
        
        m_coeffMap[m_results[i].id] = &m_results[i];
    }
    m_ret.len = m_bubs->len;
    m_ret.val = &m_results[0];

    //// do least square solve using OpenMP
    const size_t nSolid= m_samples->solid_samples[0].len +
                         m_samples->solid_samples[1].len + 
                         m_samples->solid_samples[2].len + 
                         m_samples->solid_samples[3].len + 
                         m_bottom.size(); 

    const size_t nFree       = m_samples->surf_samples[0].len;
    const size_t nRigidFluid = m_samples->surf_samples[1].len; 
    const size_t nRigidAir   = m_samples->surf_samples[2].len;

    const size_t nWidth = MAX_SRC_OUT * std::max(M_STRIDE_IN, M_STRIDE_OUT); 

    const size_t nRow   = nFree + std::max(nRigidFluid, nRigidAir) + nSolid;
    const size_t maxRow = nWidth + nRow;    // due to the tikhonov regularization

    StaticTaskSchedule<hf_bubble*> tasks;
    for(size_t i = 0;i < m_bubs->len;++ i) tasks.add_task(&m_bubs->val[i]);
    const int NUM_THREADS = tasks.num_threads();
    const int NUM_PROC    = omp_get_max_threads();

    //// allocate memory
    complex<double> *ptrUw[NUM_PROC];
    complex<double> *ptrb [NUM_PROC];
    complex<double> *orgU[NUM_PROC];    // $FITTING-ERROR
    complex<double> *orgb[NUM_PROC];    // $FITTING-ERROR
    for(int i = 0;i < NUM_THREADS;++ i)
    {
        ptrUw[i] = new complex<double>[maxRow * nWidth]; //new complex<_T>[maxRow * nWidth];
        ptrb [i] = new complex<double>[maxRow];          //new complex<_T>[maxRow];
        orgU[i]  = new complex<double>[nRow * nWidth];    // $FITTING-ERROR
        orgb[i]  = new complex<double>[nRow];             // $FITTING-ERROR
    }

    //// solve each of the bubble
    #pragma omp parallel default(none) shared(tasks, ptrUw, ptrb, orgU, orgb, cerr)
    {
        const int idx = omp_get_thread_num();
        //complex<_T> *U  = ptrU [idx];
        complex<double> *Uw = ptrUw[idx];
        complex<double> *b  = ptrb [idx];
        complex<double> *oU = orgU[idx];    // $FITTING-ERROR
        complex<double> *ob = orgb[idx];    // $FITTING-ERROR

        TikhonovQrSolver< complex<double> > solver;
        solver.epsilon() = 1E-7;

        const vector<hf_bubble*>& jobs = tasks.get_task(idx);
        for(size_t i = 0;i < jobs.size();++ i)
        {
            hf_bubble* bub = jobs[i];

            //// $FITTING-ERROR
            //complex<double>* coeff = m_coeffMap[bub->id];
            complex<double>* coeff = m_coeffMap[bub->id]->c;
            double*       residual = m_coeffMap[bub->id]->fittingError;
            //// $FITTING-ERROR

            int WIDTH_IN  = bub->regular_src[0].len * M_STRIDE_IN;
            int WIDTH_OUT = bub->regular_src[1].len * M_STRIDE_OUT;

            const int TOT_ROW_IN = nFree + nSolid + nRigidFluid + WIDTH_IN;
            const int TOT_ROW_OUT= nFree + nSolid + nRigidAir   + WIDTH_OUT;

            //////// first pass: fluid-domain solver
            // compute pressure on free surface (Aa and ba)
            compute_free_boundary_pressure(Uw, b, bub, TOT_ROW_IN, WIDTH_IN);
            // compute pressure normal derivative on solid boundary (As and bs)
            compute_solid_normal_derivative_1pass(&Uw[nFree], &b[nFree], bub, TOT_ROW_IN, WIDTH_IN);

            //// $FITTING-ERROR
            //   copy the matrix into orgU and orgb
            const int ORG_ROW_IN = TOT_ROW_IN - WIDTH_IN;
            for(int k = 0;k < WIDTH_IN;++ k)
                memcpy(&oU[k*ORG_ROW_IN], &Uw[k*TOT_ROW_IN], sizeof(complex<double>)*ORG_ROW_IN);
            memcpy(ob, b, sizeof(complex<double>)*ORG_ROW_IN);
            //// $FITTING-ERROR

            solver.set_size(TOT_ROW_IN - WIDTH_IN, WIDTH_IN);
            solver.solve(Uw, b);
            memcpy(coeff, b, sizeof(complex<double>)*WIDTH_IN);
            memcpy(m_coeffMap[bub->id]->c0, b, sizeof(complex<double>)*WIDTH_IN);

            //// $FITTING-ERROR
            // compute the first fitting error
            residual[0] = compute_residual(CblasColMajor, oU, ob, coeff, ORG_ROW_IN, WIDTH_IN);
            //// $FITTING-ERROR

            //////// second pass: air-domain solver
            // compute normal derivative on free surface
            compute_fluid_normal_derivative(Uw, b, coeff, bub, TOT_ROW_OUT, WIDTH_OUT);
            // compute normal derivative on solid boudnary
            compute_solid_normal_derivative_2pass(&Uw[nFree], bub, TOT_ROW_OUT, WIDTH_OUT);
            memset(&b[nFree], 0, sizeof(complex<double>)*(nSolid + nRigidAir));

            //// $FITTING-ERROR
            //   copy the matrix into orgU and orgb
            const int ORG_ROW_OUT = TOT_ROW_OUT - WIDTH_OUT;
            for(int k = 0;k < WIDTH_OUT;++ k)
                memcpy(&oU[k*ORG_ROW_OUT], &Uw[k*TOT_ROW_OUT], sizeof(complex<double>)*ORG_ROW_OUT);
            memcpy(ob, b, sizeof(complex<double>)*ORG_ROW_OUT);
            //// $FITTING-ERROR

            solver.set_size(TOT_ROW_OUT - WIDTH_OUT, WIDTH_OUT);
            solver.solve(Uw, b);
            memcpy(coeff, b, sizeof(complex<double>)*WIDTH_OUT);

            //// $FITTING-ERROR
            residual[1] = compute_residual(CblasColMajor, oU, ob, coeff, ORG_ROW_OUT, WIDTH_OUT);
            //// $FITTING-ERROR
        }
    }

    for(int i = 0;i < NUM_THREADS;++ i)
    {
        delete [](ptrUw[i]);
        delete [](ptrb [i]);
        delete [](orgU[i]);    // $FITTING-ERROR
        delete [](orgb[i]);    // $FITTING-ERROR
    }
    
    return &m_ret;
}

/*
 * P(f) = 0
 */
void MultipoleSolver::compute_free_boundary_pressure(
        std::complex<double> Uw[], std::complex<double> b[],
        const hf_bubble* bub, const int nrow, const int width)
{
    using namespace std;

    const int DD = M_STRIDE_IN * nrow;
    Point3d scoord;
    double kr;
    const double WAVE_NUM_FLUID = bub->freq / SOUND_SPEED_WATER;
    const hf_surf_sample* samples = &(m_samples->surf_samples[0]);

    for(int i = 0;i < samples->len;++ i)
    {
        //// compute the spherical coordinate, centered at the bubble's center
        cartesian_to_spherical(bub->pos, samples->vtx[i], scoord);
        kr = WAVE_NUM_FLUID * scoord.r;
        const double weight = sqrt(samples->w[i]);

        // first term: monopole:  s = -rho_f*omega^2*r_b^2*exp(-ikr)/r
        b[i] = exp(complex<double>(0, -kr)) * M_SQR(bub->freq * bub->rad) * 
                WATER_DENSITY * weight / scoord.r;

        // iterate on each regular source
        for(int j = 0, base = i;j < (int)bub->regular_src[0].len;++ j, base += DD)
            m_sbh.eval(WAVE_NUM_FLUID, bub->regular_src[0].val[j], samples->vtx[i], &Uw[base], nrow);
        
        cblas_zdscal(width, weight, &Uw[i], nrow);
    }
}

void MultipoleSolver::solid_1pass(
        const hf_solid_sample* sample,
        const Vector3d& normal, const hf_bubble* bub,
        const int nrow, const int width, const double wn,
        std::complex<double> Uw[], std::complex<double> b[])
{
    const int DD = M_STRIDE_IN * nrow;

    for(int i = 0;i < sample->len;++ i)
    {
        const double weight = sqrt(sample->w[i]);

        monopole_derivative(bub, wn, sample->vtx[i], normal, b[i]);
        b[i] *= -weight;

        for(int j = 0, offset = i;j < (int)bub->regular_src[0].len;++ j, offset += DD)
            m_sbh.dir_deriv(wn, 
                    bub->regular_src[0].val[j], 
                    sample->vtx[i], 
                    normal, 
                    &Uw[offset], nrow);
        cblas_zdscal(width, weight, &Uw[i], nrow);
    }
}

/*!
 * \partial P / \partial n = 0
 */
void MultipoleSolver::compute_solid_normal_derivative_1pass(
        std::complex<double> Uw[], std::complex<double> b[],
        const hf_bubble* bub, const int nrow, const int width)
{
    using namespace std;

    const double WAVE_NUM_FLUID = bub->freq / SOUND_SPEED_WATER;    // wave number
    const int DD = M_STRIDE_IN * nrow;

    int ret = 0;
    Vector3d normal;

    //// fluid domain samples on rigid body
    // iterate on each samples on rigid body
    for(int i = 0;i < m_samples->surf_samples[1].len;++ i, ++ ret)
    {
        const double weight = sqrt(m_samples->surf_samples[1].w[i]);
        monopole_derivative(bub, WAVE_NUM_FLUID, m_samples->surf_samples[1].vtx[i],
                        m_samples->surf_samples[1].nml[i], b[ret]);
        b[ret] *= -weight;

        for(int j = 0, offset = ret;j < (int)bub->regular_src[0].len;++ j, offset += DD)
            m_sbh.dir_deriv(WAVE_NUM_FLUID, 
                    bub->regular_src[0].val[j],
                    m_samples->surf_samples[1].vtx[i],
                    m_samples->surf_samples[1].nml[i],
                    &Uw[offset], nrow);
        cblas_zdscal(width, weight, &Uw[ret], nrow);
    }

    //// left/right
    normal.set(-1, 0, 0);
    for(int k = 0;k < 2;++ k, normal.x *= -1)
    {
        solid_1pass(&m_samples->solid_samples[k], normal, bub, nrow, width,
                WAVE_NUM_FLUID, &Uw[ret], &b[ret]);
        ret += m_samples->solid_samples[k].len;
    }
    //// front/back
    normal.set(0, 0, -1);
    for(int k = 2;k < 4;++ k, normal.z *= -1)
    {
        solid_1pass(&m_samples->solid_samples[k], normal, bub, nrow, width,
                WAVE_NUM_FLUID, &Uw[ret], &b[ret]);
        ret += m_samples->solid_samples[k].len;
    }

    //// bottom
    normal.set(0, -1, 0);
    for(size_t i = 0;i < m_bottom.size();++ i, ++ ret)
    {
        monopole_derivative(bub, WAVE_NUM_FLUID, m_bottom[i], normal, b[ret]);
        b[ret] *= -2*CELL_SIZE;

        for(int j = 0, offset = ret;j < (int)bub->regular_src[0].len;++ j, offset += DD)
            m_sbh.dir_deriv(WAVE_NUM_FLUID,
                    bub->regular_src[0].val[j],
                    m_bottom[i],
                    normal, &Uw[offset], nrow);
        cblas_zdscal(width, 2*CELL_SIZE, &Uw[ret], nrow);
    }
}

/*!
 * \partial P / \partial n = U'c
 */
void MultipoleSolver::compute_fluid_normal_derivative(
        std::complex<double> Uw[], std::complex<double> b[],
        const std::complex<double> coeff[], const hf_bubble* bub, 
        const int nrow, const int width)
{
    using namespace std;

    const int DD = M_STRIDE_OUT * nrow;
    const hf_surf_sample* samples = &(m_samples->surf_samples[0]);

    const int WIDTH_1PASS = bub->regular_src[0].len * M_STRIDE_IN;
    const double WAVE_NUM_FLUID = bub->freq / SOUND_SPEED_WATER;
    const double WAVE_NUM_AIR   = bub->freq / SOUND_SPEED_AIR;
    complex<double> tv, basis[WIDTH_1PASS];

    for(int i = 0;i < samples->len;++ i)
    {
        const double weight = sqrt(samples->w[i]);
        //// compute RHS
        // get the derivative on each basis
        for(int j = 0, offset = 0;j < (int)bub->regular_src[0].len;++ j, offset += M_STRIDE_IN)
            m_sbh.dir_deriv(WAVE_NUM_FLUID,
                    bub->regular_src[0].val[j],
                    samples->vtx[i],
                    samples->nml[i],
                    &basis[offset], 1);

        cblas_zdotu_sub(WIDTH_1PASS, coeff, 1, basis, 1, &tv);
        monopole_derivative(bub, WAVE_NUM_FLUID, 
                samples->vtx[i], samples->nml[i], b[i]);
        b[i] += tv;
        b[i] *= weight;

        //// compute i-th row on LHS and update weights
        for(int j = 0, offset = i;j < (int)bub->regular_src[1].len;++ j, offset += DD)
            m_shh.dir_deriv_2nd(WAVE_NUM_AIR,
                    bub->regular_src[1].val[j],
                    samples->vtx[i],
                    samples->nml[i],
                    &Uw[offset], nrow);
        //// pre-multiply a weight
        cblas_zdscal(width, weight, &Uw[i], nrow);
    }
}

void MultipoleSolver::solid_2pass(
        const hf_solid_sample* sample,
        const Vector3d& normal, const hf_bubble* bub,
        const int nrow, const int width, const double wn, 
        std::complex<double> Uw[])
{
    const int DD = M_STRIDE_OUT * nrow;

    for(int i = 0;i < sample->len;++ i)
    {
        const double weight = sqrt(sample->w[i]);

        for(int j = 0, offset = i;j < (int)bub->regular_src[1].len;++ j, offset += DD)
            m_shh.dir_deriv_2nd(wn, 
                    bub->regular_src[1].val[j],
                    sample->vtx[i],
                    normal,
                    &Uw[offset], nrow);
        cblas_zdscal(width, weight, &Uw[i], nrow);
    }
}

void MultipoleSolver::compute_solid_normal_derivative_2pass(
        std::complex<double> Uw[], const hf_bubble* bub, 
        const int nrow, const int width)
{
    using namespace std;

    const double WAVE_NUM_AIR   = bub->freq / SOUND_SPEED_AIR;
    const int DD = M_STRIDE_OUT * nrow;

    int ret = 0;
    Vector3d normal;

    //// air domain samples on rigid body
    for(int i = 0;i < m_samples->surf_samples[2].len;++ i, ++ ret)
    {
        const double weight = sqrt(m_samples->surf_samples[2].w[i]);
        for(int j = 0, offset = ret;j < (int)bub->regular_src[1].len;++ j, offset += DD)
            m_shh.dir_deriv_2nd(WAVE_NUM_AIR,
                    bub->regular_src[1].val[j],
                    m_samples->surf_samples[2].vtx[i],
                    m_samples->surf_samples[2].nml[i],
                    &Uw[offset], nrow);
        cblas_zdscal(width, weight, &Uw[ret], nrow);
    }

    //// left/right
    normal.set(-1, 0, 0);
    for(int k = 0;k < 2;++ k, normal.x *= -1)
    {
        solid_2pass(&m_samples->solid_samples[k], normal, bub, nrow,
                width, WAVE_NUM_AIR, &Uw[ret]);
        ret += m_samples->solid_samples[k].len;
    }

    //// front/back 
    normal.set(0, 0, -1);
    for(int k = 2;k < 4;++ k, normal.z *= -1)
    {
        solid_2pass(&m_samples->solid_samples[k], normal, bub, nrow, 
                width, WAVE_NUM_AIR, &Uw[ret]);
        ret += m_samples->solid_samples[k].len;
    }

    //// bottom 
    normal.set(0, -1, 0);
    for(size_t i = 0;i < m_bottom.size();++ i, ++ ret)
    {
        for(int j = 0, offset = ret;j < (int)bub->regular_src[1].len;++ j, offset += DD)
            m_shh.dir_deriv_2nd(WAVE_NUM_AIR,
                    bub->regular_src[1].val[j],
                    m_bottom[i], normal,
                    &Uw[offset], nrow);
        cblas_zdscal(width, CELL_SIZE, &Uw[ret], nrow); 
    }
}


void MultipoleSolver::monopole_derivative(const hf_bubble* bub, const double wn,
        const Point3d& pt, const Vector3d& normal, std::complex<double>& out) const
{   /* $TESTED$ */
    using namespace std;

    Vector3d dir = pt - bub->pos;
    const double r      = dir.length();
    const double invr   = 1. / r;
    const double kr     = wn * r;
    const double omegaR = bub->freq * bub->rad;

    dir *= invr;       //// normalize direction
    // compute the velocity
    out = exp(complex<double>(0, -kr)) * WATER_DENSITY *
          complex<double>(M_SQR(omegaR * invr), M_SQR(omegaR)*wn*invr);

    // project on normal direction
    out *= dir.dotProduct(normal);
}

void MultipoleSolver::sample_boundary()
{
    m_bottom.clear();
    for(int j = 2;j < (int)LATTICE_NZ-1;++ j)
    for(int i = 2;i < (int)LATTICE_NX-1;++ i)
        m_bottom.push_back(Point3d(i*CELL_SIZE, 0, j*CELL_SIZE));
}

