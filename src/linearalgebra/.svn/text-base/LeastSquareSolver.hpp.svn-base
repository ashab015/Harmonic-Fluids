/******************************************************************************
 *  File: LeastSquareSolver.hpp
 *  Implement different least-square solve methods
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
#ifndef LEAST_SQUARE_SOLVER
#   define LEAST_SQUARE_SOLVER

#include <vector>

#ifdef USE_MKL
#include <mkl.h>
#include <mkl_lapack.h>
#elif defined(__APPLE__) || defined(MACOSX)
#   undef USE_MACOS_VECLIB
#   define USE_MACOS_VECLIB
#   include <vecLib/vecLib.h>
#else
#   error ERROR: Cannot locate the LAPACK lib
#endif

#include "generic/trivial_type.hpp"

/*!
 * Solve least-square using Tikhonov regularization and QR.
 * [  A ][x] = [b]
 * [ aI ]      [0]
 * Because of the regularization, we assure the final matrix is full rank
 *
 * NOTE: the size of A should be at least (m+n)*n, and b should be at
 *       least (m + n)
 */
template <typename T>
class TikhonovQrSolver
{
    public:
        typedef T        value_type;
        typedef typename carbine::TrivialType<T>::type trivial_type;

        TikhonovQrSolver():m_m(0), m_n(0), m_lwork(0),
                m_eps((trivial_type)1E-8)
        {
            m_workspace.resize(2);
        }
        void set_size(int m, int n, int nrhs = 1);
        trivial_type& epsilon()    { return m_eps; }

        void solve(T* A, T* b);
    private:
        int             m_m;
        int             m_n;
        int             m_nrow;
        int             m_nrhs;
        int             m_lwork;
        trivial_type    m_eps;

        std::vector<value_type>     m_workspace;
};

template<>
void TikhonovQrSolver< std::complex<double> >::
set_size(int m, int n, int nrhs)
{
    m_m = m + n;
    m_n = n;
    m_nrhs = nrhs;
    m_nrow = m;

    // call the routine in query mode to get the optimal size of workspace
    std::complex<double> *a = NULL;
    std::complex<double> *b = NULL;

    int MINUS_ONE = -1;
    int ret;
    char TRANS = 'N';

    zgels(&TRANS, &m_m, &m_n, &m_nrhs,
          (MKL_Complex16*)a, &m_m,
          (MKL_Complex16*)b, &m_m,
          (MKL_Complex16*)&m_workspace[0], &MINUS_ONE,
          &ret);
    if ( ret )
    {
        fprintf(stderr, "ERROR: on zgels (query mode) %d\n", ret);
        exit(1);
    }

    // allocate workspace
    m_workspace.resize((int)m_workspace[0].real());
    m_lwork = m_workspace.size();
}

template<>
void TikhonovQrSolver< std::complex<double> >::
solve(std::complex<double>* A,
      std::complex<double>* b)
{
    using namespace std; 

    //// estimate the maximum singular value by forbenius norm
    double nrm = 0;
    complex<double> tv;
    for(int i = 0, rr = m_nrow;i < m_n;++ i, rr += m_m)
    {
        cblas_zdotc_sub(m_nrow, &A[i*m_m], 1, &A[i*m_m], 1, &tv);
        nrm += tv.real();
        memset(&A[rr], 0, sizeof(complex<double>)*m_n);
    }
    //// construct the new left-hand matrix A
    nrm = sqrt(nrm) * m_eps;
    for(int i = 0, rr = m_nrow;i < m_n;++ i, rr += (m_m+1))
        A[rr].real() = nrm;
    memset(&b[m_nrow], 0, sizeof(complex<double>)*m_n);

    int ret;
    char TRANS = 'N';
    zgels(&TRANS, &m_m, &m_n, &m_nrhs,
          (MKL_Complex16*)A, &m_m,
          (MKL_Complex16*)b, &m_m,
          (MKL_Complex16*)&m_workspace[0], &m_lwork,
          &ret);
    if ( ret )
    {
        fprintf(stderr, "ERROR: on zgels (query mode) %d\n", ret);
        exit(1);
    }
}

///////////////////////////////////////////////////////////////////////////////

template<typename T>
typename carbine::TrivialType<T>::type compute_residual(
        const int order, 
        const T *U, T *b, const T *x, 
        const int nrow, const int width);

/*
 * NOTE: after returning from this method, the vector b is overwritten by 
 *       the residual. i.e. b = Ax - b
 */
template<>
double compute_residual< std::complex<double> >(
        const int order,
        const std::complex<double> *A, 
              std::complex<double> *b,
        const std::complex<double> *x, 
        const int nrow, const int width)
{
    using namespace std;

    const complex<double> beta(-1, 0);
    complex<double>       alpha(1, 0);

    double nrm_b = cblas_dznrm2(nrow, b, 1);

    switch (order)
    {
        case CblasRowMajor:
            cblas_zgemv(CblasRowMajor, CblasNoTrans, nrow, width, &alpha, 
                        A, width, x, 1, &beta, b, 1);
            break;
        case CblasColMajor:
            cblas_zgemv(CblasColMajor, CblasNoTrans, nrow, width, &alpha, 
                        A, nrow, x, 1, &beta, b, 1);
            break;
        default:
            fprintf(stderr, "Unknown matrix layout order!\n");
            exit(1);
    }

    return cblas_dznrm2(nrow, b, 1) / nrm_b;
}

#endif
