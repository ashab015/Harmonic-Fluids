/*!
 * ParallelIter.hpp
 * author: Changxi Zheng (cxzheng@cs.cornell.edu)
 */
#ifndef PARALLEL_ITERATION_H
#   define PARALLEL_ITERATION_H

#include <omp.h>
#include <memory>
#include <algorithm>

template<typename Iter>
class ParallelIter
{
    public:
        /*!
         * Constructor should be called out side of the parallel 
         * block
         */
        ParallelIter(Iter first, size_t len):m_cnt(0), m_curIter(first)
        {
            const int s = omp_get_max_threads();
            m_load = new int[s];
            memset(m_load, 0, sizeof(int)*s);
            // do round-robine to allocate load
            for(size_t i = 0;i < len;++ m_load[i%s], ++ i);
            m_numThreads = std::min(s, (int)len);
        }

        ~ParallelIter()
        {
            delete[] m_load;
        }

        int num_threads() const 
        {
            return m_numThreads;
        }

        int next_iteration_range(Iter& start, Iter& end)
        {
            int ret;
            #pragma omp critical(updateiterrange)
            {
                start = m_curIter;
                end = m_curIter;
                ret = m_cnt;
                for(size_t i = 0;i < m_load[m_cnt];++ i, ++ end);
                ++ m_cnt;
                m_curIter = end;
            } // end pragma
            return ret;
        }

        /* for test */
        void dump(std::ostream& out) const
        {
            const int num = omp_get_max_threads();
            for(int i = 0;i < num;++ i)
                out << ' ' << m_load[i];
            out << std::endl;
        }

    private:
        int*    m_load;
        int     m_cnt;
        Iter    m_curIter;
        int     m_numThreads;
};

#endif
