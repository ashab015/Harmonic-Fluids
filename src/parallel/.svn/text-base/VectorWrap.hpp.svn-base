#ifndef PARALLEL_VECTOR_WRAP_HPP
#   define PARALLEL_VECTOR_WRAP_HPP

#include <valarray>

//! Wrap an iterator such that it could be use in the naive 
//  OpenMP parallelization.
/*!
 * \author Changxi Zheng
 * \date 2008
 */
template <typename Iter>
class VectorWrap
{
    public:
        VectorWrap(Iter first, size_t len):m_array(len)
        {
            for(size_t i = 0;i < len;++ i)
                m_array[i] = first ++;
        }

        inline Iter operator [] (int n)
        {
            return m_array[n];
        }

        inline int size() const { return m_array.size(); }

    private:
        std::valarray<Iter>     m_array;
};
#endif
