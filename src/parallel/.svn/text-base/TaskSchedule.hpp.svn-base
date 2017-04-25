#ifndef PARALLEL_TASK_SCHEDULE_HPP
#   define PARALLEL_TASK_SCHEDULE_HPP

#include <vector>
#include <algorithm>
#include <omp.h>

template <typename Task>
class StaticTaskSchedule
{
    public:
        StaticTaskSchedule():m_len(0), m_cnt(0),
            m_nthr(omp_get_max_threads()),
            m_tasks(m_nthr)
        {}

        //! call it from out of the parallel region
        void add_task(Task t)
        {
            m_tasks[m_cnt].push_back(t);
            m_cnt = (m_cnt + 1) % m_nthr;
            ++ m_len;
        }

        inline int num_threads() const
        {
            return std::min(m_nthr, m_len);
        }

        inline const std::vector<Task>& get_task(int id) const
        {
            return m_tasks[id];
        }

    private:
        int             m_len;
        int             m_cnt;
        const int       m_nthr;
        std::vector< std::vector<Task> > m_tasks;
};

#endif
