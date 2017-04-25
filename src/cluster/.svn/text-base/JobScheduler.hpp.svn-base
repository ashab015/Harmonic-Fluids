/******************************************************************************
 *  File: JobScheduler.hpp
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
#ifndef CLUSTER_JOB_SCHEDULER_HPP
#   define CLUSTER_JOB_SCHEDULER_HPP

#include <pthread.h>
#include <vector>
#include <string>
#include <string.h>
#include "Bubble.hpp"

struct ConfigCondVar;
class JobRunner;
class JobDispatcher;

class JobScheduler
{
    friend class JobDispatcher;

    public:
        typedef Bubble<double>  TBubble;

        JobScheduler(const char* conf, std::ofstream& fout);
        ~JobScheduler()
        {
            delete []m_cltLoad;
            delete []m_cltUsed;

            pthread_mutex_destroy(&m_bufLock);
            pthread_cond_destroy(&m_bufCond);
            pthread_mutex_destroy(&m_fileLock);
            pthread_mutex_destroy(&m_confLock);
            pthread_cond_destroy(&m_confCond);
            pthread_attr_destroy(&m_thdAttr);
        }

        /*
         * this method should be called at each timestep before scheduling
         * and before recomputing boundary samples
         */
        void next_step()
        {
            pthread_mutex_lock(&m_confLock);
            //// wait until all the boundary sample data have been sent
            while ( m_confFinishedNum < m_confNum )
                pthread_cond_wait(&m_confCond, &m_confLock);
            pthread_mutex_unlock(&m_confLock);

            memset(m_cltUsed, NULL, sizeof(ConfigCondVar*)*m_numClts);
            m_confNum = m_confFinishedNum = 0;
        }

        /* schedule the current bubbles */
        JobRunner* schedule(TBubble** bub_ptr, size_t nbub, double ts);

        /* this method is called when results are returned from some cluster */
        void job_finished(const JobRunner* runner);      
        void configure_finished();

        std::ofstream& require_file_writer()
        {
            pthread_mutex_lock(&m_fileLock);
            return m_fout;
        }
        void release_file_writer()
        {
            pthread_mutex_unlock(&m_fileLock);
        }

        const pthread_attr_t* thread_attributes() const
        {
            return &m_thdAttr;
        }

        const std::vector<Point3d>& fixed_sources() const
        {   return m_fixedSrcs; }
        const Point3d* get_source_ptr(int a, int b) const
        {   return &m_srcs[a][b]; }

    private:
        void generate_sources();

        std::vector<Point3d> m_srcs[2];
        std::vector<Point3d> m_fixedSrcs;

    private:
        void print_load();

        std::ofstream&                  m_fout;
        pthread_mutex_t                 m_bufLock;      // lock for buffer
        pthread_cond_t                  m_bufCond;      // if no buffer availble, wait on this
                                                        // condition variable
        pthread_mutex_t                 m_fileLock;
        pthread_mutex_t                 m_confLock;
        pthread_cond_t                  m_confCond;
        pthread_attr_t                  m_thdAttr;

        int                             m_confNum;    // how many config calls
        int                             m_confFinishedNum; // how many config calls have been finished
        size_t                          m_numClts;
        std::vector<std::string>        m_cltNames;
        std::vector<size_t>             m_cltCores;
        std::vector<size_t>             m_cltBufLen;

        int*                            m_cltLoad;      // current load on cluster
        int*                            m_cltJobNum;    // number of jobs
        ConfigCondVar**                 m_cltUsed;
};

#endif
