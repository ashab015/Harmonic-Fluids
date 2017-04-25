/******************************************************************************
 *  File: JobDispatcher.hpp
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
#ifndef CLUSTER_JOB_DISPATCHER_HPP
#   define CLUSTER_JOB_DISPATCHER_HPP

#include <pthread.h>
#include <set>
#include "JobScheduler.hpp"
#include "JobRunner.hpp"
#include "geometric/TriangleMesh.hpp"

/*
 * Dispatch the multipole solving to clusters.
 *
 * ================================================
 * Configuration file:
 * <hostname>    <# of cores>   <buffer length>
 * <hostname>    <# of cores>   <buffer length>
 */
class JobDispatcher
{
    public:
        typedef Bubble<double>      TBubble;

        JobDispatcher(const char* conf, std::ofstream& fout):
                m_scheduler(conf, fout)
        { 
            if ( pthread_attr_init(&m_thdAttr) )
            {
                fprintf(stderr, "Cannot initialize thread attributes\n");
                exit(1);
            }
            if ( pthread_attr_setdetachstate(&m_thdAttr, 
                        PTHREAD_CREATE_DETACHED) )
            {
                fprintf(stderr, "Cannot set thread attributes\n");
                exit(1);
            }
        }
        ~JobDispatcher()
        {
            pthread_attr_destroy(&m_thdAttr);
        }

        void clean_job_threads();
        /* dispatch the bubble solving at this timestep to clusters */
        void dispatch(const TriangleMesh<double>& freeBd, 
                const TriangleMesh<double>& leftBd,
                const TriangleMesh<double>& rightBd,
                const TriangleMesh<double>& frontBd,
                const TriangleMesh<double>& backBd,
                const std::set<TBubble*>& bubs, double ts);
        JobScheduler& scheduler() { return m_scheduler; }

    private:
        void write_to_file(double ts);

        std::vector<TBubble*>       m_jobs;
        std::vector<TBubble*>       m_nosolve;  // the bubbles that does't need to be solved at this timestep
        std::set<JobRunner*>        m_jobThds;  // all the current jobRunner threads
        JobScheduler                m_scheduler;
        pthread_attr_t              m_thdAttr;
};

#endif
