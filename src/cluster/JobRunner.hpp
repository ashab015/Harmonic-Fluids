/******************************************************************************
 *  File: JobRunner.hpp
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
#ifndef CLUSTER_JOB_RUNNER_HPP
#   define CLUSTER_JOB_RUNNER_HPP

#include <pthread.h>
#include "harmonic_fluids.hpp"
#include "Bubble.hpp"

struct ConfigCondVar
{
    pthread_mutex_t     mutex;
    pthread_cond_t      cond;
    bool                finished;

    ConfigCondVar():finished(false)
    {
        pthread_mutex_init(&mutex, NULL);
        pthread_cond_init(&cond, NULL);
    }

    ~ConfigCondVar()
    {
        pthread_mutex_destroy(&mutex);
        pthread_cond_destroy(&cond);
    }
};

class JobScheduler;
/*
 * This class running as a separated thread is responsible for RPC invoke
 * to solve the least-square problems on clusters.
 */
class JobRunner
{
    friend void* cluster_configure_thread(void *);
    friend void* thread_entrance(void *);
    friend class JobScheduler;

    public:
        typedef Bubble<double>  TBubble;

        JobRunner(int cltId, double ts, const char* cltName, 
                  TBubble** bptr, size_t nbub, bool needConfig, 
                  JobScheduler* js);

        ~JobRunner();

        void start();   /* start to run the current job */

        inline double timestamp() const     { return m_bubs.ts; }
        inline bool need_configure() const  { return m_needConfig; }
        inline int  num_bubbles() const     { return m_bubs.len; }
        inline bool is_finished() const     { return m_isFinished; }
        CLIENT*     get_rpc_client();

    private:
        void run();
        void write_to_file();

    private:
        int                 m_cltId;
        const char*         m_cltName;
        /* indicate if this job has been finished */
        bool                m_isFinished;
        const bool          m_needConfig;
        JobScheduler*       m_scheduler;

        ConfigCondVar*      m_configCond;
        CLIENT*             m_clnt;     // RPC client
        std::vector<double> m_iso;
        hf_bubbles          m_bubs;
        hf_results          m_ret;
        pthread_t           m_threadId;
};

struct cluster_config_arg
{
    JobRunner*      job;
    const TriangleMesh<double>* samples[5];

    cluster_config_arg(JobRunner* j, const TriangleMesh<double>* s[]):
            job(j)
    {
        memcpy(samples, s, sizeof(const TriangleMesh<double>*)*5);
    }
};

// thread function to send out the boundary 
// samples to clusters
void * cluster_configure_thread(void *);
void * thread_entrance(void *);

#endif
