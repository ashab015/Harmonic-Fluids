#include "JobDispatcher.hpp"
#include "TerminalColor.h"

void JobDispatcher::dispatch(
        const TriangleMesh<double>& freeBd,
        const TriangleMesh<double>& leftBd,
        const TriangleMesh<double>& rightBd,
        const TriangleMesh<double>& frontBd,
        const TriangleMesh<double>& backBd,
        const std::set<TBubble*>& bubs,
        double ts)
{
    using namespace std;

    if ( bubs.empty() ) 
    {
        fprintf(stderr, "JOBDISP %lf %d %d\n", ts, 0, 0);
        return;
    }
    const TriangleMesh<double>* bds[] = { 
            &freeBd, &leftBd, &rightBd, &frontBd, &backBd };

    m_jobs.clear();
    m_nosolve.clear();
    set<TBubble*>::const_iterator end = bubs.end();
    for(set<TBubble*>::const_iterator it = bubs.begin();
            it != end;++ it)
        if ( (*it)->should_sample(ts) ) 
            m_jobs.push_back(*it);
        else
            m_nosolve.push_back(*it);

    fprintf(stderr, "JOBDISP %lf %d %d\n", ts, (int)bubs.size(), (int)m_jobs.size());

    int err, conf_cnt = 0;
    pthread_t tid;
    int nbub = m_jobs.size(), bub_ptr = 0;
    while ( nbub > 0 ) 
    {
        // allocate a cluster
        JobRunner* jrunner = m_scheduler.schedule(&m_jobs[bub_ptr], nbub, ts);
        if ( jrunner->need_configure() )
        {
            ++ conf_cnt;
            // start a new thread to send geometric data if it needs to
            // send sampling data to cluster
            if ( (err = pthread_create(&tid, 
                        &m_thdAttr, 
                        cluster_configure_thread,
                        (void *)(new cluster_config_arg(jrunner, bds)))) != 0 )
            {
                fprintf(stderr, "Cannot create thread: %s\n", strerror(err));
                exit(1);
            }
        }
        nbub    -= jrunner->num_bubbles();
        bub_ptr += jrunner->num_bubbles();
        // start a new thread (JobRunner.start) to solve for bubbles on cluster
        jrunner->start();
        m_jobThds.insert(jrunner);
    }
    m_scheduler.m_confNum = conf_cnt;

    // write un-sampled bubbles to file
    write_to_file(ts);
}

void JobDispatcher::write_to_file(double ts)
{
    using namespace std;

    static const int NOSAMPLE = -99;
    if ( m_nosolve.empty() ) return;

    ofstream& fout = m_scheduler.require_file_writer();

    fout.write((const char*)&ts, sizeof(double));
    for(size_t i = 0;i < m_nosolve.size();++ i)
    {
        fout.write((const char*)&(m_nosolve[i]->id), sizeof(int));
        fout.write((const char*)&(m_nosolve[i]->iso), sizeof(double));
        fout.write((const char*)&NOSAMPLE, sizeof(int));
    }
    fout.write((const char*)&NOSAMPLE, sizeof(int));

    m_scheduler.release_file_writer();
}

/*
 * The method should be called at each timestep before calling
 * JobDispatcher::dispatch method
 *
 * It clean all the finished job threads from job pool
 */
void JobDispatcher::clean_job_threads()
{
    using namespace std;

    vector<JobRunner*> finished;
    set<JobRunner*>::iterator end = m_jobThds.end();
    for(set<JobRunner*>::iterator it = m_jobThds.begin();
            it != end;++ it)
        if ( (*it)->is_finished() )
            finished.push_back(*it);

    for(size_t i = 0;i < finished.size();++ i)
    {
        m_jobThds.erase(finished[i]);
        delete finished[i];
    }
printf("Clean job threads: %d : %d\n",
        (int)m_jobThds.size(),
        (int)finished.size());
}
