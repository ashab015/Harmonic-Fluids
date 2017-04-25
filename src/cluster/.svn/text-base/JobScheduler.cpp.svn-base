#include "JobScheduler.hpp"
#include "JobRunner.hpp"
#include <fstream>

JobScheduler::JobScheduler(const char* conf, std::ofstream& fout):
        m_confFinishedNum(0), m_confNum(0), m_fout(fout)
{
    using namespace std;

    ///// load configuration file
    char name[512];
    ifstream fin(conf);
    if ( fin.fail() )
    {
        fprintf(stderr, "Cannot open configure file: %s\n", conf);
        exit(1);
    }

    size_t nc, mc;
    CLIENT* clnt;
    while ( fin >> name >> nc >> mc)
    {
        // test if the cluster is providing service
        if ( !(clnt = clnt_create(name,
                        HARMONIC_FLUIDS_PROGRAM,
                        HARMONIC_FLUIDS_VERSION_1,
                        "tcp")) )
        {
            clnt_pcreateerror(name);
            exit(1);
        }

        m_cltNames.push_back(name);
        m_cltCores.push_back(nc);
        m_cltBufLen.push_back(mc);

        printf("CLUSTER: %s  [%d cores]\n", name, (int)nc);
        clnt_destroy(clnt);
    }
    m_numClts = m_cltCores.size();

    /////// allocate memory
    m_cltLoad   = new int[m_numClts];
    m_cltJobNum = new int[m_numClts];
    m_cltUsed   = new ConfigCondVar*[m_numClts];
    if ( !m_cltLoad || !m_cltJobNum || !m_cltUsed )
    {
        fprintf(stderr, "Cannot allocate memory for m_cltLoad or m_cltUsed\n");
        exit(1);
    }

    memset(m_cltLoad, 0, sizeof(int)*m_numClts);
    memset(m_cltJobNum, 0, sizeof(int)*m_numClts);

    ////// initialize mutex and condition variable
    pthread_mutex_init(&m_bufLock, NULL);
    pthread_cond_init(&m_bufCond, NULL);
    pthread_mutex_init(&m_fileLock, NULL);
    pthread_mutex_init(&m_confLock, NULL);
    pthread_cond_init(&m_confCond, NULL);
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

/*
 * Schedule Alg:
 * - find the cluster with lowest load/bufLen
 * - if there is a tie, find the one with minimum number of jobs
 *   (NOTE: A job may contain multiple LS solves)
 */
JobRunner* JobScheduler::schedule(TBubble** bptr, size_t nbub, double ts)
{
    int cltId = -1;     // the selected cluster ID to run the next job
    size_t num_bub;

    pthread_mutex_lock(&m_bufLock);
    do 
    {
        for(size_t i = 0;i < m_numClts;++ i)
            if ( m_cltLoad[i] < m_cltBufLen[i] )    // have buffer available on that cluster
            {
                if ( cltId < 0 )
                    cltId = (int)i;
                else
                {
                    int t = m_cltLoad[i]*m_cltBufLen[cltId] - 
                            m_cltLoad[cltId]*m_cltBufLen[i];
                    if ( t < 0 || (t == 0 && 
                         m_cltJobNum[i]*m_cltCores[cltId] < m_cltJobNum[cltId]*m_cltCores[i]) )
                        cltId = (int)i;
                } // end else
            } // end if

        // if all the clusters run out of resource, 
        // we have to wait on the m_bufCond condition variable
        if ( cltId < 0 ) 
            pthread_cond_wait(&m_bufCond, &m_bufLock);

    } while (cltId < 0);
    
    num_bub = std::min(m_cltCores[cltId], nbub);    // how many bubble to solve

    //// update the load
    ++ m_cltLoad[cltId];
    m_cltJobNum[cltId] += num_bub;
print_load();
    pthread_mutex_unlock(&m_bufLock);

    //// create JobRunner
    JobRunner* ret;
    // if this cluster is the first time to be used at the current timestep,
    // we have to send configuration data there
    if ( m_cltUsed[cltId] )
    {
        ret = new JobRunner(cltId, ts, m_cltNames[cltId].c_str(), 
                bptr, num_bub, false, this);
        ret->m_configCond = m_cltUsed[cltId];
    }
    else
    {
        ret = new JobRunner(cltId, ts, m_cltNames[cltId].c_str(), 
                bptr, num_bub, true, this);
        m_cltUsed[cltId] = ret->m_configCond;
    }

    return ret;
}

void JobScheduler::print_load()
{
    printf("LOAD:: ");
    for(int i = 0;i < m_numClts;++ i)
    {
        if ( i ) printf(" | ");
        printf("%d:%d", m_cltLoad[i], m_cltJobNum[i]);
    }
    printf("\n");
}
/*
 * Change the cltLoad and numJob record, and signal to buffer 
 * condition variable
 */
void JobScheduler::job_finished(const JobRunner* runner)
{
    pthread_mutex_lock(&m_bufLock);
    -- m_cltLoad[runner->m_cltId];
    m_cltJobNum[runner->m_cltId] -= runner->num_bubbles();
    pthread_mutex_unlock(&m_bufLock);

    // more slots available
    pthread_cond_signal(&m_bufCond);
}

void JobScheduler::configure_finished()
{
    pthread_mutex_lock(&m_confLock);
    ++ m_confFinishedNum;
    if ( m_confFinishedNum == m_confNum )
        pthread_cond_signal(&m_confCond);
    pthread_mutex_unlock(&m_confLock);
}

