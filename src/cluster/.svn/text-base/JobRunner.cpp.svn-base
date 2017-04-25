#include "JobRunner.hpp"
#include "JobScheduler.hpp"
//#include <mkl.h>

void * cluster_configure_thread(void * a)
{
    cluster_config_arg* arg = (cluster_config_arg*)a;

    int ret;
    if ( clnt_call(arg->job->m_clnt,
                HF_CONFIG,
                (xdrproc_t)xdr_samples_enc, (caddr_t)arg,
                (xdrproc_t)xdr_int, (caddr_t)&ret,
                RPC_CONFIG_TIMEOUT) != RPC_SUCCESS )
    {
        fprintf(stderr, "clnt_call Failed for child %s [Unicast]\n", 
                arg->job->m_cltName);
        clnt_perror(arg->job->m_clnt, "config call failed");
        exit(1);
    }

    pthread_mutex_lock(&(arg->job->m_configCond->mutex));
    arg->job->m_configCond->finished = true;
    pthread_mutex_unlock(&(arg->job->m_configCond->mutex));

    // notify JobRunner to start
    pthread_cond_broadcast(&(arg->job->m_configCond->cond));
    arg->job->m_scheduler->configure_finished();

    delete arg;
    return NULL;
}

void * thread_entrance(void * a)
{
    JobRunner* arg = (JobRunner*)a;
    arg->run();
    return NULL;
}

JobRunner::JobRunner(int cltId, double ts, const char* cltName, 
        TBubble** bptr, size_t nbub, bool needConfig, JobScheduler* js):
        m_cltId(cltId), m_cltName(cltName), 
        m_isFinished(false), m_needConfig(needConfig), 
        m_scheduler(js), m_configCond(NULL)
{
    if ( needConfig ) m_configCond = new ConfigCondVar;
    if ( !(m_clnt = clnt_create(m_cltName,
                    HARMONIC_FLUIDS_PROGRAM,
                    HARMONIC_FLUIDS_VERSION_1,
                    "tcp")) )
    {
        clnt_pcreateerror(m_cltName);
        exit(1);
    }
    memset(&m_ret, 0, sizeof(hf_results));

    m_bubs.ts = ts;
    // cache bubbles
    m_iso.resize(nbub);
    m_bubs.len = nbub;
    m_bubs.val = new hf_bubble[nbub];
    if ( !m_bubs.val )
    {
        fprintf(stderr, "Cannot allocate memory!\n");
        exit(1);
    }
    for(size_t i = 0;i < nbub;++ i)
    {
        m_iso[i] = bptr[i]->iso;
        m_bubs.val[i].id   = bptr[i]->id;
        m_bubs.val[i].freq = bptr[i]->omegad;
        m_bubs.val[i].rad  = bptr[i]->rad;
        m_bubs.val[i].pos  = bptr[i]->pos;

        m_bubs.val[i].regular_src[0].len = bptr[i]->numSrc[0];
        m_bubs.val[i].regular_src[0].val = new Point3d[bptr[i]->numSrc[0]];

        m_bubs.val[i].regular_src[1].len = bptr[i]->numSrc[1];
        m_bubs.val[i].regular_src[1].val = new Point3d[bptr[i]->numSrc[1]];

        //// make sure memory allocated
        if ( !m_bubs.val[i].regular_src[0].val || !m_bubs.val[i].regular_src[1].val )
        {
            fprintf(stderr, "Cannot allocate enough memory!\n");
            exit(1);
        }

        //// copy source positions
        for(int j = 0;j < bptr[i]->numSrc[0];++ j)
            m_bubs.val[i].regular_src[0].val[j] = *bptr[i]->regularSrc[0][j];
        for(int j = 0;j < bptr[i]->numSrc[1];++ j)
            m_bubs.val[i].regular_src[1].val[j] = *bptr[i]->regularSrc[1][j];
    } // end for
}

JobRunner::~JobRunner()
{
    if ( m_needConfig ) delete m_configCond;
    clnt_destroy(m_clnt);
    xdr_free((xdrproc_t)xdr_results, (char *)&m_ret);

    for(int i = 0;i < m_bubs.len;++ i)
    {
        delete [](m_bubs.val[i].regular_src[0].val);
        delete [](m_bubs.val[i].regular_src[1].val);
    }
    delete []m_bubs.val;
}

void JobRunner::start()
{
    int err;
    if ( (err = pthread_create(&m_threadId,
                    m_scheduler->thread_attributes(),
                    thread_entrance,
                    (void*)this)) != 0 )
    {
        fprintf(stderr, "Cannot create thread for JobRunner::start [%s]\n",
                strerror(err));
        exit(1);
    }
}

/*
 * call RPC on cluters to solve the least-square problem
 */
void JobRunner::run()
{
    // wait until configuration finished
    pthread_mutex_lock(&m_configCond->mutex);
    if ( !m_configCond->finished ) 
        pthread_cond_wait(&m_configCond->cond, &m_configCond->mutex);
    pthread_mutex_unlock(&m_configCond->mutex);

    if ( clnt_call(m_clnt,
                HF_SOLVE,
                (xdrproc_t)xdr_bubbles, (caddr_t)&m_bubs,
                (xdrproc_t)xdr_results, (caddr_t)&m_ret,
                RPC_SOLVE_TIMEOUT) != RPC_SUCCESS )
    {
        fprintf(stderr, "clnt_call Failed for child: %s\n",
                m_cltName);
        clnt_perror(m_clnt, "solve call failed\n");
        exit(1);
    }
    // tell scheduler to scheduler next job
    m_scheduler->job_finished(this);

    // write the results into file
    write_to_file();
    // change the state variable
    m_isFinished = true;
}

void JobRunner::write_to_file()
{
    using namespace std;
    static const int NOSAMPLE = -99;
    
    ofstream& fout = m_scheduler->require_file_writer();

    fout.write((const char*)&m_bubs.ts, sizeof(double));
    for(int i = 0;i < m_ret.len;++ i)
    {
        fout.write((const char*)&(m_bubs.val[i].id), sizeof(int));                  // bubble ID
        fout.write((const char*)&(m_iso[i]), sizeof(double));                       // bubble ISO

        fout.write((const char*)&(m_bubs.val[i].regular_src[1].len), sizeof(int));  // # of reg. src
        fout.write((const char*)m_bubs.val[i].regular_src[1].val, 
                   sizeof(Point3d) * m_bubs.val[i].regular_src[1].len);             // write fixed sources
        if ( m_ret.val[i].id != m_bubs.val[i].id )
        {
            fprintf(stderr, "ERROR: inconsistent bubble ID\n");
            exit(1);
        }
        fout.write((const char*)m_ret.val[i].c,
                sizeof(complex<double>) * m_ret.val[i].len);

        /*
        ///// REMOVE it for release
        fprintf(stderr, "JR %.16lf  %d  %d  %.16lf %.6lf %.6lf\n", 
                m_bubs.ts, 
                m_ret.val[i].id,
                m_ret.val[i].len,
                cblas_dznrm2(m_ret.val[i].len, m_ret.val[i].c, 1),
                m_ret.val[i].fittingError[0],       // $FITTING-ERROR
                m_ret.val[i].fittingError[1]);      // $FITTING-ERROR
        */
    }
    fout.write((const char*)&NOSAMPLE, sizeof(int));
    m_scheduler->release_file_writer();
}
