#include <rpc/types.h>
#include <rpc/xdr.h>
#include <rpc/auth.h>
#include <rpc/clnt.h>
#include <rpc/pmap_clnt.h>
#include <map>
#include <pthread.h>
#include <queue>
#include "TerminalColor.h"
#include "harmonic_fluids.hpp"
#include "utils/IdGenerator.h"
#include "MultipoleSolver.hpp"

#ifdef USE_MKL
#   include <mkl.h>
#endif

/*
 * It provides two methods: configure and solve.
 *
 * For configuration, it receives the boundary samples, associates it with 
 * a job ID.
 */
static int                           buf_len = 0;
static hf_samples*                   sample_buf = NULL;
static std::map<double, hf_samples*> sample_map;
static int                           sample_ptr;

static pthread_t        p_thread;
static pthread_attr_t   thread_attr;
static pthread_mutex_t  proc_lock;
static pthread_cond_t   proc_cond;

static bool as_daemon   = false;
static char* logfile    = NULL;
static FILE* logfd      = stderr;
static carbine::CircularIdGenerator id_gen(100000000);
static std::queue<int>  job_queue;

static MultipoleSolver solver;

struct thread_arg
{
    int         tid;
    hf_bubbles  bubs;
    SVCXPRT *   transp;

    thread_arg(SVCXPRT* t):transp(t) 
    {
        memset(&bubs, 0, sizeof(hf_bubbles));
    }
};

static void usage()
{
    printf("%s%dmUsage: hfsvcd -d -c <buffer length> -l <log file>%s\n", 
            TERM_COLOR_BEGIN, TERM_COLOR_RED, TERM_COLOR_END);
    printf("%s%dm", TERM_COLOR_BEGIN, TERM_COLOR_YELLOW);
    //printf("       -m <int>   maximum memory allowed to be allocated by MKL (in megabytes)\n");
    printf("       -d         start as a daemon process\n");
    printf("       -c <int>   specify the local buffer length\n");
    printf("       -l <file>  log file\n");
    printf("       -h         display help information\n");
    printf(TERM_COLOR_END);
}

static void parse_cmd(int argc, char* argv[])
{
    int opt;

    while ( (opt = getopt(argc, argv, "hdl:c:m:")) != -1 )
    {
        switch (opt)
        {
            /*
            case 'm':
                max_mkl_mem = atoll(optarg) * 1024L * 1024L;
                printf("memory limitation: %lld bytes\n", max_mkl_mem);
                break;
            */
            case 'd':
                as_daemon = true;
                break;
            case 'c':
                buf_len = atoi(optarg);
                break;
            case 'l':
                logfile = optarg;
                break;
            case 'h':
                usage();
                exit(0);
            default:
                usage();
                exit(1);
        }
    }
}

static hf_samples* next_sample_ptr()
{
    hf_samples* ret;

    sample_ptr = (sample_ptr + 1) % buf_len;
    ret = &sample_buf[sample_ptr];

    if ( ret->surf_samples[0].len )
    {
        sample_map.erase(ret->ts);
        //xdr_free((xdrproc_t)xdr_samples, (char *)ret); 
        XDR x;
        x.x_op = XDR_FREE;
        if ( !xdr_samples(&x, ret) )
        {
            fprintf(logfd, "Fail to xdr_free for xdr_samples\n");
            exit(1);
        }
    }
    memset(ret, 0, sizeof(hf_samples));
    return ret;
}

static void* solve_least_square(void* a)
{
    thread_arg* arg = (thread_arg*)a;

    pthread_mutex_lock(&proc_lock);
    job_queue.push(arg->tid);
    while ( job_queue.front() != arg->tid ) 
        pthread_cond_wait(&proc_cond, &proc_lock);
    pthread_mutex_unlock(&proc_lock);

#ifdef USE_MKL
    //// ask MKL to deallocate memory if necessary
    int mkl_mem_buf;
    MKL_INT64 mkl_mem = MKL_MemStat(&mkl_mem_buf);

    fprintf(logfd, "# MKL allocated %lld bytes in %d buffers\n", mkl_mem, mkl_mem_buf);

    if ( mkl_mem > 0 ) MKL_FreeBuffers();
#endif

    //// solve the least-square problem
    if ( !sample_map.count(arg->bubs.ts) )
    {
        fprintf(logfd, "Cannot find the samples for job [ts:%lf] in buffer\n", 
                arg->bubs.ts);
        fflush(logfd);
        exit(1);
    }
    const hf_results* ret = solver.solve(sample_map[arg->bubs.ts], &arg->bubs);

    fprintf(logfd, "MSG: finish solver.solve(...)\n");
    fflush(logfd);

    if ( !svc_sendreply(arg->transp, (xdrproc_t)xdr_results, (char *)ret) )
    {
        svcerr_systemerr(arg->transp);
        fprintf(logfd, "Fail to send reply for SOLVE!\n");
        fflush(logfd);
        exit(1);
    }

    //// free memory
    //xdr_free((xdrproc_t)xdr_bubbles, (caddr_t)&(arg->bubs));
    XDR x;
    x.x_op = XDR_FREE;
    if ( !xdr_bubbles(&x, &(arg->bubs)) )
    {
        fprintf(logfd, "Fail to xdr_free for xdr_bubbles\n");
        exit(1);
    }
    delete arg;

    pthread_mutex_lock(&proc_lock);
    job_queue.pop();
    pthread_mutex_unlock(&proc_lock);
    pthread_cond_broadcast(&proc_cond);
    return NULL;
}

/*
 * This method is called sequentially
 */
static void harmonic_fluids_program_1(struct svc_req *rqstp, register SVCXPRT *transp)
{
    thread_arg* arg;
    hf_samples* samples;
    int config_ret;

    switch (rqstp->rq_proc)
    {
        case HF_CONFIG:
            samples = next_sample_ptr();

            if ( !svc_getargs(transp, (xdrproc_t)xdr_samples, (caddr_t)samples) )
            {
                svcerr_decode(transp);
                fprintf(logfd, "ERROR: cannot get arguments\n");
                fflush(logfd);
                exit(1);
            }
            if ( sample_map.count(samples->ts) )
            {
                fprintf(logfd, "ERROR: job [ts:%lf] has alread existed in the local buffer\n",
                        samples->ts);
                fflush(logfd);
                exit(1);
            }
            // associate this sample with the job ID
            sample_map[samples->ts] = samples;

            config_ret = 0;
            if ( !svc_sendreply(transp, (xdrproc_t)xdr_int, (char *)&config_ret) )
            {
                svcerr_systemerr(transp);
                fprintf(logfd, "Fail to send reply for CONFIG!\n");
                fflush(logfd);
                exit(1);
            }
            return;

        case HF_SOLVE:
            arg = new thread_arg(transp);
            arg->tid = id_gen.next_id();
            if ( !svc_getargs(transp, (xdrproc_t)xdr_bubbles, (caddr_t)&(arg->bubs)) )
            {
                svcerr_decode(transp);
                return;
            }
            if ( pthread_create(&p_thread, 
                        &thread_attr, 
                        solve_least_square, 
                        (void *)arg) )
            {
                fprintf(logfd, "Fail to create thread!");
                fflush(logfd);
                exit(1);
            }
            return;

        default:
            svcerr_noproc(transp);
            return;
    }
}

int main (int argc, char **argv)
{
    parse_cmd(argc, argv);

    if ( buf_len == 0 )
    {
        fprintf(stderr, "Please specify buffer length with a positive value\n");
        exit(1);
    }

    //// allocate buffer
    sample_buf = new hf_samples[buf_len];
    if ( !sample_buf )
    {
        fprintf(stderr, "Cannot allocate memory for hf_samples\n");
        exit(1);
    }
    memset(sample_buf, 0, sizeof(hf_samples)*buf_len);
    sample_ptr = -1;

    //// prepare for multi-thread
    if ( pthread_attr_init(&thread_attr) )
    {
        fprintf(stderr, "Cannot initialize thread attributes\n");
        exit(1);
    }
    if ( pthread_attr_setdetachstate(&thread_attr, PTHREAD_CREATE_DETACHED) )
    {
        fprintf(stderr, "Fail to set thread attribute\n");
        exit(1);
    }

    pthread_mutex_init(&proc_lock, NULL);
    pthread_cond_init(&proc_cond, NULL);

    if ( as_daemon && !logfile )
    {
        fprintf(stderr, "Specify the log file before running as daemon\n");
        usage();
        exit(1);
    }
    if ( as_daemon )
    {
        if ( (logfd = fopen(logfile, "w")) == NULL )
        {
            fprintf(stderr, "Cannot open log file: %s\n", logfile);
            exit(1);
        }
        fprintf(logfd, "====== LOG: ======\n");
        if ( daemon(0, 0) < 0 )
        {
            fprintf(logfd, "ERROR: cannot do daemon\n");
            exit(1);
        }
        solver.logfd = logfd;
        fprintf(logfd, "MSG: Begin as a daemon\n");
        fprintf(solver.logfd, "MSG: Change log file descriptor\n");
        fprintf(logfd, "MSG: buffer length = %d\n", buf_len);
        fflush(solver.logfd);
    }

    register SVCXPRT *transp;
    pmap_unset(HARMONIC_FLUIDS_PROGRAM, HARMONIC_FLUIDS_VERSION_1);
    transp = svctcp_create(RPC_ANYSOCK, 0, 0);

    if (transp == NULL) 
    {
        fprintf(logfd, "cannot create tcp service.");
        fflush(logfd);
        exit(1);
    }
    if (!svc_register(transp, 
                HARMONIC_FLUIDS_PROGRAM, 
                HARMONIC_FLUIDS_VERSION_1, 
                harmonic_fluids_program_1, 
                IPPROTO_TCP)) 
    {
        fprintf(logfd, "unable to register (HARMONIC_FLUIDS_PROGRAM, HARMONIC_FLUIDS_VERSION_1, tcp).");
        fflush(logfd);
        exit(1);
    }

    svc_run();
    fprintf(logfd, "harmonic fluid service exited");
    fflush(logfd);
    exit(1);
    /* NOTREACHED */
}
