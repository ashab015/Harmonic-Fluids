#include "daemon.h"
#include <fcntl.h>
#include <sys/resource.h>

void daemonize(const char* logfile)
{
    int     i, fd0, fd1, fd2;
    pid_t   pid;
    struct  rlimit rl;
    struct  sigaction sa;

    umask(S_IRUSR|S_IWUSR|S_IRGRP|S_IWGRP);
    /* get number of file descriptors */
    if ( getrlimit(RLIMIT_NOFILE, &rl) < 0 )
    {
        fprintf(stderr, "Cannot get file limit\n");
        exit(1);
    }

    if ( (pid = fork()) < 0 )
    {
        fprintf(stderr, "Cannot fork\n");
        exit(1);
    }
    else if ( pid != 0 ) /* parent process */
        exit(0);
    setsid();

    sa.sa_handler = SIG_IGN;
    sigemptyset(&sa.sa_mask);
    sa.sa_flags = 0;
    if ( sigaction(SIGHUP, &sa, NULL) < 0 )
    {
        fprintf(stderr, "Cannot ignore SIGHUP\n");
        exit(1);
    }
    if ( (pid = fork()) < 0 )
    {
        fprintf(stderr, "Cannot fork\n");
        exit(1);
    }
    else if ( pid != 0 ) /* parent process */
        exit(0);

    if ( chdir("~") < 0 )
    {
        fprintf(stderr, "Cannot change to home directory\n");
        exit(0);
    }

    if ( rl.rlim_max == RLIM_INFINITY )
        rl.rlim_max = 1024;
    for(i = 0;i < rl.rlim_max;++ i)
        close(i);

    fd0 = open(logfile, O_RDWR);
    fd1 = dup(0);
    fd2 = dup(0);

    if ( fd0 != 0 || fd1 != 0 || fd2 != 0 )
    {
        fprintf(stderr, "unexpected file descriptors %d %d %d\n",
                fd0, fd1, fd2);
        exit(1);
    }
}
