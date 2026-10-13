/* rusage: run a command and print its user CPU seconds and peak RSS (kB)
 * on stderr - a stand-in for /usr/bin/time -f "%U %M" in images without it.
 * SPDX-License-Identifier: GPL-3.0-or-later */
#include <stdio.h>
#include <stdlib.h>
#include <sys/resource.h>
#include <sys/wait.h>
#include <unistd.h>

int main(int argc, char **argv)
{
    pid_t pid; int status; struct rusage ru;
    if (argc < 2) { fprintf(stderr, "usage: rusage cmd [args]\n"); return 2; }
    pid = fork();
    if (pid == 0) { execvp(argv[1], argv + 1); perror(argv[1]); _exit(127); }
    wait4(pid, &status, 0, &ru);
    fprintf(stderr, "RUSAGE user %.2f s sys %.2f s maxrss %ld kB\n",
            ru.ru_utime.tv_sec + ru.ru_utime.tv_usec / 1e6,
            ru.ru_stime.tv_sec + ru.ru_stime.tv_usec / 1e6, ru.ru_maxrss);
    return WIFEXITED(status) ? WEXITSTATUS(status) : 1;
}
