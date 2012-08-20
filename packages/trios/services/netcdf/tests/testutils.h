#ifndef _UTILS_H
#define _UTILS_H

#include <sys/param.h>
#include <limits.h>

typedef struct {
    char infname[PATH_MAX];
    char outfname[PATH_MAX];
} params;

void parse_read_args(int argc, char **argv, int rank, params *p);
void parse_write_args(int argc, char **argv, int rank, params *p);
#endif
