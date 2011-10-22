#include <stdio.h>
#include <stdarg.h>

#include <mpi.h>

#include "commsplitter.h"


void commsplitter_log(char *fmt, ...)
{
    va_list args;
    FILE *fp = stdout;

    if (commsplitter_data.debug_level <= 0) {
        return;
    }

    va_start (args, fmt);
    fprintf (fp, "COMMSPLITTER (grank=%d): ", commsplitter_data.grank);
    vfprintf (fp, fmt, args);
    va_end (args);
    fflush (fp);
}
