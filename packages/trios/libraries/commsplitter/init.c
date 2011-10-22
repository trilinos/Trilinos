#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>

#include <mpi.h>

#include "commsplitter.h"



commsplitter_t commsplitter_data;


void commsplitter_getenv(void)
{
    char *envstr=getenv("COMMSPLITTER_DEBUG_LEVEL");
    if (envstr != NULL) {
        commsplitter_data.debug_level=atoi(envstr);
    }
}

void commsplitter_init(char *app_name)
{
    int   hostname_len;
    char  hostname[MPI_MAX_PROCESSOR_NAME];

    commsplitter_data.enabled     = 1;
    commsplitter_data.debug_level = 0;
    commsplitter_data.split_comm  = MPI_COMM_NULL;

    PMPI_Comm_rank(MPI_COMM_WORLD, &commsplitter_data.grank);
    PMPI_Comm_size(MPI_COMM_WORLD, &commsplitter_data.gsize);
    PMPI_Get_processor_name(hostname, &hostname_len);

    commsplitter_getenv();

    if (commsplitter_data.grank == 0) {
        commsplitter_log("\n");
        commsplitter_log("commsplitter initialized\n");
        commsplitter_log("\n");
    }

    commsplitter_log("app_name is %s\n", app_name);
    commsplitter_log("successful init of grank=%d on %s\n", commsplitter_data.grank, hostname);

    return;
}

void commsplitter_finalize()
{
    return;
}
