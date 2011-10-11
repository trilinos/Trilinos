#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>

#include <mpi.h>

#include "commsplitter.h"



mpics_t mpics_data;


void mpics_getenv(void)
{
    char *envstr=getenv("MPICS_DEBUG_LEVEL");
    if (envstr != NULL) {
        mpics_data.debug_level=atoi(envstr);
    }
}

void mpics_init(char *app_name)
{
    int   hostname_len;
    char  hostname[MPI_MAX_PROCESSOR_NAME];

    mpics_data.enabled     = 1;
    mpics_data.debug_level = 0;
    mpics_data.split_comm  = MPI_COMM_NULL;

    PMPI_Comm_rank(MPI_COMM_WORLD, &mpics_data.grank);
    PMPI_Comm_size(MPI_COMM_WORLD, &mpics_data.gsize);
    PMPI_Get_processor_name(hostname, &hostname_len);

    mpics_getenv();

    if (mpics_data.grank == 0) {
        mpics_log("\n");
        mpics_log("mpics (comm splitter) initialized\n");
        mpics_log("\n");
    }

    mpics_log("app_name is %s\n", app_name);
    mpics_log("successful init of grank=%d on %s\n", mpics_data.grank, hostname);

    return;
}

void mpics_finalize()
{
    return;
}
