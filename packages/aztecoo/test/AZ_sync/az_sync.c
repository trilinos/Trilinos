void AZ_sync1(int proc_config[])
{
#include <stdio.h>
#include <stdlib.h>
#include "az_aztec.h"

MPI_Comm * comm = (MPI_Comm *) AZ_get_comm(proc_config);

MPI_Barrier(*comm);

} /* AZ_sync */
