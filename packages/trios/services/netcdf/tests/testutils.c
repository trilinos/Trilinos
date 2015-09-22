#include <stdio.h>
#include <limits.h>
#include <string.h>
#include "mpi.h"
#include "testutils.h"

void parse_read_args(int argc, char **argv, int rank, params *p)
{
	int inlen, outlen;
	if ( rank == 0 ) {
		if (argc == 3 ) {
			strncpy(p->infname, argv[1], PATH_MAX);
			strncpy(p->outfname, argv[2], PATH_MAX);
		} else if (argc == 0) {
			strncpy(p->infname, "pvfs:../data/test_double.nc", 
					PATH_MAX);
			strncpy(p->outfname, "pvfs:testread.nc", PATH_MAX);
		} else {
			fprintf(stderr, "Usage: %s: <source> <destination>\n", 
					argv[0]);
			MPI_Abort(MPI_COMM_WORLD, 1);
		}
		inlen = strlen(p->infname);
		outlen = strlen(p->outfname);
	}
	MPI_Bcast(&inlen, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(p->infname, inlen+1, MPI_CHAR, 0, MPI_COMM_WORLD);
	MPI_Bcast(&outlen, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(p->outfname, outlen+1, MPI_CHAR, 0, MPI_COMM_WORLD);
}

void parse_write_args(int argc, char **argv, int rank, params *p)
{
	int outlen;
	if ( rank == 0 ) {
		if (argc == 2 ) {
			strncpy(p->outfname, argv[1], PATH_MAX);
		} else if (argc == 0) {
			strncpy(p->outfname, "pvfs:testwrite.nc", PATH_MAX);
		} else {
			fprintf(stderr, "Usage: %s: <destination>\n", argv[0]);
			MPI_Abort(MPI_COMM_WORLD, 1);
		}
		outlen = strlen(p->outfname);
	}
	MPI_Bcast(&outlen, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(p->outfname, outlen+1, MPI_CHAR, 0, MPI_COMM_WORLD);
}
