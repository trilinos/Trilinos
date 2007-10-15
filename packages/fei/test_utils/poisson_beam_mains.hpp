/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#ifndef _poisson_cube_mains_hpp_
#define _poisson_cube_mains_hpp_


int poisson_main(int argc, char** argv,
		 MPI_Comm comm, int numProcs, int localProc);

int poisson3_main(int argc, char** argv,
                  MPI_Comm comm, int numProcs, int localProc);

int beam_oldfei_main(int argc, char** argv,
	      MPI_Comm comm, int numProcs, int localProc);

int beam_main(int argc, char** argv,
	       MPI_Comm comm, int numProcs, int localProc);

int feiDriver_main(int argc, char** argv,
		   MPI_Comm comm, int numProcs, int localProc);

int cFeiTester_main(int argc, char** argv,
                    MPI_Comm comm, int numProcs, int localProc);

#endif

