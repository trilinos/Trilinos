/*@HEADER
/*
************************************************************************

              Epetra: Linear Algebra Services Package 
                Copyright (2001) Sandia Corporation

Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
license for use of this work by or on behalf of the U.S. Government.

This library is free software; you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as
published by the Free Software Foundation; either version 2.1 of the
License, or (at your option) any later version.
 
This library is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.
 
You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
USA
Questions? Contact Michael A. Heroux (maherou@sandia.gov) 

************************************************************************
*/
//@HEADER
*/

/*
  file main.c
*/
#include <stdio.h>
#include "mpi.h"
#include "Petra_c_wrappers.h"

int main(int argc, char *argv[]) {
  int mypid;
  PETRA_COMM petra_comm;
  /*
     Start executable statements.
  */
  MPI_Init(&argc,&argv);
  petra_comm = petra_comm_create( MPI_COMM_WORLD );
  /*
    MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD,&myid);
  */
  mypid = petra_comm_getmypid(petra_comm);
  printf("MyPID = %d\n",mypid);
  MPI_Finalize();

  return 0;
}

/*
  end of file main.c
*/
