
C @HEADER
C ************************************************************************
C 
C               Epetra: Linear Algebra Services Package 
C                 Copyright 2011 Sandia Corporation
C 
C Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
C the U.S. Government retains certain rights in this software.
C
C Redistribution and use in source and binary forms, with or without
C modification, are permitted provided that the following conditions are
C met:
C
C 1. Redistributions of source code must retain the above copyright
C notice, this list of conditions and the following disclaimer.
C
C 2. Redistributions in binary form must reproduce the above copyright
C notice, this list of conditions and the following disclaimer in the
C documentation and/or other materials provided with the distribution.
C
C 3. Neither the name of the Corporation nor the names of the
C contributors may be used to endorse or promote products derived from
C this software without specific prior written permission.
C
C THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
C EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
C IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
C PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
C CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
C EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
C PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
C PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
C LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
C NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
C SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
C
C Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
C 
C ************************************************************************
C @HEADER

       program main
C
C---------------------------------------------------------------
C
       implicit none
       include "mpif.h"
#include "Epetra_Fortran_wrappers.h"
C
       integer i
       integer   nlocal, ierror, izero, nglobal
       parameter(izero=0)
       double precision norm1a, norm2a, norm1b, norm1apb, mynorm1b
       double precision local_mynorm1b
       double precision zero, one, two
       parameter(zero=0.0, one=1.0, two=2.0)
       EPETRA_OBJECT_PTR comm, map, vec_a, vec_b, vec_apb
       parameter (nlocal=10)
       double precision b_values(nlocal)
C
C 
       call MPI_INIT(ierror)

C     Create Epetra Comm Object
       comm = epetra_mpicomm_create1()

C     Create Epetra Map with "n" element spread evenly across the machine
C     and index base equal to 1.  Set nglobal to -1 to let map compute true nglobal value
       nglobal = -1 
       map = epetra_map_create2(nglobal, nlocal, 1, comm)
       nglobal = epetra_map_numglobalelements(map)

C     Create Epetra Vector a and fill it with the value 2.0 in each entry
       vec_a = epetra_vector_create1(map)
       ierror = epetra_vector_putscalar(vec_a, two)

C     Create vector b using storage I create and fill it with random numbers
       vec_b = epetra_vector_create2(izero, map, b_values)
       ierror = epetra_vector_random(vec_b)

C     Create vector apb and define it to be a+b
       vec_apb = epetra_vector_create1(map)
       ierror = epetra_vector_update(vec_apb, one, vec_a, one, 
     $      vec_b, zero)
       

C     Compute one and two norm of vector a
       ierror = epetra_vector_norm1(vec_a, norm1a)
       ierror = epetra_vector_norm2(vec_a, norm2a)
       
C     Compute one norms of b and apb
       ierror = epetra_vector_norm1(vec_b, norm1b)
       ierror = epetra_vector_norm1(vec_apb, norm1apb)

       if (nlocal.lt.20) then
          if (epetra_comm_mypid(comm).eq.0) print*,'Vector A = '
          call epetra_vector_print(vec_a)
          if (epetra_comm_mypid(comm).eq.0) print*,'Vector B = '
          call epetra_vector_print(vec_b)
          if (epetra_comm_mypid(comm).eq.0) print*,'Vector A + B = '
          call epetra_vector_print(vec_apb)
       endif
       
       local_mynorm1b = 0.0
       do 10 i=1,nlocal
          local_mynorm1b = local_mynorm1b + abs(b_values(i))
 10    continue

       call MPI_ALLREDUCE(local_mynorm1b, mynorm1b, 1, 
     $                    MPI_DOUBLE_PRECISION, MPI_SUM, 
     $                    MPI_COMM_WORLD, ierror)

       
C     write results
       if (epetra_comm_mypid(comm).eq.0) then
          write(*,99) 'One Norm of Vector A of 2s  = ', norm1a
          write(*,99) '                    Should  = ', two*nglobal
          write(*,99) 'Two Norm of Vector A of 2s  = ', norm2a
          write(*,99) 'One Norm of Random Vector B = ', norm1b
          write(*,99) '                    Should  = ', mynorm1b
          write(*,99) 'One Norm of A plus B        = ', norm1apb
       endif
 99    format(' ',a,f15.7)

C     Delete all object we created
       call epetra_vector_destroy(vec_apb)
       call epetra_vector_destroy(vec_b)
       call epetra_vector_destroy(vec_a)
       call epetra_map_destroy(map)
       call epetra_comm_destroy(comm)

       call MPI_FINALIZE(ierror)
C
       stop
       end

