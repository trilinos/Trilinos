!! 
!! @HEADER
!! *****************************************************************************
!!  Zoltan Toolkit for Load-balancing, Partitioning, Ordering and Coloring
!!
!! Copyright 2012 NTESS and the Zoltan contributors.
!! SPDX-License-Identifier: BSD-3-Clause
!! *****************************************************************************
!! @HEADER
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                                          //
!! File:      driver.cc                                                     //
!! Project:   Local HSFC Ordering                                           //
!! Author:    Michael Wolf                                                  //
!! Date Started:  11/02/2009                                                //
!!                                                                          //
!! Description:                                                             //
!!              File tests local HSFC ordering for simple test problem      //
!!                                                                          //
!! $Id: driver.cc 11 2009-11-10 00:15:18Z mmwolf $                          //
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program exampleRCB
   use mpi_h
   use zoltan
   use zoltanRCBex

   implicit none

   integer(Zoltan_INT) :: error
   real(Zoltan_FLOAT) :: version

!!!  numGlobObjs, numLocObjs, GIDs, xcoords, ycoords defined in zoltanRCBex module

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   interface

     subroutine readInputObjects(fname,numGlobObjs,numLocObjs,GIDs,xs,ys)
        character (len=*) :: fname
        integer :: numGlobObjs, numLocObjs
        integer, allocatable :: GIDs(:)  
        real, allocatable :: xs(:), ys(:)  
     end subroutine readInputObjects

   end interface
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Body of program
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   call MPI_Init(error)

   error = Zoltan_Initialize(version)

   call readInputObjects("mesh.txt",numGlobObjs,numLocObjs,GIDs,xcoords,ycoords)

   call partitionMeshWithRCB()

   call visualizePartition()

   deallocate(GIDs)
   deallocate(xcoords,ycoords)

   !! function in zoltanRCBex module that cleans up Zoltan data structures
   call zoltanCleanUp()

   call MPI_Finalize(error)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end program exampleRCB
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine readInputObjects(fname, numGlobObjs, numLocObjs,GIDs,xcoords,ycoords)
  use mpi_h
  implicit none

  character (len=*) :: fname

  integer :: numGlobObjs, numLocObjs
  integer, allocatable :: GIDs(:)  
  real, allocatable :: xcoords(:), ycoords(:)  

  ! Local declarations
  integer :: fnum = 2, i, currIndx
  integer :: myRank, numProcs, mpi_ierr
  integer :: tmpGID
  real :: tmpX, tmpY


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  call MPI_Comm_rank(MPI_COMM_WORLD, myrank, mpi_ierr)
  call MPI_Comm_size(MPI_COMM_WORLD, numProcs, mpi_ierr)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Scan data to determine # local objects
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  open(unit = fnum, file = fname)
  read(fnum,*) numGlobObjs

  numLocObjs=0

  do i = 1, numGlobObjs, 1
    read(fnum,*) tmpGID, tmpX, tmpY

    !! assumes gids start at 1, gives round robin initial distribution
    if ( MOD(tmpGID-1,numProcs) == myRank) then
      numLocObjs = numLocObjs + 1
    end if
  end do

  close(fnum)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Allocate data for my part of mesh
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  allocate(GIDs(numLocObjs))
  allocate(xcoords(numLocObjs))
  allocate(ycoords(numLocObjs))
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Fill data for my part of mesh
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  open(unit = fnum, file = fname)

  read(fnum,*) numGlobObjs

  currIndx = 1
  do i = 1, numGlobObjs, 1
    read(fnum,*) tmpGID, tmpX, tmpY

    !! assumes gids start at 1, gives round robin initial distribution
    if ( MOD(tmpGID-1,numProcs) == myRank) then
      GIDs(currIndx) = tmpGID
      xcoords(currIndx) = tmpX
      ycoords(currIndx) = tmpY
      currIndx = currIndx + 1
    end if
  end do

  close(fnum)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end subroutine readInputObjects
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Works for this specific 6x6 mesh
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine showSimpleMeshPartitions(myProc, numIDs, GIDs, parts)
  use mpi_h
  implicit none

  integer :: myProc, numIDs
  integer :: GIDs(*), parts(*)  

  !! Local variables
  integer :: partAssign(36)
  integer :: allPartAssign(36)
  integer :: i, j, part, mpi_ierr
  integer :: partRow(6)

  data partAssign/ 36 * 0/

  do i = 1, numIDs, 1
    partAssign(GIDs(i)) = parts(i);
  end do

  call MPI_Reduce(partAssign, allPartAssign, 36, MPI_INTEGER, MPI_MAX, 0, MPI_COMM_WORLD, mpi_ierr);

  if (myProc == 0) then
    do i=30, 0, -6
      do j = 1, 6, 1
        partRow(j) = allPartAssign(i+j);
      end do

      write(*,'(I1,A,I1,A,I1,A,I1,A,I1,A,I1)') partRow(1), '-----', partRow(2), '-----', partRow(3), '-----', &
                 partRow(4), '-----', partRow(5), '-----', partRow(6)

      if (i > 0) then
        write(*,'(A)') '|     |     |     |     |     |'
      end if
    end do
    write(*,*) 
  end if

end subroutine showSimpleMeshPartitions
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine visualizePartition()
  use mpi_h
  use zoltanRCBex
  implicit none

  !! Variables defined in zoltanRCBex module:
  !!      numLocObjs, GIDs, numExport, exportLocalGids, exportToPart 

  !! local variables
  integer :: parts(numLocObjs)
  integer :: myrank, i, error

  call MPI_Comm_rank(MPI_COMM_WORLD, myrank, error)


  do i=1, numLocObjs, 1
    parts(i) = myRank;
  end do

  if (myRank== 0) then
    write (*,*) 'Mesh part assignments before calling Zoltan'
  end if

  call showSimpleMeshPartitions(myRank, numLocObjs, GIDs, parts);

  do i=1, numExport, 1
    parts(exportLocalGids(i)) = exportToPart(i)
  end do


  if (myRank == 0) then
    write (*,*) 'Mesh part assignments after calling Zoltan'
  end if

  call showSimpleMeshPartitions(myRank, numLocObjs, GIDs, parts)

end subroutine visualizePartition
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
