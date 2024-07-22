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
! ***************************************************************************
!
!   Code imported to Zoltan zdrive from
!
!   zoltanParams_read_file.c
!
!   Read Zoltan parameters from a file, call Zoltan to set the parameters
!
!   zoltanParams library
!
!   Jim Teresco
!
!   Department of Computer Science
!   Williams College
!
!   and
!
!   Computer Science Research Institute
!   Sandia National Laboratories
!
!

!  Translated to Fortran by Bill Mitchell, NIST, March 2007

module dr_param_file

!#include <stdio.h>
!#include <mpi.h>
!#include <stdlib.h>
!#include "zoltan.h"

use zoltan
use mpi_h

implicit none
private

! Are other routines public?  If so, add them separated by commas. Use
! an amperstand at the end of the line to continue on another line.

public ztnPrm_read_file

! standard error and standard out are not standardized in Fortran 90, but most
! compilers use units 0 (most) and 6 (universal).  Change here if necessary.

integer, parameter :: stderr = 0, stdout = 6

!#define DEBUG 1

logical, parameter :: DEBUG = .false.

!struct zoltanParams_list_entry {
!  char *param;
!  char *value;
!  struct zoltanParams_list_entry *next;
!};

! unknown length character strings are tricky in Fortran, just pick a
! length that is hopefully big enough.  Change it here if necessary.

integer, parameter :: MAX_CHAR_LEN = 128

type ztnPrm_list_entry
   character(len=MAX_CHAR_LEN) :: param
   character(len=MAX_CHAR_LEN) :: value
   type (ztnPrm_list_entry), pointer :: next
end type ztnPrm_list_entry

!struct zoltanParams_hier_struct {
!  int partition;
!  struct zoltanParams_list_entry *first;
!};

type ztnPrm_hier_struct
   integer :: partition
   type (ztnPrm_list_entry), pointer :: first
end type ztnPrm_hier_struct

!static struct zoltanParams_hier_struct **zph = NULL;
!static int num_levels = 0;
!static MPI_Comm comm;

! I think the use of **zph is as an allocatable array

type(ztnPrm_hier_struct), save, allocatable :: zph(:)
integer :: num_levels = 0
integer :: comm

contains


!static void check_level(int level) {
!
!  if (!zph) {
!    fprintf(stderr,"check_level: must set number of levels first\n");
!    return;
!  }
!
!  if (level >= num_levels) {
!    fprintf(stderr,"check_level: invalid level\n");
!  }
!}

subroutine check_level(level)
integer :: level

  if (.not. allocated(zph)) then
    write(stderr,*) "check_level: must set number of levels first"
    return
  endif

  if (level >= num_levels) then
    write(stderr,*) "check_level: invalid level"
  endif
end subroutine check_level

!void zoltanParams_hier_free() {
!  int i;
!
!  if (!zph) {
!    fprintf(stderr, "zoltanParams_hier_free warning: not allocated\n");
!    return;
!  }
!
!  for (i=0; i<num_levels; i++) {
!    free(zph[i]);
!  }
!
!  free(zph);
!}

subroutine ztnPrm_hier_free()
  integer :: i

  if (.not. allocated(zph)) then
    write(stderr,*) "ztnPrm_hier_free warning: not allocated"
    return
  endif

  deallocate(zph)
end subroutine ztnPrm_hier_free

!void zoltanParams_hier_set_num_levels(int levels) {
!  int i;

subroutine ztnPrm_hier_set_num_levels(levels)
  integer :: levels
  integer :: i, astat

!#ifdef DEBUG
!  printf("(zoltanParams_hier_set_num_levels) setting to %d\n", levels);  
!#endif

  if (DEBUG) then
    write(stdout,*) "(ztnPrm_hier_set_num_levels) setting to ",levels
  endif

!  if (zph) {
!    fprintf(stderr,"zoltanParams_hier_set_num_levels warning: already initialized, reinitializing\n");
!    zoltanParams_hier_free();
!  }

  if (allocated(zph)) then
    write(stderr,*) "ztnPrm_hier_set_num_levels warning: already initialized, reinitializing"
    call ztnPrm_hier_free()
  endif

!  if (levels <= 0) {
!    fprintf(stderr, "(zoltanParams_hier_set_num_levels) num levels must be positive\n");
!    return;
!  }

  if (levels <= 0) then
    write(stderr,*) "(ztnPrm_hier_set_num_levels) num levels must be positive"
    return
  endif

!  num_levels = levels;
!
!  SAFE_MALLOC(zph, struct zoltanParams_hier_struct **, 
!	      sizeof(struct zoltanParams_hier_struct *) * levels);
!
!  for (i=0; i<levels; i++) {
!    SAFE_MALLOC(zph[i],  struct zoltanParams_hier_struct *, 
!		sizeof (struct zoltanParams_hier_struct));
!    zph[i]->partition = 0;
!    zph[i]->first = NULL;
!  }

  num_levels = levels

  allocate(zph(0:levels-1),stat=astat)
  if (astat /= 0) then
    write(stderr,*) "allocation failed in ztnPrm_hier_set_num_level"
    stop
  endif

  do i=0,levels-1
    zph(i)%partition = 0
    nullify(zph(i)%first)
  end do

!}

end subroutine ztnPrm_hier_set_num_levels

!void zoltanParams_hier_set_partition(int level, int partition) {
!
!#ifdef DEBUG
!  int mypid;
!  MPI_Comm_rank(comm, &mypid);
!
!  printf("[%d] will compute partition %d at level %d\n", 
!	 mypid, partition, level);
!#endif
!
!  check_level(level);
!
!  zph[level]->partition = partition;
!}

subroutine ztnPrm_hier_set_partition(level,partition)
integer :: level, partition
integer :: mypid, ierr

  if (DEBUG) then
    call MPI_Comm_rank(comm,mypid,ierr)

    write(stdout,*) "[",mypid,"] will compute partition ",partition," at level ",level
  endif

  call check_level(level)

  zph(level)%partition = partition
end subroutine ztnPrm_hier_set_partition

!void zoltanParams_hier_set_param(int level, char *param, char *value) {
!  struct zoltanParams_list_entry *newparam, *nextparam;

subroutine ztnPrm_hier_set_param(level,param,value)
integer :: level
character(len=*) :: param, value
type(ztnPrm_list_entry), pointer :: newparam, nextparam
integer :: mypid, ierr, astat

!#ifdef DEBUG
!  int mypid;
!  MPI_Comm_rank(comm, &mypid);
!  printf("[%d] will set param <%s> to <%s> at level %d\n", 
!	 mypid, param, value, level); 
!#endif

  if (DEBUG) then
    call MPI_Comm_rank(comm,mypid,ierr)
    write(stdout,*) "[",mypid,"] will set param ",trim(param)," to ",trim(value)," at level ",level
  endif

!  check_level(level);
!
!  SAFE_MALLOC(newparam, struct zoltanParams_list_entry *,
!	      sizeof(struct zoltanParams_list_entry));

  call check_level(level)

  allocate(newparam,stat=astat)
  if (astat /= 0) then
    write(stderr,*) "allocation failed in ztnPrm_hier_set_param"
    stop
  endif

!  newparam->param = strdup(param);
!  newparam->value = strdup(value);
!  newparam->next = NULL;

  newparam%param = param
  newparam%value = value
  nullify(newparam%next)
  
!  if (!zph[level]->first) {
!    zph[level]->first = newparam;
!    return;
!  }

  if (.not. associated(zph(level)%first)) then
    zph(level)%first => newparam
    return
  endif

!  nextparam = zph[level]->first;
!  while (nextparam->next) nextparam=nextparam->next;
!  nextparam->next = newparam;    

  nextparam => zph(level)%first
  do while (associated(nextparam%next))
    nextparam => nextparam%next
  end do
  nextparam%next => newparam

!}

end subroutine ztnPrm_hier_set_param

!int zoltanParams_hier_get_num_levels() {
!
!  return num_levels;
!}

function ztnPrm_hier_get_num_levels()
integer :: ztnPrm_hier_get_num_levels

  ztnPrm_hier_get_num_levels = num_levels
end function ztnPrm_hier_get_num_levels

!int zoltanParams_hier_get_part(int level) {
!
!  check_level(level);
!
!  return zph[level]->partition;
!}

function ztnPrm_hier_get_part(level)
integer :: level
integer :: ztnPrm_hier_get_part

  call check_level(level)

  ztnPrm_hier_get_part = zph(level)%partition
end function ztnPrm_hier_get_part

!void zoltanParams_hier_use_params(int level, struct Zoltan_Struct *zz, int *ierr) {
!  struct zoltanParams_list_entry *nextparam;
!
!  *ierr = ZOLTAN_OK;
!  check_level(level);
!  
!  nextparam = zph[level]->first;
!
!  while (nextparam) {
!    *ierr = Zoltan_Set_Param(zz, nextparam->param, nextparam->value);
!    if (*ierr != ZOLTAN_OK) return;
!    nextparam = nextparam->next;
!  }
!  
!}

subroutine ztnPrm_hier_use_params(level,zz,ierr)
integer :: level
type(Zoltan_Struct), pointer :: zz
integer :: ierr
type(ztnPrm_list_entry), pointer :: nextparam

  ierr = ZOLTAN_OK
  call check_level(level)

  nextparam => zph(level)%first

  do while (associated(nextparam))
    ierr = Zoltan_Set_Param(zz, nextparam%param, nextparam%value)
    if (ierr /= ZOLTAN_OK) return
    nextparam => nextparam%next
  end do
end subroutine ztnPrm_hier_use_params

!static int get_num_levels(void *data, int *ierr) {
!
!  *ierr = ZOLTAN_OK;
!  return zoltanParams_hier_get_num_levels();
!}

function get_num_levels(data, ierr)
integer(Zoltan_INT), intent(in) :: data(*)
integer(Zoltan_INT), intent(out) :: ierr
integer(Zoltan_INT) :: get_num_levels

  ierr = ZOLTAN_OK
  get_num_levels = ztnPrm_hier_get_num_levels()
end function get_num_levels

!static int get_part(void *data, int level, int *ierr) {
!
!  *ierr = ZOLTAN_OK;
!
!  return ztnPrm_hier_get_part(level);
!}

function get_part(data, level, ierr)
integer(Zoltan_INT), intent(in) :: data(*)
integer(Zoltan_INT), intent(in) :: level
integer(Zoltan_INT), intent(out) :: ierr
integer(Zoltan_INT) :: get_part

  ierr = ZOLTAN_OK

  get_part = ztnPrm_hier_get_part(level)
end function get_part

!static void get_method(void *data, int level, struct Zoltan_Struct *zz,
!		       int *ierr) {
!
!  zoltanParams_hier_use_params(level, zz, ierr);
!}

subroutine get_method(data,level,azz,ierr)
integer(Zoltan_INT), intent(in) :: data(*)
integer(Zoltan_INT), intent(in) :: level
type(Zoltan_Struct), intent(in), target :: azz
integer(Zoltan_INT), intent(out) :: ierr
type(Zoltan_Struct), pointer :: zz

  zz => azz

  call ztnPrm_hier_use_params(level, zz, ierr)
end subroutine get_method

!void zoltanParams_set_comm(MPI_Comm thecomm) {
!
!   remember the comm passed in 
!  MPI_Comm_dup(thecomm, &comm);
!}

subroutine ztnPrm_set_comm(thecomm)
integer :: thecomm
integer :: ierr

! remember the comm passed in

  call MPI_Comm_dup(thecomm, comm, ierr)
end subroutine ztnPrm_set_comm

!void zoltanParams_hier_setup(struct Zoltan_Struct *zz) {
!
!   make sure the hierarchical balancing callbacks are in place 
!  if (Zoltan_Set_Fn(zz, ZOLTAN_HIER_NUM_LEVELS_FN_TYPE, 
!		    (void (*)()) get_num_levels, NULL) == ZOLTAN_FATAL) {
!    fprintf(stderr,"zoltanParams_hier_setup: set NUM_LEVELS callback failed\n");
!  }
!
!  if (Zoltan_Set_Fn(zz, ZOLTAN_HIER_PARTITION_FN_TYPE, 
!		    (void (*)()) get_part, NULL) == ZOLTAN_FATAL) {
!    fprintf(stderr,"zoltanParams_hier_setup: set PARTITION callback failed\n");
!  }
!
!  if (Zoltan_Set_Fn(zz, ZOLTAN_HIER_METHOD_FN_TYPE, 
!		    (void (*)()) get_method, NULL) == ZOLTAN_FATAL) {
!    fprintf(stderr,"zoltanParams_hier_setup: set METHOD callback failed\n");
!  }   
!}

subroutine ztnPrm_hier_setup(zz)
type(Zoltan_Struct), pointer :: zz
integer(Zoltan_INT) :: dummy(1) = (/0/)

! make sure the hierarchical balancing callbacks are in place

  if (Zoltan_Set_Hier_Num_Levels_Fn(zz, get_num_levels, dummy) == &
      ZOLTAN_FATAL) then
    write(stderr,*) "ztnPrm_hier_setup: set NUM_LEVELS callback failed"
  endif

  if (Zoltan_Set_Hier_Part_Fn(zz, get_part, dummy) == &
      ZOLTAN_FATAL) then
    write(stderr,*) "ztnPrm_hier_setup: set PARTITION callback failed"
  endif

  if (Zoltan_Set_Hier_Method_Fn(zz, get_method, dummy) == &
      ZOLTAN_FATAL) then
    write(stderr,*) "ztnPrm_hier_setup: set METHOD callback failed"
  endif

end subroutine ztnPrm_hier_setup

!
!
!  zoltanParams_read_file
!
!  Set up the given Zoltan_Struct with parameters as specified
!  in the given file.
!
!  File format:
!
!  Lines of the format:
!  ZOLTAN_PARAM PARAM_VALUE
!
!  If the parameter is LB_METHOD set to HIER, the next part of the file
!  is interpreted as hierarchical balancing parameters:
!
!  num_levels
!  level 0 partitions for each proc
!  level 0 parameters
!  end with LEVEL END
!  level 1 partitions for each proc
!  level 1 parameters
!  end with LEVEL END
!  ...
!
!  End file with EOF
!
!

!void zoltanParams_read_file(struct Zoltan_Struct *lb, char *file, 
!			    MPI_Comm thecomm) {
!  FILE *fp;
!  char str1[500], str2[500];
!  int numlevels, level, partition, proc;
!  int ierr;
!  int mypid, numprocs;

subroutine ztnPrm_read_file(lb, file, thecomm)
type(Zoltan_Struct), pointer :: lb
character(len=*) :: file
integer :: thecomm
integer :: fp
character(len=500) :: str1, str2
integer :: numlevels, level, proc
integer :: ierr
integer :: mypid, numprocs
logical :: not2
integer, allocatable :: partition(:)

!   remember the comm passed in 
!  MPI_Comm_dup(thecomm, &comm);
!
!  MPI_Comm_rank(comm, &mypid);
!  MPI_Comm_size(comm, &numprocs);

! remember the comm passed in
  call MPI_Comm_dup(thecomm, comm, ierr)

  call MPI_Comm_rank(comm, mypid, ierr)
  call MPI_Comm_size(comm, numprocs, ierr)

!  fp = fopen(file, "r");
!  if (!fp) {
!    fprintf(stderr,"Cannot open file %s for reading", file);
!    return;
!  } 

! Assume unit 9 is available.  If it isn't, an error will be reported and
! you can change it to some other positive integer, not too big.

  fp = 9
  open(unit=fp,file=trim(file),action="read",iostat=ierr)
  if (ierr /= 0) then
    write(stderr,*) "cannot open file ",trim(file)," for reading"
    return
  endif

!#ifdef DEBUG
!  if (mypid == 0) {
!    printf("Reading Zoltan parameters from file %s\n", file);
!  }
!#endif

  if (DEBUG) then
    if (mypid == 0) then
      write(stdout,*) "Reading Zoltan parameters from file ",trim(file)
    endif
  endif

!  while (fscanf(fp, "%s %s\n", str1, str2) == 2) {

  do
    call myread(fp, str1, str2, not2)
    if (not2) exit

!    ierr = Zoltan_Set_Param(lb, str1, str2);
!    if (ierr != ZOLTAN_OK) {
!      fprintf(stderr,"Zoltan_Set_Param failed to set param <%s> to <%s>",str1,str2);
!    }
!#ifdef DEBUG
!    else {
!      if (mypid == 0) {
!	printf("Set Zoltan parameter <%s> to <%s>\n", str1, str2);
!      }
!    }
!#endif

    ! get rid of the leading space left on str2
    str2 = adjustl(str2)
    ierr = Zoltan_Set_Param(lb, trim(str1), trim(str2))
    if (ierr /= ZOLTAN_OK) then
      write(stderr,*) "Zoltan_Set_Param failed to set param ",trim(str1)," to ",trim(str2)
    endif
    if (DEBUG) then
      if (ierr == ZOLTAN_OK) then
        if (mypid == 0) then
          write(stdout,*) "Set Zoltan parameter ",trim(str1)," to ",trim(str2)
        endif
      endif
    endif

!    if (strcmp(str1,"LB_METHOD") == 0 && strcmp(str2,"HIER") == 0) {

    if (trim(str1) == "LB_METHOD" .and. trim(str2) == "HIER") then

!      zoltanParams_hier_setup(lb);

      call ztnPrm_hier_setup(lb)
 
!       the rest of the file contains hierarchical balancing parameters 
!      fscanf(fp, "%d", &numlevels);

! the rest of the file contains hierarchical balancing parameters
! The line containing numlevels is already in str1 (NO - it's next in the file)

      read(fp,*) numlevels

!#ifdef DEBUG
!      printf("[%d] read in numlevels=%d\n", mypid, numlevels);
!#endif

      if (DEBUG) then
        write(stdout,*) "[",mypid,"] read in numlevels=",numlevels
      endif

!      zoltanParams_hier_set_num_levels(numlevels);

      call ztnPrm_hier_set_num_levels(numlevels)

!      for (level=0; level<numlevels; level++) {
!	 first, a list of partitions for each proc should be in the file 
!	for (proc=0; proc<numprocs; proc++) {
!	  fscanf(fp, "%d", &partition);
!	  if (proc == mypid) zoltanParams_hier_set_partition(level, partition);
!	}

      allocate(partition(0:numprocs-1))
      ! probably should check that allocate succeeded
      do level=0,numlevels-1
         read(fp,*) partition ! assumes the line has exactly numprocs numbers
         call ztnPrm_hier_set_partition(level,partition(mypid))

!	 then parameters until we get LEVEL END 
!	while ((fscanf(fp, "%s %s\n", str1, str2) == 2) &&
!	       (strcmp(str1, "LEVEL") != 0) &&
!	       (strcmp(str2, "END") != 0)) {
!	  
!	  zoltanParams_hier_set_param(level, str1, str2);
!	}
!      }

! then parameters until we get LEVEL END
        do
          read(fp,*) str1, str2
          str2 = adjustl(str2)
          if (trim(str1) == "LEVEL" .and. trim(str2) == "END") exit
          call ztnPrm_hier_set_param(level, str1, str2)
        end do
      end do
      deallocate(partition)

!    }

    endif

!  }

  end do

!  fclose(fp);

  close(fp)

!}

end subroutine ztnPrm_read_file

! Fortran will generate an error if we try to read 2 strings and there is
! only 1 there.  So we have to read the whole line into a string and
! see if there are 1 or 2 strings in there.  Then read the individual
! strings and return them with a flag indicating if there are 1 or 2.

subroutine myread(runit,str1,str2,not2)
integer :: runit
character(len=*) :: str1, str2
logical :: not2
integer :: iostat

! assume 1000 is plenty long for an input line.

character(len=1000) :: line

! read the whole input line

  read(runit,"(A)",iostat=iostat) line

! end of file?
  if (iostat /= 0) then
     not2 = .true.
  else

! remove leading blanks

  line = adjustl(line)

! read the first string

  read(line,*) str1

! if the length of the whole line with leading and trailing blanks removed
! is the same as the length of the first string, then there is only 1 string

  if (len_trim(line) == len_trim(str1)) then
    not2 = .true.

! otherwise, read the second line

  else
    not2 = .false.
    read(line(len_trim(str1)+1:),"(A)") str2

  endif
  endif

end subroutine myread

end module dr_param_file
