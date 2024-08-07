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
module dr_input
use zoltan
use dr_const
use mpi_h
implicit none
private

public :: read_cmd_file, check_inp, brdcst_cmd_info, gen_par_filename, &
          PARIO_INFO, NEMESIS_FILE, CHACO_FILE, MM_FILE

!--------------------------------------------------------------------------
! Purpose: Determine file types for command files and read in the parallel 
!          ExodusII command file.                                          
!          Taken from nemesis utilites nem_spread and nem_join.            
!--------------------------------------------------------------------------
! Author(s):  Matthew M. St.John (9226)                                    
!   Translated to Fortran by William F. Mitchell
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
! Revision History:                                                        
!    14 April 1999:       Date of creation.                                
!      02 September 1999:   Fortran translation
!--------------------------------------------------------------------------


integer(Zoltan_INT), parameter :: NEMESIS_FILE = 0, &
                             CHACO_FILE   = 1, &
                             MM_FILE = 2

integer(Zoltan_INT), parameter :: MAX_INPUT_STR_LN = 4096 ! maximum string length
                                                     ! for read_string()

! Structure used to store the information necessary for parallel I/O. 

type PARIO_INFO
   
  integer(Zoltan_INT) :: init_dist_pins
  integer(Zoltan_INT) :: dsk_list_cnt
  integer(Zoltan_INT), pointer :: dsk_list(:)
  integer(Zoltan_INT) :: rdisk
  integer(Zoltan_INT) :: num_dsk_ctrlrs   ! The number of disk controllers.     
  integer(Zoltan_INT) :: pdsk_add_fact    ! The offset from zero used by the    
                                     ! the target machine.                 
  integer(Zoltan_INT) :: zeros
                       ! 1 - if the target machine uses leading zeros when 
                       !     designating the disk number (eg - the paragon 
                       !     uses /pfs/io_01)                              
                       ! 0 - if it does not (eg - the tflop uses           
                       !     /pfs/tmp_1)                                   

  integer(Zoltan_INT) :: file_type        ! input file type 

  ! The root location of the parallel disks 
  character(len=FILENAME_MAX+1) :: pdsk_root

  ! The subdirectory to write files to 
  character(len=FILENAME_MAX+1) :: pdsk_subdir

  ! The base name of the input file. 
  character(len=FILENAME_MAX+1) :: pexo_fname

end type PARIO_INFO

contains

function lowercase(string)
character(len=*), intent(in) :: string
character(len=len(string)) :: lowercase

! returns the string converted to lower case

integer, parameter :: int_A = iachar("A"), &
                      int_Z = iachar("Z"), &
                      a_A_diff = iachar("a") - iachar("A")
integer :: i, int_char

do i=1,len(string)
   int_char = iachar(string(i:i))
   if (int_char >= int_A .and. int_char <= int_Z) then
      lowercase(i:i) = achar(int_char + a_A_diff)
   else
      lowercase(i:i) = string(i:i)
   endif
end do

end function lowercase

!***************************************************************************
!***************************************************************************
logical function read_cmd_file(filename, prob, pio_info)
character(len=*) :: filename
type(PROB_INFO) :: prob
type(PARIO_INFO) :: pio_info

!
! *          This function reads the ASCII parallel-exodus command file.
! *
! *   Input
! *   -----
! *   filename - The name of the command file.
! *   pio_info - parallel I/O information.
! 

! I'm really coping out here.  This was written for debugging and initial
! testing of the Fortran interface, and doesn't need the full capability
! of reading the command files.  So currently it assumes the command file
! looks like the Chaco test files that existed at the time it was written.
! WFM 9/1/99

! local declarations 
  integer, parameter :: file_cmd = 11
  character(len=MAX_INPUT_STR_LN + 1) :: inp_line, command, temp_string
  integer :: iostat
  logical :: more_params

!**************************** BEGIN EXECUTION *****************************

!   Open the file 
  open(unit=file_cmd,file=filename,action='read',iostat=iostat)
  if (iostat /= 0) then
    read_cmd_file = .false.
    return
  endif

! assume no more than 15 parameters 
  allocate(prob%params(0:15))
  prob%num_params = 1
  prob%params(0)%str(0) = "DEBUG_MEMORY"
  prob%params(0)%str(1) = "1"

!   Begin parsing the input file 
  do ! while not end of data
    read(unit=file_cmd,fmt="(a)",iostat=iostat) inp_line
    if (iostat /= 0) exit ! end of data

!     skip any line that is a comment 
    if (inp_line == '') cycle
    if (inp_line(1:1) == '#') cycle

! find what is before the equal sign

    command = inp_line(1:index(inp_line,"=")-1)

! if there is a tab, take what is before the tab

    if (index(command,"	") /= 0) then ! "tab"
       command = command(1:index(command,"	")-1) ! "tab"
    endif

! ****** File Name ******

    if (lowercase(trim(command)) == "file name") then
! assumes there is one blank between "=" and the file name
       pio_info%pexo_fname = trim(inp_line(index(inp_line,"=")+2:))
    endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! hacks to allow more of input file to be read (KDD, 10/2000) !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if (lowercase(trim(command)) == "decomposition method") then
! assumes there is one blank between "=" and the method
       prob%method = trim(inp_line(index(inp_line,"=")+2:))
    endif

    if (lowercase(trim(command)) == "file type") then
! assumes there is one blank between "=" and the file type
       if (lowercase(trim(inp_line(index(inp_line,"=")+2:))) == "chaco") then
          pio_info%file_type = CHACO_FILE
       else if (lowercase(trim(inp_line(index(inp_line,"=")+2:))) == "matrixmarket") then
          pio_info%file_type = MM_FILE
       else
          print *, "Error:  zfdrive can read only Chaco or MatrixMarket format files."  
          read_cmd_file = .false.
       endif
    endif

    if ((lowercase(trim(command)) == "zoltan parameters").or. &
        (lowercase(trim(command)) == "zoltan parameter")) then
! assumes there is one blank between "=" and the parameter name
       temp_string = lowercase(trim(inp_line(index(inp_line,"=")+2:)))
! assumes no blanks between second "=" and the parameter name
! skip the input line if there are no parameters specified on it.
       if (index(temp_string,"=").gt.0) then
          more_params = .true.
          do while (more_params)
             if (index(temp_string, ",").gt.0) then
                more_params = .true.
                prob%params(prob%num_params)%str(1) = &
                    temp_string(index(temp_string,"=")+1:index(temp_string,",")-1)
             else
                more_params = .false.
                prob%params(prob%num_params)%str(1) = &
                    temp_string(index(temp_string,"=")+1:)
             endif
             prob%params(prob%num_params)%str(0) = &
                 temp_string(1:index(temp_string,"=")-1)
             prob%num_params = prob%num_params+1
             if (more_params) then
                temp_string = temp_string(index(temp_string,",")+1:)
                do while (temp_string(1:1).eq." ")  !skip white space
                   temp_string = temp_string(2:)
                enddo
             endif
          enddo
       endif
    endif

    if (lowercase(trim(command)) == "test multi callbacks") then
!       assumes there is one blank between "=" and the input value
        Test_Multi_Callbacks = iachar(trim(inp_line(index(inp_line,"= ")+2:))) - iachar('0')
    endif

    if (lowercase(trim(command)) == "test graph callbacks") then
        Test_Graph_Callbacks = iachar(trim(inp_line(index(inp_line,"= ")+2:))) - iachar('0')
    endif

    if (lowercase(trim(command)) == "test hypergraph callbacks") then
        Test_Hypergraph_Callbacks = iachar(trim(inp_line(index(inp_line,"= ")+2:))) - iachar('0')
    endif

    if (lowercase(trim(command)) == "test local partitions") then
!       assumes there is one blank between "=" and the input value
        Test_Local_Partitions = iachar(trim(inp_line(index(inp_line,"=")+2:))) - iachar('0')
    endif

    if (lowercase(trim(command)) == "test generate files") then
!       assumes there is one blank between "=" and the input value
        Test_Gen_Files = iachar(trim(inp_line(index(inp_line,"= ")+2:))) - iachar('0')
    endif

    if (lowercase(trim(command)) == "test drops") then
!       assumes there is one blank between "=" and the input value
        Test_Drops = iachar(trim(inp_line(index(inp_line,"= ")+2:))) - iachar('0')
    endif

    if (lowercase(trim(command)) == "zdrive action") then
!       assumes there is one blank between "=" and the input value
        Driver_Action = iachar(trim(inp_line(index(inp_line,"= ")+2:))) - iachar('0')
    endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! end of hacks to allow more of input file to be read (KDD, 10/2000) !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Added JDT 3/2007 to save the zoltanparams file name
! Note: we assume a good format here:
! zoltanparams file = thefile.name
    if (lowercase(trim(command)) == "zoltanparams file") then
       prob%ztnPrm_file = lowercase(trim(inp_line(index(inp_line,"=")+2:)))
    endif
       
!

! The other commands are not processed.  In the initial tests they either
! always have the same value (in which case they are set after this loop)
! or they do not appear.

  end do ! while not end of data

! Assume parallel disk info is number=0

  pio_info%num_dsk_ctrlrs = 0

!   Close the command file 
  close(file_cmd)

  read_cmd_file = .true.

end function read_cmd_file

!***************************************************************************
!***************************************************************************
!***************************************************************************
logical function check_inp(prob, pio_info)
type(PROB_INFO) :: prob
type(PARIO_INFO) :: pio_info

!**************************** BEGIN EXECUTION *****************************

!   check for the parallel Nemesis file for proc 0 
  if (len_trim(pio_info%pexo_fname) <= 0) then
    print *, "fatal: must specify file base name"
    check_inp = .false.
    return
  endif

! Not supporting NEMESIS
!   default file type is nemesis 
!  if (pio_info->file_type < 0) pio_info->file_type = NEMESIS_FILE;
!
!#ifndef ZOLTAN_NEMESIS
!   
!   * if compiling without the ZOLTAN_NEMESIS flag (i.e., not linking with 
!   * Nemesis library), can't use NEMESIS_FILE file type.
!   
!
!  if (pio_info->file_type == NEMESIS_FILE) {
!    Gen_Error(0, "fatal: must link with Nemesis libraries for Nemesis "
!                 "file types");
!    return 0;
!  }
!#endif  !ZOLTAN_NEMESIS 

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                 Check the parallel IO specifications                      
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!   check that there is a list of disks, or a number of raids 
  if ((pio_info%dsk_list_cnt <= 0) .and. (pio_info%num_dsk_ctrlrs < 0)) then
    pio_info%num_dsk_ctrlrs = 0 !  default to single directory 
  endif

!   default is not to have preceeding 0's in the disk names 
  if (pio_info%zeros < 0) pio_info%zeros = 0

!   most systems that we deal with start their files systems with 1 not 0 
  if (pio_info%pdsk_add_fact < 0) pio_info%pdsk_add_fact = 1

!  
!   * if there are parallel disks, then the root and subdir locations must
!   * be specified
!   
  if (pio_info%num_dsk_ctrlrs > 0 .or. pio_info%dsk_list_cnt > 0) then
    if (len_trim(pio_info%pdsk_root) == 0) then
      print *, "fatal: must specify parallel disk root name"
      check_inp = .false.
      return
    endif
    if (len_trim(pio_info%pdsk_subdir) == 0) then
      print *, "fatal: must specify parallel disk subdirectory"
      check_inp = .false.
      return
    endif
  else
    if (len_trim(pio_info%pdsk_root) == 0) then
      pio_info%pdsk_root = "." !  default is execution directory 
    endif
  endif

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                 Check the Zoltan specifications                           
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!  
!   * Make sure a load-balancing method was provided.
!   
  if (len_trim(prob%method) == 0) then
    print *, "fatal: load balance method must be specified"
    check_inp = .false.
    return
  endif

  check_inp = .true.
end function check_inp

!***************************************************************************
!***************************************************************************
!***************************************************************************
subroutine brdcst_cmd_info(Proc, prob, pio_info)
integer(Zoltan_INT) :: Proc
type(PROB_INFO) :: prob
type(PARIO_INFO) :: pio_info

! local declarations 
  integer(Zoltan_INT) :: ctrl_id
  integer(Zoltan_INT) :: size
  integer(Zoltan_INT) :: ierr, i
!**************************** BEGIN EXECUTION *****************************

  call MPI_Bcast(Test_Multi_Callbacks, 1, MPI_INTEGER, 0, zoltan_get_global_comm(), ierr)
  call MPI_Bcast(Test_Local_Partitions, 1, MPI_INTEGER, 0, zoltan_get_global_comm(), ierr)
  call MPI_Bcast(Test_Drops, 1, MPI_INTEGER, 0, zoltan_get_global_comm(), ierr)
  call MPI_Bcast(Test_Gen_Files, 1, MPI_INTEGER, 0, zoltan_get_global_comm(), ierr)
  call MPI_Bcast(Driver_Action, 1, MPI_INTEGER, 0, zoltan_get_global_comm(), ierr)
  call MPI_Bcast(pio_info%dsk_list_cnt, 1, MPI_INTEGER, 0, zoltan_get_global_comm(), ierr)
  call MPI_Bcast(pio_info%rdisk, 1, MPI_INTEGER, 0, zoltan_get_global_comm(), ierr)
  call MPI_Bcast(pio_info%num_dsk_ctrlrs, 1, MPI_INTEGER, 0, zoltan_get_global_comm(), ierr)
  call MPI_Bcast(pio_info%pdsk_add_fact, 1, MPI_INTEGER, 0, zoltan_get_global_comm(), ierr)
  call MPI_Bcast(pio_info%zeros, 1, MPI_INTEGER, 0, zoltan_get_global_comm(), ierr)
  call MPI_Bcast(pio_info%file_type, 1, MPI_INTEGER, 0, zoltan_get_global_comm(), ierr)
  call MPI_Bcast(pio_info%pdsk_root, len(pio_info%pdsk_root), MPI_CHARACTER, &
            0, zoltan_get_global_comm(), ierr)
  call MPI_Bcast(pio_info%pdsk_subdir, len(pio_info%pdsk_root), MPI_CHARACTER, &
            0, zoltan_get_global_comm(), ierr)
  call MPI_Bcast(pio_info%pexo_fname, len(pio_info%pdsk_root), MPI_CHARACTER, &
            0, zoltan_get_global_comm(), ierr)

  if(pio_info%dsk_list_cnt > 0) then
    if(Proc /= 0) then
      allocate(pio_info%dsk_list(0:pio_info%dsk_list_cnt-1))
    endif

    call MPI_Bcast(pio_info%dsk_list, pio_info%dsk_list_cnt, MPI_INTEGER, &
                   0, zoltan_get_global_comm(), ierr)
  endif

!   broadcast the param file name 
  call MPI_Bcast(prob%ztnPrm_file, len(prob%ztnPrm_file), MPI_CHARACTER, &
       0, zoltan_get_global_comm(), ierr)

!   and broadcast the problem specifications 
  call MPI_Bcast(prob%method, len(prob%method), MPI_CHARACTER, 0,zoltan_get_global_comm(),ierr)
  call MPI_Bcast(prob%num_params, 1, MPI_INTEGER, 0, zoltan_get_global_comm(), ierr)
  if (prob%num_params > 0) then
    size = len(prob%params(0)%str(0))
    if (Proc /= 0) then
      allocate(prob%params(0:prob%num_params-1))
    endif
    do i=0,prob%num_params-1
      call MPI_Bcast(prob%params(i)%str(0), size, MPI_CHARACTER, 0, &
                     zoltan_get_global_comm(), ierr)
      call MPI_Bcast(prob%params(i)%str(1), size, MPI_CHARACTER, 0, &
                     zoltan_get_global_comm(), ierr)
    end do
  endif

!   now calculate where the file for this processor is 
  if(pio_info%dsk_list_cnt <= 0) then
    if (pio_info%num_dsk_ctrlrs > 0) then
      ctrl_id = mod(Proc,pio_info%num_dsk_ctrlrs)
      pio_info%rdisk = ctrl_id + pio_info%pdsk_add_fact
    endif
  else
    ctrl_id = mod(Proc,pio_info%dsk_list_cnt)
    pio_info%rdisk = pio_info%dsk_list(ctrl_id)
  endif

end subroutine brdcst_cmd_info

!***************************************************************************
!***************************************************************************
!***************************************************************************
subroutine gen_par_filename(scalar_fname, par_fname, pio_info, proc_for, nprocs)
character(len=*) :: scalar_fname, par_fname
type(PARIO_INFO) :: pio_info
integer(Zoltan_INT) :: proc_for, nprocs

!----------------------------------------------------------------------------
! *
! *      Author(s):     Gary Hennigan (1421)
!        Translated to Fortran by William F. Mitchell
! *----------------------------------------------------------------------------
! *      Function which generates the name of a parallel file for a
! *      particular processor. The function does this by appending
! *      "N.p" to the end of the input parameter "scalar_fname", where:
! *
! *              N - The number of processors utilized
! *              p - The processor ID.
! *
! *      In addition, the location of the parallel disk system is prepended
! *      to each file name.
! *---------------------------------------------------------------------------
! *      Example:
! *
! *        scalar_fname = "Parallel-exoII-"   (Input)
! *        par_fname    = "/raid/io_01/tmp/rf_crew/Parallel-exoII-8.0" (Output)
! *
! *      where, for this example:
! *
! *              N = 8 processors
! *              p = 0 particular processor ID
! *---------------------------------------------------------------------------
! *      Revision History:
! *
! *              05 November 1993:    Date of Creation
!                02 September 1999:   Fortran translation
! *---------------------------------------------------------------------------
! 


!        Local variables      

  integer(Zoltan_INT) :: iTemp1
  integer(Zoltan_INT) :: iMaxDigit, iMyDigit
  character(len=FILENAME_MAX) :: cTemp
  character(len=6) :: frmat
  character(len=32) :: nproc_str, myproc_str
  character(len=2) :: rdisk_str

!************************ EXECUTION BEGINS ******************************

!  
!   * Find out the number of digits needed to specify the processor ID.
!   * This allows numbers like 01-99, i.e., prepending zeros to the
!   * name to preserve proper alphabetic sorting of the files.
!   

  iMaxDigit = 0
  iTemp1 = nprocs
  do while (iTemp1 >= 1)
    iTemp1 = iTemp1/10
    iMaxDigit = iMaxDigit + 1
  end do

  iMyDigit = 0
  iTemp1 = proc_for
  do while (iTemp1 >= 1)
    iTemp1 = iTemp1/10
    iMyDigit = iMyDigit + 1
  end do

! create the character strings containing the numbers

  frmat=""
  write(frmat,"(a2,i1,a1)") "(I",iMaxDigit,")"
  write(nproc_str,frmat) nprocs
  frmat=""
  write(frmat,"(a2,i1,a1,i1,a1)") "(I",iMaxDigit,".",iMaxDigit,")"
  write(myproc_str,frmat) proc_for

! create the filename with the digit suffixes

  cTemp = trim(scalar_fname)//"."//trim(nproc_str)//"."//trim(myproc_str)

!  
!   * Finally, generate the complete file specification for the parallel
!   * file used by this processor.
!   
  if (pio_info%num_dsk_ctrlrs > 0) then
    if(pio_info%zeros /= 0) then
      write(rdisk_str,"(I2.2)") pio_info%rdisk
    else
      write(rdisk_str,"(I2)") pio_info%rdisk
      rdisk_str = adjustl(rdisk_str)
    endif
    par_fname = trim(pio_info%pdsk_root)//trim(rdisk_str)//"/"//&
                trim(pio_info%pdsk_subdir)//trim(cTemp)
  else
    par_fname = trim(pio_info%pdsk_root)//"/"//trim(cTemp)
  endif

end subroutine gen_par_filename

end module dr_input
