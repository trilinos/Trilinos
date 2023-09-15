!! 
!! @HEADER
!!
!!!!**********************************************************************
!!
!!  Zoltan Toolkit for Load-balancing, Partitioning, Ordering and Coloring
!!                  Copyright 2012 Sandia Corporation
!!
!! Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
!! the U.S. Government retains certain rights in this software.
!!
!! Redistribution and use in source and binary forms, with or without
!! modification, are permitted provided that the following conditions are
!! met:
!!
!! 1. Redistributions of source code must retain the above copyright
!! notice, this list of conditions and the following disclaimer.
!!
!! 2. Redistributions in binary form must reproduce the above copyright
!! notice, this list of conditions and the following disclaimer in the
!! documentation and/or other materials provided with the distribution.
!!
!! 3. Neither the name of the Corporation nor the names of the
!! contributors may be used to endorse or promote products derived from
!! this software without specific prior written permission.
!!
!! THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
!! EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
!! IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
!! PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
!! CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
!! EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
!! PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
!! PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
!! LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
!! NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
!! SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!!
!! Questions? Contact Karen Devine	kddevin@sandia.gov
!!                    Erik Boman	egboman@sandia.gov
!!
!!!!**********************************************************************
!!
!! @HEADER
 !!

!--------------------------------------------------------------------------
! Purpose: Driver for dynamic load-balance library, ZOLTAN.                
!                                                                          
!--------------------------------------------------------------------------
! Author(s):  Matthew M. St.John (9226)                                    
!   Translated to Fortran by William F. Mitchell
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
! Revision History:                                                        
!                                                                          
!    30 March 1999:    Date of creation                                    
!       1 September 1999: Fortran translation
!--------------------------------------------------------------------------

!**************************************************************************
!**************************************************************************
!**************************************************************************


program fdriver
use zoltan
use mpi_h
use zoltan_user_data
use dr_const
use dr_input
use dr_chaco_io
use dr_loadbal
use dr_mm_io
use dr_sort
implicit none

! Local declarations. 
  character(len=64)  :: cmd_file

  real(Zoltan_FLOAT) :: version

  integer(Zoltan_INT) :: Proc, Num_Proc
  integer(Zoltan_INT) :: error, i

  type(PARIO_INFO) :: pio_info
  type(PROB_INFO) :: prob

  character(len=MPI_MAX_PROCESSOR_NAME) :: procname
  integer(Zoltan_INT) :: namelen
  integer(Zoltan_INT) :: alloc_stat

! KDDKDD
!  character(len=32) :: filename
!  character(len=10) :: mystring
!  integer :: fp
! KDDKDD
! interface blocks for external procedures

interface

   logical function read_mesh(Proc, Num_Proc, prob, pio_info)
   use zoltan
   use zoltan_user_data
   use dr_const
   use dr_input
   integer(Zoltan_INT) :: Proc, Num_Proc
   type(PROB_INFO) :: prob
   type(PARIO_INFO) :: pio_info
   end function read_mesh

   subroutine print_input_info(fp, Num_Proc, prob)
   use zoltan
   use dr_const
   integer(Zoltan_INT) :: fp
   integer(Zoltan_INT) :: Num_Proc
   type(PROB_INFO) :: prob
   end subroutine print_input_info

   logical function output_results(cmd_file, Proc, Num_Proc, prob, pio_info, &
                                   elements)
   use zoltan
   use dr_const
   use dr_input
   use zoltan_user_data
   character(len=*) :: cmd_file
   integer(Zoltan_INT) :: Proc, Num_Proc
   type(PROB_INFO) :: prob
   type(PARIO_INFO) :: pio_info
   type(ELEM_INFO), pointer :: elements(:)
   end function output_results

end interface

!**************************** BEGIN EXECUTION *****************************

!   initialize MPI 
  call MPI_Init(error)

!   get some machine information 
  call MPI_Comm_rank(zoltan_get_global_comm(), Proc, error)
  call MPI_Comm_size(zoltan_get_global_comm(), Num_Proc, error)

  call MPI_Get_processor_name(procname, namelen, error)
  print *,"Processor ",Proc," of ",Num_Proc," on host ",procname(1:namelen)

! Set the input file

  cmd_file = "zdrive.inp"

!   initialize Zoltan 
  error = Zoltan_Initialize(version)
  if (error /= ZOLTAN_OK) then
    print *, "fatal: Zoltan_Initialize returned error code, ", error
    goto 9999
  endif

!   initialize some variables 

  allocate(Mesh, stat=alloc_stat)
  if (alloc_stat /= 0) then
    print *, "fatal: insufficient memory"
    goto 9999
  endif

  nullify(Mesh%eb_names,Mesh%eb_ids,Mesh%eb_cnts,Mesh%eb_nnodes, &
               Mesh%eb_nattrs,Mesh%ecmap_id,Mesh%ecmap_cnt,Mesh%ecmap_elemids,&
               Mesh%ecmap_sideids,Mesh%ecmap_neighids,Mesh%elements, &
               Mesh%hgid,Mesh%hindex,Mesh%hvertex)
  Mesh%necmap = 0
  Mesh%nhedges = 0

  pio_info%init_dist_pins =  1  !INITIAL_LINEAR
  pio_info%dsk_list_cnt   = -1
  pio_info%num_dsk_ctrlrs = -1
  pio_info%pdsk_add_fact  = -1
  pio_info%zeros          = -1
  pio_info%file_type      = -1
  pio_info%pdsk_root      = ''
  pio_info%pdsk_subdir    = ''
  pio_info%pexo_fname     = ''

  prob%method             = ''
  prob%num_params         = 0
  prob%ztnPrm_file        = ''
  nullify(prob%params)

!   Read in the ascii input file 
  if(Proc == 0) then
    print *
    print *
    print *,"Reading the command file, ", cmd_file
    if(.not. read_cmd_file(cmd_file, prob, pio_info)) then
      print *, 'fatal: Could not read in the command file "',cmd_file,'"!'
      goto 9999
    endif

    if (.not. check_inp(prob, pio_info)) then
      print *, "fatal: Error in user specified parameters."
      goto 9999
    endif

    call print_input_info(6, Num_Proc, prob)
  endif

!   broadcast the command info to all of the processor 
  call brdcst_cmd_info(Proc, prob, pio_info)

!  
!   * now read in the mesh and element information.
!   * This is the only function call to do this. Upon return,
!   * the mesh struct and the elements array should be filled.
!   
  if (.not. read_mesh(Proc, Num_Proc, prob, pio_info)) then
      print *, "fatal: Error returned from read_mesh"
      goto 9999
  endif

! KDDKDD  TEMPORARY OUTPUT
!  if (Mesh%nhedges > 0) then
!    write (mystring, "(i1)" ) Proc
!    filename = "helpme."//mystring
!    open(unit=fp,file=filename,action="write")
!    write(fp,*) "Hyperedges:"
!    do i = 0, Mesh%nhedges-1
!      write(fp,*) "Edge ", Mesh%hgid(i)
!      do j = Mesh%hindex(i), Mesh%hindex(i+1)-1
!        write(fp, *) "        ", Mesh%hvertex(j)
!     enddo
!   enddo
!   close(unit=fp)
! endif
! KDDKDD  END TEMPORARY OUTPUT

!  
!   * now run zoltan to get a new load balance and perform
!   * the migration
!   
  if (.not. run_zoltan(Proc, prob, pio_info)) then
      print *, "fatal: Error returned from run_zoltan"
      goto 9999
  endif

!  
!   * output the results
!   
  if (.not. output_results(cmd_file, Proc, Num_Proc, prob, pio_info, Mesh%elements)) then
      print *, "fatal: Error returned from output_results"
      goto 9999
  endif

9999 continue
  if (associated(Mesh%elements)) then
    do i = 0, Mesh%elem_array_len-1
      call free_element_arrays(Mesh%elements(i))
    end do
    deallocate(Mesh%elements)
  endif
  if (associated(Mesh%hgid)) deallocate(Mesh%hgid)
  if (associated(Mesh%hindex)) deallocate(Mesh%hindex)
  if (associated(Mesh%hvertex)) deallocate(Mesh%hvertex)
  if (associated(Mesh)) deallocate(Mesh)
  if (associated(prob%params)) deallocate(prob%params)
  call Zoltan_Memory_Stats()
  call MPI_Finalize(error)

end program fdriver

!***************************************************************************
!***************************************************************************
!***************************************************************************
! This function determines which input file type is being used,
! * and calls the appropriate read function. If a new type of input
! * file is added to the driver, then a section needs to be added for
! * it here.
! *---------------------------------------------------------------------------
logical function read_mesh(Proc, Num_Proc, prob, pio_info)
use zoltan
use zoltan_user_data
use dr_const
use dr_input
use dr_chaco_io
use dr_mm_io
implicit none
  integer(Zoltan_INT) :: Proc
  integer(Zoltan_INT) :: Num_Proc
  type(PROB_INFO) :: prob
  type(PARIO_INFO) :: pio_info

! local declarations 
!-----------------------------Execution Begins------------------------------
  if (pio_info%file_type == CHACO_FILE) then
    if (.not. read_chaco_mesh(Proc, Num_Proc, prob, pio_info, Mesh%elements)) then
        print *, "fatal: Error returned from read_chaco_mesh"
        read_mesh = .false.
        return
    endif
  else if (pio_info%file_type == MM_FILE) then
    if (.not. read_mm_file(Proc, Num_Proc, prob, pio_info)) then
        print *, "fatal: Error returned from read_mm_file"
        read_mesh = .false.
        return
    endif
! not supporting NEMESIS yet
!  else if (pio_info->file_type == NEMESIS_FILE) {
!    if (!read_exoII_mesh(Proc, Num_Proc, prob, pio_info, elements)) {
!        Gen_Error(0, "fatal: Error returned from read_exoII_mesh\n");
!        return 0;
!    }
!  }
  else
    print *, "fatal: Input file type not supported."
    read_mesh = .false.
    return
  endif
  read_mesh = .true.
  return
end function read_mesh

!***************************************************************************
!***************************************************************************
subroutine print_input_info(fp, Num_Proc, prob)
use zoltan
use dr_const
implicit none
integer(Zoltan_INT) :: fp
integer(Zoltan_INT) :: Num_Proc
type(PROB_INFO) :: prob

integer :: i

  write(fp,*) "Input values:"
  write(fp,*) "  ",DRIVER_NAME," version ", VER_STR
  write(fp,*) "  Total number of Processors = ", Num_Proc
  write(fp,*)

  write(fp,*)
  write(fp,*) "  Performing load balance using ", prob%method
  write(fp,*) "  Parameters:"
  do i = 0, prob%num_params-1
    write(fp,*) "    ",trim(prob%params(i)%str(0)),"  ",trim(prob%params(i)%str(1))
  end do

  write(fp,*) "##########################################################"
end subroutine print_input_info

!************************************************************************

logical function output_results(cmd_file, Proc, Num_Proc, prob, pio_info, &
                                elements)
use zoltan
use dr_const
use dr_input
use zoltan_user_data
use dr_sort
character(len=*) :: cmd_file
integer(Zoltan_INT) :: Proc, Num_Proc
type(PROB_INFO) :: prob
type(PARIO_INFO) :: pio_info
type(ELEM_INFO), pointer :: elements(:)

!
! * For the first swipe at this, don't try to create a new
! * exodus/nemesis file or anything. Just get the global ids,
! * sort them, and print them to a new ascii file.
! 

!   Local declarations. 
  character(len=FILENAME_MAX+1) :: par_out_fname, ctemp

  integer(Zoltan_INT), allocatable :: global_ids(:), parts(:), index(:)
  integer(Zoltan_INT), allocatable :: orders(:), iperms(:)
  integer(Zoltan_INT) ::    i, j, alloc_stat

  integer ::  fp=21

  interface
   subroutine echo_cmd_file(fp, cmd_file)
   character(len=*) :: cmd_file
   integer :: fp
   end subroutine echo_cmd_file

  end interface

!**************************** BEGIN EXECUTION *****************************

  allocate(global_ids(0:Mesh%num_elems),stat=alloc_stat)
  if (alloc_stat /= 0) then
    print *, "fatal: insufficient memory"
    output_results = .false.
    return
  endif
  allocate(parts(0:Mesh%num_elems),stat=alloc_stat)
  if (alloc_stat /= 0) then
    print *, "fatal: insufficient memory"
    output_results = .false.
    return
  endif
  allocate(index(0:Mesh%num_elems),stat=alloc_stat)
  if (alloc_stat /= 0) then
    print *, "fatal: insufficient memory"
    output_results = .false.
    return
  endif
  allocate(orders(0:Mesh%num_elems),stat=alloc_stat)
  if (alloc_stat /= 0) then
    print *, "fatal: insufficient memory"
    output_results = .false.
    return
  endif
  allocate(iperms(0:Mesh%num_elems),stat=alloc_stat)
  if (alloc_stat /= 0) then
    print *, "fatal: insufficient memory"
    output_results = .false.
    return
  endif

  j = 0
  do i = 0, Mesh%elem_array_len-1
    if (elements(i)%globalID >= 0) then
      global_ids(j) = elements(i)%globalID
      parts(j) = elements(i)%my_part
      orders(j) = elements(i)%perm_value;
      iperms(j) = elements(i)%invperm_value;
      index(j) = j
      j = j+1
    endif
  end do

  call dr_sort_index(0, Mesh%num_elems-1, global_ids, index)

!   generate the parallel filename for this processor 
  ctemp = pio_info%pexo_fname(1:len_trim(pio_info%pexo_fname))//".out"
  call gen_par_filename(ctemp, par_out_fname, pio_info, Proc, Num_Proc)

  open(unit=fp,file=par_out_fname,action="write")
  if (Proc == 0) then
    call echo_cmd_file(fp, cmd_file)
  endif

  write(fp,*) "Global element ids assigned to processor ", Proc
  write(fp,*) "GID	Part	Perm	IPerm"
  do i = 0, Mesh%num_elems-1
    j = index(i)
    write(fp,*) global_ids(j),"	", parts(j), "	", orders(j), "	", iperms(j) 
  end do

  close(fp)
  deallocate(global_ids)
  deallocate(parts)
  deallocate(index)
  deallocate(orders)
  deallocate(iperms)

  output_results = .true.
end function output_results

!************************************************************************
subroutine echo_cmd_file(fp, cmd_file)
character(len=*) :: cmd_file
integer :: fp
integer, parameter :: file_cmd = 11
character(len=4096+1) :: inp_line

! Routine to echo the input file into the output results (so that
! we know what conditions were used to produce a given result).


!   Open the file 
  open(unit=file_cmd,file=cmd_file,action='read',iostat=iostat)
  if (iostat /= 0) then
    print *, "Error:  Could not find command file ", cmd_file
    return
  endif

  do
    read(unit=file_cmd,fmt="(a)",iostat=iostat) inp_line
    if (iostat /= 0) exit ! end of data
    write(fp, *) trim(inp_line)
  end do

  close(file_cmd)
end subroutine echo_cmd_file
