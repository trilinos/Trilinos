!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Zoltan Library for Parallel Applications                                   !
! For more info, see the README file in the top-level Zoltan directory.      ! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  CVS File Information :
!     $RCSfile$
!     $Author$
!     $Date$
!     $Revision$
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!--------------------------------------------------------------------
! File dr_const.h translated to Fortran by William F. Mitchell
!
! Revision History:
!
!    1 September 1999:   Translated to Fortran
!--------------------------------------------------------------------

module dr_const
use zoltan
use zoltan_user_data
implicit none
private

public :: DRIVER_NAME, VER_STR, PROB_INFO, MESH_INFO, Mesh, &
          FILENAME_MAX, MAX_PARAMETER_LEN, Parameter_Pair

public :: Test_Multi_Callbacks
public :: Test_Local_Partitions
public :: Test_Drops
public :: Test_Gen_Files
public :: Driver_Action

!/*****************************************************************************
! *  Definitions for the Zoltan library driver program.
! *****************************************************************************/

character(len=7), parameter :: DRIVER_NAME = "zfdrive"
character(len=3), parameter :: VER_STR = "1.0"

! A global variable indicating whether list-based (multi) query functions
! should be registered.  Default is 0.
integer(Zoltan_INT) :: Test_Multi_Callbacks = 0
integer(Zoltan_INT) :: Test_Local_Partitions = 0
integer(Zoltan_INT) :: Test_Drops = 0
integer(Zoltan_INT) :: Test_Gen_Files = 0
integer(Zoltan_INT) :: Driver_Action = 1

!/* If it doesn't get defined in stdio.h then use this as a default */
integer(Zoltan_INT), parameter :: FILENAME_MAX = 1024

integer(Zoltan_INT), parameter :: MAX_NP_ELEM = 27     ! max nodes per element
integer(Zoltan_INT), parameter :: MAX_DIM     =  3     ! max number of dimensions
!integer(Zoltan_INT), parameter :: MAX_EB_NAME_LEN = 32 ! chars for element block
integer(Zoltan_INT), parameter :: MAX_PARAMETER_LEN = 128 ! chars in parameter

!/*
! * Structure used to describe an element. Each processor will
! * allocate an array of these structures.
!   Actual definition has been moved to zoltan_user_data so that User_Data
!   can use it.
! */
!type ELEM_INFO
!  integer(Zoltan_INT) :: border   ! set to 1 if this element is a border element
!  integer(Zoltan_INT) :: globalID ! Global ID of this element, the local ID is the
!                             ! position in the array of elements
!  integer(Zoltan_INT) :: elem_blk ! element block number which this element is in
!  real(Zoltan_FLOAT)  :: cpu_wgt  ! the computational weight associated with the elem
!  real(Zoltan_FLOAT)  :: mem_wgt  ! the memory weight associated with the elem
!  real(Zoltan_FLOAT), pointer ::  coord(:,:) ! array for the coordinates of the
!                                        ! element. For Nemesis meshes, nodal
!                                        ! coordinates are stored; for Chaco
!                                        ! graphs with geometry, one set of
!                                        ! coords is stored.
!  integer(Zoltan_INT), pointer :: connect(:) ! list of nodes that make up this
!                                        ! element, the node numbers in this
!                                        ! list are global and not local
!  integer(Zoltan_INT), pointer :: adj(:)  ! list of adjacent elements .
!                         ! For Nemesis input, the list is ordered by
!                         ! side number, to encode side-number info needed to
!                         ! rebuild communication maps.  Value -1 represents 
!                         ! sides with no neighboring element (e.g., along mesh
!                         ! boundaries).  Chaco doesn't have "sides," so the 
!                         ! ordering is irrelevent for Chaco input.
!  integer(Zoltan_INT), pointer :: adj_proc(:) ! list of processors for adjacent
!                                         ! elements
!  real(Zoltan_FLOAT), pointer :: edge_wgt(:)  ! edge weights for adjacent elements
!  integer(Zoltan_INT) :: nadj    ! number of entries in adj
!  integer(Zoltan_INT) :: adj_len ! allocated length of adj/adj_proc/edge_wgt arrays
!end type

!/*
! * structure for general mesh information
! */
!/* Structure used to store information about the mesh */
! Moved to zoltan_user_data.f90 for use in callback functions.
!type MESH_INFO
!  integer(Zoltan_INT) :: num_nodes     ! number of nodes on this processor
!  integer(Zoltan_INT) :: num_elems     ! number of elements on this processor
!  integer(Zoltan_INT) :: num_dims      ! number of dimensions for the mesh
!  integer(Zoltan_INT) :: num_el_blks   ! number of element blocks in the mesh
!  integer(Zoltan_INT) :: num_node_sets ! number of node sets in the mesh
!  integer(Zoltan_INT) :: num_side_sets ! number of side sets in the mesh
!  character(len=MAX_EB_NAME_LEN), pointer :: eb_names(:) ! element block element
!                                                         ! names
!  integer(Zoltan_INT), pointer :: eb_ids(:)    ! element block ids
!  integer(Zoltan_INT), pointer :: eb_cnts(:)   ! number of elements in each element
!                                          ! block
!  integer(Zoltan_INT), pointer :: eb_nnodes(:) ! number of nodes per element in each
!                                          ! element block
!                                          ! for Nemesis meshes, this value
!                                          ! depends on element type;
!                                          ! for Chaco graphs, only one "node"
!                                          ! per element.
!  integer(Zoltan_INT), pointer :: eb_nattrs(:) ! number of attributes per element in
!                                          ! each element block
!  integer(Zoltan_INT) :: elem_array_len ! length that the ELEM_INFO array is
!                                   ! allocated for. Need to know this when array
!                                   ! is not completely filled during migration
!  integer(Zoltan_INT) :: necmap         ! number of elemental communication maps.
!  integer(Zoltan_INT), pointer :: ecmap_id(:)       ! IDs of each elemental
!                                               ! communication map.
!  integer(Zoltan_INT), pointer :: ecmap_cnt(:)      ! number of elements in each
!                                               ! elemental communication map.
!  integer(Zoltan_INT), pointer :: ecmap_elemids(:)  ! element ids of elements for
!                                               ! all elemental communication
!                                               ! maps. (local numbering)
!  integer(Zoltan_INT), pointer :: ecmap_sideids(:)  ! side ids of elements for all
!                                               ! elemental communication maps.
!  integer(Zoltan_INT), pointer :: ecmap_neighids(:) ! elements ids of neighboring
!                                               ! elements for all elemental
!                                               ! communication maps. 
!                                               ! (global numbering)
!  type(ELEM_INFO), pointer :: elements(:)      ! array of elements in the mesh.
!end type
!
!/*
! * global struct for mesh description
! * The Zoltan callback functions need both the element information struct
! * array and the mesh information struct. It is a lot easier to just pass
! * the element struct array as a data pointer, and have the mesh information
! * as a global variable.
! */

type(MESH_INFO),pointer :: Mesh


!/* typedef for parameter strings. 
!   Parameters are specified as pairs
!   of strings:
!   param_str = value_str              */

type Parameter_Pair
   character(len=MAX_PARAMETER_LEN) :: str(0:1)
end type Parameter_Pair

!/* Structure for the problem description. */

type PROB_INFO
  character(len=32) :: method      ! this is the method string that will
                                   ! be passed unchanged to Zoltan
  integer(Zoltan_INT) :: num_params     ! number of parameters read.
  type(Parameter_Pair), pointer :: params(:) !/* parameter array to be passed to
                                             ! Zoltan.  Parameters are specified
                                             ! as pairs of strings:
                                             ! param_str = value_str
end type

end module dr_const
