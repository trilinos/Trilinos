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

module zoltan_user_data
use zoltan_types
implicit none

! User defined data types for passing data to the query functions.  These can
! be used any way you want, but one suggestion is to use them as "wrapper"
! types for your own user defined types, e.g.
! type Zoltan_User_Data_1 
!    type(your_type), pointer :: ptr
! end type
! Exactly four data types must be defined, but you don't have to use them.

integer(Zoltan_INT), parameter :: MAX_EB_NAME_LEN = 32 ! chars for element block

!/*
! * Structure used to describe an element. Each processor will
! * allocate an array of these structures.
!   Moved here from dr_consts.f90 so that User_Data can use it
! */
type ELEM_INFO
  integer(Zoltan_INT) :: border   ! set to 1 if this element is a border element
  integer(Zoltan_INT) :: globalID ! Global ID of this element; local ID is the
                                  ! position in the array of elements
  integer(Zoltan_INT) :: elem_blk ! elem block number which this element is in
  integer(Zoltan_INT) :: my_part  ! partition to which this element is assigned
  integer(Zoltan_INT) :: perm_value  ! permutation value
  integer(Zoltan_INT) :: invperm_value  ! inverse permutation value
  real(Zoltan_FLOAT)  :: cpu_wgt  ! computational weight associated with elem
  real(Zoltan_FLOAT)  :: mem_wgt  ! the memory weight associated with the elem
  real(Zoltan_FLOAT), pointer ::  coord(:,:) ! array for the coordinates of the
                                        ! element. For Nemesis meshes, nodal
                                        ! coordinates are stored; for Chaco
                                        ! graphs with geometry, one set of
                                        ! coords is stored.
  integer(Zoltan_INT), pointer :: connect(:) ! list of nodes that make up this
                                        ! element, the node numbers in this
                                        ! list are global and not local
  integer(Zoltan_INT), pointer :: adj(:)  ! list of adjacent elements .
                         ! For Nemesis input, the list is ordered by
                         ! side number, to encode side-number info needed to
                         ! rebuild communication maps.  Value -1 represents
                         ! sides with no neighboring element (e.g., along mesh
                         ! boundaries).  Chaco doesn't have "sides," so the
                         ! ordering is irrelevent for Chaco input.
  integer(Zoltan_INT), pointer :: adj_proc(:) ! list of processors for adjacent
                                         ! elements
  real(Zoltan_FLOAT), pointer :: edge_wgt(:)  ! edge weights for adjacent elements
  integer(Zoltan_INT) :: nadj    ! number of entries in adj
  integer(Zoltan_INT) :: adj_len ! allocated length of adj/adj_proc/edge_wgt arrays
end type

!/*
! * structure for general mesh information
! */
!/* Structure used to store information about the mesh */
type MESH_INFO
  integer(Zoltan_INT) :: num_nodes     ! number of nodes on this processor
  integer(Zoltan_INT) :: num_elems     ! number of elements on this processor
  integer(Zoltan_INT) :: num_dims      ! number of dimensions for the mesh
  integer(Zoltan_INT) :: num_el_blks   ! number of element blocks in the mesh
  integer(Zoltan_INT) :: num_node_sets ! number of node sets in the mesh
  integer(Zoltan_INT) :: num_side_sets ! number of side sets in the mesh
  character(len=MAX_EB_NAME_LEN), pointer :: eb_names(:) ! element block element
                                                         ! names
  integer(Zoltan_INT), pointer :: eb_ids(:)    ! element block ids
  integer(Zoltan_INT), pointer :: eb_cnts(:)   ! number of elements in each element
                                          ! block
  integer(Zoltan_INT), pointer :: eb_nnodes(:) ! number of nodes per element in each
                                          ! element block
                                          ! for Nemesis meshes, this value
                                          ! depends on element type;
                                          ! for Chaco graphs, only one "node"
                                          ! per element.
  integer(Zoltan_INT), pointer :: eb_nattrs(:) ! number of attributes per element in
                                          ! each element block
  integer(Zoltan_INT) :: elem_array_len ! length that the ELEM_INFO array is
                                   ! allocated for. Need to know this when array
                                   ! is not completely filled during migration
  integer(Zoltan_INT) :: necmap         ! number of elemental communication maps.
  integer(Zoltan_INT), pointer :: ecmap_id(:)       ! IDs of each elemental
                                              ! communication map.
  integer(Zoltan_INT), pointer :: ecmap_cnt(:)      ! number of elements in each
                                               ! elemental communication map.
  integer(Zoltan_INT), pointer :: ecmap_elemids(:)  ! element ids of elements for
                                               ! all elemental communication
                                               ! maps. (local numbering)
  integer(Zoltan_INT), pointer :: ecmap_sideids(:)  ! side ids of elements for all
                                               ! elemental communication maps.
  integer(Zoltan_INT), pointer :: ecmap_neighids(:) ! elements ids of neighboring
                                               ! elements for all elemental
                                               ! communication maps.
                                               ! (global numbering)
  type(ELEM_INFO), pointer :: elements(:)      ! array of elements in the mesh.
end type

type Zoltan_User_Data_1
   type(ELEM_INFO), pointer :: ptr(:)
end type Zoltan_User_Data_1

type Zoltan_User_Data_2
   type(MESH_INFO), pointer :: ptr
end type Zoltan_User_Data_2

type Zoltan_User_Data_3
   integer :: dummy
end type Zoltan_User_Data_3

type Zoltan_User_Data_4
   integer :: dummy
end type Zoltan_User_Data_4

end module zoltan_user_data
