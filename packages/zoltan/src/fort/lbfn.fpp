
subroutine LB_to_ZZ(lb, zz)
type(LB_Struct) INTENT_IN lb
type(Zoltan_Struct) INTENT_OUT zz
zz%addr = lb%addr
#ifdef ABSOFT
! workaround for a bug in the Absoft compiler
zz%dummy = lb%dummy
#endif
end subroutine LB_to_ZZ

subroutine ZZ_to_LB(zz, lb)
type(Zoltan_Struct) INTENT_IN zz
type(LB_Struct) INTENT_OUT lb
lb%addr = zz%addr
#ifdef ABSOFT
! workaround for a bug in the Absoft compiler
lb%dummy = zz%dummy
#endif
end subroutine ZZ_to_LB

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function LBf90_Create(communicator)
type(LB_Struct), pointer :: LBf90_Create
type(Zoltan_Struct), pointer :: zz
integer INTENT_IN communicator
allocate(LBf90_Create)
zz => Zf90_Create(communicator)
if (associated(zz)) then
  call ZZ_to_LB(zz, LBf90_Create)
  deallocate(zz)
  nullify(zz)
else
  deallocate(LBf90_Create)
  nullify(LBf90_Create)
endif
end function LBf90_Create

subroutine LBf90_Destroy(lb)
type(LB_Struct), pointer :: lb
type(Zoltan_Struct), pointer :: zz
if (associated(lb)) then
  allocate(zz)
  call LB_to_ZZ(lb,zz)
  call Zf90_Destroy(zz)
  deallocate(lb)
  nullify(lb)
endif
end subroutine LBf90_Destroy

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function LBf90_Set_Fn0f(lb,fn_type,fn_ptr)
integer(Zoltan_INT) :: LBf90_Set_Fn0f
type(LB_Struct) INTENT_IN lb
type(ZOLTAN_FN_TYPEF) INTENT_IN fn_type
integer(Zoltan_INT), external :: fn_ptr
integer(Zoltan_INT), dimension(Zoltan_PTR_LENGTH) :: lb_addr
integer(Zoltan_INT) :: nbytes, i
nbytes = Zoltan_PTR_LENGTH
do i=1,nbytes
   lb_addr(i) = ichar(lb%addr%addr(i:i))
end do
#ifdef NASOFTWARE
LBf90_Set_Fn0f = Zfw_Set_Fn0f(lb_addr,nbytes,fn_type%choice,loc(fn_ptr))
#else
LBf90_Set_Fn0f = Zfw_Set_Fn0f(lb_addr,nbytes,fn_type%choice,fn_ptr)
#endif
end function LBf90_Set_Fn0f

function LBf90_Set_Fn0s(lb,fn_type,fn_ptr)
integer(Zoltan_INT) :: LBf90_Set_Fn0s
type(LB_Struct) INTENT_IN lb
type(ZOLTAN_FN_TYPES) INTENT_IN fn_type
external fn_ptr
integer(Zoltan_INT), dimension(Zoltan_PTR_LENGTH) :: lb_addr
integer(Zoltan_INT) :: nbytes, i
nbytes = Zoltan_PTR_LENGTH
do i=1,nbytes
   lb_addr(i) = ichar(lb%addr%addr(i:i))
end do
#ifdef NASOFTWARE
LBf90_Set_Fn0s = Zfw_Set_Fn0s(lb_addr,nbytes,fn_type%choice,loc(fn_ptr))
#else
LBf90_Set_Fn0s = Zfw_Set_Fn0s(lb_addr,nbytes,fn_type%choice,fn_ptr)
#endif
end function LBf90_Set_Fn0s

function LBf90_Set_Fn1f(lb,fn_type,fn_ptr,data)
integer(Zoltan_INT) :: LBf90_Set_Fn1f
type(LB_Struct) INTENT_IN lb
type(ZOLTAN_FN_TYPEF) INTENT_IN fn_type
integer(Zoltan_INT), external :: fn_ptr
integer(Zoltan_INT) INTENT_IN data(*)
integer(Zoltan_INT), dimension(Zoltan_PTR_LENGTH) :: lb_addr
integer(Zoltan_INT) :: nbytes, i
nbytes = Zoltan_PTR_LENGTH
do i=1,nbytes
   lb_addr(i) = ichar(lb%addr%addr(i:i))
end do
#ifdef NASOFTWARE
LBf90_Set_Fn1f = Zfw_Set_Fn1f(lb_addr,nbytes,fn_type%choice,loc(fn_ptr),data)
#else
LBf90_Set_Fn1f = Zfw_Set_Fn1f(lb_addr,nbytes,fn_type%choice,fn_ptr,data)
#endif
end function LBf90_Set_Fn1f

function LBf90_Set_Fn1s(lb,fn_type,fn_ptr,data)
integer(Zoltan_INT) :: LBf90_Set_Fn1s
type(LB_Struct) INTENT_IN lb
type(ZOLTAN_FN_TYPES) INTENT_IN fn_type
external fn_ptr
integer(Zoltan_INT) INTENT_IN data(*)
integer(Zoltan_INT), dimension(Zoltan_PTR_LENGTH) :: lb_addr
integer(Zoltan_INT) :: nbytes, i
nbytes = Zoltan_PTR_LENGTH
do i=1,nbytes
   lb_addr(i) = ichar(lb%addr%addr(i:i))
end do
#ifdef NASOFTWARE
LBf90_Set_Fn1s = Zfw_Set_Fn1s(lb_addr,nbytes,fn_type%choice,loc(fn_ptr),data)
#else
LBf90_Set_Fn1s = Zfw_Set_Fn1s(lb_addr,nbytes,fn_type%choice,fn_ptr,data)
#endif
end function LBf90_Set_Fn1s

function LBf90_Set_Fn2f(lb,fn_type,fn_ptr,data)
integer(Zoltan_INT) :: LBf90_Set_Fn2f
type(LB_Struct) INTENT_IN lb
type(ZOLTAN_FN_TYPEF) INTENT_IN fn_type
integer(Zoltan_INT), external :: fn_ptr
real(Zoltan_FLOAT) INTENT_IN data(*)
integer(Zoltan_INT), dimension(Zoltan_PTR_LENGTH) :: lb_addr
integer(Zoltan_INT) :: nbytes, i
nbytes = Zoltan_PTR_LENGTH
do i=1,nbytes
   lb_addr(i) = ichar(lb%addr%addr(i:i))
end do
#ifdef NASOFTWARE
LBf90_Set_Fn2f = Zfw_Set_Fn2f(lb_addr,nbytes,fn_type%choice,loc(fn_ptr),data)
#else
LBf90_Set_Fn2f = Zfw_Set_Fn2f(lb_addr,nbytes,fn_type%choice,fn_ptr,data)
#endif
end function LBf90_Set_Fn2f

function LBf90_Set_Fn2s(lb,fn_type,fn_ptr,data)
integer(Zoltan_INT) :: LBf90_Set_Fn2s
type(LB_Struct) INTENT_IN lb
type(ZOLTAN_FN_TYPES) INTENT_IN fn_type
external fn_ptr
real(Zoltan_FLOAT) INTENT_IN data(*)
integer(Zoltan_INT), dimension(Zoltan_PTR_LENGTH) :: lb_addr
integer(Zoltan_INT) :: nbytes, i
nbytes = Zoltan_PTR_LENGTH
do i=1,nbytes
   lb_addr(i) = ichar(lb%addr%addr(i:i))
end do
#ifdef NASOFTWARE
LBf90_Set_Fn2s = Zfw_Set_Fn2s(lb_addr,nbytes,fn_type%choice,loc(fn_ptr),data)
#else
LBf90_Set_Fn2s = Zfw_Set_Fn2s(lb_addr,nbytes,fn_type%choice,fn_ptr,data)
#endif
end function LBf90_Set_Fn2s

function LBf90_Set_Fn3f(lb,fn_type,fn_ptr,data)
integer(Zoltan_INT) :: LBf90_Set_Fn3f
type(LB_Struct) INTENT_IN lb
type(ZOLTAN_FN_TYPEF) INTENT_IN fn_type
integer(Zoltan_INT), external :: fn_ptr
real(Zoltan_DOUBLE) INTENT_IN data(*)
integer(Zoltan_INT), dimension(Zoltan_PTR_LENGTH) :: lb_addr
integer(Zoltan_INT) :: nbytes, i
nbytes = Zoltan_PTR_LENGTH
do i=1,nbytes
   lb_addr(i) = ichar(lb%addr%addr(i:i))
end do
#ifdef NASOFTWARE
LBf90_Set_Fn3f = Zfw_Set_Fn3f(lb_addr,nbytes,fn_type%choice,loc(fn_ptr),data)
#else
LBf90_Set_Fn3f = Zfw_Set_Fn3f(lb_addr,nbytes,fn_type%choice,fn_ptr,data)
#endif
end function LBf90_Set_Fn3f

function LBf90_Set_Fn3s(lb,fn_type,fn_ptr,data)
integer(Zoltan_INT) :: LBf90_Set_Fn3s
type(LB_Struct) INTENT_IN lb
type(ZOLTAN_FN_TYPES) INTENT_IN fn_type
external fn_ptr
real(Zoltan_DOUBLE) INTENT_IN data(*)
integer(Zoltan_INT), dimension(Zoltan_PTR_LENGTH) :: lb_addr
integer(Zoltan_INT) :: nbytes, i
nbytes = Zoltan_PTR_LENGTH
do i=1,nbytes
   lb_addr(i) = ichar(lb%addr%addr(i:i))
end do
#ifdef NASOFTWARE
LBf90_Set_Fn3s = Zfw_Set_Fn3s(lb_addr,nbytes,fn_type%choice,loc(fn_ptr),data)
#else
LBf90_Set_Fn3s = Zfw_Set_Fn3s(lb_addr,nbytes,fn_type%choice,fn_ptr,data)
#endif
end function LBf90_Set_Fn3s

function LBf90_Set_Fn8f(lb,fn_type,fn_ptr,data)
integer(Zoltan_INT) :: LBf90_Set_Fn8f
type(LB_Struct) INTENT_IN lb
type(ZOLTAN_FN_TYPEF) INTENT_IN fn_type
integer(Zoltan_INT), external :: fn_ptr
type(LB_User_Data_1) INTENT_IN data
integer(Zoltan_INT), dimension(Zoltan_PTR_LENGTH) :: lb_addr
integer(Zoltan_INT) :: nbytes, i
nbytes = Zoltan_PTR_LENGTH
do i=1,nbytes
   lb_addr(i) = ichar(lb%addr%addr(i:i))
end do
#ifdef NASOFTWARE
LBf90_Set_Fn8f = Zfw_Set_Fn8f(lb_addr,nbytes,fn_type%choice,loc(fn_ptr),data)
#else
LBf90_Set_Fn8f = Zfw_Set_Fn8f(lb_addr,nbytes,fn_type%choice,fn_ptr,data)
#endif
end function LBf90_Set_Fn8f

function LBf90_Set_Fn8s(lb,fn_type,fn_ptr,data)
integer(Zoltan_INT) :: LBf90_Set_Fn8s
type(LB_Struct) INTENT_IN lb
type(ZOLTAN_FN_TYPES) INTENT_IN fn_type
external fn_ptr
type(LB_User_Data_1) INTENT_IN data
integer(Zoltan_INT), dimension(Zoltan_PTR_LENGTH) :: lb_addr
integer(Zoltan_INT) :: nbytes, i
nbytes = Zoltan_PTR_LENGTH
do i=1,nbytes
   lb_addr(i) = ichar(lb%addr%addr(i:i))
end do
#ifdef NASOFTWARE
LBf90_Set_Fn8s = Zfw_Set_Fn8s(lb_addr,nbytes,fn_type%choice,loc(fn_ptr),data)
#else
LBf90_Set_Fn8s = Zfw_Set_Fn8s(lb_addr,nbytes,fn_type%choice,fn_ptr,data)
#endif
end function LBf90_Set_Fn8s

function LBf90_Set_Fn9f(lb,fn_type,fn_ptr,data)
integer(Zoltan_INT) :: LBf90_Set_Fn9f
type(LB_Struct) INTENT_IN lb
type(ZOLTAN_FN_TYPEF) INTENT_IN fn_type
integer(Zoltan_INT), external :: fn_ptr
type(LB_User_Data_2) INTENT_IN data
integer(Zoltan_INT), dimension(Zoltan_PTR_LENGTH) :: lb_addr
integer(Zoltan_INT) :: nbytes, i
nbytes = Zoltan_PTR_LENGTH
do i=1,nbytes
   lb_addr(i) = ichar(lb%addr%addr(i:i))
end do
#ifdef NASOFTWARE
LBf90_Set_Fn9f = Zfw_Set_Fn9f(lb_addr,nbytes,fn_type%choice,loc(fn_ptr),data)
#else
LBf90_Set_Fn9f = Zfw_Set_Fn9f(lb_addr,nbytes,fn_type%choice,fn_ptr,data)
#endif
end function LBf90_Set_Fn9f

function LBf90_Set_Fn9s(lb,fn_type,fn_ptr,data)
integer(Zoltan_INT) :: LBf90_Set_Fn9s
type(LB_Struct) INTENT_IN lb
type(ZOLTAN_FN_TYPES) INTENT_IN fn_type
external fn_ptr
type(LB_User_Data_2) INTENT_IN data
integer(Zoltan_INT), dimension(Zoltan_PTR_LENGTH) :: lb_addr
integer(Zoltan_INT) :: nbytes, i
nbytes = Zoltan_PTR_LENGTH
do i=1,nbytes
   lb_addr(i) = ichar(lb%addr%addr(i:i))
end do
#ifdef NASOFTWARE
LBf90_Set_Fn9s = Zfw_Set_Fn9s(lb_addr,nbytes,fn_type%choice,loc(fn_ptr),data)
#else
LBf90_Set_Fn9s = Zfw_Set_Fn9s(lb_addr,nbytes,fn_type%choice,fn_ptr,data)
#endif
end function LBf90_Set_Fn9s

function LBf90_Set_FnAf(lb,fn_type,fn_ptr,data)
integer(Zoltan_INT) :: LBf90_Set_FnAf
type(LB_Struct) INTENT_IN lb
type(ZOLTAN_FN_TYPEF) INTENT_IN fn_type
integer(Zoltan_INT), external :: fn_ptr
type(LB_User_Data_3) INTENT_IN data
integer(Zoltan_INT), dimension(Zoltan_PTR_LENGTH) :: lb_addr
integer(Zoltan_INT) :: nbytes, i
nbytes = Zoltan_PTR_LENGTH
do i=1,nbytes
   lb_addr(i) = ichar(lb%addr%addr(i:i))
end do
#ifdef NASOFTWARE
LBf90_Set_FnAf = Zfw_Set_FnAf(lb_addr,nbytes,fn_type%choice,loc(fn_ptr),data)
#else
LBf90_Set_FnAf = Zfw_Set_FnAf(lb_addr,nbytes,fn_type%choice,fn_ptr,data)
#endif
end function LBf90_Set_FnAf

function LBf90_Set_FnAs(lb,fn_type,fn_ptr,data)
integer(Zoltan_INT) :: LBf90_Set_FnAs
type(LB_Struct) INTENT_IN lb
type(ZOLTAN_FN_TYPES) INTENT_IN fn_type
external fn_ptr
type(LB_User_Data_3) INTENT_IN data
integer(Zoltan_INT), dimension(Zoltan_PTR_LENGTH) :: lb_addr
integer(Zoltan_INT) :: nbytes, i
nbytes = Zoltan_PTR_LENGTH
do i=1,nbytes
   lb_addr(i) = ichar(lb%addr%addr(i:i))
end do
#ifdef NASOFTWARE
LBf90_Set_FnAs = Zfw_Set_FnAs(lb_addr,nbytes,fn_type%choice,loc(fn_ptr),data)
#else
LBf90_Set_FnAs = Zfw_Set_FnAs(lb_addr,nbytes,fn_type%choice,fn_ptr,data)
#endif
end function LBf90_Set_FnAs

function LBf90_Set_FnBf(lb,fn_type,fn_ptr,data)
integer(Zoltan_INT) :: LBf90_Set_FnBf
type(LB_Struct) INTENT_IN lb
type(ZOLTAN_FN_TYPEF) INTENT_IN fn_type
integer(Zoltan_INT), external :: fn_ptr
type(LB_User_Data_4) INTENT_IN data
integer(Zoltan_INT), dimension(Zoltan_PTR_LENGTH) :: lb_addr
integer(Zoltan_INT) :: nbytes, i
nbytes = Zoltan_PTR_LENGTH
do i=1,nbytes
   lb_addr(i) = ichar(lb%addr%addr(i:i))
end do
#ifdef NASOFTWARE
LBf90_Set_FnBf = Zfw_Set_FnBf(lb_addr,nbytes,fn_type%choice,loc(fn_ptr),data)
#else
LBf90_Set_FnBf = Zfw_Set_FnBf(lb_addr,nbytes,fn_type%choice,fn_ptr,data)
#endif
end function LBf90_Set_FnBf

function LBf90_Set_FnBs(lb,fn_type,fn_ptr,data)
integer(Zoltan_INT) :: LBf90_Set_FnBs
type(LB_Struct) INTENT_IN lb
type(ZOLTAN_FN_TYPES) INTENT_IN fn_type
external fn_ptr
type(LB_User_Data_4) INTENT_IN data
integer(Zoltan_INT), dimension(Zoltan_PTR_LENGTH) :: lb_addr
integer(Zoltan_INT) :: nbytes, i
nbytes = Zoltan_PTR_LENGTH
do i=1,nbytes
   lb_addr(i) = ichar(lb%addr%addr(i:i))
end do
#ifdef NASOFTWARE
LBf90_Set_FnBs = Zfw_Set_FnBs(lb_addr,nbytes,fn_type%choice,loc(fn_ptr),data)
#else
LBf90_Set_FnBs = Zfw_Set_FnBs(lb_addr,nbytes,fn_type%choice,fn_ptr,data)
#endif
end function LBf90_Set_FnBs

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function LBf90_LB_Set_Method(lb,string)
integer(Zoltan_INT) :: LBf90_LB_Set_Method
type(LB_Struct) INTENT_IN lb
character(len=*) INTENT_IN string
type(Zoltan_Struct) :: zz
call LB_to_ZZ(lb,zz)
LBf90_LB_Set_Method = Zf90_Set_Param(zz,"LB_METHOD",string)
end function LBf90_LB_Set_Method

function LBf90_Set_Param(lb,param_name,new_value)
integer(Zoltan_INT) :: LBf90_Set_Param
type(LB_Struct) INTENT_IN lb
character(len=*) INTENT_IN param_name, new_value
type(Zoltan_Struct) :: zz
call LB_to_ZZ(lb,zz)
LBf90_Set_Param = Zf90_Set_Param(zz,param_name,new_value)
end function LBf90_Set_Param

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function LBf90_LB_Balance(lb,changes,num_gid_entries,num_lid_entries, &
                       num_import,import_global_ids, &
                       import_local_ids,import_procs,num_export, &
                       export_global_ids,export_local_ids,export_procs)
integer(Zoltan_INT) :: LBf90_LB_Balance
type(LB_Struct) INTENT_IN lb
logical, intent(out) :: changes
integer(Zoltan_INT), intent(out) :: num_gid_entries, num_lid_entries
integer(Zoltan_INT), intent(out) :: num_import, num_export
integer(Zoltan_INT), pointer, dimension(:) :: import_global_ids, export_global_ids
integer(Zoltan_INT), pointer, dimension(:) :: import_local_ids, export_local_ids
integer(Zoltan_INT), pointer, dimension(:) :: import_procs, export_procs
type(Zoltan_Struct) :: zz
call LB_to_ZZ(lb,zz)
LBf90_LB_Balance = Zf90_LB_Balance(zz,changes, &
                             num_gid_entries, num_lid_entries, &
                             num_import,import_global_ids,import_local_ids, &
                             import_procs,num_export,export_global_ids, &
                             export_local_ids,export_procs)
end function LBf90_LB_Balance


function LBf90_LB_Eval(lb,print_stats,nobj,obj_wgt, &
                    ncuts,cut_wgt,nboundary,nadj)
integer(Zoltan_INT) :: LBf90_LB_Eval
type(LB_Struct) INTENT_IN lb
logical INTENT_IN print_stats
integer(Zoltan_INT), intent(out), optional :: nobj, ncuts, nboundary, nadj
real(Zoltan_FLOAT), intent(out), optional :: obj_wgt(*), cut_wgt(*)
type(Zoltan_Struct) :: zz
call LB_to_ZZ(lb,zz)
LBf90_LB_Eval = Zf90_LB_Eval(zz,print_stats,nobj,obj_wgt, &
                              ncuts,cut_wgt,nboundary,nadj)
end function LBf90_LB_Eval

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function LBf90_LB_Point_Assign(lb,coords,proc)
integer(Zoltan_INT) :: LBf90_LB_Point_Assign
type(LB_Struct) INTENT_IN lb
real(Zoltan_DOUBLE), dimension(*) INTENT_IN coords
integer(Zoltan_INT), intent(out) :: proc
type(Zoltan_Struct) :: zz
call LB_to_ZZ(lb,zz)
LBf90_LB_Point_Assign = Zf90_LB_Point_Assign(zz,coords,proc)
end function LBf90_LB_Point_Assign

function LBf90_LB_Box_Assign(lb,xmin,ymin,zmin,xmax,ymax,zmax,procs,numprocs)
integer(Zoltan_INT) :: LBf90_LB_Box_Assign
type(LB_Struct) INTENT_IN lb
real(Zoltan_DOUBLE) INTENT_IN xmin,ymin,zmin,xmax,ymax,zmax
integer(Zoltan_INT), intent(out), dimension(*) :: procs
integer(Zoltan_INT), intent(out) :: numprocs
type(Zoltan_Struct) :: zz
call LB_to_ZZ(lb,zz)
LBf90_LB_Box_Assign = Zf90_LB_Box_Assign(zz,xmin,ymin,zmin,xmax,ymax,zmax, &
                                          procs,numprocs)
end function LBf90_LB_Box_Assign

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function LBf90_Compute_Destinations(lb, &
                       num_import,import_global_ids, &
                       import_local_ids,import_procs,num_export, &
                       export_global_ids,export_local_ids,export_procs)
integer(Zoltan_INT) :: LBf90_Compute_Destinations
type(LB_Struct) INTENT_IN lb
integer(Zoltan_INT) INTENT_IN num_import
integer(Zoltan_INT), intent(out) :: num_export
integer(Zoltan_INT), pointer, dimension(:) :: import_global_ids, export_global_ids
integer(Zoltan_INT), pointer, dimension(:) :: import_local_ids, export_local_ids
integer(Zoltan_INT), pointer, dimension(:) :: import_procs, export_procs
type(Zoltan_Struct) :: zz
call LB_to_ZZ(lb,zz)
LBf90_Compute_Destinations = Zf90_Compute_Destinations(zz, &
                             num_import,import_global_ids,import_local_ids, &
                             import_procs,num_export,export_global_ids, &
                             export_local_ids,export_procs)
end function LBf90_Compute_Destinations


function LBf90_Help_Migrate(lb, &
                       num_import,import_global_ids, &
                       import_local_ids,import_procs,num_export, &
                       export_global_ids,export_local_ids,export_procs)
integer(Zoltan_INT) :: LBf90_Help_Migrate
type(LB_Struct) INTENT_IN lb
integer(Zoltan_INT) INTENT_IN num_import, num_export
integer(Zoltan_INT), pointer, dimension(:) :: import_global_ids, export_global_ids
integer(Zoltan_INT), pointer, dimension(:) :: import_local_ids, export_local_ids
integer(Zoltan_INT), pointer, dimension(:) :: import_procs, export_procs
type(Zoltan_Struct) :: zz
call LB_to_ZZ(lb,zz)
LBf90_Help_Migrate = Zf90_Help_Migrate(zz, &
                             num_import,import_global_ids,import_local_ids, &
                             import_procs,num_export,export_global_ids, &
                             export_local_ids,export_procs)
end function LBf90_Help_Migrate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine LBf90_Reftree_Get_Child_Order(lb,order,ierr)
type(LB_Struct) INTENT_IN lb
integer(Zoltan_INT), intent(inout), dimension(*) :: order
integer(Zoltan_INT), intent(out) :: ierr
type(Zoltan_Struct) :: zz
call LB_to_ZZ(lb,zz)
call Zf90_Reftree_Get_Child_Order(zz,order,ierr)
end subroutine LBf90_Reftree_Get_Child_Order

