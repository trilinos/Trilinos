Program TestExoRead_Part
   include 'exodusII.inc'
   !Use exodusIIF90

   Integer                      :: exoid,cpu_ws,io_ws,mod_sz,ierr,i,nparts
   Integer                      :: num_dim,num_nodes,num_elem,num_elem_blk,num_node_sets,num_side_sets
   Integer                      :: time_step = 1
   Integer                      :: var_index = 2
   Integer                      :: istart,iend,len
   Real,dimension(:),Pointer    :: var_values
   character*(MXLNLN)           :: titl
   Real                         :: vers


   cpu_ws = 0
   io_ws = 0

   exoid = exopen ("test.exo", EXWRIT, cpu_ws, io_ws, vers, ierr)
   write (*, '("after exopen, error = ",i3)') ierr

   write (*, '("test.exo is an EXODUSII file; version ", f4.2)') vers
   write (*, '("  I/O word size",i2)') io_ws
   write (*, '("  CPU word size",i2)') cpu_ws
   mod_sz = exlgmd(exoid)
   write (*, '("  Model Size",i2)') mod_sz

   call exgini (exoid, titl, num_dim, num_nodes, num_elem,num_elem_blk, num_node_sets, num_side_sets, ierr)
   write (*, '("after exgini, error = ", i3)' ) ierr

   write (*, '("database parameters")')
   write (*, '("   title = ",a81)') titl
   write (*, '("   num_dim = ", i3 )') num_dim
   write (*, '("   num_nodes = ", i3 )') num_nodes
   write (*, '("   num_elem = ", i3 )') num_elem
   write (*, '("   num_elem_blk = ", i3 )') num_elem_blk
   write (*, '("   num_node_sets = ", i3 )') num_node_sets
   write (*, '("   num_side_sets = ", i3)') num_side_sets

   Allocate(var_values(num_nodes))
   call exgnv (exoid, time_step, var_index, num_nodes, var_values,ierr)
   write (*, '("after exgnv, error = ", i3)' ) ierr
   write(*,*) var_values

!   call expnv (exoid, time_step, var_index, num_nodes, var_values,ierr)
   DeAllocate(var_values)

   nparts=4
   Do i = 1, nparts
      istart = (i-1)*num_nodes/nparts+1
      iend   = i*num_nodes/nparts
      len    = iend - istart + 1
      write(*,'("   chunk ",i3," -- ",i3, ":")') istart,iend
      Allocate(var_values(len))
      call exgpv(exoid, time_step, EX_NODAL,var_index,0,istart,len,var_values,ierr)
      write (*, '("after exgpv, error = ", i3)' ) ierr
      write(*,*) '   *', var_values
      var_values = var_values + 1
      call exppv(exoid,time_step,EX_NODAL,var_index,0,istart,len,var_values,ierr)
      write (*, '("after exppv, error = ", i3)' ) ierr
      DeAllocate(var_values)
   End Do

   Do i = 1, nparts
      istart = (i-1)*num_nodes/nparts+1
      iend   = i*num_nodes/nparts
      len    = iend - istart + 1
      write(*,'("   chunk ",i3," -- ",i3, ":")') istart,iend
      Allocate(var_values(len))
      call exgpcc(exoid,istart,len,1,var_values,ierr)
      write (*, '("after exgpcc, error = ", i3)' ) ierr
      write(*,*) '   coord', var_values
      var_values = var_values + 1
      call exppcc(exoid,istart,len,1,var_values,ierr)
      write (*, '("after expcc, error = ", i3)' ) ierr
      DeAllocate(var_values)
   End Do

   call exclos (exoid, ierr)
   write (*, '("after exclos, error = ", i3)' ) ierr
End Program TestExoRead_Part
