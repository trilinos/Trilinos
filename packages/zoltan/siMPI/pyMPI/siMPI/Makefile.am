#**************************************************************************#
#* FILE   **************        Makefile.am        ************************#
#**************************************************************************#
#* Author: Patrick Miller July 19 2002					  *#
#* Copyright (C) 2002 University of California Regents			  *#
#**************************************************************************#
#*                                                                        *#
#**************************************************************************#

#We don't want to generally add this flag
#AM_CFLAGS = -fPIC

#bin_PROGRAMS = mpicc mpif77 mpirun
bin_PROGRAMS = mpirun
#noinst_PROGRAMS = smoketest #fortran_grinder fortran_parameterize

lib_LIBRARIES = libsimpi.a
#include_HEADERS = mpi.h mpif.h mpi_implementation.h mpi_prototypes.h mpi_profile.h mpi_config.h _MPI_UTILITY.h
include_HEADERS = mpi.h mpi_implementation.h mpi_prototypes.h mpi_profile.h mpi_config.h _MPI_UTILITY.h

#fortran_grinder_SOURCES = \
#    fortran_grinder.c

#fortran_parameterize_SOURCES = \
#    fortran_parameterize.c

#mpicc_SOURCES = \
#    mpicc.c

#mpif77_SOURCES = \
#    mpif77.c

mpirun_SOURCES = \
    mpirun.c

#smoketest_SOURCES = \
#    smoketest.c

#Since we don't need Fortran, I removed the following files from the list below:
# F_MPI_Finalize.c, F_MPI_Init.c, F_MPI_Reduce.c

libsimpi_a_SOURCES = \
    MPI_Comm_create.c \
    MPI_Comm_free.c \
    MPI_Comm_group.c \
    MPI_Comm_rank.c \
    MPI_Comm_size.c \
    MPI_ERRORS_ARE_FATAL.c \
    MPI_ERRORS_RETURN.c \
    PMPI_Abort.c \
    MPI_Finalize.c \
    MPI_Init.c \
    PMPI_Comm_create.c \
    PMPI_Comm_free.c \
    PMPI_Comm_group.c \
    PMPI_Comm_rank.c \
    PMPI_Comm_size.c \
    PMPI_Finalize.c \
    PMPI_Init.c \
    _MPI_CHECK_STATUS.c \
    _MPI_ERR_ROUTINE.c \
    PMPI_Send.c \
    PMPI_Msend.c \
    MPI_Send.c \
    _MPI_UTILITY.c \
    _MPI_SENDREC_UTIL.c \
    _MPI_TYPE_UTIL.c \
    _MPI_variables.c \
    PMPI_Recv.c \
    MPI_Recv.c \
    PMPI_Isend.c \
    MPI_Isend.c \
    PMPI_Irecv.c \
    MPI_Irecv.c \
    PMPI_Bsend.c \
    MPI_Bsend.c \
    PMPI_Reduce.c \
    MPI_Reduce.c \
    PMPI_Bcast.c \
    _MPI_BCAST.c \
    _MPI_OP_UTIL.c \
    PMPI_Op_create.c \
    MPI_Op_create.c \
    PMPI_Op_free.c \
    MPI_Op_free.c \
    PMPI_Sendrecv.c \
    MPI_Sendrecv.c \
    PMPI_Wait.c \
    MPI_Wait.c \
    _MPI_ABORT.c \
    mpi_abort.c \
    _MPI_COMM_UTIL.c \
    _MPI_REQ_UTIL.c \
    PMPI_Ssend.c \
    MPI_Ssend.c \
    PMPI_Rsend.c \
    MPI_Rsend.c \
    PMPI_Waitany.c \
    PMPI_Scatter.c \
    MPI_Scatter.c \
    PMPI_Gather.c \
    MPI_Gather.c \
    PMPI_Issend.c \
    MPI_Issend.c \
    PMPI_Ibsend.c \
    MPI_Ibsend.c \
    PMPI_Irsend.c \
    MPI_Irsend.c \
    MPIO_Test.c \
    MPIO_Wait.c \
    PMPIO_Test.c \
    PMPIO_Wait.c \
    MPI_Address.c \
    MPI_Allgather.c \
    MPI_Allgatherv.c \
    MPI_Allreduce.c \
    MPI_Alltoall.c \
    MPI_Alltoallv.c \
    MPI_Attr_delete.c \
    MPI_Attr_get.c \
    MPI_Attr_put.c \
    MPI_Barrier.c \
    MPI_Bsend_init.c \
    MPI_Buffer_attach.c \
    MPI_Buffer_detach.c \
    MPI_Cancel.c \
    MPI_Cartdim_get.c \
    MPI_Errhandler_get.c \
    MPI_Info_dup.c \
    MPI_Cart_coords.c \
    MPI_Cart_create.c \
    MPI_Cart_get.c \
    MPI_Cart_map.c \
    MPI_Cart_rank.c \
    MPI_Cart_shift.c \
    MPI_Cart_sub.c \
    MPI_Comm_compare.c \
    MPI_Comm_dup.c \
    MPI_Comm_get_name.c \
    MPI_Comm_remote_group.c \
    MPI_Comm_remote_size.c \
    MPI_Comm_set_name.c \
    MPI_Comm_split.c \
    MPI_Comm_test_inter.c \
    PMPI_Address.c \
    PMPI_Allgather.c \
    PMPI_Allgatherv.c \
    PMPI_Allreduce.c \
    PMPI_Alltoall.c \
    PMPI_Alltoallv.c \
    PMPI_Attr_delete.c \
    PMPI_Attr_get.c \
    PMPI_Attr_put.c \
    PMPI_Barrier.c \
    PMPI_Bsend_init.c \
    PMPI_Buffer_attach.c \
    PMPI_Buffer_detach.c \
    PMPI_Cartdim_get.c \
    PMPI_Errhandler_get.c \
    PMPI_Info_dup.c \
    PMPI_Cancel.c \
    PMPI_Cart_coords.c \
    PMPI_Cart_create.c \
    PMPI_Cart_get.c \
    PMPI_Cart_map.c \
    PMPI_Cart_rank.c \
    PMPI_Cart_shift.c \
    PMPI_Cart_sub.c \
    PMPI_Comm_compare.c \
    PMPI_Comm_dup.c \
    PMPI_Comm_get_name.c \
    PMPI_Comm_remote_group.c \
    PMPI_Comm_remote_size.c \
    PMPI_Comm_set_name.c \
    PMPI_Comm_split.c \
    PMPI_Comm_test_inter.c \
    MPI_File_c2f.c \
    MPI_File_close.c \
    MPI_File_delete.c \
    MPI_File_f2c.c \
    MPI_File_get_amode.c \
    MPI_File_get_atomicity.c \
    MPI_File_get_byte_offset.c \
    MPI_File_get_errhandler.c \
    MPI_File_get_group.c \
    MPI_File_get_info.c \
    MPI_File_get_position.c \
    MPI_File_get_position_shared.c \
    MPI_File_get_size.c \
    MPI_File_get_type_extent.c \
    MPI_File_get_view.c \
    MPI_File_iread.c \
    MPI_File_iread_at.c \
    MPI_File_iread_shared.c \
    MPI_File_iwrite.c \
    MPI_File_iwrite_at.c \
    MPI_File_iwrite_shared.c \
    MPI_File_open.c \
    MPI_File_preallocate.c \
    MPI_File_read.c \
    MPI_File_read_all.c \
    MPI_File_read_all_begin.c \
    MPI_File_read_all_end.c \
    MPI_File_read_at.c \
    MPI_File_read_at_all.c \
    MPI_File_read_at_all_begin.c \
    MPI_File_read_at_all_end.c \
    MPI_File_read_ordered.c \
    MPI_File_read_ordered_begin.c \
    MPI_File_read_ordered_end.c \
    MPI_File_read_shared.c \
    MPI_File_seek.c \
    MPI_File_seek_shared.c \
    MPI_File_set_atomicity.c \
    MPI_File_set_errhandler.c \
    MPI_File_set_info.c \
    MPI_File_set_size.c \
    MPI_File_set_view.c \
    MPI_Dims_create.c \
    MPI_Errhandler_create.c \
    MPI_Errhandler_free.c \
    MPI_Errhandler_set.c \
    MPI_Error_class.c \
    MPI_Error_string.c \
    PMPI_File_c2f.c \
    PMPI_File_close.c \
    PMPI_File_delete.c \
    PMPI_File_f2c.c \
    PMPI_File_get_amode.c \
    PMPI_File_get_atomicity.c \
    PMPI_File_get_byte_offset.c \
    PMPI_File_get_errhandler.c \
    PMPI_File_get_group.c \
    PMPI_File_get_info.c \
    PMPI_File_get_position.c \
    PMPI_File_get_position_shared.c \
    PMPI_File_get_size.c \
    PMPI_File_get_type_extent.c \
    PMPI_File_get_view.c \
    PMPI_File_iread.c \
    PMPI_File_iread_at.c \
    PMPI_File_iread_shared.c \
    PMPI_File_iwrite.c \
    PMPI_File_iwrite_at.c \
    PMPI_File_iwrite_shared.c \
    PMPI_File_open.c \
    PMPI_File_preallocate.c \
    PMPI_File_read.c \
    PMPI_File_read_all.c \
    PMPI_File_read_all_begin.c \
    PMPI_File_read_all_end.c \
    PMPI_File_read_at.c \
    PMPI_File_read_at_all.c \
    PMPI_File_read_at_all_begin.c \
    PMPI_File_read_at_all_end.c \
    PMPI_File_read_ordered.c \
    PMPI_File_read_ordered_begin.c \
    PMPI_File_read_ordered_end.c \
    PMPI_File_read_shared.c \
    PMPI_File_seek.c \
    PMPI_File_seek_shared.c \
    PMPI_File_set_atomicity.c \
    PMPI_File_set_errhandler.c \
    PMPI_File_set_info.c \
    PMPI_File_set_size.c \
    PMPI_File_set_view.c \
    PMPI_Dims_create.c \
    PMPI_Errhandler_create.c \
    PMPI_Errhandler_free.c \
    PMPI_Errhandler_set.c \
    PMPI_Error_class.c \
    PMPI_Error_string.c \
    MPI_File_sync.c \
    MPI_File_write.c \
    MPI_File_write_all.c \
    MPI_File_write_all_begin.c \
    MPI_File_write_all_end.c \
    MPI_File_write_at.c \
    MPI_File_write_at_all.c \
    MPI_File_write_at_all_begin.c \
    MPI_File_write_at_all_end.c \
    MPI_File_write_ordered.c \
    MPI_File_write_ordered_begin.c \
    MPI_File_write_ordered_end.c \
    MPI_File_write_shared.c \
    MPI_Finalized.c \
    MPI_Gatherv.c \
    MPI_Get_count.c \
    MPI_Get_elements.c \
    MPI_Get_processor_name.c \
    MPI_Get_version.c \
    MPI_Graph_create.c \
    MPI_Graph_get.c \
    MPI_Graph_map.c \
    MPI_Graph_neighbors.c \
    MPI_Graph_neighbors_count.c \
    MPI_Graphdims_get.c \
    MPI_Group_compare.c \
    MPI_Group_difference.c \
    MPI_Group_excl.c \
    MPI_Group_free.c \
    MPI_Group_incl.c \
    MPI_Group_intersection.c \
    MPI_Group_range_excl.c \
    MPI_Group_range_incl.c \
    MPI_Group_rank.c \
    MPI_Group_size.c \
    MPI_Group_translate_ranks.c \
    MPI_Group_union.c \
    MPI_Info_c2f.c \
    MPI_Info_create.c \
    MPI_Info_delete.c \
    MPI_Info_f2c.c \
    MPI_Info_free.c \
    MPI_Info_get.c \
    MPI_Info_get_nkeys.c \
    MPI_Info_get_nthkey.c \
    MPI_Info_get_valuelen.c \
    MPI_Info_set.c \
    MPI_Init_thread.c \
    MPI_Initialized.c \
    MPI_Int2handle.c \
    MPI_Intercomm_create.c \
    MPI_Intercomm_merge.c \
    MPI_Iprobe.c \
    MPI_Keyval_create.c \
    MPI_Keyval_free.c \
    MPI_Pack.c \
    MPI_Pack_size.c \
    MPI_Pcontrol.c \
    MPI_Probe.c \
    MPI_Recv_init.c \
    MPI_Reduce_scatter.c \
    MPI_Request_c2f.c \
    MPI_Request_free.c \
    MPI_Rsend_init.c \
    MPI_Scan.c \
    MPI_Scatterv.c \
    MPI_Send_init.c \
    MPI_Sendrecv_replace.c \
    MPI_Ssend_init.c \
    MPI_Start.c \
    MPI_Startall.c \
    MPI_Status_c2f.c \
    MPI_Status_set_cancelled.c \
    MPI_Status_set_elements.c \
    MPI_Test.c \
    MPI_Test_cancelled.c \
    MPI_Testall.c \
    MPI_Testany.c \
    MPI_Testsome.c \
    MPI_Topo_test.c \
    MPI_Type_commit.c \
    MPI_Type_contiguous.c \
    MPI_Type_create_darray.c \
    MPI_Type_create_indexed_block.c \
    MPI_Type_create_subarray.c \
    MPI_Type_extent.c \
    MPI_Type_free.c \
    MPI_Type_get_contents.c \
    MPI_Type_get_envelope.c \
    MPI_Type_hindexed.c \
    MPI_Type_hvector.c \
    MPI_Type_indexed.c \
    MPI_Type_lb.c \
    MPI_Type_size.c \
    MPI_Type_struct.c \
    MPI_Type_ub.c \
    MPI_Type_vector.c \
    MPI_Unpack.c \
    MPI_Waitsome.c \
    MPI_Waitany.c \
    MPI_Waitall.c \
    PMPI_Waitall.c \
    MPI_Wtick.c \
    MPI_Wtime.c \
    PMPI_File_sync.c \
    PMPI_File_write.c \
    PMPI_File_write_all.c \
    PMPI_File_write_all_begin.c \
    PMPI_File_write_all_end.c \
    PMPI_File_write_at.c \
    PMPI_File_write_at_all.c \
    PMPI_File_write_at_all_begin.c \
    PMPI_File_write_at_all_end.c \
    PMPI_File_write_ordered.c \
    PMPI_File_write_ordered_begin.c \
    PMPI_File_write_ordered_end.c \
    PMPI_File_write_shared.c \
    PMPI_Finalized.c \
    PMPI_Gatherv.c \
    PMPI_Get_count.c \
    PMPI_Get_elements.c \
    PMPI_Get_processor_name.c \
    PMPI_Get_version.c \
    PMPI_Graph_create.c \
    PMPI_Graph_get.c \
    PMPI_Graph_map.c \
    PMPI_Graph_neighbors.c \
    PMPI_Graph_neighbors_count.c \
    PMPI_Graphdims_get.c \
    PMPI_Group_compare.c \
    PMPI_Group_difference.c \
    PMPI_Group_excl.c \
    PMPI_Group_free.c \
    PMPI_Group_incl.c \
    PMPI_Group_intersection.c \
    PMPI_Group_range_excl.c \
    PMPI_Group_range_incl.c \
    PMPI_Group_rank.c \
    PMPI_Group_size.c \
    PMPI_Group_translate_ranks.c \
    PMPI_Group_union.c \
    PMPI_Info_c2f.c \
    PMPI_Info_create.c \
    PMPI_Info_delete.c \
    PMPI_Info_f2c.c \
    PMPI_Info_free.c \
    PMPI_Info_get.c \
    PMPI_Info_get_nkeys.c \
    PMPI_Info_get_nthkey.c \
    PMPI_Info_get_valuelen.c \
    PMPI_Info_set.c \
    PMPI_Init_thread.c \
    PMPI_Initialized.c \
    PMPI_Int2handle.c \
    PMPI_Intercomm_create.c \
    PMPI_Intercomm_merge.c \
    PMPI_Iprobe.c \
    PMPI_Keyval_create.c \
    PMPI_Keyval_free.c \
    PMPI_Pack.c \
    PMPI_Pack_size.c \
    PMPI_Pcontrol.c \
    PMPI_Probe.c \
    PMPI_Recv_init.c \
    PMPI_Reduce_scatter.c \
    PMPI_Request_c2f.c \
    PMPI_Request_free.c \
    PMPI_Rsend_init.c \
    PMPI_Scan.c \
    PMPI_Scatterv.c \
    PMPI_Send_init.c \
    PMPI_Sendrecv_replace.c \
    PMPI_Ssend_init.c \
    PMPI_Start.c \
    PMPI_Startall.c \
    PMPI_Status_c2f.c \
    PMPI_Status_set_cancelled.c \
    PMPI_Status_set_elements.c \
    PMPI_Test.c \
    PMPI_Test_cancelled.c \
    PMPI_Testall.c \
    PMPI_Testany.c \
    PMPI_Testsome.c \
    PMPI_Topo_test.c \
    PMPI_Type_commit.c \
    PMPI_Type_contiguous.c \
    PMPI_Type_create_darray.c \
    PMPI_Type_create_indexed_block.c \
    PMPI_Type_create_subarray.c \
    PMPI_Type_extent.c \
    PMPI_Type_free.c \
    PMPI_Type_get_contents.c \
    PMPI_Type_get_envelope.c \
    PMPI_Type_hindexed.c \
    PMPI_Type_hvector.c \
    PMPI_Type_indexed.c \
    PMPI_Type_lb.c \
    PMPI_Type_size.c \
    PMPI_Type_struct.c \
    PMPI_Type_ub.c \
    PMPI_Type_vector.c \
    PMPI_Unpack.c \
    PMPI_Waitsome.c \
    PMPI_Wtime.c \
    PMPI_Wtick.c


#We don't need any Fortran-related code
#mpif.h: mpi.h fortran_parameterize
#	$(RM) -f mpif.h
#	./fortran_parameterize mpi.h > mpif.h

#F_MPI_Init.c: MPI_Init.c fortran_grinder
#	$(RM) -f F_MPI_Init.c
#	./fortran_grinder $(srcdir)/MPI_Init.c '(MPI_FORTRAN_INTEGER* status) { int istatus = MPI_Init(0,0); *status = istatus; }' > F_MPI_Init.c

#F_MPI_Finalize.c: MPI_Finalize.c fortran_grinder
#	$(RM) -f F_MPI_Finalize.c
#	./fortran_grinder $(srcdir)/MPI_Finalize.c '(MPI_FORTRAN_INTEGER* status) { int istatus = MPI_Finalize(); *status = istatus; }' > F_MPI_Finalize.c

#F_MPI_Reduce.c: MPI_Reduce.c fortran_grinder
#	$(RM) -f F_MPI_Reduce.c
#	./fortran_grinder $(srcdir)/MPI_Reduce.c > F_MPI_Reduce.c

$(libsimpi_a_OBJECTS): mpi.h mpi_implementation.h _MPI_UTILITY.h mpi_prototypes.h mpi_profile.h
#$(mpicc_OBJECTS): mpi_path.h
#$(mpif77_OBJECTS): mpi_path.h mpicc.c
#$(fortran_grinder_OBJECTS): mpi_config.h
#$(fortran_parameterize_OBJECTS): mpi_config.h

TESTS = smoketest \
	regressiontest

#EXTRA_DIST = boot mpi_implementation.h mpi_path.h _MPI_UTILITY.h mpi_prototypes.h mpi_profile.h
EXTRA_DIST = boot mpi_implementation.h _MPI_UTILITY.h mpi_prototypes.h mpi_profile.h

#smoketest: smoketest.c ./mpicc
#	./mpicc -o smoketest smoketest.c
