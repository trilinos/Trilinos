ml_aggregate.o : ../src/ml_aggregate.c ../src/ml_aggregate.h \
        ../src/ml_comm.h ../src/ml_defs.h ../src/ml_lapack.h ../src/ml_memory.h ../src/ml_operator.h
	$(CC) -c $(CFLAGS) ../src/ml_aggregate.c -o $@

ml_aztec_utils.o : ../src/ml_aztec_utils.c ../src/ml_aggregate.h \
        ../src/ml_aztec_utils.h ../src/ml_bicgstabl.h ../src/ml_cg.h ../src/ml_comm.h ../src/ml_comminfoagx.h ../src/ml_defs.h ../src/ml_elementagx.h ../src/ml_ggraph.h ../src/ml_gmres.h ../src/ml_gridagx.h ../src/ml_gridfunc.h ../src/ml_include.h ../src/ml_intlist.h ../src/ml_krylov.h ../src/ml_mat_formats.h ../src/ml_memory.h ../src/ml_op_utils.h ../src/ml_operatoragx.h ../src/ml_rap.h ../src/ml_rbm.h ../src/ml_smooth_prolongator.h ../src/ml_solver.h ../src/ml_struct.h ../src/ml_utils.h ../src/ml_vec.h ../src/mli_solver.h
	$(CC) -c $(CFLAGS) ../src/ml_aztec_utils.c -o $@

ml_bdrypts.o : ../src/ml_bdrypts.c ../src/ml_bdrypts.h ../src/ml_defs.h \
        ../src/ml_memory.h
	$(CC) -c $(CFLAGS) ../src/ml_bdrypts.c -o $@

ml_check.o : ../src/ml_check.c ../src/ml_1level.h ../src/ml_bdrypts.h \
        ../src/ml_check.h ../src/ml_comm.h ../src/ml_comminfoagx.h ../src/ml_comminfoop.h ../src/ml_csolve.h ../src/ml_defs.h ../src/ml_grid.h ../src/ml_gridfunc.h ../src/ml_krylov.h ../src/ml_mapper.h ../src/ml_mat_formats.h ../src/ml_memory.h ../src/ml_operator.h ../src/ml_operatoragx.h ../src/ml_rap.h ../src/ml_smoother.h ../src/ml_solver.h ../src/ml_struct.h ../src/ml_utils.h ../src/ml_vec.h
	$(CC) -c $(CFLAGS) ../src/ml_check.c -o $@

ml_comm.o : ../src/ml_comm.c ../src/ml_comm.h ../src/ml_defs.h \
        ../src/ml_memory.h
	$(CC) -c $(CFLAGS) ../src/ml_comm.c -o $@

ml_comminfoagx.o : ../src/ml_comminfoagx.c ../src/ml_comminfoagx.h \
        ../src/ml_defs.h ../src/ml_memory.h
	$(CC) -c $(CFLAGS) ../src/ml_comminfoagx.c -o $@

ml_comminfoop.o : ../src/ml_comminfoop.c ../src/ml_aggregate.h \
        ../src/ml_aztec_utils.h ../src/ml_bicgstabl.h ../src/ml_cg.h ../src/ml_comm.h ../src/ml_comminfoagx.h ../src/ml_defs.h ../src/ml_elementagx.h ../src/ml_ggraph.h ../src/ml_gmres.h ../src/ml_gridagx.h ../src/ml_gridfunc.h ../src/ml_include.h ../src/ml_intlist.h ../src/ml_krylov.h ../src/ml_mat_formats.h ../src/ml_memory.h ../src/ml_op_utils.h ../src/ml_operatoragx.h ../src/ml_rap.h ../src/ml_rbm.h ../src/ml_smooth_prolongator.h ../src/ml_solver.h ../src/ml_struct.h ../src/ml_utils.h ../src/ml_vec.h ../src/mli_solver.h
	$(CC) -c $(CFLAGS) ../src/ml_comminfoop.c -o $@

ml_csolve.o : ../src/ml_csolve.c ../src/ml_1level.h ../src/ml_csolve.h \
        ../src/ml_defs.h ../src/ml_memory.h
	$(CC) -c $(CFLAGS) ../src/ml_csolve.c -o $@

ml_elementagx.o : ../src/ml_elementagx.c ../src/ml_elementagx.h \
        ../src/ml_memory.h
	$(CC) -c $(CFLAGS) ../src/ml_elementagx.c -o $@

ml_exch_row.o : ../src/ml_exch_row.c ../src/ml_1level.h ../src/ml_bdrypts.h \
        ../src/ml_comm.h ../src/ml_comminfoop.h ../src/ml_csolve.h ../src/ml_defs.h ../src/ml_grid.h ../src/ml_gridfunc.h ../src/ml_krylov.h ../src/ml_mapper.h ../src/ml_mat_formats.h ../src/ml_operator.h ../src/ml_operatoragx.h ../src/ml_rap.h ../src/ml_smoother.h ../src/ml_solver.h ../src/ml_struct.h ../src/ml_utils.h ../src/ml_vec.h
	$(CC) -c $(CFLAGS) ../src/ml_exch_row.c -o $@

ml_get_basis.o : ../src/ml_get_basis.c ../src/ml_defs.h ../src/ml_gridfunc.h \
        ../src/ml_memory.h
	$(CC) -c $(CFLAGS) ../src/ml_get_basis.c -o $@

ml_agg_MIS.o : ../src/ml_agg_MIS.c ../src/ml_aggregate.h \
        ../src/ml_comm.h ../src/ml_comminfoop.h ../src/ml_defs.h ../src/ml_ggraph.h ../src/ml_gridfunc.h ../src/ml_lapack.h ../src/ml_mat_formats.h ../src/ml_memory.h ../src/ml_operator.h
	$(CC) -c $(CFLAGS) ../src/ml_agg_MIS.c -o $@

ml_grid.o : ../src/ml_grid.c ../src/ml_defs.h ../src/ml_grid.h \
        ../src/ml_gridfunc.h ../src/ml_memory.h
	$(CC) -c $(CFLAGS) ../src/ml_grid.c -o $@

ml_gridagx.o : ../src/ml_gridagx.c ../src/ml_defs.h ../src/ml_elementagx.h \
        ../src/ml_gridagx.h ../src/ml_intlist.h ../src/ml_memory.h
	$(CC) -c $(CFLAGS) ../src/ml_gridagx.c -o $@

ml_gridfunc.o : ../src/ml_gridfunc.c ../src/ml_defs.h ../src/ml_gridfunc.h \
        ../src/ml_memory.h
	$(CC) -c $(CFLAGS) ../src/ml_gridfunc.c -o $@

ml_intlist.o : ../src/ml_intlist.c ../src/ml_defs.h ../src/ml_intlist.h \
        ../src/ml_memory.h
	$(CC) -c $(CFLAGS) ../src/ml_intlist.c -o $@

ml_lapack.o : ../src/ml_lapack.c ../src/ml_defs.h
	$(CC) -c $(CFLAGS) ../src/ml_lapack.c -o $@

ml_mapper.o : ../src/ml_mapper.c ../src/ml_defs.h ../src/ml_mapper.h \
        ../src/ml_memory.h
	$(CC) -c $(CFLAGS) ../src/ml_mapper.c -o $@

ml_mat_formats.o : ../src/ml_mat_formats.c ../src/ml_aggregate.h \
        ../src/ml_aztec_utils.h ../src/ml_bicgstabl.h ../src/ml_cg.h ../src/ml_comm.h ../src/ml_comminfoagx.h ../src/ml_defs.h ../src/ml_elementagx.h ../src/ml_ggraph.h ../src/ml_gmres.h ../src/ml_gridagx.h ../src/ml_gridfunc.h ../src/ml_include.h ../src/ml_intlist.h ../src/ml_krylov.h ../src/ml_mat_formats.h ../src/ml_memory.h ../src/ml_op_utils.h ../src/ml_operatoragx.h ../src/ml_rap.h ../src/ml_rbm.h ../src/ml_smooth_prolongator.h ../src/ml_solver.h ../src/ml_struct.h ../src/ml_utils.h ../src/ml_vec.h ../src/mli_solver.h
	$(CC) -c $(CFLAGS) ../src/ml_mat_formats.c -o $@

ml_matmat_mult.o : ../src/ml_matmat_mult.c ../src/ml_1level.h \
        ../src/ml_bdrypts.h ../src/ml_comm.h ../src/ml_comminfoagx.h ../src/ml_comminfoop.h ../src/ml_csolve.h ../src/ml_defs.h ../src/ml_grid.h ../src/ml_gridfunc.h ../src/ml_krylov.h ../src/ml_mapper.h ../src/ml_mat_formats.h ../src/ml_operator.h ../src/ml_operatoragx.h ../src/ml_rap.h ../src/ml_smoother.h ../src/ml_solver.h ../src/ml_struct.h ../src/ml_utils.h ../src/ml_vec.h
	$(CC) -c $(CFLAGS) ../src/ml_matmat_mult.c -o $@

ml_memory.o : ../src/ml_memory.c ../src/ml_comm.h ../src/ml_defs.h \
        ../src/ml_memory.h
	$(CC) -c $(CFLAGS) ../src/ml_memory.c -o $@

ml_op_utils.o : ../src/ml_op_utils.c ../src/ml_aggregate.h \
        ../src/ml_aztec_utils.h ../src/ml_bicgstabl.h ../src/ml_cg.h ../src/ml_comm.h ../src/ml_comminfoagx.h ../src/ml_defs.h ../src/ml_elementagx.h ../src/ml_ggraph.h ../src/ml_gmres.h ../src/ml_gridagx.h ../src/ml_gridfunc.h ../src/ml_include.h ../src/ml_intlist.h ../src/ml_krylov.h ../src/ml_mat_formats.h ../src/ml_memory.h ../src/ml_op_utils.h ../src/ml_operatoragx.h ../src/ml_rap.h ../src/ml_rbm.h ../src/ml_smooth_prolongator.h ../src/ml_solver.h ../src/ml_struct.h ../src/ml_utils.h ../src/ml_vec.h ../src/mli_solver.h
	$(CC) -c $(CFLAGS) ../src/ml_op_utils.c -o $@

ml_operator.o : ../src/ml_operator.c ../src/ml_1level.h ../src/ml_bdrypts.h \
        ../src/ml_defs.h ../src/ml_memory.h ../src/ml_operator.h ../src/ml_operatoragx.h ../src/ml_vec.h
	$(CC) -c $(CFLAGS) ../src/ml_operator.c -o $@

ml_operatoragx.o : ../src/ml_operatoragx.c ../src/ml_comm.h \
        ../src/ml_comminfoagx.h ../src/ml_memory.h ../src/ml_operatoragx.h ../src/ml_struct.h
	$(CC) -c $(CFLAGS) ../src/ml_operatoragx.c -o $@

ml_rap.o : ../src/ml_rap.c ../src/ml_1level.h ../src/ml_bdrypts.h \
        ../src/ml_comm.h ../src/ml_comminfoagx.h ../src/ml_comminfoop.h ../src/ml_csolve.h ../src/ml_defs.h ../src/ml_grid.h ../src/ml_gridfunc.h ../src/ml_krylov.h ../src/ml_mapper.h ../src/ml_mat_formats.h ../src/ml_operator.h ../src/ml_operatoragx.h ../src/ml_rap.h ../src/ml_smoother.h ../src/ml_solver.h ../src/ml_struct.h ../src/ml_utils.h ../src/ml_vec.h
	$(CC) -c $(CFLAGS) ../src/ml_rap.c -o $@

ml_rap_utils.o : ../src/ml_rap_utils.c ../src/ml_1level.h \
        ../src/ml_bdrypts.h ../src/ml_comm.h ../src/ml_comminfoop.h ../src/ml_csolve.h ../src/ml_defs.h ../src/ml_grid.h ../src/ml_gridfunc.h ../src/ml_krylov.h ../src/ml_mapper.h ../src/ml_mat_formats.h ../src/ml_operator.h ../src/ml_operatoragx.h ../src/ml_rap.h ../src/ml_smoother.h ../src/ml_solver.h ../src/ml_struct.h ../src/ml_utils.h ../src/ml_vec.h
	$(CC) -c $(CFLAGS) ../src/ml_rap_utils.c -o $@

ml_setup.o : ../src/ml_setup.c ../src/ml_1level.h ../src/ml_bdrypts.h \
        ../src/ml_comm.h ../src/ml_comminfoagx.h ../src/ml_comminfoop.h ../src/ml_csolve.h ../src/ml_defs.h ../src/ml_elementagx.h ../src/ml_grid.h ../src/ml_gridagx.h ../src/ml_gridfunc.h ../src/ml_intlist.h ../src/ml_krylov.h ../src/ml_mapper.h ../src/ml_mat_formats.h ../src/ml_memory.h ../src/ml_operator.h ../src/ml_operatoragx.h ../src/ml_rap.h ../src/ml_smoother.h ../src/ml_solver.h ../src/ml_struct.h ../src/ml_utils.h ../src/ml_vec.h
	$(CC) -c $(CFLAGS) ../src/ml_setup.c -o $@

ml_smooth_prolongator.o : ../src/ml_smooth_prolongator.c ../src/ml_aggregate.h \
        ../src/ml_aztec_utils.h ../src/ml_bicgstabl.h ../src/ml_cg.h ../src/ml_comm.h ../src/ml_comminfoagx.h ../src/ml_defs.h ../src/ml_elementagx.h ../src/ml_ggraph.h ../src/ml_gmres.h ../src/ml_gridagx.h ../src/ml_gridfunc.h ../src/ml_include.h ../src/ml_intlist.h ../src/ml_krylov.h ../src/ml_mat_formats.h ../src/ml_memory.h ../src/ml_op_utils.h ../src/ml_operator.h ../src/ml_operatoragx.h ../src/ml_rap.h ../src/ml_rbm.h ../src/ml_smooth_prolongator.h ../src/ml_solver.h ../src/ml_struct.h ../src/ml_utils.h ../src/ml_vec.h ../src/mli_solver.h
	$(CC) -c $(CFLAGS) ../src/ml_smooth_prolongator.c -o $@

ml_smoother.o : ../src/ml_smoother.c ../src/ml_1level.h ../src/ml_defs.h \
        ../src/ml_include.h ../src/ml_lapack.h ../src/ml_memory.h ../src/ml_smoother.h
	$(CC) -c $(CFLAGS) ../src/ml_smoother.c -o $@

ml_solver.o : ../src/ml_solver.c ../src/ml_aggregate.h \
        ../src/ml_aztec_utils.h ../src/ml_bicgstabl.h ../src/ml_cg.h ../src/ml_comm.h ../src/ml_comminfoagx.h ../src/ml_defs.h ../src/ml_elementagx.h ../src/ml_ggraph.h ../src/ml_gmres.h ../src/ml_gridagx.h ../src/ml_gridfunc.h ../src/ml_include.h ../src/ml_intlist.h ../src/ml_krylov.h ../src/ml_mat_formats.h ../src/ml_memory.h ../src/ml_op_utils.h ../src/ml_operatoragx.h ../src/ml_rap.h ../src/ml_rbm.h ../src/ml_smooth_prolongator.h ../src/ml_solver.h ../src/ml_struct.h ../src/ml_utils.h ../src/ml_vec.h ../src/mli_solver.h
	$(CC) -c $(CFLAGS) ../src/ml_solver.c -o $@

ml_struct.o : ../src/ml_struct.c ../src/ml_1level.h ../src/ml_aggregate.h \
        ../src/ml_bdrypts.h ../src/ml_comm.h ../src/ml_comminfoop.h ../src/ml_csolve.h ../src/ml_defs.h ../src/ml_grid.h ../src/ml_gridfunc.h ../src/ml_krylov.h ../src/ml_lapack.h ../src/ml_mapper.h ../src/ml_mat_formats.h ../src/ml_operator.h ../src/ml_operatoragx.h ../src/ml_rap.h ../src/ml_smooth_prolongator.h ../src/ml_smoother.h ../src/ml_solver.h ../src/ml_struct.h ../src/ml_utils.h ../src/ml_vec.h
	$(CC) -c $(CFLAGS) ../src/ml_struct.c -o $@

ml_superlu.o : ../src/ml_superlu.c ../src/ml_1level.h ../src/ml_bdrypts.h \
        ../src/ml_comm.h ../src/ml_comminfoop.h ../src/ml_csolve.h ../src/ml_defs.h ../src/ml_grid.h ../src/ml_gridfunc.h ../src/ml_krylov.h ../src/ml_mapper.h ../src/ml_mat_formats.h ../src/ml_memory.h ../src/ml_operator.h ../src/ml_operatoragx.h ../src/ml_rap.h ../src/ml_smoother.h ../src/ml_solver.h ../src/ml_struct.h ../src/ml_utils.h ../src/ml_vec.h
	$(CC) -c $(CFLAGS) ../src/ml_superlu.c -o $@

ml_utils.o : ../src/ml_utils.c ../src/ml_comm.h ../src/ml_defs.h \
        ../src/ml_memory.h ../src/ml_struct.h ../src/ml_utils.h
	$(CC) -c $(CFLAGS) ../src/ml_utils.c -o $@

ml_vec.o : ../src/ml_vec.c ../src/ml_comm.h ../src/ml_defs.h \
        ../src/ml_memory.h ../src/ml_vec.h
	$(CC) -c $(CFLAGS) ../src/ml_vec.c -o $@

seg_precond.o : ../src/seg_precond.c ../src/ml_1level.h ../src/ml_aggregate.h \
        ../src/ml_aztec_utils.h ../src/ml_bdrypts.h ../src/ml_bicgstabl.h ../src/ml_cg.h ../src/ml_comm.h ../src/ml_comminfoagx.h ../src/ml_comminfoop.h ../src/ml_csolve.h ../src/ml_defs.h ../src/ml_elementagx.h ../src/ml_ggraph.h ../src/ml_gmres.h ../src/ml_grid.h ../src/ml_gridagx.h ../src/ml_gridfunc.h ../src/ml_include.h ../src/ml_intlist.h ../src/ml_krylov.h ../src/ml_mapper.h ../src/ml_mat_formats.h ../src/ml_memory.h ../src/ml_op_utils.h ../src/ml_operator.h ../src/ml_operatoragx.h ../src/ml_rap.h ../src/ml_rbm.h ../src/ml_smooth_prolongator.h ../src/ml_smoother.h ../src/ml_solver.h ../src/ml_struct.h ../src/ml_utils.h ../src/ml_vec.h ../src/mli_solver.h ../src/seg_precond.h
	$(CC) -c $(CFLAGS) ../src/seg_precond.c -o $@

ml_ggraph.o : ../src/ml_ggraph.c ../src/ml_comm.h ../src/ml_comminfoop.h \
        ../src/ml_defs.h ../src/ml_ggraph.h ../src/ml_gridfunc.h ../src/ml_mat_formats.h ../src/ml_memory.h ../src/ml_operator.h
	$(CC) -c $(CFLAGS) ../src/ml_ggraph.c -o $@

ml_krylov.o : ../src/ml_krylov.c ../src/ml_comm.h ../src/ml_krylov.h \
        ../src/ml_operator.h
	$(CC) -c $(CFLAGS) ../src/ml_krylov.c -o $@

ml_cg.o : ../src/ml_cg.c ../src/ml_cg.h ../src/ml_krylov.h
	$(CC) -c $(CFLAGS) ../src/ml_cg.c -o $@

driver.o : ../src/driver.c ../src/ml_aggregate.h ../src/ml_aztec_utils.h \
        ../src/ml_bicgstabl.h ../src/ml_cg.h ../src/ml_comm.h ../src/ml_comminfoagx.h ../src/ml_defs.h ../src/ml_elementagx.h ../src/ml_ggraph.h ../src/ml_gmres.h ../src/ml_gridagx.h ../src/ml_gridfunc.h ../src/ml_include.h ../src/ml_intlist.h ../src/ml_krylov.h ../src/ml_mat_formats.h ../src/ml_memory.h ../src/ml_op_utils.h ../src/ml_operatoragx.h ../src/ml_rap.h ../src/ml_rbm.h ../src/ml_smooth_prolongator.h ../src/ml_solver.h ../src/ml_struct.h ../src/ml_utils.h ../src/ml_vec.h ../src/mli_solver.h
	$(CC) -c $(CFLAGS) ../src/driver.c -o $@

ml_xxt.o : ../src/ml_xxt.c ../src/ml_aggregate.h ../src/ml_aztec_utils.h \
        ../src/ml_bicgstabl.h ../src/ml_cg.h ../src/ml_comm.h ../src/ml_comminfoagx.h ../src/ml_defs.h ../src/ml_elementagx.h ../src/ml_ggraph.h ../src/ml_gmres.h ../src/ml_gridagx.h ../src/ml_gridfunc.h ../src/ml_include.h ../src/ml_intlist.h ../src/ml_krylov.h ../src/ml_mat_formats.h ../src/ml_memory.h ../src/ml_op_utils.h ../src/ml_operatoragx.h ../src/ml_rap.h ../src/ml_rbm.h ../src/ml_smooth_prolongator.h ../src/ml_solver.h ../src/ml_struct.h ../src/ml_utils.h ../src/ml_vec.h ../src/mli_solver.h
	$(CC) -c $(CFLAGS) ../src/ml_xxt.c -o $@

ml_bicgstabl.o : ../src/ml_bicgstabl.c ../src/ml_bicgstabl.h \
        ../src/ml_krylov.h
	$(CC) -c $(CFLAGS) ../src/ml_bicgstabl.c -o $@

ml_gmres.o : ../src/ml_gmres.c ../src/ml_gmres.h ../src/ml_krylov.h
	$(CC) -c $(CFLAGS) ../src/ml_gmres.c -o $@

ml_pde.o : ../src/ml_pde.c ../src/ml_comm.h ../src/ml_defs.h \
        ../src/ml_include.h ../src/ml_memory.h ../src/ml_pde.h ../src/mli_solver.h
	$(CC) -c $(CFLAGS) ../src/ml_pde.c -o $@

mli_solver.o : ../src/mli_solver.c ../src/ml_include.h ../src/mli_solver.h
	$(CC) -c $(CFLAGS) ../src/mli_solver.c -o $@

ml_rbm.o : ../src/ml_rbm.c ../src/ml_rbm.h
	$(CC) -c $(CFLAGS) ../src/ml_rbm.c -o $@

convertSund2AZdatafile.o : ../examples/convertSund2AZdatafile.c
	$(CC) -c $(CFLAGS) ../examples/convertSund2AZdatafile.c -o $@

ml_example1d.o : ../examples/ml_example1d.c ../src/ml_aggregate.h \
        ../src/ml_aztec_utils.h ../src/ml_bicgstabl.h ../src/ml_cg.h ../src/ml_comm.h ../src/ml_comminfoagx.h ../src/ml_defs.h ../src/ml_elementagx.h ../src/ml_ggraph.h ../src/ml_gmres.h ../src/ml_gridagx.h ../src/ml_gridfunc.h ../src/ml_include.h ../src/ml_intlist.h ../src/ml_krylov.h ../src/ml_mat_formats.h ../src/ml_memory.h ../src/ml_op_utils.h ../src/ml_operatoragx.h ../src/ml_rap.h ../src/ml_rbm.h ../src/ml_smooth_prolongator.h ../src/ml_solver.h ../src/ml_struct.h ../src/ml_utils.h ../src/ml_vec.h ../src/mli_solver.h
	$(CC) -c $(CFLAGS) ../examples/ml_example1d.c -o $@

ml_example1dGS.o : ../examples/ml_example1dGS.c ../src/ml_aggregate.h \
        ../src/ml_aztec_utils.h ../src/ml_bicgstabl.h ../src/ml_cg.h ../src/ml_comm.h ../src/ml_comminfoagx.h ../src/ml_defs.h ../src/ml_elementagx.h ../src/ml_ggraph.h ../src/ml_gmres.h ../src/ml_gridagx.h ../src/ml_gridfunc.h ../src/ml_include.h ../src/ml_intlist.h ../src/ml_krylov.h ../src/ml_mat_formats.h ../src/ml_memory.h ../src/ml_op_utils.h ../src/ml_operatoragx.h ../src/ml_rap.h ../src/ml_rbm.h ../src/ml_smooth_prolongator.h ../src/ml_solver.h ../src/ml_struct.h ../src/ml_utils.h ../src/ml_vec.h ../src/mli_solver.h
	$(CC) -c $(CFLAGS) ../examples/ml_example1dGS.c -o $@

ml_example2d.o : ../examples/ml_example2d.c ../src/ml_aggregate.h \
        ../src/ml_aztec_utils.h ../src/ml_bicgstabl.h ../src/ml_cg.h ../src/ml_comm.h ../src/ml_comminfoagx.h ../src/ml_defs.h ../src/ml_elementagx.h ../src/ml_ggraph.h ../src/ml_gmres.h ../src/ml_gridagx.h ../src/ml_gridfunc.h ../src/ml_include.h ../src/ml_intlist.h ../src/ml_krylov.h ../src/ml_mat_formats.h ../src/ml_memory.h ../src/ml_op_utils.h ../src/ml_operatoragx.h ../src/ml_rap.h ../src/ml_rbm.h ../src/ml_smooth_prolongator.h ../src/ml_solver.h ../src/ml_struct.h ../src/ml_utils.h ../src/ml_vec.h ../src/mli_solver.h
	$(CC) -c $(CFLAGS) ../examples/ml_example2d.c -o $@

ml_example3d.o : ../examples/ml_example3d.c ../src/ml_aggregate.h \
        ../src/ml_aztec_utils.h ../src/ml_bicgstabl.h ../src/ml_cg.h ../src/ml_comm.h ../src/ml_comminfoagx.h ../src/ml_defs.h ../src/ml_elementagx.h ../src/ml_ggraph.h ../src/ml_gmres.h ../src/ml_gridagx.h ../src/ml_gridfunc.h ../src/ml_include.h ../src/ml_intlist.h ../src/ml_krylov.h ../src/ml_mat_formats.h ../src/ml_memory.h ../src/ml_op_utils.h ../src/ml_operatoragx.h ../src/ml_rap.h ../src/ml_rbm.h ../src/ml_smooth_prolongator.h ../src/ml_solver.h ../src/ml_struct.h ../src/ml_utils.h ../src/ml_vec.h ../src/mli_solver.h
	$(CC) -c $(CFLAGS) ../examples/ml_example3d.c -o $@

ml_readex.o : ../examples/ml_readex.c ../src/ml_aggregate.h \
        ../src/ml_aztec_utils.h ../src/ml_bicgstabl.h ../src/ml_cg.h ../src/ml_comm.h ../src/ml_comminfoagx.h ../src/ml_defs.h ../src/ml_elementagx.h ../src/ml_ggraph.h ../src/ml_gmres.h ../src/ml_gridagx.h ../src/ml_gridfunc.h ../src/ml_include.h ../src/ml_intlist.h ../src/ml_krylov.h ../src/ml_mat_formats.h ../src/ml_memory.h ../src/ml_op_utils.h ../src/ml_operatoragx.h ../src/ml_rap.h ../src/ml_rbm.h ../src/ml_smooth_prolongator.h ../src/ml_solver.h ../src/ml_struct.h ../src/ml_utils.h ../src/ml_vec.h ../src/mli_solver.h
	$(CC) -c $(CFLAGS) ../examples/ml_readex.c -o $@

seg_readex.o : ../examples/seg_readex.c ../src/ml_aggregate.h \
        ../src/ml_aztec_utils.h ../src/ml_bicgstabl.h ../src/ml_cg.h ../src/ml_comm.h ../src/ml_comminfoagx.h ../src/ml_defs.h ../src/ml_elementagx.h ../src/ml_ggraph.h ../src/ml_gmres.h ../src/ml_gridagx.h ../src/ml_gridfunc.h ../src/ml_include.h ../src/ml_intlist.h ../src/ml_krylov.h ../src/ml_mat_formats.h ../src/ml_memory.h ../src/ml_op_utils.h ../src/ml_operatoragx.h ../src/ml_rap.h ../src/ml_rbm.h ../src/ml_smooth_prolongator.h ../src/ml_solver.h ../src/ml_struct.h ../src/ml_utils.h ../src/ml_vec.h ../src/mli_solver.h
	$(CC) -c $(CFLAGS) ../examples/seg_readex.c -o $@

ml_ex1d.o : ../examples/ml_ex1d.c ../src/ml_aggregate.h \
        ../src/ml_aztec_utils.h ../src/ml_bicgstabl.h ../src/ml_cg.h ../src/ml_comm.h ../src/ml_comminfoagx.h ../src/ml_defs.h ../src/ml_elementagx.h ../src/ml_ggraph.h ../src/ml_gmres.h ../src/ml_gridagx.h ../src/ml_gridfunc.h ../src/ml_include.h ../src/ml_intlist.h ../src/ml_krylov.h ../src/ml_mat_formats.h ../src/ml_memory.h ../src/ml_op_utils.h ../src/ml_operatoragx.h ../src/ml_rap.h ../src/ml_rbm.h ../src/ml_smooth_prolongator.h ../src/ml_solver.h ../src/ml_struct.h ../src/ml_utils.h ../src/ml_vec.h ../src/mli_solver.h
	$(CC) -c $(CFLAGS) ../examples/ml_ex1d.c -o $@

ml_star2d.o : ../examples/ml_star2d.c ../src/ml_aggregate.h \
        ../src/ml_aztec_utils.h ../src/ml_bicgstabl.h ../src/ml_cg.h ../src/ml_comm.h ../src/ml_comminfoagx.h ../src/ml_defs.h ../src/ml_elementagx.h ../src/ml_ggraph.h ../src/ml_gmres.h ../src/ml_gridagx.h ../src/ml_gridfunc.h ../src/ml_include.h ../src/ml_intlist.h ../src/ml_krylov.h ../src/ml_mat_formats.h ../src/ml_memory.h ../src/ml_op_utils.h ../src/ml_operatoragx.h ../src/ml_rap.h ../src/ml_rbm.h ../src/ml_smooth_prolongator.h ../src/ml_solver.h ../src/ml_struct.h ../src/ml_utils.h ../src/ml_vec.h ../src/mli_solver.h
	$(CC) -c $(CFLAGS) ../examples/ml_star2d.c -o $@

oldml_readex.o : ../examples/oldml_readex.c ../src/ml_aggregate.h \
        ../src/ml_aztec_utils.h ../src/ml_bicgstabl.h ../src/ml_cg.h ../src/ml_comm.h ../src/ml_comminfoagx.h ../src/ml_defs.h ../src/ml_elementagx.h ../src/ml_ggraph.h ../src/ml_gmres.h ../src/ml_gridagx.h ../src/ml_gridfunc.h ../src/ml_include.h ../src/ml_intlist.h ../src/ml_krylov.h ../src/ml_mat_formats.h ../src/ml_memory.h ../src/ml_op_utils.h ../src/ml_operatoragx.h ../src/ml_rap.h ../src/ml_rbm.h ../src/ml_smooth_prolongator.h ../src/ml_solver.h ../src/ml_struct.h ../src/ml_utils.h ../src/ml_vec.h ../src/mli_solver.h
	$(CC) -c $(CFLAGS) ../examples/oldml_readex.c -o $@

new_readex.o : ../examples/new_readex.c ../src/ml_aggregate.h \
        ../src/ml_aztec_utils.h ../src/ml_bicgstabl.h ../src/ml_cg.h ../src/ml_comm.h ../src/ml_comminfoagx.h ../src/ml_defs.h ../src/ml_elementagx.h ../src/ml_ggraph.h ../src/ml_gmres.h ../src/ml_gridagx.h ../src/ml_gridfunc.h ../src/ml_include.h ../src/ml_intlist.h ../src/ml_krylov.h ../src/ml_mat_formats.h ../src/ml_memory.h ../src/ml_op_utils.h ../src/ml_operatoragx.h ../src/ml_rap.h ../src/ml_rbm.h ../src/ml_smooth_prolongator.h ../src/ml_solver.h ../src/ml_struct.h ../src/ml_utils.h ../src/ml_vec.h ../src/mli_solver.h
	$(CC) -c $(CFLAGS) ../examples/new_readex.c -o $@

