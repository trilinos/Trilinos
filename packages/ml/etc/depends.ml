ml_agg_MIS.o : ../flat/ml_agg_MIS.c ../flat/ml_aggregate.h \
        ../flat/ml_comm.h ../flat/ml_comminfoop.h ../flat/ml_defs.h ../flat/ml_ggraph.h ../flat/ml_gridfunc.h ../flat/ml_lapack.h ../flat/ml_mat_formats.h ../flat/ml_memory.h ../flat/ml_operator.h
	$(CC) -c $(CFLAGS) ../flat/ml_agg_MIS.c -o $@

ml_agg_coupled.o : ../flat/ml_agg_coupled.c ../flat/ml_aggregate.h \
        ../flat/ml_comm.h ../flat/ml_defs.h ../flat/ml_lapack.h ../flat/ml_memory.h ../flat/ml_operator.h
	$(CC) -c $(CFLAGS) ../flat/ml_agg_coupled.c -o $@

ml_agg_dd.o : ../flat/ml_agg_dd.c ../flat/ml_aggregate.h \
        ../flat/ml_comm.h ../flat/ml_defs.h ../flat/ml_lapack.h ../flat/ml_memory.h ../flat/ml_operator.h
	$(CC) -c $(CFLAGS) ../flat/ml_agg_dd.c -o $@

ml_agg_genP.o : ../flat/ml_agg_genP.c ../flat/ml_agg_genP.h \
        ../flat/ml_aggregate.h ../flat/ml_aztec_utils.h ../flat/ml_bicgstabl.h ../flat/ml_cg.h ../flat/ml_comm.h ../flat/ml_defs.h ../flat/ml_elementagx.h ../flat/ml_gmres.h ../flat/ml_grid.h ../flat/ml_gridagx.h ../flat/ml_gridfunc.h ../flat/ml_include.h ../flat/ml_intlist.h ../flat/ml_krylov.h ../flat/ml_op_utils.h ../flat/ml_operator.h ../flat/ml_operatoragx.h ../flat/ml_pde.h ../flat/ml_solver.h ../flat/ml_struct.h ../flat/ml_vec.h ../flat/mli_solver.h
	$(CC) -c $(CFLAGS) ../flat/ml_agg_genP.c -o $@

ml_agg_uncoupled.o : ../flat/ml_agg_uncoupled.c ../flat/ml_aggregate.h \
        ../flat/ml_comm.h ../flat/ml_defs.h ../flat/ml_lapack.h ../flat/ml_memory.h ../flat/ml_operator.h
	$(CC) -c $(CFLAGS) ../flat/ml_agg_uncoupled.c -o $@

ml_aggregate.o : ../flat/ml_aggregate.c ../flat/ml_aggregate.h \
        ../flat/ml_comm.h ../flat/ml_defs.h ../flat/ml_lapack.h ../flat/ml_memory.h ../flat/ml_operator.h
	$(CC) -c $(CFLAGS) ../flat/ml_aggregate.c -o $@

ml_check.o : ../flat/ml_check.c ../flat/ml_1level.h ../flat/ml_bdrypts.h \
        ../flat/ml_check.h ../flat/ml_comm.h ../flat/ml_comminfoagx.h ../flat/ml_comminfoop.h ../flat/ml_csolve.h ../flat/ml_defs.h ../flat/ml_grid.h ../flat/ml_gridfunc.h ../flat/ml_krylov.h ../flat/ml_mapper.h ../flat/ml_mat_formats.h ../flat/ml_memory.h ../flat/ml_operator.h ../flat/ml_operatoragx.h ../flat/ml_rap.h ../flat/ml_smoother.h ../flat/ml_solver.h ../flat/ml_struct.h ../flat/ml_utils.h ../flat/ml_vec.h
	$(CC) -c $(CFLAGS) ../flat/ml_check.c -o $@

ml_ggraph.o : ../flat/ml_ggraph.c ../flat/ml_comm.h ../flat/ml_comminfoop.h \
        ../flat/ml_defs.h ../flat/ml_ggraph.h ../flat/ml_gridfunc.h ../flat/ml_mat_formats.h ../flat/ml_memory.h ../flat/ml_operator.h
	$(CC) -c $(CFLAGS) ../flat/ml_ggraph.c -o $@

ml_comm.o : ../flat/ml_comm.c ../flat/ml_comm.h ../flat/ml_defs.h \
        ../flat/ml_memory.h
	$(CC) -c $(CFLAGS) ../flat/ml_comm.c -o $@

ml_comminfoagx.o : ../flat/ml_comminfoagx.c ../flat/ml_comminfoagx.h \
        ../flat/ml_defs.h ../flat/ml_memory.h
	$(CC) -c $(CFLAGS) ../flat/ml_comminfoagx.c -o $@

ml_comminfoop.o : ../flat/ml_comminfoop.c ../flat/ml_agg_genP.h \
        ../flat/ml_aggregate.h ../flat/ml_aztec_utils.h ../flat/ml_bicgstabl.h ../flat/ml_cg.h ../flat/ml_comm.h ../flat/ml_defs.h ../flat/ml_elementagx.h ../flat/ml_gmres.h ../flat/ml_grid.h ../flat/ml_gridagx.h ../flat/ml_gridfunc.h ../flat/ml_include.h ../flat/ml_intlist.h ../flat/ml_krylov.h ../flat/ml_operator.h ../flat/ml_operatoragx.h ../flat/ml_pde.h ../flat/ml_solver.h ../flat/ml_struct.h ../flat/ml_vec.h ../flat/mli_solver.h
	$(CC) -c $(CFLAGS) ../flat/ml_comminfoop.c -o $@

ml_exch_row.o : ../flat/ml_exch_row.c ../flat/ml_1level.h \
        ../flat/ml_bdrypts.h ../flat/ml_comm.h ../flat/ml_comminfoop.h ../flat/ml_csolve.h ../flat/ml_defs.h ../flat/ml_grid.h ../flat/ml_gridfunc.h ../flat/ml_krylov.h ../flat/ml_mapper.h ../flat/ml_mat_formats.h ../flat/ml_operator.h ../flat/ml_operatoragx.h ../flat/ml_rap.h ../flat/ml_smoother.h ../flat/ml_solver.h ../flat/ml_struct.h ../flat/ml_utils.h ../flat/ml_vec.h
	$(CC) -c $(CFLAGS) ../flat/ml_exch_row.c -o $@

ml_bdrypts.o : ../flat/ml_bdrypts.c ../flat/ml_bdrypts.h \
        ../flat/ml_defs.h ../flat/ml_memory.h
	$(CC) -c $(CFLAGS) ../flat/ml_bdrypts.c -o $@

ml_elementagx.o : ../flat/ml_elementagx.c ../flat/ml_elementagx.h \
        ../flat/ml_memory.h
	$(CC) -c $(CFLAGS) ../flat/ml_elementagx.c -o $@

ml_get_basis.o : ../flat/ml_get_basis.c ../flat/ml_defs.h \
        ../flat/ml_gridfunc.h ../flat/ml_memory.h
	$(CC) -c $(CFLAGS) ../flat/ml_get_basis.c -o $@

ml_grid.o : ../flat/ml_grid.c ../flat/ml_defs.h ../flat/ml_grid.h \
        ../flat/ml_gridfunc.h ../flat/ml_memory.h
	$(CC) -c $(CFLAGS) ../flat/ml_grid.c -o $@

ml_gridagx.o : ../flat/ml_gridagx.c ../flat/ml_defs.h ../flat/ml_elementagx.h \
        ../flat/ml_gridagx.h ../flat/ml_intlist.h ../flat/ml_memory.h
	$(CC) -c $(CFLAGS) ../flat/ml_gridagx.c -o $@

ml_gridfunc.o : ../flat/ml_gridfunc.c ../flat/ml_defs.h ../flat/ml_gridfunc.h \
        ../flat/ml_memory.h
	$(CC) -c $(CFLAGS) ../flat/ml_gridfunc.c -o $@

ml_mapper.o : ../flat/ml_mapper.c ../flat/ml_defs.h ../flat/ml_mapper.h \
        ../flat/ml_memory.h
	$(CC) -c $(CFLAGS) ../flat/ml_mapper.c -o $@

ml_pde.o : ../flat/ml_pde.c ../flat/ml_1level.h ../flat/ml_bdrypts.h \
        ../flat/ml_comm.h ../flat/ml_comminfoop.h ../flat/ml_csolve.h ../flat/ml_defs.h ../flat/ml_grid.h ../flat/ml_gridfunc.h ../flat/ml_include.h ../flat/ml_krylov.h ../flat/ml_mapper.h ../flat/ml_mat_formats.h ../flat/ml_memory.h ../flat/ml_operator.h ../flat/ml_operatoragx.h ../flat/ml_pde.h ../flat/ml_rap.h ../flat/ml_smoother.h ../flat/ml_solver.h ../flat/ml_struct.h ../flat/ml_utils.h ../flat/ml_vec.h ../flat/mli_solver.h
	$(CC) -c $(CFLAGS) ../flat/ml_pde.c -o $@

ml_setup.o : ../flat/ml_setup.c ../flat/ml_1level.h ../flat/ml_bdrypts.h \
        ../flat/ml_comm.h ../flat/ml_comminfoagx.h ../flat/ml_comminfoop.h ../flat/ml_csolve.h ../flat/ml_defs.h ../flat/ml_elementagx.h ../flat/ml_grid.h ../flat/ml_gridagx.h ../flat/ml_gridfunc.h ../flat/ml_intlist.h ../flat/ml_krylov.h ../flat/ml_mapper.h ../flat/ml_mat_formats.h ../flat/ml_memory.h ../flat/ml_operator.h ../flat/ml_operatoragx.h ../flat/ml_rap.h ../flat/ml_smoother.h ../flat/ml_solver.h ../flat/ml_struct.h ../flat/ml_utils.h ../flat/ml_vec.h
	$(CC) -c $(CFLAGS) ../flat/ml_setup.c -o $@

ml_bicgstabl.o : ../flat/ml_bicgstabl.c ../flat/ml_bicgstabl.h \
        ../flat/ml_krylov.h
	$(CC) -c $(CFLAGS) ../flat/ml_bicgstabl.c -o $@

ml_cg.o : ../flat/ml_cg.c ../flat/ml_cg.h ../flat/ml_krylov.h
	$(CC) -c $(CFLAGS) ../flat/ml_cg.c -o $@

ml_gmres.o : ../flat/ml_gmres.c ../flat/ml_gmres.h ../flat/ml_krylov.h
	$(CC) -c $(CFLAGS) ../flat/ml_gmres.c -o $@

ml_krylov.o : ../flat/ml_krylov.c ../flat/ml_comm.h ../flat/ml_krylov.h \
        ../flat/ml_operator.h
	$(CC) -c $(CFLAGS) ../flat/ml_krylov.c -o $@

driver.o : ../flat/driver.c ../flat/ml_agg_genP.h ../flat/ml_aggregate.h \
        ../flat/ml_aztec_utils.h ../flat/ml_bicgstabl.h ../flat/ml_cg.h ../flat/ml_comm.h ../flat/ml_defs.h ../flat/ml_elementagx.h ../flat/ml_gmres.h ../flat/ml_grid.h ../flat/ml_gridagx.h ../flat/ml_gridfunc.h ../flat/ml_include.h ../flat/ml_intlist.h ../flat/ml_krylov.h ../flat/ml_operator.h ../flat/ml_operatoragx.h ../flat/ml_pde.h ../flat/ml_solver.h ../flat/ml_struct.h ../flat/ml_vec.h ../flat/mli_solver.h
	$(CC) -c $(CFLAGS) ../flat/driver.c -o $@

ml_edmond.o : ../flat/ml_edmond.c ../flat/ml_1level.h ../flat/ml_agg_genP.h \
        ../flat/ml_aggregate.h ../flat/ml_aztec_utils.h ../flat/ml_bdrypts.h ../flat/ml_comm.h ../flat/ml_comminfoop.h ../flat/ml_csolve.h ../flat/ml_defs.h ../flat/ml_grid.h ../flat/ml_gridfunc.h ../flat/ml_include.h ../flat/ml_krylov.h ../flat/ml_lapack.h ../flat/ml_mapper.h ../flat/ml_mat_formats.h ../flat/ml_operator.h ../flat/ml_operatoragx.h ../flat/ml_rap.h ../flat/ml_smoother.h ../flat/ml_solver.h ../flat/ml_struct.h ../flat/ml_utils.h ../flat/ml_vec.h
	$(CC) -c $(CFLAGS) ../flat/ml_edmond.c -o $@

ml_seg_precond.o : ../flat/ml_seg_precond.c ../flat/ml_1level.h \
        ../flat/ml_agg_genP.h ../flat/ml_aggregate.h ../flat/ml_aztec_utils.h ../flat/ml_bdrypts.h ../flat/ml_bicgstabl.h ../flat/ml_cg.h ../flat/ml_comm.h ../flat/ml_comminfoop.h ../flat/ml_csolve.h ../flat/ml_defs.h ../flat/ml_elementagx.h ../flat/ml_gmres.h ../flat/ml_grid.h ../flat/ml_gridagx.h ../flat/ml_gridfunc.h ../flat/ml_include.h ../flat/ml_intlist.h ../flat/ml_krylov.h ../flat/ml_mapper.h ../flat/ml_mat_formats.h ../flat/ml_operator.h ../flat/ml_operatoragx.h ../flat/ml_pde.h ../flat/ml_rap.h ../flat/ml_seg_precond.h ../flat/ml_smoother.h ../flat/ml_solver.h ../flat/ml_struct.h ../flat/ml_utils.h ../flat/ml_vec.h ../flat/mli_solver.h
	$(CC) -c $(CFLAGS) ../flat/ml_seg_precond.c -o $@

ml_smoother_edmond.o : ../flat/ml_smoother_edmond.c ../flat/ml_1level.h \
        ../flat/ml_aztec_utils.h ../flat/ml_defs.h ../flat/ml_include.h ../flat/ml_lapack.h ../flat/ml_memory.h ../flat/ml_smoother.h
	$(CC) -c $(CFLAGS) ../flat/ml_smoother_edmond.c -o $@

ml_struct.o : ../flat/ml_struct.c ../flat/ml_1level.h ../flat/ml_agg_genP.h \
        ../flat/ml_aggregate.h ../flat/ml_aztec_utils.h ../flat/ml_bdrypts.h ../flat/ml_comm.h ../flat/ml_comminfoop.h ../flat/ml_csolve.h ../flat/ml_defs.h ../flat/ml_grid.h ../flat/ml_gridfunc.h ../flat/ml_include.h ../flat/ml_krylov.h ../flat/ml_lapack.h ../flat/ml_mapper.h ../flat/ml_mat_formats.h ../flat/ml_memory.h ../flat/ml_operator.h ../flat/ml_operatoragx.h ../flat/ml_rap.h ../flat/ml_smoother.h ../flat/ml_solver.h ../flat/ml_struct.h ../flat/ml_utils.h ../flat/ml_vec.h
	$(CC) -c $(CFLAGS) ../flat/ml_struct.c -o $@

mli_solver.o : ../flat/mli_solver.c ../flat/ml_1level.h \
        ../flat/ml_bdrypts.h ../flat/ml_comm.h ../flat/ml_comminfoop.h ../flat/ml_csolve.h ../flat/ml_defs.h ../flat/ml_grid.h ../flat/ml_gridfunc.h ../flat/ml_include.h ../flat/ml_krylov.h ../flat/ml_mapper.h ../flat/ml_mat_formats.h ../flat/ml_operator.h ../flat/ml_operatoragx.h ../flat/ml_rap.h ../flat/ml_smoother.h ../flat/ml_solver.h ../flat/ml_struct.h ../flat/ml_utils.h ../flat/ml_vec.h ../flat/mli_solver.h
	$(CC) -c $(CFLAGS) ../flat/mli_solver.c -o $@

ml_mat_formats.o : ../flat/ml_mat_formats.c ../flat/ml_agg_genP.h \
        ../flat/ml_aggregate.h ../flat/ml_aztec_utils.h ../flat/ml_bicgstabl.h ../flat/ml_cg.h ../flat/ml_comm.h ../flat/ml_defs.h ../flat/ml_elementagx.h ../flat/ml_gmres.h ../flat/ml_grid.h ../flat/ml_gridagx.h ../flat/ml_gridfunc.h ../flat/ml_include.h ../flat/ml_intlist.h ../flat/ml_krylov.h ../flat/ml_operator.h ../flat/ml_operatoragx.h ../flat/ml_pde.h ../flat/ml_solver.h ../flat/ml_struct.h ../flat/ml_vec.h ../flat/mli_solver.h
	$(CC) -c $(CFLAGS) ../flat/ml_mat_formats.c -o $@

ml_matmat_mult.o : ../flat/ml_matmat_mult.c ../flat/ml_1level.h \
        ../flat/ml_bdrypts.h ../flat/ml_comm.h ../flat/ml_comminfoagx.h ../flat/ml_comminfoop.h ../flat/ml_csolve.h ../flat/ml_defs.h ../flat/ml_grid.h ../flat/ml_gridfunc.h ../flat/ml_krylov.h ../flat/ml_mapper.h ../flat/ml_mat_formats.h ../flat/ml_operator.h ../flat/ml_operatoragx.h ../flat/ml_rap.h ../flat/ml_smoother.h ../flat/ml_solver.h ../flat/ml_struct.h ../flat/ml_utils.h ../flat/ml_vec.h
	$(CC) -c $(CFLAGS) ../flat/ml_matmat_mult.c -o $@

ml_op_utils.o : ../flat/ml_op_utils.c ../flat/ml_agg_genP.h \
        ../flat/ml_aggregate.h ../flat/ml_aztec_utils.h ../flat/ml_bicgstabl.h ../flat/ml_cg.h ../flat/ml_comm.h ../flat/ml_defs.h ../flat/ml_elementagx.h ../flat/ml_gmres.h ../flat/ml_grid.h ../flat/ml_gridagx.h ../flat/ml_gridfunc.h ../flat/ml_include.h ../flat/ml_intlist.h ../flat/ml_krylov.h ../flat/ml_operator.h ../flat/ml_operatoragx.h ../flat/ml_pde.h ../flat/ml_solver.h ../flat/ml_struct.h ../flat/ml_vec.h ../flat/mli_solver.h
	$(CC) -c $(CFLAGS) ../flat/ml_op_utils.c -o $@

ml_operator.o : ../flat/ml_operator.c ../flat/ml_1level.h \
        ../flat/ml_bdrypts.h ../flat/ml_defs.h ../flat/ml_memory.h ../flat/ml_operator.h ../flat/ml_operatoragx.h ../flat/ml_vec.h
	$(CC) -c $(CFLAGS) ../flat/ml_operator.c -o $@

ml_operatoragx.o : ../flat/ml_operatoragx.c ../flat/ml_comm.h \
        ../flat/ml_comminfoagx.h ../flat/ml_memory.h ../flat/ml_operatoragx.h ../flat/ml_struct.h
	$(CC) -c $(CFLAGS) ../flat/ml_operatoragx.c -o $@

ml_rap.o : ../flat/ml_rap.c ../flat/ml_1level.h ../flat/ml_bdrypts.h \
        ../flat/ml_comm.h ../flat/ml_comminfoagx.h ../flat/ml_comminfoop.h ../flat/ml_csolve.h ../flat/ml_defs.h ../flat/ml_grid.h ../flat/ml_gridfunc.h ../flat/ml_krylov.h ../flat/ml_mapper.h ../flat/ml_mat_formats.h ../flat/ml_operator.h ../flat/ml_operatoragx.h ../flat/ml_rap.h ../flat/ml_smoother.h ../flat/ml_solver.h ../flat/ml_struct.h ../flat/ml_utils.h ../flat/ml_vec.h
	$(CC) -c $(CFLAGS) ../flat/ml_rap.c -o $@

ml_rap_utils.o : ../flat/ml_rap_utils.c ../flat/ml_1level.h \
        ../flat/ml_bdrypts.h ../flat/ml_comm.h ../flat/ml_comminfoop.h ../flat/ml_csolve.h ../flat/ml_defs.h ../flat/ml_grid.h ../flat/ml_gridfunc.h ../flat/ml_krylov.h ../flat/ml_mapper.h ../flat/ml_mat_formats.h ../flat/ml_operator.h ../flat/ml_operatoragx.h ../flat/ml_rap.h ../flat/ml_smoother.h ../flat/ml_solver.h ../flat/ml_struct.h ../flat/ml_utils.h ../flat/ml_vec.h
	$(CC) -c $(CFLAGS) ../flat/ml_rap_utils.c -o $@

ml_csolve.o : ../flat/ml_csolve.c ../flat/ml_1level.h ../flat/ml_csolve.h \
        ../flat/ml_defs.h ../flat/ml_memory.h
	$(CC) -c $(CFLAGS) ../flat/ml_csolve.c -o $@

ml_smoother.o : ../flat/ml_smoother.c ../flat/ml_1level.h \
        ../flat/ml_aztec_utils.h ../flat/ml_defs.h ../flat/ml_include.h ../flat/ml_lapack.h ../flat/ml_memory.h ../flat/ml_smoother.h
	$(CC) -c $(CFLAGS) ../flat/ml_smoother.c -o $@

ml_solver.o : ../flat/ml_solver.c ../flat/ml_agg_genP.h \
        ../flat/ml_aggregate.h ../flat/ml_aztec_utils.h ../flat/ml_bicgstabl.h ../flat/ml_cg.h ../flat/ml_comm.h ../flat/ml_defs.h ../flat/ml_elementagx.h ../flat/ml_gmres.h ../flat/ml_grid.h ../flat/ml_gridagx.h ../flat/ml_gridfunc.h ../flat/ml_include.h ../flat/ml_intlist.h ../flat/ml_krylov.h ../flat/ml_operator.h ../flat/ml_operatoragx.h ../flat/ml_pde.h ../flat/ml_solver.h ../flat/ml_struct.h ../flat/ml_vec.h ../flat/mli_solver.h
	$(CC) -c $(CFLAGS) ../flat/ml_solver.c -o $@

ml_superlu.o : ../flat/ml_superlu.c ../flat/ml_1level.h \
        ../flat/ml_bdrypts.h ../flat/ml_comm.h ../flat/ml_comminfoop.h ../flat/ml_csolve.h ../flat/ml_defs.h ../flat/ml_grid.h ../flat/ml_gridfunc.h ../flat/ml_krylov.h ../flat/ml_mapper.h ../flat/ml_mat_formats.h ../flat/ml_memory.h ../flat/ml_operator.h ../flat/ml_operatoragx.h ../flat/ml_rap.h ../flat/ml_smoother.h ../flat/ml_solver.h ../flat/ml_struct.h ../flat/ml_utils.h ../flat/ml_vec.h
	$(CC) -c $(CFLAGS) ../flat/ml_superlu.c -o $@

ml_xxt.o : ../flat/ml_xxt.c ../flat/ml_agg_genP.h ../flat/ml_aggregate.h \
        ../flat/ml_aztec_utils.h ../flat/ml_bicgstabl.h ../flat/ml_cg.h ../flat/ml_comm.h ../flat/ml_defs.h ../flat/ml_elementagx.h ../flat/ml_gmres.h ../flat/ml_grid.h ../flat/ml_gridagx.h ../flat/ml_gridfunc.h ../flat/ml_include.h ../flat/ml_intlist.h ../flat/ml_krylov.h ../flat/ml_operator.h ../flat/ml_operatoragx.h ../flat/ml_pde.h ../flat/ml_solver.h ../flat/ml_struct.h ../flat/ml_vec.h ../flat/mli_solver.h
	$(CC) -c $(CFLAGS) ../flat/ml_xxt.c -o $@

ml_aztec_utils.o : ../flat/ml_aztec_utils.c ../flat/ml_agg_genP.h \
        ../flat/ml_aggregate.h ../flat/ml_aztec_utils.h ../flat/ml_bicgstabl.h ../flat/ml_cg.h ../flat/ml_comm.h ../flat/ml_defs.h ../flat/ml_elementagx.h ../flat/ml_gmres.h ../flat/ml_grid.h ../flat/ml_gridagx.h ../flat/ml_gridfunc.h ../flat/ml_include.h ../flat/ml_intlist.h ../flat/ml_krylov.h ../flat/ml_operator.h ../flat/ml_operatoragx.h ../flat/ml_pde.h ../flat/ml_solver.h ../flat/ml_struct.h ../flat/ml_vec.h ../flat/mli_solver.h
	$(CC) -c $(CFLAGS) ../flat/ml_aztec_utils.c -o $@

ml_intlist.o : ../flat/ml_intlist.c ../flat/ml_defs.h ../flat/ml_intlist.h \
        ../flat/ml_memory.h
	$(CC) -c $(CFLAGS) ../flat/ml_intlist.c -o $@

ml_lapack.o : ../flat/ml_lapack.c ../flat/ml_aztec_lapack.h \
        ../flat/ml_defs.h ../flat/ml_superlu_lapack.h ../flat/ml_vendor_lapack.h
	$(CC) -c $(CFLAGS) ../flat/ml_lapack.c -o $@

ml_memory.o : ../flat/ml_memory.c ../flat/ml_comm.h ../flat/ml_defs.h \
        ../flat/ml_memory.h
	$(CC) -c $(CFLAGS) ../flat/ml_memory.c -o $@

ml_rbm.o : ../flat/ml_rbm.c ../flat/ml_rbm.h
	$(CC) -c $(CFLAGS) ../flat/ml_rbm.c -o $@

ml_utils.o : ../flat/ml_utils.c ../flat/ml_comm.h ../flat/ml_defs.h \
        ../flat/ml_memory.h ../flat/ml_utils.h
	$(CC) -c $(CFLAGS) ../flat/ml_utils.c -o $@

ml_vec.o : ../flat/ml_vec.c ../flat/ml_comm.h ../flat/ml_defs.h \
        ../flat/ml_memory.h ../flat/ml_vec.h
	$(CC) -c $(CFLAGS) ../flat/ml_vec.c -o $@

ml_read_elas.o : ../examples/ml_read_elas.c ../flat/ml_agg_genP.h \
        ../flat/ml_aggregate.h ../flat/ml_aztec_utils.h ../flat/ml_bicgstabl.h ../flat/ml_cg.h ../flat/ml_comm.h ../flat/ml_defs.h ../flat/ml_elementagx.h ../flat/ml_gmres.h ../flat/ml_grid.h ../flat/ml_gridagx.h ../flat/ml_gridfunc.h ../flat/ml_include.h ../flat/ml_intlist.h ../flat/ml_krylov.h ../flat/ml_operator.h ../flat/ml_operatoragx.h ../flat/ml_pde.h ../flat/ml_solver.h ../flat/ml_struct.h ../flat/ml_vec.h ../flat/mli_solver.h
	$(CC) -c $(CFLAGS) ../examples/ml_read_elas.c -o $@

convertSund2AZdatafile.o : ../examples/convertSund2AZdatafile.c
	$(CC) -c $(CFLAGS) ../examples/convertSund2AZdatafile.c -o $@

ml_ex1d.o : ../examples/ml_ex1d.c ../flat/ml_agg_genP.h \
        ../flat/ml_aggregate.h ../flat/ml_aztec_utils.h ../flat/ml_bicgstabl.h ../flat/ml_cg.h ../flat/ml_comm.h ../flat/ml_defs.h ../flat/ml_elementagx.h ../flat/ml_gmres.h ../flat/ml_grid.h ../flat/ml_gridagx.h ../flat/ml_gridfunc.h ../flat/ml_include.h ../flat/ml_intlist.h ../flat/ml_krylov.h ../flat/ml_operator.h ../flat/ml_operatoragx.h ../flat/ml_pde.h ../flat/ml_solver.h ../flat/ml_struct.h ../flat/ml_vec.h ../flat/mli_solver.h
	$(CC) -c $(CFLAGS) ../examples/ml_ex1d.c -o $@

ml_example1d.o : ../examples/ml_example1d.c ../flat/ml_agg_genP.h \
        ../flat/ml_aggregate.h ../flat/ml_aztec_utils.h ../flat/ml_bicgstabl.h ../flat/ml_cg.h ../flat/ml_comm.h ../flat/ml_defs.h ../flat/ml_elementagx.h ../flat/ml_gmres.h ../flat/ml_grid.h ../flat/ml_gridagx.h ../flat/ml_gridfunc.h ../flat/ml_include.h ../flat/ml_intlist.h ../flat/ml_krylov.h ../flat/ml_operator.h ../flat/ml_operatoragx.h ../flat/ml_pde.h ../flat/ml_solver.h ../flat/ml_struct.h ../flat/ml_vec.h ../flat/mli_solver.h
	$(CC) -c $(CFLAGS) ../examples/ml_example1d.c -o $@

ml_example1dGS.o : ../examples/ml_example1dGS.c ../flat/ml_agg_genP.h \
        ../flat/ml_aggregate.h ../flat/ml_aztec_utils.h ../flat/ml_bicgstabl.h ../flat/ml_cg.h ../flat/ml_comm.h ../flat/ml_defs.h ../flat/ml_elementagx.h ../flat/ml_gmres.h ../flat/ml_grid.h ../flat/ml_gridagx.h ../flat/ml_gridfunc.h ../flat/ml_include.h ../flat/ml_intlist.h ../flat/ml_krylov.h ../flat/ml_operator.h ../flat/ml_operatoragx.h ../flat/ml_pde.h ../flat/ml_solver.h ../flat/ml_struct.h ../flat/ml_vec.h ../flat/mli_solver.h
	$(CC) -c $(CFLAGS) ../examples/ml_example1dGS.c -o $@

ml_example2d.o : ../examples/ml_example2d.c ../flat/ml_agg_genP.h \
        ../flat/ml_aggregate.h ../flat/ml_aztec_utils.h ../flat/ml_bicgstabl.h ../flat/ml_cg.h ../flat/ml_comm.h ../flat/ml_defs.h ../flat/ml_elementagx.h ../flat/ml_gmres.h ../flat/ml_grid.h ../flat/ml_gridagx.h ../flat/ml_gridfunc.h ../flat/ml_include.h ../flat/ml_intlist.h ../flat/ml_krylov.h ../flat/ml_operator.h ../flat/ml_operatoragx.h ../flat/ml_pde.h ../flat/ml_solver.h ../flat/ml_struct.h ../flat/ml_vec.h ../flat/mli_solver.h
	$(CC) -c $(CFLAGS) ../examples/ml_example2d.c -o $@

ml_example3d.o : ../examples/ml_example3d.c ../flat/ml_agg_genP.h \
        ../flat/ml_aggregate.h ../flat/ml_aztec_utils.h ../flat/ml_bicgstabl.h ../flat/ml_cg.h ../flat/ml_comm.h ../flat/ml_defs.h ../flat/ml_elementagx.h ../flat/ml_gmres.h ../flat/ml_grid.h ../flat/ml_gridagx.h ../flat/ml_gridfunc.h ../flat/ml_include.h ../flat/ml_intlist.h ../flat/ml_krylov.h ../flat/ml_operator.h ../flat/ml_operatoragx.h ../flat/ml_pde.h ../flat/ml_solver.h ../flat/ml_struct.h ../flat/ml_vec.h ../flat/mli_solver.h
	$(CC) -c $(CFLAGS) ../examples/ml_example3d.c -o $@

ml_readex.o : ../examples/ml_readex.c ../flat/ml_agg_genP.h \
        ../flat/ml_aggregate.h ../flat/ml_aztec_utils.h ../flat/ml_bicgstabl.h ../flat/ml_cg.h ../flat/ml_comm.h ../flat/ml_defs.h ../flat/ml_elementagx.h ../flat/ml_gmres.h ../flat/ml_grid.h ../flat/ml_gridagx.h ../flat/ml_gridfunc.h ../flat/ml_include.h ../flat/ml_intlist.h ../flat/ml_krylov.h ../flat/ml_operator.h ../flat/ml_operatoragx.h ../flat/ml_pde.h ../flat/ml_solver.h ../flat/ml_struct.h ../flat/ml_vec.h ../flat/mli_solver.h
	$(CC) -c $(CFLAGS) ../examples/ml_readex.c -o $@

ml_star2d.o : ../examples/ml_star2d.c ../flat/ml_agg_genP.h \
        ../flat/ml_aggregate.h ../flat/ml_aztec_utils.h ../flat/ml_bicgstabl.h ../flat/ml_cg.h ../flat/ml_comm.h ../flat/ml_defs.h ../flat/ml_elementagx.h ../flat/ml_gmres.h ../flat/ml_grid.h ../flat/ml_gridagx.h ../flat/ml_gridfunc.h ../flat/ml_include.h ../flat/ml_intlist.h ../flat/ml_krylov.h ../flat/ml_operator.h ../flat/ml_operatoragx.h ../flat/ml_pde.h ../flat/ml_solver.h ../flat/ml_struct.h ../flat/ml_vec.h ../flat/mli_solver.h
	$(CC) -c $(CFLAGS) ../examples/ml_star2d.c -o $@

mlguide.o : ../examples/mlguide.c ../flat/ml_agg_genP.h \
        ../flat/ml_aggregate.h ../flat/ml_aztec_utils.h ../flat/ml_bicgstabl.h ../flat/ml_cg.h ../flat/ml_comm.h ../flat/ml_defs.h ../flat/ml_elementagx.h ../flat/ml_gmres.h ../flat/ml_grid.h ../flat/ml_gridagx.h ../flat/ml_gridfunc.h ../flat/ml_include.h ../flat/ml_intlist.h ../flat/ml_krylov.h ../flat/ml_operator.h ../flat/ml_operatoragx.h ../flat/ml_pde.h ../flat/ml_solver.h ../flat/ml_struct.h ../flat/ml_vec.h ../flat/mli_solver.h
	$(CC) -c $(CFLAGS) ../examples/mlguide.c -o $@

mlguide_par.o : ../examples/mlguide_par.c ../flat/ml_agg_genP.h \
        ../flat/ml_aggregate.h ../flat/ml_aztec_utils.h ../flat/ml_bicgstabl.h ../flat/ml_cg.h ../flat/ml_comm.h ../flat/ml_defs.h ../flat/ml_elementagx.h ../flat/ml_gmres.h ../flat/ml_grid.h ../flat/ml_gridagx.h ../flat/ml_gridfunc.h ../flat/ml_include.h ../flat/ml_intlist.h ../flat/ml_krylov.h ../flat/ml_operator.h ../flat/ml_operatoragx.h ../flat/ml_pde.h ../flat/ml_solver.h ../flat/ml_struct.h ../flat/ml_vec.h ../flat/mli_solver.h
	$(CC) -c $(CFLAGS) ../examples/mlguide_par.c -o $@

new_readex.o : ../examples/new_readex.c ../flat/ml_agg_genP.h \
        ../flat/ml_aggregate.h ../flat/ml_aztec_utils.h ../flat/ml_bicgstabl.h ../flat/ml_cg.h ../flat/ml_comm.h ../flat/ml_defs.h ../flat/ml_elementagx.h ../flat/ml_gmres.h ../flat/ml_grid.h ../flat/ml_gridagx.h ../flat/ml_gridfunc.h ../flat/ml_include.h ../flat/ml_intlist.h ../flat/ml_krylov.h ../flat/ml_operator.h ../flat/ml_operatoragx.h ../flat/ml_pde.h ../flat/ml_solver.h ../flat/ml_struct.h ../flat/ml_vec.h ../flat/mli_solver.h
	$(CC) -c $(CFLAGS) ../examples/new_readex.c -o $@

oldml_readex.o : ../examples/oldml_readex.c ../flat/ml_agg_genP.h \
        ../flat/ml_aggregate.h ../flat/ml_aztec_utils.h ../flat/ml_bicgstabl.h ../flat/ml_cg.h ../flat/ml_comm.h ../flat/ml_defs.h ../flat/ml_elementagx.h ../flat/ml_gmres.h ../flat/ml_grid.h ../flat/ml_gridagx.h ../flat/ml_gridfunc.h ../flat/ml_include.h ../flat/ml_intlist.h ../flat/ml_krylov.h ../flat/ml_operator.h ../flat/ml_operatoragx.h ../flat/ml_pde.h ../flat/ml_solver.h ../flat/ml_struct.h ../flat/ml_vec.h ../flat/mli_solver.h
	$(CC) -c $(CFLAGS) ../examples/oldml_readex.c -o $@

seg_readex.o : ../examples/seg_readex.c ../flat/ml_agg_genP.h \
        ../flat/ml_aggregate.h ../flat/ml_aztec_utils.h ../flat/ml_bicgstabl.h ../flat/ml_cg.h ../flat/ml_comm.h ../flat/ml_defs.h ../flat/ml_elementagx.h ../flat/ml_gmres.h ../flat/ml_grid.h ../flat/ml_gridagx.h ../flat/ml_gridfunc.h ../flat/ml_include.h ../flat/ml_intlist.h ../flat/ml_krylov.h ../flat/ml_operator.h ../flat/ml_operatoragx.h ../flat/ml_pde.h ../flat/ml_solver.h ../flat/ml_struct.h ../flat/ml_vec.h ../flat/mli_solver.h
	$(CC) -c $(CFLAGS) ../examples/seg_readex.c -o $@

ml_recirc.o : ../examples/ml_recirc.c ../flat/ml_agg_genP.h \
        ../flat/ml_aggregate.h ../flat/ml_aztec_utils.h ../flat/ml_bicgstabl.h ../flat/ml_cg.h ../flat/ml_comm.h ../flat/ml_defs.h ../flat/ml_elementagx.h ../flat/ml_gmres.h ../flat/ml_grid.h ../flat/ml_gridagx.h ../flat/ml_gridfunc.h ../flat/ml_include.h ../flat/ml_intlist.h ../flat/ml_krylov.h ../flat/ml_operator.h ../flat/ml_operatoragx.h ../flat/ml_pde.h ../flat/ml_solver.h ../flat/ml_struct.h ../flat/ml_vec.h ../flat/mli_solver.h
	$(CC) -c $(CFLAGS) ../examples/ml_recirc.c -o $@

ml_read_salsa.o : ../examples/ml_read_salsa.c ../flat/ml_agg_genP.h \
        ../flat/ml_aggregate.h ../flat/ml_aztec_utils.h ../flat/ml_bicgstabl.h ../flat/ml_cg.h ../flat/ml_comm.h ../flat/ml_defs.h ../flat/ml_elementagx.h ../flat/ml_gmres.h ../flat/ml_grid.h ../flat/ml_gridagx.h ../flat/ml_gridfunc.h ../flat/ml_include.h ../flat/ml_intlist.h ../flat/ml_krylov.h ../flat/ml_operator.h ../flat/ml_operatoragx.h ../flat/ml_pde.h ../flat/ml_solver.h ../flat/ml_struct.h ../flat/ml_vec.h ../flat/mli_solver.h
	$(CC) -c $(CFLAGS) ../examples/ml_read_salsa.c -o $@

