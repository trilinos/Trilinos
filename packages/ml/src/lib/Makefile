include ../../etc/makefile.$(ARCH)

#=======================================================================
# ML source files
#=======================================================================

ML_OBJ = ml_aztec_utils.o  ml_ggraph.o        ml_memory.o \
         ml_struct.o       ml_comm.o          ml_gridagx.o \
         ml_operatoragx.o  ml_superlu.o       ml_comminfoagx.o \
         ml_grid.o         ml_rap.o           ml_aggregate.o \
         ml_elementagx.o   ml_intlist.o       ml_rap_utils.o \
         ml_utils.o        ml_exch_row.o      ml_xxt.o \
         ml_setup.o        ml_vec.o           ml_get_basis.o \
         ml_matmat_mult.o  ml_solver.o        ml_check.o \
         ml_bdrypts.o      ml_mapper.o        ml_gridfunc.o \
         ml_smoother.o     ml_csolve.o        ml_operator.o \
         ml_op_utils.o     ml_comminfoop.o    ml_mat_formats.o \
         ml_agg_genP.o     ml_lapack.o        mli_solver.o \
         ml_krylov.o       ml_cg.o            ml_gmres.o \
         ml_bicgstabl.o    ml_pde.o           ml_agg_dd.o \
         ml_agg_coupled.o  ml_agg_uncoupled.o ml_agg_MIS.o \
         ml_amg.o          ml_amg_genP.o      ml_amg_MIS.o


#=======================================================================
# ML library file
#=======================================================================

libml.a: $(ML_OBJ)
	$(AR) rcv libml.a $(ML_OBJ); \
	$(RANLIB) libml.a

clean:
	@echo "junk" > dummy.o
	@rm -f *.o *.a core 

