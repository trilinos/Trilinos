
# The following line helps the test harness recover from build errors.

all-local:
	@echo ""
	@echo "Trilinos package Zoltan subdirectory src built successfully."
	@echo ""

ALL_H = \
	$(srcdir)/all/all_allo_const.h

ALL = \
	$(srcdir)/all/all_allo.c

COLORING_H = \
	$(srcdir)/coloring/coloring_const.h \
	$(srcdir)/coloring/coloring.h \
	$(srcdir)/coloring/bucket.h \
	$(srcdir)/coloring/g2l_hash.h

COLORING = \
	$(srcdir)/coloring/coloring.c \
	$(srcdir)/coloring/color_test.c \
	$(srcdir)/coloring/bucket.c \
	$(srcdir)/coloring/g2l_hash.c

GRAPH_H = \
	$(srcdir)/graph/graph_const.h \
	$(srcdir)/graph/graph.h

GRAPH = \
	$(srcdir)/graph/graph.c

HA_H = \
	$(srcdir)/ha/ha_const.h \
	$(srcdir)/ha/ha_ovis.h 

HA = \
	$(srcdir)/ha/divide_machine.c \
	$(srcdir)/ha/get_processor_name.c \
	$(srcdir)/ha/ha_ovis.c 

HIER_H = \
	$(srcdir)/hier/hier.h

HIER = \
	$(srcdir)/hier/hier.c \
	$(srcdir)/hier/hier_free_struct.c

HSFC_H = \
	$(srcdir)/hsfc/hsfc.h \
	$(srcdir)/hsfc/hsfc_const.h \
	$(srcdir)/hsfc/hsfc_hilbert_const.h

HSFC = \
	$(srcdir)/hsfc/hsfc_box_assign.c \
	$(srcdir)/hsfc/hsfc.c \
	$(srcdir)/hsfc/hsfc_hilbert.c \
	$(srcdir)/hsfc/hsfc_point_assign.c

INCLUDE_H = \
	$(srcdir)/include/lbi_const.h \
	$(srcdir)/include/zoltan_align.h \
	$(srcdir)/include/zoltan_cpp.h \
	$(srcdir)/include/zoltan.h \
	$(srcdir)/include/zoltan_eval.h \
	$(srcdir)/include/zoltan_types.h \
	$(srcdir)/include/zoltan_comm.h \
	$(srcdir)/include/zoltan_comm_cpp.h \
	$(srcdir)/include/zoltan_timer_cpp.h \
	$(srcdir)/include/zoltan_timer.h \
	$(srcdir)/include/zoltan_dd.h \
	$(srcdir)/include/zoltan_dd_cpp.h \
	$(srcdir)/include/zoltan_partition_tree.h \
	$(srcdir)/include/zoltan_mem.h

#       $(srcdir)/Zoltan_ConfigDefs.h

LB_H = \
	$(srcdir)/lb/lb_const.h \
	$(srcdir)/lb/lb_init_const.h

LB = \
	$(srcdir)/lb/lb_balance.c \
	$(srcdir)/lb/lb_box_assign.c \
	$(srcdir)/lb/lb_copy.c \
	$(srcdir)/lb/lb_eval.c \
	$(srcdir)/lb/lb_free.c \
	$(srcdir)/lb/lb_init.c \
	$(srcdir)/lb/lb_invert.c \
	$(srcdir)/lb/lb_migrate.c \
	$(srcdir)/lb/lb_part2proc.c \
	$(srcdir)/lb/lb_point_assign.c \
	$(srcdir)/lb/lb_remap.c \
	$(srcdir)/lb/lb_set_fn.c \
	$(srcdir)/lb/lb_set_method.c \
	$(srcdir)/lb/lb_set_part_sizes.c

MATRIX_H = \
	$(srcdir)/matrix/matrix.h

MATRIX = \
	$(srcdir)/matrix/matrix_build.c \
	$(srcdir)/matrix/matrix_distribute.c \
	$(srcdir)/matrix/matrix_operations.c \
	$(srcdir)/matrix/matrix_sym.c \
	$(srcdir)/matrix/matrix_utils.c 

ORDER_H = \
	$(srcdir)/order/order_const.h \
	$(srcdir)/order/hsfcOrder.h 

ORDER = \
	$(srcdir)/order/order.c \
	$(srcdir)/order/order_struct.c \
	$(srcdir)/order/order_tools.c \
	$(srcdir)/order/hsfcOrder.c \
	$(srcdir)/order/perm.c

PAR_H = \
	$(srcdir)/par/par_average_const.h \
	$(srcdir)/par/par_bisect_const.h \
	$(srcdir)/par/par_const.h \
	$(srcdir)/par/par_median_const.h \
	$(srcdir)/par/par_tflops_special_const.h

PAR = \
	$(srcdir)/par/par_average.c \
	$(srcdir)/par/par_bisect.c \
	$(srcdir)/par/par_median.c \
	$(srcdir)/par/par_median_randomized.c \
	$(srcdir)/par/par_stats.c \
	$(srcdir)/par/par_sync.c \
	$(srcdir)/par/par_tflops_special.c

PARAMS_H = \
	$(srcdir)/params/key_params.h \
	$(srcdir)/params/params_const.h

PARAMS = \
	$(srcdir)/params/assign_param_vals.c \
	$(srcdir)/params/bind_param.c \
	$(srcdir)/params/check_param.c \
	$(srcdir)/params/free_params.c \
	$(srcdir)/params/key_params.c \
	$(srcdir)/params/print_params.c \
	$(srcdir)/params/set_param.c

TPLS_H = \
	$(srcdir)/tpls/third_library_const.h \
	$(srcdir)/tpls/third_library.h \
	$(srcdir)/tpls/third_library_tools.h \
	$(srcdir)/tpls/graph_util.h 

TPLS = \
	$(srcdir)/tpls/build_graph.c \
	$(srcdir)/tpls/postprocessing.c \
	$(srcdir)/tpls/preprocessing.c \
	$(srcdir)/tpls/scatter_graph.c \
	$(srcdir)/tpls/third_library.c \
	$(srcdir)/tpls/verify_graph.c

EXTRA_PARMETIS_INTERFACE_H = \
	$(srcdir)/tpls/parmetis_interface.h

EXTRA_PARMETIS_INTERFACE = \
	$(srcdir)/tpls/parmetis_interface.c

if BUILD_PARMETIS

PARMETIS_INTERFACE_H = $(EXTRA_PARMETIS_INTERFACE_H)

PARMETIS_INTERFACE = $(EXTRA_PARMETIS_INTERFACE)

else

PARMETIS_INTERFACE_H =

PARMETIS_INTERFACE =

endif

EXTRA_SCOTCH_INTERFACE_H = \
	$(srcdir)/tpls/scotch_interface.h

EXTRA_SCOTCH_INTERFACE = \
	$(srcdir)/tpls/scotch_interface.c

if BUILD_SCOTCH

SCOTCH_INTERFACE_H = $(EXTRA_SCOTCH_INTERFACE_H)

SCOTCH_INTERFACE = $(EXTRA_SCOTCH_INTERFACE)

else

SCOTCH_INTERFACE_H =

SCOTCH_INTERFACE =

endif

PHG_H = \
	$(srcdir)/phg/phg_comm.h \
	$(srcdir)/phg/phg_const.h \
	$(srcdir)/phg/phg_distrib.h \
	$(srcdir)/phg/phg.h \
	$(srcdir)/phg/phg_lookup.h \
	$(srcdir)/phg/phg_verbose.h \
	$(srcdir)/phg/phg_hypergraph.h \
	$(srcdir)/phg/phg_tree.h \
	$(srcdir)/phg/phg_util.h

PHG = \
	$(srcdir)/phg/phg_build.c \
	$(srcdir)/phg/phg_build_calls.c \
	$(srcdir)/phg/phg.c \
	$(srcdir)/phg/phg_lookup.c \
	$(srcdir)/phg/phg_verbose.c \
	$(srcdir)/phg/phg_coarse.c \
	$(srcdir)/phg/phg_comm.c \
	$(srcdir)/phg/phg_distrib.c \
	$(srcdir)/phg/phg_gather.c \
	$(srcdir)/phg/phg_hypergraph.c \
	$(srcdir)/phg/phg_match.c \
	$(srcdir)/phg/phg_order.c \
	$(srcdir)/phg/phg_parkway.c \
	$(srcdir)/phg/phg_patoh.c \
	$(srcdir)/phg/phg_plot.c \
	$(srcdir)/phg/phg_rdivide.c \
	$(srcdir)/phg/phg_refinement.c \
	$(srcdir)/phg/phg_scale.c \
	$(srcdir)/phg/phg_serialpartition.c \
	$(srcdir)/phg/phg_util.c \
	$(srcdir)/phg/phg_tree.c \
	$(srcdir)/phg/phg_partition_tree.c \
	$(srcdir)/phg/phg_Vcycle.c

RCB_H = \
	$(srcdir)/rcb/create_proc_list_const.h \
	$(srcdir)/rcb/rcb_const.h \
	$(srcdir)/rcb/rcb.h \
	$(srcdir)/rcb/rib_const.h \
	$(srcdir)/rcb/rib.h \
	$(srcdir)/rcb/shared.h

RCB = \
	$(srcdir)/rcb/box_assign.c \
	$(srcdir)/rcb/create_proc_list.c \
	$(srcdir)/rcb/inertial1d.c \
	$(srcdir)/rcb/inertial2d.c \
	$(srcdir)/rcb/inertial3d.c \
	$(srcdir)/rcb/point_assign.c \
	$(srcdir)/rcb/rcb_box.c \
	$(srcdir)/rcb/rcb_partition_tree.c \
	$(srcdir)/rcb/rcb.c \
	$(srcdir)/rcb/rcb_util.c \
	$(srcdir)/rcb/rib.c \
	$(srcdir)/rcb/rib_util.c \
	$(srcdir)/rcb/shared.c

REFTREE_H = \
	$(srcdir)/reftree/reftree_const.h \
	$(srcdir)/reftree/reftree.h

REFTREE = \
	$(srcdir)/reftree/reftree_build.c \
	$(srcdir)/reftree/reftree_coarse_path.c \
	$(srcdir)/reftree/reftree_hash.c \
	$(srcdir)/reftree/reftree_part.c

SIMPLE = \
	$(srcdir)/simple/block.c \
	$(srcdir)/simple/cyclic.c \
	$(srcdir)/simple/random.c

SIMPLE_H = \
	$(srcdir)/simple/simple_const.h

TIMER_H = \
	$(srcdir)/timer/timer_const.h

TIMER = \
	$(srcdir)/timer/timer_params.c

UTILITIES_H = \
	$(srcdir)/Utilities/Communication/comm.h \
	$(srcdir)/Utilities/DDirectory/DD.h \
	$(srcdir)/Utilities/Timer/timer.h \
	$(srcdir)/Utilities/shared/zoltan_util.h \
	$(srcdir)/Utilities/shared/zoltan_id.h

UTILITIES = \
	$(srcdir)/Utilities/Communication/comm_exchange_sizes.c \
	$(srcdir)/Utilities/Communication/comm_invert_map.c \
	$(srcdir)/Utilities/Communication/comm_do.c \
	$(srcdir)/Utilities/Communication/comm_do_reverse.c \
	$(srcdir)/Utilities/Communication/comm_info.c \
	$(srcdir)/Utilities/Communication/comm_create.c \
	$(srcdir)/Utilities/Communication/comm_resize.c \
	$(srcdir)/Utilities/Communication/comm_sort_ints.c \
	$(srcdir)/Utilities/Communication/comm_destroy.c \
	$(srcdir)/Utilities/Communication/comm_invert_plan.c \
	$(srcdir)/Utilities/Timer/zoltan_timer.c \
	$(srcdir)/Utilities/Timer/timer.c \
	$(srcdir)/Utilities/DDirectory/DD_Memory.c \
	$(srcdir)/Utilities/DDirectory/DD_Find.c \
	$(srcdir)/Utilities/DDirectory/DD_Destroy.c \
	$(srcdir)/Utilities/DDirectory/DD_Set_Neighbor_Hash_Fn3.c \
	$(srcdir)/Utilities/DDirectory/DD_Remove.c \
	$(srcdir)/Utilities/DDirectory/DD_Create.c \
	$(srcdir)/Utilities/DDirectory/DD_Update.c \
	$(srcdir)/Utilities/DDirectory/DD_Stats.c \
	$(srcdir)/Utilities/DDirectory/DD_Hash2.c \
	$(srcdir)/Utilities/DDirectory/DD_Print.c \
	$(srcdir)/Utilities/DDirectory/DD_Set_Neighbor_Hash_Fn2.c \
	$(srcdir)/Utilities/DDirectory/DD_Set_Hash_Fn.c \
	$(srcdir)/Utilities/DDirectory/DD_Set_Neighbor_Hash_Fn1.c \
	$(srcdir)/Utilities/Memory/mem.c \
	$(srcdir)/Utilities/shared/zoltan_align.c \
	$(srcdir)/Utilities/shared/zoltan_id.c

# ?? 
#	$(srcdir)/Utilities/Communication/Epetra_ZoltanMpiDistributor.cpp \
#	$(srcdir)/Utilities/DDirectory/DD_Main_2.cpp \
#	$(srcdir)/Utilities/Communication/comm_main_2.cpp

ZZ_H = \
	$(srcdir)/zz/zz_const.h \
	$(srcdir)/zz/zz_heap.h \
	$(srcdir)/zz/zz_id_const.h \
	$(srcdir)/zz/zz_rand.h \
	$(srcdir)/zz/zz_sort.h \
	$(srcdir)/zz/murmur3.h \
	$(srcdir)/zz/zz_util_const.h

ZZ = \
	$(srcdir)/zz/zz_coord.c \
	$(srcdir)/zz/zz_gen_files.c \
	$(srcdir)/zz/zz_hash.c \
	$(srcdir)/zz/murmur3.c \
	$(srcdir)/zz/zz_map.c \
	$(srcdir)/zz/zz_heap.c \
	$(srcdir)/zz/zz_init.c \
	$(srcdir)/zz/zz_obj_list.c \
	$(srcdir)/zz/zz_rand.c \
	$(srcdir)/zz/zz_set_fn.c \
	$(srcdir)/zz/zz_sort.c \
	$(srcdir)/zz/zz_struct.c \
	$(srcdir)/zz/zz_back_trace.c \
	$(srcdir)/zz/zz_util.c

EXTRA_FORTRAN_H = \
    $(srcdir)/fort/cwrap_fmangle.h \
    $(srcdir)/fort/set_borderobjlist.fn \
    $(srcdir)/fort/set_borderobjlist.if \
    $(srcdir)/fort/set_childlist.fn \
    $(srcdir)/fort/set_childlist.if \
    $(srcdir)/fort/set_childweight.fn \
    $(srcdir)/fort/set_childweight.if \
    $(srcdir)/fort/set_coarseobjlist.fn \
    $(srcdir)/fort/set_coarseobjlist.if \
    $(srcdir)/fort/set_edgelist.fn \
    $(srcdir)/fort/set_edgelist.if \
    $(srcdir)/fort/set_edgelistmulti.fn \
    $(srcdir)/fort/set_edgelistmulti.if \
    $(srcdir)/fort/set_firstborderobj.fn \
    $(srcdir)/fort/set_firstborderobj.if \
    $(srcdir)/fort/set_firstcoarseobj.fn \
    $(srcdir)/fort/set_firstcoarseobj.if \
    $(srcdir)/fort/set_firstobj.fn \
    $(srcdir)/fort/set_firstobj.if \
    $(srcdir)/fort/set_fixedobjlist.fn \
    $(srcdir)/fort/set_fixedobjlist.if \
    $(srcdir)/fort/set_geom.fn \
    $(srcdir)/fort/set_geom.if \
    $(srcdir)/fort/set_geommulti.fn \
    $(srcdir)/fort/set_geommulti.if \
    $(srcdir)/fort/set_hgcs.fn \
    $(srcdir)/fort/set_hgcs.if \
    $(srcdir)/fort/set_hgedgeweights.fn \
    $(srcdir)/fort/set_hgedgeweights.if \
    $(srcdir)/fort/set_hgsizecs.fn \
    $(srcdir)/fort/set_hgsizecs.if \
    $(srcdir)/fort/set_hgsizeedgeweights.fn \
    $(srcdir)/fort/set_hgsizeedgeweights.if \
    $(srcdir)/fort/set_hiermethod.fn \
    $(srcdir)/fort/set_hiermethod.if \
    $(srcdir)/fort/set_hiernumlevels.fn \
    $(srcdir)/fort/set_hiernumlevels.if \
    $(srcdir)/fort/set_hierpartition.fn \
    $(srcdir)/fort/set_hierpartition.if \
    $(srcdir)/fort/set_midmigrate.fn \
    $(srcdir)/fort/set_midmigrate.if \
    $(srcdir)/fort/set_midmigratepp.fn \
    $(srcdir)/fort/set_midmigratepp.if \
    $(srcdir)/fort/set_nextborderobj.fn \
    $(srcdir)/fort/set_nextborderobj.if \
    $(srcdir)/fort/set_nextcoarseobj.fn \
    $(srcdir)/fort/set_nextcoarseobj.if \
    $(srcdir)/fort/set_nextobj.fn \
    $(srcdir)/fort/set_nextobj.if \
    $(srcdir)/fort/set_numborderobj.fn \
    $(srcdir)/fort/set_numborderobj.if \
    $(srcdir)/fort/set_numchild.fn \
    $(srcdir)/fort/set_numchild.if \
    $(srcdir)/fort/set_numcoarseobj.fn \
    $(srcdir)/fort/set_numcoarseobj.if \
    $(srcdir)/fort/set_numedges.fn \
    $(srcdir)/fort/set_numedges.if \
    $(srcdir)/fort/set_numedgesmulti.fn \
    $(srcdir)/fort/set_numedgesmulti.if \
    $(srcdir)/fort/set_numfixedobj.fn \
    $(srcdir)/fort/set_numfixedobj.if \
    $(srcdir)/fort/set_numgeom.fn \
    $(srcdir)/fort/set_numgeom.if \
    $(srcdir)/fort/set_numobj.fn \
    $(srcdir)/fort/set_numobj.if \
    $(srcdir)/fort/set_objlist.fn \
    $(srcdir)/fort/set_objlist.if \
    $(srcdir)/fort/set_objsize.fn \
    $(srcdir)/fort/set_objsize.if \
    $(srcdir)/fort/set_objsizemulti.fn \
    $(srcdir)/fort/set_objsizemulti.if \
    $(srcdir)/fort/set_packobj.fn \
    $(srcdir)/fort/set_packobj.if \
    $(srcdir)/fort/set_packobjmulti.fn \
    $(srcdir)/fort/set_packobjmulti.if \
    $(srcdir)/fort/set_partition.fn \
    $(srcdir)/fort/set_partition.if \
    $(srcdir)/fort/set_partitionmulti.fn \
    $(srcdir)/fort/set_partitionmulti.if \
    $(srcdir)/fort/set_postmigrate.fn \
    $(srcdir)/fort/set_postmigrate.if \
    $(srcdir)/fort/set_postmigratepp.fn \
    $(srcdir)/fort/set_postmigratepp.if \
    $(srcdir)/fort/set_premigrate.fn \
    $(srcdir)/fort/set_premigrate.if \
    $(srcdir)/fort/set_premigratepp.fn \
    $(srcdir)/fort/set_premigratepp.if \
    $(srcdir)/fort/set_unpackobj.fn \
    $(srcdir)/fort/set_unpackobj.if \
    $(srcdir)/fort/set_unpackobjmulti.fn \
    $(srcdir)/fort/set_unpackobjmulti.if 

EXTRA_FORTRAN = \
    $(srcdir)/fort/ztypes.f90 \
    $(srcdir)/fort/cwrap.c \
    $(srcdir)/fort/zoltan_user_data.f90 \
    $(srcdir)/fort/fwrap.f90 

if BUILD_ZOLTAN_F90_INTERFACE
FORTRAN_H =$(EXTRA_FORTRAN_H)
FORTRAN = $(EXTRA_FORTRAN)
else
FORTRAN_H =
FORTRAN =
endif

# ------------------------------------------------------------------------
# Zoltan library specifications
# ------------------------------------------------------------------------
ZOLTAN_LIB = libzoltan.a

ZOLTAN_H = \
	$(ALL_H) \
	$(COLORING_H) \
	$(GRAPH_H) \
	$(HA_H) \
	$(HIER_H) \
	$(HSFC_H) \
	$(INCLUDE_H) \
	$(LB_H) \
	$(MATRIX_H) \
	$(ORDER_H) \
	$(PAR_H) \
	$(PARAMS_H) \
	$(TPLS_H) \
	$(PARMETIS_INTERFACE_H) \
	$(PHG_H) \
	$(RCB_H) \
	$(REFTREE_H) \
	$(SCOTCH_INTERFACE_H) \
	$(TIMER_H) \
	$(SIMPLE_H) \
	$(UTILITIES_H) \
	$(ZZ_H) \
	$(FORTRAN_H)

libzoltan_a_SOURCES = \
	$(ALL) \
	$(COLORING) \
	$(GRAPH) \
	$(HA) \
	$(HIER) \
	$(HSFC) \
	$(LB) \
	$(MATRIX) \
	$(ORDER) \
	$(PAR) \
	$(PARAMS) \
	$(TPLS) \
	$(PARMETIS_INTERFACE) \
	$(PHG) \
	$(RCB) \
	$(REFTREE) \
	$(SCOTCH_INTERFACE) \
	$(SIMPLE) \
	$(TIMER) \
	$(UTILITIES) \
	$(ZZ) \
	$(FORTRAN)

EXTRA_libzoltan_a_SOURCES = \
	$(EXTRA_PARMETIS_INTERFACE_H) \
	$(EXTRA_PARMETIS_INTERFACE) \
	$(EXTRA_SCOTCH_INTERFACE_H) \
	$(EXTRA_SCOTCH_INTERFACE) \
	$(EXTRA_FORTRAN_H) \
	$(EXTRA_FORTRAN)

if F90_MODULE_LOWERCASE
FORTRAN_MODULES = \
	zoltan_types.mod \
	zoltan.mod \
	zoltan_user_data.mod
else
FORTRAN_MODULES = \
	ZOLTAN_TYPES.mod \
	ZOLTAN.mod \
	ZOLTAN_USER_DATA.mod
endif

include $(top_builddir)/Makefile.export.zoltan

EXPORT_INCLUDES = $(ZOLTAN_LIBRARY_INCLUDES)

AM_CPPFLAGS = $(EXPORT_INCLUDES) @MPI_RECV_LIMIT_FLAG@

# ------------------------------------------------------------------------
# For using a special archiver
# ------------------------------------------------------------------------

# SPECIAL NOTE: Why is the "else" there?  For some reason, automake
# fails to define <lib>_AR of USE_ALTERNATE_AR is false, so we have to
# define it manually.

if USE_ALTERNATE_AR
libzoltan_a_AR = $(ALTERNATE_AR)
else
libzoltan_a_AR = $(AR) cru
endif

# ------------------------------------------------------------------------
# Some C++ compilers create extra .o-files for templates. We need to
# be sure to include these, and this is the hack to do it.
# ------------------------------------------------------------------------

libzoltan_a_LIBADD = $(XTRALDADD)

# ------------------------------------------------------------------------
# List of all libraries to install in $(libexecdir)
# ------------------------------------------------------------------------

#np# replace new_package with the name of the package being autotool'ed here
lib_LIBRARIES = $(ZOLTAN_LIB)

# ------------------------------------------------------------------------
# List of all headers to install in $(includedir)
# ------------------------------------------------------------------------

# Eventually we will try to install fewer headers
include_HEADERS = $(INCLUDE_H) 

# ------------------------------------------------------------------------
# List of all Fortran Module files to install in $(includedir)
# ------------------------------------------------------------------------
if BUILD_ZOLTAN_F90_INTERFACE
include_HEADERS += $(FORTRAN_MODULES)
endif

# ------------------------------------------------------------------------
# Other files we want in the tarball
# ------------------------------------------------------------------------

EXTRA_FORTRAN_INTERFACE = \
    $(srcdir)/fort/README \
    $(srcdir)/fort/makefile 
    
ZOLTAN_MATLAB = \
    $(srcdir)/matlab/README \
    $(srcdir)/matlab/mmread.m \
    $(srcdir)/matlab/mmwrite.m \
    $(srcdir)/matlab/plotcolors.m \
    $(srcdir)/matlab/zdrive.matlab \
    $(srcdir)/matlab/zoltPartSpy.m \
    $(srcdir)/matlab/zoltan.m

ZOLTAN_POSTPROCESS = \
    $(srcdir)/postprocessing/bin/toscotchperm.py \
    $(srcdir)/postprocessing/ordering/README \
    $(srcdir)/postprocessing/ordering/mm_input.c \
    $(srcdir)/postprocessing/ordering/mm_input.h \
    $(srcdir)/postprocessing/ordering/order_eval.c \
    $(srcdir)/postprocessing/ordering/order_eval.h \
    $(srcdir)/postprocessing/ordering/ordering.py \
    $(srcdir)/postprocessing/ordering/toscotchperm.py \
    $(srcdir)/postprocessing/ordering/zdrive.inp \
    $(srcdir)/postprocessing/ordering/zdrive.inp.parmetis \
    $(srcdir)/postprocessing/ordering/zdrive.inp.scotch

NATIVE_ZOLTAN_BUILD =  \
    $(srcdir)/makefile \
    $(srcdir)/makefile_sub \
    $(srcdir)/makefile.zscript \
    $(srcdir)/Utilities/makefile \
    $(srcdir)/Utilities/makefile_sub \
    $(srcdir)/Utilities/Config/Config.RedStorm \
    $(srcdir)/Utilities/Config/Config.aix \
    $(srcdir)/Utilities/Config/Config.octopiLamPureTest \
    $(srcdir)/Utilities/Config/Config.octopiLamTest \
    $(srcdir)/Utilities/Config/Config.octopiPureTest \
    $(srcdir)/Utilities/Config/Config.bgl \
    $(srcdir)/Utilities/Config/Config.octopiTest \
    $(srcdir)/Utilities/Config/Config.odin \
    $(srcdir)/Utilities/Config/Config.osf1 \
    $(srcdir)/Utilities/Config/Config.pure \
    $(srcdir)/Utilities/Config/Config.pure_siMPI \
    $(srcdir)/Utilities/Config/Config.g++ \
    $(srcdir)/Utilities/Config/Config.purecov \
    $(srcdir)/Utilities/Config/Config.gcc  \
    $(srcdir)/Utilities/Config/Config.purecov2 \
    $(srcdir)/Utilities/Config/Config.gcc64 \
    $(srcdir)/Utilities/Config/Config.qed \
    $(srcdir)/Utilities/Config/Config.gcc_siMPI \
    $(srcdir)/Utilities/Config/Config.quantify \
    $(srcdir)/Utilities/Config/Config.sadl14088 \
    $(srcdir)/Utilities/Config/Config.generic \
    $(srcdir)/Utilities/Config/Config.sgi32 \
    $(srcdir)/Utilities/Config/Config.ibm \
    $(srcdir)/Utilities/Config/Config.sgi64 \
    $(srcdir)/Utilities/Config/Config.solaris \
    $(srcdir)/Utilities/Config/Config.sun \
    $(srcdir)/Utilities/Config/Config.linux \
    $(srcdir)/Utilities/Config/Config.linux32 \
    $(srcdir)/Utilities/Config/Config.linux64 \
    $(srcdir)/Utilities/Config/Config.sunos_forte7 \
    $(srcdir)/Utilities/Config/Config.linux64lam \
    $(srcdir)/Utilities/Config/Config.tbird \
    $(srcdir)/Utilities/Config/Config.linux64lamcov \
    $(srcdir)/Utilities/Config/Config.linux64lampure \
    $(srcdir)/Utilities/Config/Config.linux64pure \
    $(srcdir)/Utilities/Config/Config.linux64scotch \
    $(srcdir)/Utilities/Config/Config.linux_cca \
    $(srcdir)/Utilities/Config/Config.ubgl \
    $(srcdir)/Utilities/Config/Config.vampir \
    $(srcdir)/Utilities/Config/Config.mac32 \
    $(srcdir)/Utilities/Config/Config.x1 \
    $(srcdir)/Utilities/Config/Config.mac64 \
    $(srcdir)/Utilities/Config/Config.mac64_siMPI \
    $(srcdir)/Utilities/Config/Config.xt4scotch \
    $(srcdir)/Utilities/Config/Config.zscript

EXTRA_DIST = $(ZOLTAN_H) \
    $(srcdir)/Utilities/Communication/README \
    $(srcdir)/Utilities/Memory/README \
    $(srcdir)/Utilities/Timer/README \
    $(srcdir)/Utilities/DDirectory/README \
    $(srcdir)/Utilities/README \
    $(NATIVE_ZOLTAN_BUILD) \
    $(ZOLTAN_MATLAB) $(ZOLTAN_POSTPROCESS)

#EXTRA_DIST = $(ZOLTAN_H) \
#    $(srcdir)/Utilities/Communication/README \
#    $(srcdir)/Utilities/Memory/README \
#    $(srcdir)/Utilities/Timer/README \
#    $(srcdir)/Utilities/DDirectory/README \
#    $(srcdir)/Utilities/README \
#    $(NATIVE_ZOLTAN_BUILD) \
#    $(EXTRA_FORTRAN_INTERFACE) \
#    $(srcdir)/fdriver_old/README \
#    $(srcdir)/fdriver_old/makefile \
#    $(ZOLTAN_MATLAB) $(ZOLTAN_POSTPROCESS)

# ------------------------------------------------------------------------
# Special stuff to install in our special $(execincludedir)
# ------------------------------------------------------------------------

#np# replace new_package with the name of the package being autotool'ed here
# SPECIAL NOTE: New_Package_config.h is a machine-dependent file, so we need
# to install it in the machine-dependent directory. However, that is
# not a default installation directory, so we had to create it
# special.

# All Trilinos headers are now installed in the same directory
execincludedir = $(includedir)
#np# replace new_package with the name of the package being autotool'ed here
nodist_execinclude_HEADERS = include/Zoltan_config.h

# Add the absolute path to src/fort to F90 build line, so F90 can find
# the files included in fwrap.f90.  
AM_FCFLAGS = -I@abs_top_srcdir@/src/fort 
AM_FFLAGS = -I@abs_top_srcdir@/src/fort 
# Build zfdrive only if HAVE_MPI; I haven't figured out a proper build of
# siMPI's Fortran interface yet (and serial zfdrive is not an important
# use case).

FDRIVER_SUBDIR =
if BUILD_TESTS
if BUILD_ZOLTAN_F90_INTERFACE
if HAVE_MPI
FDRIVER_SUBDIR += fdriver
endif
endif
endif

#SUBDIRS = . 
if BUILD_TESTS
SUBDIRS = . driver $(FDRIVER_SUBDIR)
else
SUBDIRS = .
endif


# ------------------------------------------------------------------------
# Files to be deleted by 'make maintainer-clean'
# ------------------------------------------------------------------------

MAINTAINERCLEANFILES = Makefile.in

MOSTLYCLEANFILES = *.mod

# ------------------------------------------------------------------------
# Explicit FORTRAN dependencies to allow parallel builds (make -j)
# KDD:  4/28/15  Specifying the dependencies does not seem to help 
# KDD:           the parallel F90 builds (although articles online 
# KDD:           indicate it should).  I'll comment them out
# KDD:           and use .NOTPARALLEL: instead.
# ------------------------------------------------------------------------

#if BUILD_ZOLTAN_F90_INTERFACE
#if F90_MODULE_LOWERCASE
#zoltan.mod:  zoltan_types.mod zoltan_user_data.mod $(srcdir)/fort/fwrap.f90 \
#             $(EXTRA_FORTRAN_H)
#zoltan_user_data.mod:  $(srcdir)/fort/zoltan_user_data.f90 zoltan_types.mod
#zoltan_types.mod:  $(srcdir)/fort/ztypes.f90
#else
#ZOLTAN.mod:  ZOLTAN_TYPES.mod ZOLTAN_USER_DATA.mod $(srcdir)/fort/fwrap.f90 \
#             $(EXTRA_FORTRAN_H)
#ZOLTAN_USER_DATA.mod:  $(srcdir)/fort/zoltan_user_data.f90 ZOLTAN_TYPES.mod
#ZOLTAN_TYPES.mod:  $(srcdir)/fort/ztypes.f90
#endif
#endif

if BUILD_ZOLTAN_F90_INTERFACE
.NOTPARALLEL:
endif
