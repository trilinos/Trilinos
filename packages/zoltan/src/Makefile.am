# @HEADER
##############################################################################
# Zoltan Library for Parallel Applications                                   #
# Copyright (c) 2000,2001,2002, Sandia National Laboratories.                #
# This software is distributed under the GNU Lesser General Public License.  #
# For more info, see the README file in the top-level Zoltan directory.      #
##############################################################################
# @HEADER

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
	$(srcdir)/coloring/g2l_hash.h

COLORING = \
	$(srcdir)/coloring/coloring.c \
	$(srcdir)/coloring/color_test.c \
	$(srcdir)/coloring/g2l_hash.c

HA_H = \
	$(srcdir)/ha/ha_const.h \
	$(srcdir)/ha/ha_drum.h

HA = \
	$(srcdir)/ha/build_machine_desc.c \
	$(srcdir)/ha/divide_machine.c \
	$(srcdir)/ha/get_processor_name.c \
	$(srcdir)/ha/ha_drum.c

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
	$(srcdir)/include/zoltan_types.h \
	$(srcdir)/include/zoltan_comm.h \
	$(srcdir)/include/zoltan_comm_cpp.h \
	$(srcdir)/include/zoltan_timer_cpp.h \
	$(srcdir)/include/zoltan_timer.h \
	$(srcdir)/include/zoltan_dd.h \
	$(srcdir)/include/zoltan_dd_cpp.h \
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

EXTRA_OCT_H = \
	$(srcdir)/oct/costs_const.h \
	$(srcdir)/oct/costs.h \
	$(srcdir)/oct/dfs_const.h \
	$(srcdir)/oct/dfs.h \
	$(srcdir)/oct/migoct_const.h \
	$(srcdir)/oct/migreg_const.h \
	$(srcdir)/oct/migreg.h \
	$(srcdir)/oct/migtags_const.h \
	$(srcdir)/oct/msg_const.h \
	$(srcdir)/oct/octant_const.h \
	$(srcdir)/oct/octree_const.h \
	$(srcdir)/oct/octupdate_const.h \
	$(srcdir)/oct/octupdate.h \
	$(srcdir)/oct/oct_util_const.h \
	$(srcdir)/oct/rootlist_const.h \
	$(srcdir)/oct/SFC.h

EXTRA_OCT = \
	$(srcdir)/oct/costs.c \
	$(srcdir)/oct/dfs.c \
	$(srcdir)/oct/migoct.c \
	$(srcdir)/oct/migreg.c \
	$(srcdir)/oct/migtags.c \
	$(srcdir)/oct/msg.c \
	$(srcdir)/oct/octant.c \
	$(srcdir)/oct/oct_plot.c \
	$(srcdir)/oct/octree.c \
	$(srcdir)/oct/octupdate.c \
	$(srcdir)/oct/oct_util.c \
	$(srcdir)/oct/output.c \
	$(srcdir)/oct/rootlist.c

if BUILD_OCTREEPARTITIONING

OCT_H = $(EXTRA_OCT_H)

OCT = $(EXTRA_OCT)

else

OCT_H =

OCT =

endif

ORDER_H = \
	$(srcdir)/order/order_const.h

ORDER = \
	$(srcdir)/order/order.c \
	$(srcdir)/order/order_struct.c \
	$(srcdir)/order/order_tools.c \
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

PARMETIS_H = \
	$(srcdir)/parmetis/third_library_const.h \
	$(srcdir)/parmetis/third_library.h \
	$(srcdir)/parmetis/third_library_tools.h

PARMETIS = \
	$(srcdir)/parmetis/build_graph.c \
	$(srcdir)/parmetis/postprocessing.c \
	$(srcdir)/parmetis/preprocessing.c \
	$(srcdir)/parmetis/scatter_graph.c \
	$(srcdir)/parmetis/third_library.c \
	$(srcdir)/parmetis/verify_graph.c

EXTRA_PARMETIS_INTERFACE_H = \
	$(srcdir)/parmetis/parmetis_interface.h

EXTRA_PARMETIS_INTERFACE = \
	$(srcdir)/parmetis/parmetis_interface.c

if BUILD_PARMETIS

PARMETIS_INTERFACE_H = $(EXTRA_PARMETIS_INTERFACE_H)

PARMETIS_INTERFACE = $(EXTRA_PARMETIS_INTERFACE)

else

PARMETIS_INTERFACE_H =

PARMETIS_INTERFACE =

endif

EXTRA_SCOTCH_INTERFACE_H = \
	$(srcdir)/parmetis/scotch_interface.h

EXTRA_SCOTCH_INTERFACE = \
	$(srcdir)/parmetis/scotch_interface.c

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
	$(srcdir)/phg/phg_hypergraph.h \
	$(srcdir)/phg/phg_util.h

PHG = \
	$(srcdir)/phg/phg_build.c \
	$(srcdir)/phg/phg_build_calls.c \
	$(srcdir)/phg/phg.c \
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
	$(srcdir)/simple/random.c

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
	$(srcdir)/zz/zz_util_const.h

ZZ = \
	$(srcdir)/zz/zz_coord.c \
	$(srcdir)/zz/zz_gen_files.c \
	$(srcdir)/zz/zz_hash.c \
	$(srcdir)/zz/zz_heap.c \
	$(srcdir)/zz/zz_init.c \
	$(srcdir)/zz/zz_obj_list.c \
	$(srcdir)/zz/zz_rand.c \
	$(srcdir)/zz/zz_set_fn.c \
	$(srcdir)/zz/zz_sort.c \
	$(srcdir)/zz/zz_struct.c \
	$(srcdir)/zz/zz_util.c

# ------------------------------------------------------------------------
# Zoltan library specifications
# ------------------------------------------------------------------------
ZOLTAN_LIB = libzoltan.a

ZOLTAN_H = \
	$(ALL_H) \
	$(COLORING_H) \
	$(HA_H) \
	$(HIER_H) \
	$(HSFC_H) \
	$(INCLUDE_H) \
	$(LB_H) \
	$(OCT_H) \
	$(ORDER_H) \
	$(PAR_H) \
	$(PARAMS_H) \
	$(PARMETIS_H) \
	$(PARMETIS_INTERFACE_H) \
	$(PHG_H) \
	$(RCB_H) \
	$(REFTREE_H) \
	$(SCOTCH_INTERFACE_H) \
	$(TIMER_H) \
	$(UTILITIES_H) \
	$(ZZ_H)

libzoltan_a_SOURCES = \
	$(ALL) \
	$(COLORING) \
	$(HA) \
	$(HIER) \
	$(HSFC) \
	$(LB) \
	$(OCT) \
	$(ORDER) \
	$(PAR) \
	$(PARAMS) \
	$(PARMETIS) \
	$(PARMETIS_INTERFACE) \
	$(PHG) \
	$(RCB) \
	$(REFTREE) \
	$(SCOTCH_INTERFACE) \
	$(SIMPLE) \
	$(TIMER) \
	$(UTILITIES) \
	$(ZZ)

EXTRA_libzoltan_a_SOURCES = \
	$(EXTRA_OCT_H) \
	$(EXTRA_OCT) \
	$(EXTRA_PARMETIS_INTERFACE_H) \
	$(EXTRA_PARMETIS_INTERFACE) \
	$(EXTRA_SCOTCH_INTERFACE_H) \
	$(EXTRA_SCOTCH_INTERFACE)

include $(top_builddir)/Makefile.export.zoltan

if USING_GNUMAKE
EXPORT_INCLUDES = $(shell $(PERL_EXE) $(top_srcdir)/config/strip_dup_incl_paths.pl $(ZOLTAN_LIBRARY_INCLUDES))
else
EXPORT_INCLUDES = $(ZOLTAN_LIBRARY_INCLUDES)
endif

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
# Other files we want in the tarball
# ------------------------------------------------------------------------

FORTRAN_INTERFACE = \
    $(srcdir)/fort/README \
    $(srcdir)/fort/cwrap.c \
    $(srcdir)/fort/fwrap.fpp \
    $(srcdir)/fort/lb_user_const.f90 \
    $(srcdir)/fort/lbfn.fpp \
    $(srcdir)/fort/lbif.fpp \
    $(srcdir)/fort/makefile \
    $(srcdir)/fort/set_borderobjlist.fn \
    $(srcdir)/fort/set_borderobjlist.fn.lbfn \
    $(srcdir)/fort/set_borderobjlist.if \
    $(srcdir)/fort/set_borderobjlist.if.lbif \
    $(srcdir)/fort/set_childlist.fn \
    $(srcdir)/fort/set_childlist.fn.lbfn \
    $(srcdir)/fort/set_childlist.if \
    $(srcdir)/fort/set_childlist.if.lbif \
    $(srcdir)/fort/set_childweight.fn \
    $(srcdir)/fort/set_childweight.fn.lbfn \
    $(srcdir)/fort/set_childweight.if \
    $(srcdir)/fort/set_childweight.if.lbif \
    $(srcdir)/fort/set_coarseobjlist.fn \
    $(srcdir)/fort/set_coarseobjlist.fn.lbfn \
    $(srcdir)/fort/set_coarseobjlist.if \
    $(srcdir)/fort/set_coarseobjlist.if.lbif \
    $(srcdir)/fort/set_edgelist.fn \
    $(srcdir)/fort/set_edgelist.fn.lbfn \
    $(srcdir)/fort/set_edgelist.if \
    $(srcdir)/fort/set_edgelist.if.lbif \
    $(srcdir)/fort/set_edgelistmulti.fn \
    $(srcdir)/fort/set_edgelistmulti.if \
    $(srcdir)/fort/set_firstborderobj.fn \
    $(srcdir)/fort/set_firstborderobj.fn.lbfn \
    $(srcdir)/fort/set_firstborderobj.if \
    $(srcdir)/fort/set_firstborderobj.if.lbif \
    $(srcdir)/fort/set_firstcoarseobj.fn \
    $(srcdir)/fort/set_firstcoarseobj.fn.lbfn \
    $(srcdir)/fort/set_firstcoarseobj.if \
    $(srcdir)/fort/set_firstcoarseobj.if.lbif \
    $(srcdir)/fort/set_firstobj.fn \
    $(srcdir)/fort/set_firstobj.fn.lbfn \
    $(srcdir)/fort/set_firstobj.if \
    $(srcdir)/fort/set_firstobj.if.lbif \
    $(srcdir)/fort/set_fixedobjlist.fn \
    $(srcdir)/fort/set_fixedobjlist.if \
    $(srcdir)/fort/set_geom.fn \
    $(srcdir)/fort/set_geom.fn.lbfn \
    $(srcdir)/fort/set_geom.if \
    $(srcdir)/fort/set_geom.if.lbif \
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
    $(srcdir)/fort/set_midmigrate.fn.lbfn \
    $(srcdir)/fort/set_midmigrate.if \
    $(srcdir)/fort/set_midmigrate.if.lbif \
    $(srcdir)/fort/set_midmigratepp.fn \
    $(srcdir)/fort/set_midmigratepp.if \
    $(srcdir)/fort/set_nextborderobj.fn \
    $(srcdir)/fort/set_nextborderobj.fn.lbfn \
    $(srcdir)/fort/set_nextborderobj.if \
    $(srcdir)/fort/set_nextborderobj.if.lbif \
    $(srcdir)/fort/set_nextcoarseobj.fn \
    $(srcdir)/fort/set_nextcoarseobj.fn.lbfn \
    $(srcdir)/fort/set_nextcoarseobj.if \
    $(srcdir)/fort/set_nextcoarseobj.if.lbif \
    $(srcdir)/fort/set_nextobj.fn \
    $(srcdir)/fort/set_nextobj.fn.lbfn \
    $(srcdir)/fort/set_nextobj.if \
    $(srcdir)/fort/set_nextobj.if.lbif \
    $(srcdir)/fort/set_numborderobj.fn \
    $(srcdir)/fort/set_numborderobj.fn.lbfn \
    $(srcdir)/fort/set_numborderobj.if \
    $(srcdir)/fort/set_numborderobj.if.lbif \
    $(srcdir)/fort/set_numchild.fn \
    $(srcdir)/fort/set_numchild.fn.lbfn \
    $(srcdir)/fort/set_numchild.if \
    $(srcdir)/fort/set_numchild.if.lbif \
    $(srcdir)/fort/set_numcoarseobj.fn \
    $(srcdir)/fort/set_numcoarseobj.fn.lbfn \
    $(srcdir)/fort/set_numcoarseobj.if \
    $(srcdir)/fort/set_numcoarseobj.if.lbif \
    $(srcdir)/fort/set_numedges.fn \
    $(srcdir)/fort/set_numedges.fn.lbfn \
    $(srcdir)/fort/set_numedges.if \
    $(srcdir)/fort/set_numedges.if.lbif \
    $(srcdir)/fort/set_numedgesmulti.fn \
    $(srcdir)/fort/set_numedgesmulti.if \
    $(srcdir)/fort/set_numfixedobj.fn \
    $(srcdir)/fort/set_numfixedobj.if \
    $(srcdir)/fort/set_numgeom.fn \
    $(srcdir)/fort/set_numgeom.fn.lbfn \
    $(srcdir)/fort/set_numgeom.if \
    $(srcdir)/fort/set_numgeom.if.lbif \
    $(srcdir)/fort/set_numobj.fn \
    $(srcdir)/fort/set_numobj.fn.lbfn \
    $(srcdir)/fort/set_numobj.if \
    $(srcdir)/fort/set_numobj.if.lbif \
    $(srcdir)/fort/set_objlist.fn \
    $(srcdir)/fort/set_objlist.fn.lbfn \
    $(srcdir)/fort/set_objlist.if \
    $(srcdir)/fort/set_objlist.if.lbif \
    $(srcdir)/fort/set_objsize.fn \
    $(srcdir)/fort/set_objsize.fn.lbfn \
    $(srcdir)/fort/set_objsize.if \
    $(srcdir)/fort/set_objsize.if.lbif \
    $(srcdir)/fort/set_objsizemulti.fn \
    $(srcdir)/fort/set_objsizemulti.if \
    $(srcdir)/fort/set_packobj.fn \
    $(srcdir)/fort/set_packobj.fn.lbfn \
    $(srcdir)/fort/set_packobj.if \
    $(srcdir)/fort/set_packobj.if.lbif \
    $(srcdir)/fort/set_packobjmulti.fn \
    $(srcdir)/fort/set_packobjmulti.if \
    $(srcdir)/fort/set_partition.fn \
    $(srcdir)/fort/set_partition.if \
    $(srcdir)/fort/set_partitionmulti.fn \
    $(srcdir)/fort/set_partitionmulti.if \
    $(srcdir)/fort/set_postmigrate.fn \
    $(srcdir)/fort/set_postmigrate.fn.lbfn \
    $(srcdir)/fort/set_postmigrate.if \
    $(srcdir)/fort/set_postmigrate.if.lbif \
    $(srcdir)/fort/set_postmigratepp.fn \
    $(srcdir)/fort/set_postmigratepp.if \
    $(srcdir)/fort/set_premigrate.fn \
    $(srcdir)/fort/set_premigrate.fn.lbfn \
    $(srcdir)/fort/set_premigrate.if \
    $(srcdir)/fort/set_premigrate.if.lbif \
    $(srcdir)/fort/set_premigratepp.fn \
    $(srcdir)/fort/set_premigratepp.if \
    $(srcdir)/fort/set_unpackobj.fn \
    $(srcdir)/fort/set_unpackobj.fn.lbfn \
    $(srcdir)/fort/set_unpackobj.if \
    $(srcdir)/fort/set_unpackobj.if.lbif \
    $(srcdir)/fort/set_unpackobjmulti.fn \
    $(srcdir)/fort/set_unpackobjmulti.if \
    $(srcdir)/fort/sppr.f90 \
    $(srcdir)/fort/spprinc.64bit \
    $(srcdir)/fort/spprinc.absoft \
    $(srcdir)/fort/spprinc.most \
    $(srcdir)/fort/spprinc.nas \
    $(srcdir)/fort/spprinc.sun \
    $(srcdir)/fort/zoltan_user_data.f90 \
    $(srcdir)/fort/ztypes.fpp \
    $(srcdir)/fort/ztypes.lbif
    
FORTRAN_DRIVER = \
    $(srcdir)/fdriver/Makefile.am \
    $(srcdir)/fdriver/Makefile.in \
    $(srcdir)/fdriver/README.mpich \
    $(srcdir)/fdriver/farg_nagf90.f \
    $(srcdir)/fdriver/farg_nagf95.f \
    $(srcdir)/fdriver/farg_nas.f \
    $(srcdir)/fdriver/farg_typical.f \
    $(srcdir)/fdriver/fdr_chaco_io.f90 \
    $(srcdir)/fdriver/fdr_const.f90 \
    $(srcdir)/fdriver/fdr_input.f90 \
    $(srcdir)/fdriver/fdr_loadbal.f90 \
    $(srcdir)/fdriver/fdr_main.f90 \
    $(srcdir)/fdriver/fdr_migrate.f90 \
    $(srcdir)/fdriver/fdr_mm_io.f90 \
    $(srcdir)/fdriver/fdr_param_file.f90 \
    $(srcdir)/fdriver/fdr_sort.f90 \
    $(srcdir)/fdriver/lampmpi.f90 \
    $(srcdir)/fdriver/makefile \
    $(srcdir)/fdriver/mmio.f \
    $(srcdir)/fdriver/mpi_h.f \
    $(srcdir)/fdriver/mpistring.c \
    $(srcdir)/fdriver/zoltan_user_data.f90

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
    $(srcdir)/Utilities/Config/Config.nistlinux \
    $(srcdir)/Utilities/Config/Config.nistsun \
    $(srcdir)/Utilities/Config/Config.octopiLamPureTest \
    $(srcdir)/Utilities/Config/Config.octopiLamTest \
    $(srcdir)/Utilities/Config/Config.octopiPureTest \
    $(srcdir)/Utilities/Config/Config.bgl \
    $(srcdir)/Utilities/Config/Config.octopiTest \
    $(srcdir)/Utilities/Config/Config.odin \
    $(srcdir)/Utilities/Config/Config.osf1 \
    $(srcdir)/Utilities/Config/Config.dec \
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
    $(srcdir)/Utilities/Config/Config.icc \
    $(srcdir)/Utilities/Config/Config.icc-valor \
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
    $(srcdir)/Utilities/Config/Config.x86_freebsd-mpich \
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
    $(FORTRAN_INTERFACE) $(FORTRAN_DRIVER) \
    $(ZOLTAN_MATLAB) $(ZOLTAN_POSTPROCESS)

# ------------------------------------------------------------------------
# Special stuff to install in our special $(execincludedir)
# ------------------------------------------------------------------------

#np# replace new_package with the name of the package being autotool'ed here
# SPECIAL NOTE: New_Package_config.h is a machine-dependent file, so we need
# to install it in the machine-dependent directory. However, that is
# not a default installation directory, so we had to create it
# special.

# All Trilinos headers are now installed in the same directory
execincludedir = $(exec_prefix)/include
#np# replace new_package with the name of the package being autotool'ed here
#nodist_execinclude_HEADERS = Zoltan_config.h

if BUILD_FORTRAN90DRIVER
FDRIVER_SUBDIR = #fdriver
endif

SUBDIRS = . driver $(FDRIVER_SUBDIR)

# ------------------------------------------------------------------------
# Files to be deleted by 'make maintainer-clean'
# ------------------------------------------------------------------------

MAINTAINERCLEANFILES = Makefile.in 

