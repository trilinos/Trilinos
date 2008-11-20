
INCLUDE(TrilinosPackageLibraryHelpers)

#
# A) Package-specific configuration options
#

TRILINOS_PACKAGE_CONFIGURE_FILE(${PACKAGE_NAME}_config.h)

#
# B) Define the header and source files (and directories)
#

#
# src
#

SET(HEADERS "")
SET(NONINSTALLHEADERS "")
SET(SOURCES "")

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})

SET(HEADERS ${HEADERS}
  ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}_config.h
  )

#
# src/include
#

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/include)

SET(HEADERS ${HEADERS}
  include/lbi_const.h
  include/zoltan_align.h
  include/zoltan_cpp.h
  include/zoltan.h
  include/zoltan_types.h
  include/zoltan_comm.h
  include/zoltan_comm_cpp.h
  include/zoltan_timer_cpp.h
  include/zoltan_timer.h
  include/zoltan_dd.h
  include/zoltan_dd_cpp.h
  include/zoltan_mem.h  
  )

#
# src/all
#

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/all)

SET(NOINSTALLHEADERS ${NOINSTALLHEADERS}
  all/all_allo_const.h
  )

SET(SOURCES ${SOURCES}
  all/all_allo.c
  )

#
# src/coloring
#

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/coloring)

SET(NOINSTALLHEADERS ${NOINSTALLHEADERS}
  coloring/coloring_const.h
  coloring/coloring.h
  coloring/g2l_hash.h
  )

SET(SOURCES ${SOURCES}
  coloring/coloring.c
  coloring/color_test.c
  coloring/g2l_hash.c
  )

#
# src/ha
#

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/ha)

SET(NOINSTALLHEADERS ${NOINSTALLHEADERS}
  ha/ha_const.h
  ha/ha_drum.h
  )

SET(SOURCES ${SOURCES}
  ha/build_machine_desc.c
  ha/divide_machine.c
  ha/get_processor_name.c
  ha/ha_drum.c
  )

#
# src/hier
#

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/hier)

SET(NOINSTALLHEADERS ${NOINSTALLHEADERS}
  hier/hier.h
  )

SET(SOURCES ${SOURCES}
  hier/hier.c
  hier/hier_free_struct.c
  )

#
# src/hsfc
#

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/hsfc)

SET(NOINSTALLHEADERS ${NOINSTALLHEADERS}
  hsfc/hsfc.h
  hsfc/hsfc_const.h
  hsfc/hsfc_hilbert_const.h
  )

SET(SOURCES ${SOURCES}
  hsfc/hsfc_box_assign.c
  hsfc/hsfc.c
  hsfc/hsfc_hilbert.c
  hsfc/hsfc_point_assign.c
  )

#
# src/lb
#

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/lb)

SET(NOINSTALLHEADERS ${NOINSTALLHEADERS}
  lb/lb_const.h
  lb/lb_init_const.h
  )

SET(SOURCES ${SOURCES}
  lb/lb_balance.c
  lb/lb_box_assign.c
  lb/lb_copy.c
  lb/lb_eval.c
  lb/lb_free.c
  lb/lb_init.c
  lb/lb_invert.c
  lb/lb_migrate.c
  lb/lb_part2proc.c
  lb/lb_point_assign.c
  lb/lb_remap.c
  lb/lb_set_fn.c
  lb/lb_set_method.c
  lb/lb_set_part_sizes.c
  )

#
# src/oct
#

IF (${PACKAGE_NAME}_ENABLE_OCT)

  INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/oct)
  
  SET(NOINSTALLHEADERS ${NOINSTALLHEADERS}
    oct/costs_const.h
    oct/costs.h
    oct/dfs_const.h
    oct/dfs.h
    oct/migoct_const.h
    oct/migreg_const.h
    oct/migreg.h
    oct/migtags_const.h
    oct/msg_const.h
    oct/octant_const.h
    oct/octree_const.h
    oct/octupdate_const.h
    oct/octupdate.h
    oct/oct_util_const.h
    oct/rootlist_const.h
    oct/SFC.h
    )
  
  SET(SOURCES ${SOURCES}
    oct/costs.c
    oct/dfs.c
    oct/migoct.c
    oct/migreg.c
    oct/migtags.c
    oct/msg.c
    oct/octant.c
    oct/oct_plot.c
    oct/octree.c
    oct/octupdate.c
    oct/oct_util.c
    oct/output.c
    oct/rootlist.c
    )

ENDIF()

#
# src/order
#

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/order)

SET(NOINSTALLHEADERS ${NOINSTALLHEADERS}
  order/order_const.h
  )

SET(SOURCES ${SOURCES}
  order/order.c
  order/order_struct.c
  order/order_tools.c
  order/perm.c
  )

#
# src/par
#

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/par)

SET(NOINSTALLHEADERS ${NOINSTALLHEADERS}
  par/par_average_const.h
  par/par_bisect_const.h
  par/par_const.h
  par/par_median_const.h
  par/par_tflops_special_const.h
  )

SET(SOURCES ${SOURCES}
  par/par_average.c
  par/par_bisect.c
  par/par_median.c
  par/par_median_randomized.c
  par/par_stats.c
  par/par_sync.c
  par/par_tflops_special.c
  )

#
# src/params
#

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/params)

SET(NOINSTALLHEADERS ${NOINSTALLHEADERS}
  params/key_params.h
  params/params_const.h
  )

SET(SOURCES ${SOURCES}
  params/assign_param_vals.c
  params/bind_param.c
  params/check_param.c
  params/free_params.c
  params/key_params.c
  params/print_params.c
  params/set_param.c
  )



#
# src/parmetis
#

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/parmetis)

SET(NOINSTALLHEADERS ${NOINSTALLHEADERS}
  parmetis/third_library_const.h
  parmetis/third_library.h
  parmetis/third_library_tools.h
  )

SET(SOURCES ${SOURCES}
  parmetis/build_graph.c
  parmetis/postprocessing.c
  parmetis/preprocessing.c
  parmetis/scatter_graph.c
  parmetis/third_library.c
  parmetis/verify_graph.c
  )

IF (${PACKAGE_NAME}_ENABLE_ParMETIS)

  SET(NOINSTALLHEADERS ${NOINSTALLHEADERS}
    parmetis/parmetis_interface.h
    )
  
  SET(SOURCES ${SOURCES}
    parmetis/parmetis_interface.c
    )

ENDIF()

IF (${PACKAGE_NAME}_ENABLE_Scotch)

  SET(NOINSTALLHEADERS ${NOINSTALLHEADERS}
  parmetis/scotch_interface.h
    )
  
  SET(SOURCES ${SOURCES}
  parmetis/scotch_interface.c
    )

ENDIF()

#
# src/phg
#

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/phg)

SET(NOINSTALLHEADERS ${NOINSTALLHEADERS}
  phg/phg_comm.h
  phg/phg_const.h
  phg/phg_distrib.h
  phg/phg.h
  phg/phg_hypergraph.h
  phg/phg_util.h
  )

SET(SOURCES ${SOURCES}
  phg/phg_build.c
  phg/phg_build_calls.c
  phg/phg.c
  phg/phg_coarse.c
  phg/phg_comm.c
  phg/phg_distrib.c
  phg/phg_gather.c
  phg/phg_hypergraph.c
  phg/phg_match.c
  phg/phg_order.c
  phg/phg_parkway.c
  phg/phg_patoh.c
  phg/phg_plot.c
  phg/phg_rdivide.c
  phg/phg_refinement.c
  phg/phg_scale.c
  phg/phg_serialpartition.c
  phg/phg_util.c
  phg/phg_Vcycle.c
  )

#
# src/rcb
#

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/rcb)

SET(NOINSTALLHEADERS ${NOINSTALLHEADERS}
  rcb/create_proc_list_const.h
  rcb/rcb_const.h
  rcb/rcb.h
  rcb/rib_const.h
  rcb/rib.h
  rcb/shared.h
  )

SET(SOURCES ${SOURCES}
  rcb/box_assign.c
  rcb/create_proc_list.c
  rcb/inertial1d.c
  rcb/inertial2d.c
  rcb/inertial3d.c
  rcb/point_assign.c
  rcb/rcb_box.c
  rcb/rcb.c
  rcb/rcb_util.c
  rcb/rib.c
  rcb/rib_util.c
  rcb/shared.c
  )

#
# src/reftree
#

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/reftree)

SET(NOINSTALLHEADERS ${NOINSTALLHEADERS}
  reftree/reftree_const.h
  reftree/reftree.h
  )

SET(SOURCES ${SOURCES}
  reftree/reftree_build.c
  reftree/reftree_coarse_path.c
  reftree/reftree_hash.c
  reftree/reftree_part.c
  )

#
# src/simple
#

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/simple)

SET(SOURCES ${SOURCES}
  simple/block.c
  simple/random.c
  )

#
# src/timer
#

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/timer)

SET(NOINSTALLHEADERS ${NOINSTALLHEADERS}
  timer/timer_const.h
  )

SET(SOURCES ${SOURCES}
  timer/timer_params.c
  )

#
# src/Utilities
#

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/Utilities/Communication)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/Utilities/DDirectory)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/Utilities/Timer)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/Utilities/shared)

SET(NOINSTALLHEADERS ${NOINSTALLHEADERS}
  Utilities/Communication/comm.h
  Utilities/DDirectory/DD.h
  Utilities/Timer/timer.h
  Utilities/shared/zoltan_util.h
  Utilities/shared/zoltan_id.h
  )

SET(SOURCES ${SOURCES}
  Utilities/Communication/comm_exchange_sizes.c
  Utilities/Communication/comm_invert_map.c
  Utilities/Communication/comm_do.c
  Utilities/Communication/comm_do_reverse.c
  Utilities/Communication/comm_info.c
  Utilities/Communication/comm_create.c
  Utilities/Communication/comm_resize.c
  Utilities/Communication/comm_sort_ints.c
  Utilities/Communication/comm_destroy.c
  Utilities/Communication/comm_invert_plan.c
  Utilities/Timer/zoltan_timer.c
  Utilities/Timer/timer.c
  Utilities/DDirectory/DD_Find.c
  Utilities/DDirectory/DD_Destroy.c
  Utilities/DDirectory/DD_Set_Neighbor_Hash_Fn3.c
  Utilities/DDirectory/DD_Remove.c
  Utilities/DDirectory/DD_Create.c
  Utilities/DDirectory/DD_Update.c
  Utilities/DDirectory/DD_Stats.c
  Utilities/DDirectory/DD_Hash2.c
  Utilities/DDirectory/DD_Print.c
  Utilities/DDirectory/DD_Set_Neighbor_Hash_Fn2.c
  Utilities/DDirectory/DD_Set_Hash_Fn.c
  Utilities/DDirectory/DD_Set_Neighbor_Hash_Fn1.c
  Utilities/Memory/mem.c
  Utilities/shared/zoltan_align.c
  Utilities/shared/zoltan_id.c
  )

#
# src/zz
#

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/zz)

SET(NOINSTALLHEADERS ${NOINSTALLHEADERS}
  zz/zz_const.h
  zz/zz_heap.h
  zz/zz_id_const.h
  zz/zz_rand.h
  zz/zz_sort.h
  zz/zz_util_const.h
  )

SET(SOURCES ${SOURCES}
  zz/zz_coord.c
  zz/zz_gen_files.c
  zz/zz_hash.c
  zz/zz_heap.c
  zz/zz_init.c
  zz/zz_obj_list.c
  zz/zz_rand.c
  zz/zz_set_fn.c
  zz/zz_sort.c
  zz/zz_struct.c
  zz/zz_util.c
  )


#
# C) Define the targets for package's library(s)
#

IF (NOT Trilinos_ENABLE_MPI)
  SET(DEPLIBS simpi)
ENDIF()

TRILINOS_PACKAGE_ADD_LIBRARY(
  zoltan
  HEADERS ${HEADERS}
  NOINSTALLHEADERS ${NOINSTALLHEADERS}
  SOURCES ${SOURCES}
  DEPLIBS ${DEPLIBS}
  )

#
# D) Export the dependancy variables of this package for other packages to use
#

TRILINOS_PACKAGE_EXPORT_DEPENDENCY_VARIABLES()

ADD_SUBDIRECTORY(driver)
