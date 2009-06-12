
INCLUDE(PackageLibraryMacros)

#
# A) Package-specific configuration options
#

ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_METIS)
SET(HAVE_METIS ${${PACKAGE_NAME}_ENABLE_METIS})

ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_ParMETIS)
SET(HAVE_PARMETIS ${${PACKAGE_NAME}_ENABLE_ParMETIS})

ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_PaToH)
SET(HAVE_PATOH ${${PACKAGE_NAME}_ENABLE_PaToH})

ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_Scotch)
SET(HAVE_SCOTCH ${${PACKAGE_NAME}_ENABLE_Scotch})

ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_Zlib)
SET(HAVE_GZIP ${${PACKAGE_NAME}_ENABLE_Zlib})

PACKAGE_CONFIGURE_FILE(${PACKAGE_NAME}_config.h)

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

APPEND_SET(HEADERS
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

APPEND_SET(NOINSTALLHEADERS
  all/all_allo_const.h
  )

APPEND_SET(SOURCES
  all/all_allo.c
  )

#
# src/coloring
#

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/coloring)

APPEND_SET(NOINSTALLHEADERS
  coloring/coloring_const.h
  coloring/coloring.h
  coloring/g2l_hash.h
  )

APPEND_SET(SOURCES
  coloring/coloring.c
  coloring/color_test.c
  coloring/g2l_hash.c
  )

#
# src/ha
#

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/ha)

APPEND_SET(NOINSTALLHEADERS
  ha/ha_const.h
  ha/ha_drum.h
  )

APPEND_SET(SOURCES
  ha/build_machine_desc.c
  ha/divide_machine.c
  ha/get_processor_name.c
  ha/ha_drum.c
  )

#
# src/hier
#

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/hier)

APPEND_SET(NOINSTALLHEADERS
  hier/hier.h
  )

APPEND_SET(SOURCES
  hier/hier.c
  hier/hier_free_struct.c
  )

#
# src/hsfc
#

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/hsfc)

APPEND_SET(NOINSTALLHEADERS
  hsfc/hsfc.h
  hsfc/hsfc_const.h
  hsfc/hsfc_hilbert_const.h
  )

APPEND_SET(SOURCES
  hsfc/hsfc_box_assign.c
  hsfc/hsfc.c
  hsfc/hsfc_hilbert.c
  hsfc/hsfc_point_assign.c
  )

#
# src/lb
#

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/lb)

APPEND_SET(NOINSTALLHEADERS
  lb/lb_const.h
  lb/lb_init_const.h
  )

APPEND_SET(SOURCES
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

  APPEND_SET(NOINSTALLHEADERS
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

  APPEND_SET(SOURCES
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

APPEND_SET(NOINSTALLHEADERS
  order/order_const.h
  )

APPEND_SET(SOURCES
  order/order.c
  order/order_struct.c
  order/order_tools.c
  order/perm.c
  )

#
# src/par
#

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/par)

APPEND_SET(NOINSTALLHEADERS
  par/par_average_const.h
  par/par_bisect_const.h
  par/par_const.h
  par/par_median_const.h
  par/par_tflops_special_const.h
  )

APPEND_SET(SOURCES
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

APPEND_SET(NOINSTALLHEADERS
  params/key_params.h
  params/params_const.h
  )

APPEND_SET(SOURCES
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

APPEND_SET(NOINSTALLHEADERS
  parmetis/third_library_const.h
  parmetis/third_library.h
  parmetis/third_library_tools.h
  )

APPEND_SET(SOURCES
  parmetis/build_graph.c
#  parmetis/build_graph2.c
  parmetis/postprocessing.c
  parmetis/preprocessing.c
  parmetis/scatter_graph.c
  parmetis/third_library.c
  parmetis/verify_graph.c
  )

IF (${PACKAGE_NAME}_ENABLE_ParMETIS OR ${PACKAGE_NAME}_ENABLE_METIS)

  APPEND_SET(NOINSTALLHEADERS
    parmetis/parmetis_interface.h
    )

  APPEND_SET(SOURCES
    parmetis/parmetis_interface.c
    )

ENDIF()

IF (${PACKAGE_NAME}_ENABLE_Scotch)

  APPEND_SET(NOINSTALLHEADERS
  parmetis/scotch_interface.h
    )

  APPEND_SET(SOURCES
  parmetis/scotch_interface.c
    )

ENDIF()

#
# src/phg
#

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/phg)

APPEND_SET(NOINSTALLHEADERS
  phg/phg_comm.h
  phg/phg_const.h
  phg/phg_distrib.h
  phg/phg.h
  phg/phg_lookup.h
  phg/phg_verbose.h
  phg/phg_hypergraph.h
  phg/phg_tree.h
  phg/phg_util.h
  )

APPEND_SET(SOURCES
  phg/phg_build.c
  phg/phg_build_calls.c
  phg/phg.c
  phg/phg_lookup.c
  phg/phg_verbose.c
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
  phg/phg_tree.c
  phg/phg_two_ways.c
  phg/phg_util.c
  phg/phg_Vcycle.c
  )

#
# src/rcb
#

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/rcb)

APPEND_SET(NOINSTALLHEADERS
  rcb/create_proc_list_const.h
  rcb/rcb_const.h
  rcb/rcb.h
  rcb/rib_const.h
  rcb/rib.h
  rcb/shared.h
  )

APPEND_SET(SOURCES
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

APPEND_SET(NOINSTALLHEADERS
  reftree/reftree_const.h
  reftree/reftree.h
  )

APPEND_SET(SOURCES
  reftree/reftree_build.c
  reftree/reftree_coarse_path.c
  reftree/reftree_hash.c
  reftree/reftree_part.c
  )

#
# src/simple
#

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/simple)

APPEND_SET(SOURCES
  simple/block.c
  simple/random.c
  )

#
# src/timer
#

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/timer)

APPEND_SET(NOINSTALLHEADERS
  timer/timer_const.h
  )

APPEND_SET(SOURCES
  timer/timer_params.c
  )

#
# src/Utilities
#

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/Utilities/Communication)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/Utilities/DDirectory)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/Utilities/Timer)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/Utilities/shared)

APPEND_SET(NOINSTALLHEADERS
  Utilities/Communication/comm.h
  Utilities/DDirectory/DD.h
  Utilities/Timer/timer.h
  Utilities/shared/zoltan_util.h
  Utilities/shared/zoltan_id.h
  )

APPEND_SET(SOURCES
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

APPEND_SET(NOINSTALLHEADERS
  zz/zz_const.h
  zz/zz_heap.h
  zz/zz_id_const.h
  zz/zz_rand.h
  zz/zz_sort.h
  zz/zz_util_const.h
  )

APPEND_SET(SOURCES
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

IF (NOT TPL_ENABLE_MPI)
  SET(DEPLIBS simpi)
ENDIF()

PACKAGE_ADD_LIBRARY(
  zoltan
  HEADERS ${HEADERS}
  NOINSTALLHEADERS ${NOINSTALLHEADERS}
  SOURCES ${SOURCES}
  DEPLIBS ${DEPLIBS}
  )
