
INCLUDE(TribitsLibraryMacros)

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

ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_OVIS)
SET(HAVE_OVIS ${${PACKAGE_NAME}_ENABLE_OVIS})

ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_Zlib)
SET(HAVE_GZIP ${${PACKAGE_NAME}_ENABLE_Zlib})

ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_PURIFY)
SET(HAVE_PURIFY ${${PACKAGE_NAME}_ENABLE_PURIFY})

TRIBITS_CONFIGURE_FILE(${PACKAGE_NAME}_config.h)

#
# B) Define the header and source files (and directories)
#

#
# src
#

SET(HEADERS "")
SET(NONINSTALLHEADERS "")
SET(SOURCES "")
SET(MODULES "")

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
  include/zoltan_eval.h
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
  coloring/bucket.h
  )

APPEND_SET(SOURCES
  coloring/coloring.c
  coloring/color_test.c
  coloring/g2l_hash.c
  coloring/bucket.c
  )

#
# src/graph
#

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/graph)

APPEND_SET(NOINSTALLHEADERS
  graph/graph.h
  graph/graph_const.h
  )

APPEND_SET(SOURCES
  graph/graph.c
  )


#
# src/ha
#

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/ha)

APPEND_SET(NOINSTALLHEADERS
  ha/ha_const.h
  ha/ha_ovis.h
  )

APPEND_SET(SOURCES
  ha/divide_machine.c
  ha/get_processor_name.c
  ha/ha_ovis.c
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
# src/matrix
#

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/matrix)

APPEND_SET(NOINSTALLHEADERS
  matrix/zoltan_matrix.h
  )

APPEND_SET(SOURCES
  matrix/matrix_build.c
  matrix/matrix_distribute.c
  matrix/matrix_operations.c
  matrix/matrix_sym.c
  matrix/matrix_utils.c
  )

#
# src/order
#

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/order)

APPEND_SET(NOINSTALLHEADERS
  order/order_const.h
  )

APPEND_SET(SOURCES
  order/hsfcOrder.c
  order/order.c
  order/order_struct.c
  order/order_tools.c
  order/perm.c
  )

IF (${PACKAGE_NAME}_ENABLE_HUND)
  APPEND_SET(SOURCES
    order/hund.c
    )
ENDIF()
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
# src/tpls
#

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/tpls)

APPEND_SET(NOINSTALLHEADERS
  tpls/third_library_const.h
  tpls/third_library.h
  tpls/third_library_tools.h
  tpls/graph_util.h
  )

APPEND_SET(SOURCES
  tpls/build_graph.c
  tpls/postprocessing.c
  tpls/preprocessing.c
  tpls/scatter_graph.c
  tpls/third_library.c
  tpls/verify_graph.c
  )

IF (${PACKAGE_NAME}_ENABLE_ParMETIS OR ${PACKAGE_NAME}_ENABLE_METIS)

  APPEND_SET(NOINSTALLHEADERS
    tpls/parmetis_interface.h
    )

  APPEND_SET(SOURCES
    tpls/parmetis_interface.c
    )

ENDIF()

IF (${PACKAGE_NAME}_ENABLE_Scotch)

  APPEND_SET(NOINSTALLHEADERS
  tpls/scotch_interface.h
    )

  APPEND_SET(SOURCES
  tpls/scotch_interface.c
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
  simple/cyclic.c
  simple/random.c
  )

APPEND_SET(NOINSTALLHEADERS
  simple/simple_const.h
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
  Utilities/DDirectory/DD_Memory.h
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
  Utilities/DDirectory/DD_Memory.c
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
  zz/zz_hash.h
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
  zz/zz_map.c
  zz/zz_heap.c
  zz/zz_init.c
  zz/zz_obj_list.c
  zz/zz_rand.c
  zz/zz_set_fn.c
  zz/zz_sort.c
  zz/zz_struct.c
  zz/zz_back_trace.c
  zz/zz_util.c
  )


#
# src/fort
#

IF(${PROJECT_NAME}_ENABLE_Fortran AND BUILD_ZOLTAN_F90_INTERFACE)
  INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/fort)

  APPEND_SET(NOINSTALLHEADERS
    fort/cwrap_fmangle.h
    fort/set_borderobjlist.fn
    fort/set_borderobjlist.if
    fort/set_childlist.fn
    fort/set_childlist.if
    fort/set_childweight.fn
    fort/set_childweight.if
    fort/set_coarseobjlist.fn
    fort/set_coarseobjlist.if
    fort/set_edgelist.fn
    fort/set_edgelist.if
    fort/set_edgelistmulti.fn
    fort/set_edgelistmulti.if
    fort/set_firstborderobj.fn
    fort/set_firstborderobj.if
    fort/set_firstcoarseobj.fn
    fort/set_firstcoarseobj.if
    fort/set_firstobj.fn
    fort/set_firstobj.if
    fort/set_fixedobjlist.fn
    fort/set_fixedobjlist.if
    fort/set_geom.fn
    fort/set_geom.if
    fort/set_geommulti.fn
    fort/set_geommulti.if
    fort/set_hgcs.fn
    fort/set_hgcs.if
    fort/set_hgedgeweights.fn
    fort/set_hgedgeweights.if
    fort/set_hgsizecs.fn
    fort/set_hgsizecs.if
    fort/set_hgsizeedgeweights.fn
    fort/set_hgsizeedgeweights.if
    fort/set_hiermethod.fn
    fort/set_hiermethod.if
    fort/set_hiernumlevels.fn
    fort/set_hiernumlevels.if
    fort/set_hierpartition.fn
    fort/set_hierpartition.if
    fort/set_midmigrate.fn
    fort/set_midmigrate.if
    fort/set_midmigratepp.fn
    fort/set_midmigratepp.if
    fort/set_nextborderobj.fn
    fort/set_nextborderobj.if
    fort/set_nextcoarseobj.fn
    fort/set_nextcoarseobj.if
    fort/set_nextobj.fn
    fort/set_nextobj.if
    fort/set_numborderobj.fn
    fort/set_numborderobj.if
    fort/set_numchild.fn
    fort/set_numchild.if
    fort/set_numcoarseobj.fn
    fort/set_numcoarseobj.if
    fort/set_numedges.fn
    fort/set_numedges.if
    fort/set_numedgesmulti.fn
    fort/set_numedgesmulti.if
    fort/set_numfixedobj.fn
    fort/set_numfixedobj.if
    fort/set_numgeom.fn
    fort/set_numgeom.if
    fort/set_numobj.fn
    fort/set_numobj.if
    fort/set_objlist.fn
    fort/set_objlist.if
    fort/set_objsize.fn
    fort/set_objsize.if
    fort/set_objsizemulti.fn
    fort/set_objsizemulti.if
    fort/set_packobj.fn
    fort/set_packobj.if
    fort/set_packobjmulti.fn
    fort/set_packobjmulti.if
    fort/set_partition.fn
    fort/set_partition.if
    fort/set_partitionmulti.fn
    fort/set_partitionmulti.if
    fort/set_postmigrate.fn
    fort/set_postmigrate.if
    fort/set_postmigratepp.fn
    fort/set_postmigratepp.if
    fort/set_premigrate.fn
    fort/set_premigrate.if
    fort/set_premigratepp.fn
    fort/set_premigratepp.if
    fort/set_unpackobj.fn
    fort/set_unpackobj.if
    fort/set_unpackobjmulti.fn
    fort/set_unpackobjmulti.if
    )

  APPEND_SET(SOURCES 
    fort/ztypes.f90
    fort/zoltan_user_data.f90
    fort/fwrap.f90
    fort/cwrap.c
    )

#  FILE(GLOB MODULES ${CMAKE_CURRENT_BINARY_DIR}/*.mod)
#  MESSAGE("KDDKDD MODULES " ${MODULES})

  APPEND_SET(MODULES
    ${CMAKE_CURRENT_BINARY_DIR}/zoltan.mod
    ${CMAKE_CURRENT_BINARY_DIR}/zoltan_types.mod
    ${CMAKE_CURRENT_BINARY_DIR}/zoltan_user_data.mod
    )
    
ENDIF()


#
# C) Define the targets for package's library(s)
#

IF (NOT TPL_ENABLE_MPI)
  SET(DEPLIBS simpi)
ENDIF()

# Add libm linking for ceil, floor and such.
LIST(APPEND DEPLIBS m)

TRIBITS_ADD_LIBRARY(
  zoltan
  HEADERS ${HEADERS}
  NOINSTALLHEADERS ${NOINSTALLHEADERS}
  SOURCES ${SOURCES}
  DEPLIBS ${DEPLIBS}
  )

INSTALL(FILES ${MODULES} DESTINATION ${CMAKE_INSTALL_PREFIX}/include)

