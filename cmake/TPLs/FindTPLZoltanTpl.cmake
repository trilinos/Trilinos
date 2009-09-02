INCLUDE(TPLDeclareLibraries)

TPL_DECLARE_LIBRARIES( ZoltanTpl
  REQUIRED_HEADERS lbi_const.h zoltan_comm_cpp.h zoltan_cpp.h zoltan_dd.h zoltan_mem.h zoltan_timer.h zoltan_align.h zoltan_comm.h zoltan_dd_cpp.h zoltan.h zoltan_timer_cpp.h zoltan_types.h
  REQUIRED_LIBS_NAMES "zoltan"
  )
