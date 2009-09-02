INCLUDE(TPLDeclareLibraries)

TPL_DECLARE_LIBRARIES( ExodusII
  #REQUIRED_HEADERS exodus.h
  REQUIRED_HEADERS exodusII.h netcdf.h
  REQUIRED_LIBS_NAMES exoIIv2c netcdf
  )
