INCLUDE(TPLDeclareLibraries)

TPL_DECLARE_LIBRARIES( Boost
  REQUIRED_HEADERS boost/version.hpp boost/mpl/at.hpp
  )

# This broke trilinos configuration:
#  REQUIRED_LIBS_NAMES "program_options"
