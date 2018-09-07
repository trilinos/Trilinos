INCLUDE(TrilinosCreateClientTemplateHeaders)

# Function that follows the Tpetra convention for mangling C++ types
# so that they can be used as C preprocessor macro arguments.
#
# TYPE_MANGLED_OUT [out] The mangled type name.
#
# TYPE_IN [in] The type to mangle.
FUNCTION(TPETRANEW_MANGLE_TEMPLATE_PARAMETER TYPE_MANGLED_OUT TYPE_IN)
  STRING(REPLACE "<" "0" TMP0 "${TYPE_IN}")
  STRING(REPLACE ">" "0" TMP1 "${TMP0}")
  STRING(REPLACE "::" "_" TMP2 "${TMP1}")
  # Spaces (as in "long long") get squished out.
  STRING(REPLACE " " "" TMP3 "${TMP2}")
  SET(${TYPE_MANGLED_OUT} ${TMP3} PARENT_SCOPE)
ENDFUNCTION(TPETRANEW_MANGLE_TEMPLATE_PARAMETER)

# Function that turns a valid Scalar (or Packet) template parameter
# into a macro name (all caps, with no white space and no punctuation
# other than underscore).
#
# NAME_OUT [out] The mangled type name.
#
# NAME_IN [in] The type to mangle.
FUNCTION(TPETRA_SCALAR_MACRO_NAME NAME_OUT NAME_IN)
  STRING(COMPARE EQUAL "${NAME_IN}" "__float128" IS_FLOAT128)
  IF(IS_FLOAT128)
    # __float128 is a special case; we remove the __ from the macro name.
    SET(${NAME_OUT} "FLOAT128" PARENT_SCOPE)
  ELSE()
    STRING(COMPARE EQUAL "${NAME_IN}" "std::complex<float>" IS_COMPLEX_FLOAT)
    IF(IS_COMPLEX_FLOAT)
      SET(${NAME_OUT} "COMPLEX_FLOAT" PARENT_SCOPE)
    ELSE()
      STRING(COMPARE EQUAL "${NAME_IN}" "std::complex<double>" IS_COMPLEX_DOUBLE)
      IF(IS_COMPLEX_DOUBLE)
        SET(${NAME_OUT} "COMPLEX_DOUBLE" PARENT_SCOPE)
      ELSE()
        # Convert to upper case, convert double colons to underscores,
        # and hope for the best.
        #
        # It would be nice if CMake were consistent about where output
        # arguments go.  Alas, this is not to be.  TOUPPER puts the
        # output argument last; REPLACE puts it after the search and
        # substitute strings, before the input string.
        STRING(TOUPPER "${NAME_IN}" TMP0)
        STRING(REPLACE "::" "_" TMP1 "${TMP0}")
        STRING(REPLACE " " "_" TMP2 "${TMP1}")
        SET(${NAME_OUT} ${TMP2} PARENT_SCOPE)
      ENDIF()
    ENDIF()
  ENDIF()
ENDFUNCTION(TPETRA_SCALAR_MACRO_NAME)

# Function to generate one .cpp file for the given Scalar (or Packet,
# in the case of DistObject) template argument, for explicit
# instantiation over that template parameter. This is meant to be called
# by TPETRANEW_PROCESS_ALL_S_TEMPLATES. This function takes the names
# already mangled, to avoid unnecessary string processing overhead.
#
# OUTPUT_FILE [out] Name of the generated .cpp file.
#
# TEMPLATE_FILE [in] Name of the input .tmpl "template" file.  This
#   function does string substitution in that file, using the input
#   arguments of this function.  For example, @SC_MACRO_EXPR@ (Scalar
#   macro expression) gets substituted for the value of this
#   function's SC_MACRO_EXPR input argument.
#
# CLASS_NAME [in] Name of the Tpetra class (without namespace
#   qualifiers; must live in the Tpetra namespace)
#
# CLASS_MACRO_NAME [in] Name of the Tpetra class, suitably mangled for
#   use in a macro name.
#
# SC_MANGLED_NAME [in] Name of the Scalar (SC) or Packet type, mangled
#   for use as a macro argument (e.g., spaces and colons removed).
#
# SC_MACRO_EXPR [in] Expression that asks Tpetra whether the given
#   Scalar (SC) type is supported.
#
FUNCTION(TPETRANEW_PROCESS_ONE_S_TEMPLATE OUTPUT_FILE TEMPLATE_FILE CLASS_NAME CLASS_MACRO_NAME SC_MANGLED_NAME SC_MACRO_EXPR)
  STRING(REPLACE "ETI_SC.tmpl" "${CLASS_NAME}_${SC_MACRO_NAME}.cpp" OUT_FILE "${TEMPLATE_FILE}")
  CONFIGURE_FILE("${TEMPLATE_FILE}" "${OUT_FILE}")

  SET(${OUTPUT_FILE} ${OUT_FILE} PARENT_SCOPE)
ENDFUNCTION(TPETRANEW_PROCESS_ONE_S_TEMPLATE)

# Function to generate .cpp files for ETI of a Tpetra class, over all
# enabled Scalar template arguments.  We generate one .cpp file for
# each Scalar template argument for which Tpetra does ETI.
#
# OUTPUT_FILES [out] List of the generated .cpp files.
#
# TEMPLATE_FILE [in] Name of the input .tmpl "template" file.  This
#   function does string substitution in that file, using the input
#   arguments of this function.  For example, @SC_MACRO_EXPR@ (Scalar
#   macro expression) gets substituted for the value of this
#   function's SC_MACRO_EXPR input argument.
#
# CLASS_NAME [in] Name of the Tpetra class (without namespace
#   qualifiers; must live in the Tpetra namespace)
#
# CLASS_MACRO_NAME [in] Name of the Tpetra class, suitably mangled for
#   use in a macro name.
#
# SCALAR_TYPES [in] All Scalar (or Packet, in the case of DistObject)
#   types over which to do ETI for the given class.  This may include
#   integer Scalar types, if appropriate for that class.
#
# GLOBAL_ORDINAL_TYPE [in] The global ordinal type that Tpetra will
#   use for this class.
#
# MUST_HAVE_SCALAR_INT [in] (Boolean) Whether the class must be
#   instantiated with Scalar = int, even if int is not the default
#   global_ordinal_type.
#
FUNCTION(TPETRANEW_PROCESS_ALL_S_TEMPLATES OUTPUT_FILES TEMPLATE_FILE CLASS_NAME CLASS_MACRO_NAME SCALAR_TYPES GLOBAL_ORDINAL_TYPE MUST_HAVE_SCALAR_INT)
  SET(OUT_FILES "")

  FOREACH(SC ${SCALAR_TYPES})
    TPETRANEW_MANGLE_TEMPLATE_PARAMETER(SC_MANGLED "${SC}")
    TPETRA_SCALAR_MACRO_NAME(SC_MACRO_NAME "${SC}")

    # We don't need a macro expression to enable instantiation,
    # if SC is the global ordinal type, int, or char (special
    # case for DistObject).
    STRING(COMPARE EQUAL "${SC}" "${GLOBAL_ORDINAL_TYPE}" IS_GO)
    STRING(COMPARE EQUAL "${SC}" "int" IS_INT)
    STRING(COMPARE EQUAL "${SC}" "char" IS_CHAR)    

    # We don't need a macro expression if SC is either int or the
    # global ordinal type, since we always instantiate for SC = those
    # types.
    IF(IS_GO OR IS_INT OR IS_CHAR)
      SET(SC_MACRO_EXPR "")
    ELSE()
      SET(SC_MACRO_EXPR "&& defined(HAVE_TPETRA_INST_${SC_MACRO_NAME})")
    ENDIF()

    TPETRANEW_PROCESS_ONE_S_TEMPLATE(OUT_FILE "${TEMPLATE_FILE}" "${CLASS_NAME}" "${CLASS_MACRO_NAME}" "${SC_MANGLED}" "${SC_MACRO_EXPR}")
    LIST(APPEND OUT_FILES ${OUT_FILE})
  ENDFOREACH() # SC

  # This is the standard CMake idiom for setting an output variable so
  # that the caller can see the result.
  SET(${OUTPUT_FILES} ${OUT_FILES} PARENT_SCOPE)
ENDFUNCTION(TPETRANEW_PROCESS_ALL_S_TEMPLATES)

#
# A) Package-specific configuration options
#

# TpetraNew just uses TpetraCore_config.h for macros.
#
#TRIBITS_CONFIGURE_FILE(${PACKAGE_NAME}_config.h)

#
# B) Define the header and source files (and directories)
#

SET(HEADERS "")
SET(SOURCES "")

SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR})
APPEND_GLOB(HEADERS ${DIR}/*.h)
APPEND_GLOB(HEADERS ${DIR}/*.hpp)
APPEND_GLOB(SOURCES ${DIR}/*.cpp)
TRILINOS_CREATE_CLIENT_TEMPLATE_HEADERS(${DIR})

# Must glob the binary dir last to get all of the auto-generated headers
SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_BINARY_DIR})
APPEND_GLOB(HEADERS ${DIR}/*.hpp)
# APPEND_SET(HEADERS ${DIR}/${PACKAGE_NAME}_config.h )
# APPEND_SET(HEADERS ${DIR}/${PACKAGE_NAME}_ETIHelperMacros.h)

#
# TpetraNew just inherits this property from TpetraCore.
#
ASSERT_DEFINED (TpetraCore_ENABLE_EXPLICIT_INSTANTIATION)
IF (TpetraCore_ENABLE_EXPLICIT_INSTANTIATION)

  # Generate source (.cpp) files for classes that need explicit
  # template instantiation over the Scalar (or Packet, in the case of
  # DistObject) template parameter.

  SET (MultiVector_ETI_SCALARS ${TpetraCore_ETI_SCALARS})
  IF (NOT Tpetra_INST_INT_INT)
    # GO = int is disabled, so the Scalar = GlobalOrdinal case doesn't
    # cover Scalar = int.  We need Scalar = int for communication of
    # process ranks.
    LIST(APPEND MultiVector_ETI_SCALARS "int")
  ENDIF()

  # Generate ETI .cpp files for TpetraNew::MultiVector.
  #
  # mfh 07 Sep 2018: TpetraNew doesn't have MultiVector yet.
  # Once it does, uncomment the two commands below.
  #
  # TPETRANEW_PROCESS_ALL_S_TEMPLATES(NEW_MULTIVECTOR_OUTPUT_FILES "TpetraNew_ETI_SC.tmpl" "MultiVector" "MULTIVECTOR" "${MultiVector_ETI_SCALARS}" "${Tpetra_DEFAULT_GLOBAL_ORDINAL}" TRUE)
  # LIST(APPEND SOURCES ${NEW_MULTIVECTOR_OUTPUT_FILES})

  # DistObject requires the same Scalar type instantiations as
  # MultiVector, plus Scalar = char (for CrsMatrix).
  SET (DistObject_ETI_SCALARS ${MultiVector_ETI_SCALARS})
  LIST(APPEND DistObject_ETI_SCALARS "char")

  # Generate ETI .cpp files for TpetraNew::DistObject.
  TPETRANEW_PROCESS_ALL_S_TEMPLATES(NEW_DISTOBJECT_OUTPUT_FILES "TpetraNew_ETI_SC.tmpl" "DistObject" "DISTOBJECT" "${DistObject_ETI_SCALARS}" "${Tpetra_DEFAULT_GLOBAL_ORDINAL}" TRUE)
  LIST(APPEND SOURCES ${NEW_DISTOBJECT_OUTPUT_FILES})

  #MESSAGE(STATUS "SOURCES = ${SOURCES}")
ENDIF()

#
# C) Define the targets for package's library(s)
#

TRIBITS_ADD_LIBRARY(
  tpetranew
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  ADDED_LIB_TARGET_NAME_OUT TPETRANEW_LIBNAME
  )

# We need to set the linker language explicitly here for CUDA builds.
SET_PROPERTY(
  TARGET ${TPETRANEW_LIBNAME}
  APPEND PROPERTY LINKER_LANGUAGE CXX
  )

#
# Make a trivial change to this comment if you add / remove a file
# either to / from this directory.  That ensures that running "make"
# will also rerun CMake in order to regenerate Makefiles.
#
