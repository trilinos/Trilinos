MESSAGE(STATUS "${PACKAGE_NAME}: Processing ETI / test support")

# This CMake module generates the following header file, which gets
# written to the build directory (like other header files that CMake
# generates).  The file contains macros that do instantiation over a
# finite set of template parameter combinations.  We use the macros
# both for ETI (explicit template instantiation), and for tests.
# Thus, the macros need to be generated even if ETI is OFF.
SET(${PACKAGE_NAME}_ETI_FILE TpetraCore_ETIHelperMacros.h)
SET(${PACKAGE_NAME}_ETI_FILE_PATH ${Tpetra_BINARY_DIR}/core/src/${${PACKAGE_NAME}_ETI_FILE})

#
# Users have the option to generate the above header file themselves.
# We prefer that users let Trilinos generate the header file.
# However, folks who make intense use of TriBITS sometimes find that
# reusing a previously generated header file shaves a couple minutes
# off their CMake configuration time.  Thus, we give them that option.
#

ADVANCED_SET(Tpetra_USE_STATIC_ETI_MACROS_HEADER_FILE ""
  CACHE PATH
  "If set, gives the path to a static version of the file ${${PACKAGE_NAME}_ETI_FILE}.  If not set (the default), then the file is generated automatically."
  )

IF(Tpetra_USE_STATIC_ETI_MACROS_HEADER_FILE)
  # The user wants us to accept their header file and not generate one.
  MESSAGE("-- NOTE: Skipping generation and using provided static file"
     " '${Tpetra_USE_STATIC_ETI_MACROS_HEADER_FILE}'")
  CONFIGURE_FILE(
    ${Tpetra_USE_STATIC_ETI_MACROS_HEADER_FILE}
    ${${PACKAGE_NAME}_ETI_FILE_PATH}
    COPYONY
    )
  RETURN()
ENDIF()

#
# The user wants us to generate the header file.  This is the default
# behavior.
#

# Tpetra ETI type fields.  S, LO, GO, N correspond to the four
# template parameters of most Tpetra classes: Scalar, LocalOrdinal,
# GlobalOrdinal, and Node.  SIN and SOUT are for classes and functions
# that convert between two different Scalar types.  CS and DS serve a
# similar function.
SET(${PACKAGE_NAME}_ETI_FIELDS "SIN|SOUT|S|LO|GO|N|CS|DS")

# mfh 11 Feb 2015: If a given Node type doesn't support certain
# combinations of Scalar, LO, GO, etc. types, this would be the place
# to add in only those type combinations that the Node type supports.
# We used to do this for Node = KokkosClassic::ThrustGPUNode.  I much
# prefer instead that all type combinations _build_, but that
# unsupported types for a particular Node get stub implementations
# that just throw.  This simplifies the CMake configuration process.

# Set up a pattern that excludes all Scalar types that are also
# possible GlobalOrdinal types, but includes all other types.
# TriBITS' ETI system knows how to interpret this pattern.
#
# FIXME (mfh 17 Aug 2015) A better way to do this would be to subtract
# away all enabled GlobalOrdinal types.  Plus, what if someone really
# wants a CrsMatrix<int,...>?
TRIBITS_ETI_TYPE_EXPANSION(${PACKAGE_NAME}_ETI_EXCLUDE_SET_ORDINAL_SCALAR "S=int|unsigned|unsigned int|long|unsigned long|long long|unsigned long long" "LO=.*" "GO=.*" "N=.*")


# Excludes for MultiVector
TRIBITS_ETI_TYPE_EXPANSION(${PACKAGE_NAME}_ETI_EXCLUDE_SET_GORD_NOT_INT "S=.*" "LO=.*" "GO=int" "N=.*")


# TriBITS' ETI system expects a set of types to be a string, delimited
# by |.  Each template parameter (e.g., Scalar, LocalOrdinal, ...) has
# its own set.  The JOIN commands below set up those lists.  We use
# the following sets that Tpetra defines:
#
# Scalar: ${PACKAGE_NAME}_ETI_SCALARS
# LocalOrdinal: ${PACKAGE_NAME}_ETI_LORDS
# GlobalOrdinal: ${PACKAGE_NAME}_ETI_GORDS
# Node: ${PACKAGE_NAME}_ETI_NODES
#
# Note that the Scalar set from Tpetra includes the Scalar =
# GlobalOrdinal case.  We have to exclude that explicitly in what
# follows.
JOIN(${PACKAGE_NAME}_ETI_LORDS   "|" FALSE ${${PACKAGE_NAME}_ETI_LORDS}  )
JOIN(${PACKAGE_NAME}_ETI_GORDS   "|" FALSE ${${PACKAGE_NAME}_ETI_GORDS}  )
JOIN(${PACKAGE_NAME}_ETI_NODES   "|" FALSE ${${PACKAGE_NAME}_ETI_NODES}  )
JOIN(${PACKAGE_NAME}_ETI_SCALARS "|" FALSE ${${PACKAGE_NAME}_ETI_SCALARS})

MESSAGE(STATUS "Enabled Scalar types:        ${${PACKAGE_NAME}_ETI_SCALARS}")
MESSAGE(STATUS "Enabled LocalOrdinal types:  ${${PACKAGE_NAME}_ETI_LORDS}")
MESSAGE(STATUS "Enabled GlobalOrdinal types: ${${PACKAGE_NAME}_ETI_GORDS}")
MESSAGE(STATUS "Enabled Node types:          ${${PACKAGE_NAME}_ETI_NODES}")

# Construct the "type expansion" string that TriBITS' ETI system
# expects.  Even if ETI is OFF, we will use this to generate macros
# for instantiating tests.
TRIBITS_ETI_TYPE_EXPANSION(SingleScalarInsts
  "S=${${PACKAGE_NAME}_ETI_SCALARS}"
  "N=${${PACKAGE_NAME}_ETI_NODES}"
  "LO=${${PACKAGE_NAME}_ETI_LORDS}"
  "GO=${${PACKAGE_NAME}_ETI_GORDS}")

# Set up the set of enabled type combinations, in a format that
# TriBITS understands.
#
# mfh 17 Aug 2015: I don't exactly understand what's going on here,
# but it looks like if ETI is enabled, we let users modify
# TpetraCore_ETI_LIBRARYSET, and if it's not, we don't.
ASSERT_DEFINED(TpetraCore_ENABLE_EXPLICIT_INSTANTIATION)
IF(TpetraCore_ENABLE_EXPLICIT_INSTANTIATION)
  TRIBITS_ADD_ETI_INSTANTIATIONS(TpetraCore ${SingleScalarInsts})
  # mfh 17 Aug 2015: For some reason, these are empty unless ETI is
  # ON.  That seems to be OK, though, unless users want to exclude
  # enabling certain type combinations when ETI is OFF.
  MESSAGE(STATUS "Excluded type combinations: ${${PACKAGE_NAME}_ETI_EXCLUDE_SET}:${${PACKAGE_NAME}_ETI_EXCLUDE_SET_INT}")
ELSE()
  TRIBITS_ETI_TYPE_EXPANSION(${PACKAGE_NAME}_ETI_LIBRARYSET
    "S=${${PACKAGE_NAME}_ETI_SCALARS}"
    "N=${${PACKAGE_NAME}_ETI_NODES}"
    "LO=${${PACKAGE_NAME}_ETI_LORDS}"
    "GO=${${PACKAGE_NAME}_ETI_GORDS}")
ENDIF()
MESSAGE(STATUS "Set of enabled types, before exclusions: ${${PACKAGE_NAME}_ETI_LIBRARYSET}")

#
# Generate the instantiation macros.  These go into
# TpetraCore_ETIHelperMacros.h, which is generated from
# TpetraCore_ETIHelperMacros.h.in (in this directory).
#
TRIBITS_ETI_GENERATE_MACROS(
    "${${PACKAGE_NAME}_ETI_FIELDS}" "${${PACKAGE_NAME}_ETI_LIBRARYSET}"
    "${${PACKAGE_NAME}_ETI_EXCLUDE_SET};${${PACKAGE_NAME}_ETI_EXCLUDE_SET_ORDINAL_SCALAR}"
    list_of_manglings eti_typedefs
    "TPETRA_INSTANTIATE_SLGN_NO_ORDINAL_SCALAR(S,LO,GO,N)" TPETRA_ETIMACRO_SLGN_NO_ORDINAL_SCALAR
    "TPETRA_INSTANTIATE_SLG_NO_ORDINAL_SCALAR(S,LO,GO)"    TPETRA_ETIMACRO_SLG_NO_ORDINAL_SCALAR
    "TPETRA_INSTANTIATE_SL_NO_ORDINAL_SCALAR(S,LO)"        TPETRA_ETIMACRO_SL_NO_ORDINAL_SCALAR
    "TPETRA_INSTANTIATE_SN_NO_ORDINAL_SCALAR(S,N)"        TPETRA_ETIMACRO_SN_NO_ORDINAL_SCALAR
    "TPETRA_INSTANTIATE_S_NO_ORDINAL_SCALAR(S)"            TPETRA_ETIMACRO_S_NO_ORDINAL_SCALAR)
TRIBITS_ETI_GENERATE_MACROS(
    "${${PACKAGE_NAME}_ETI_FIELDS}" "${${PACKAGE_NAME}_ETI_LIBRARYSET}"
    "${${PACKAGE_NAME}_ETI_EXCLUDE_SET};${${PACKAGE_NAME}_ETI_EXCLUDE_SET_INT}"
    list_of_manglings eti_typedefs
    "TPETRA_INSTANTIATE_GLGN(GO,LO,GO,N)"           TPETRA_ETIMACRO_GLGN
    "TPETRA_INSTANTIATE_LGN(LO,GO,N)"               TPETRA_ETIMACRO_LGN
    "TPETRA_INSTANTIATE_GLG(GO,LO,GO)"              TPETRA_ETIMACRO_GLG
    "TPETRA_INSTANTIATE_LG(LO,GO)"                  TPETRA_ETIMACRO_LG
    "TPETRA_INSTANTIATE_SL(S,LO)"                   TPETRA_ETIMACRO_SL
    "TPETRA_INSTANTIATE_SN(S,N)"                    TPETRA_ETIMACRO_SN
    "TPETRA_INSTANTIATE_GN(GO,N)"                   TPETRA_ETIMACRO_GN
    "TPETRA_INSTANTIATE_S(S)"                       TPETRA_ETIMACRO_S
    "TPETRA_INSTANTIATE_L(LO)"                      TPETRA_ETIMACRO_L
    "TPETRA_INSTANTIATE_G(GO)"                      TPETRA_ETIMACRO_G
    "TPETRA_INSTANTIATE_N(N)"                       TPETRA_ETIMACRO_N
    "TPETRA_INSTANTIATE_TSLGN(CS,DS,LO,GO,N)"       TPETRA_ETIMACRO_TSLGN
    "TPETRA_INSTANTIATE_TSLG(CS,DS,LO,GO)"          TPETRA_ETIMACRO_TSLG
    "TPETRA_INSTANTIATE_CONVERT(SOUT,SIN,LO,GO,N)"  TPETRA_ETIMACRO_CONVERT
    "TPETRA_INSTANTIATE_CONVERT_SSL(SOUT,SIN,LO)"    TPETRA_ETIMACRO_CONVERT_SSL)

#SET(${PACKAGE_NAME}_SCALAR_INT_ETI_FIELDS "LO|GO|N")
#TRIBITS_ETI_TYPE_EXPANSION(${PACKAGE_NAME}_SCALAR_INT_ETI_LIBRARYSET
#  "LO=${${PACKAGE_NAME}_ETI_LORDS}"
#  "GO=${${PACKAGE_NAME}_ETI_GORDS}"
#  "N=${${PACKAGE_NAME}_ETI_NODES}")


# For Mutlivector and friends
TRIBITS_ETI_GENERATE_MACROS(
    "${${PACKAGE_NAME}_ETI_FIELDS}" "${${PACKAGE_NAME}_ETI_LIBRARYSET}"
    "${${PACKAGE_NAME}_ETI_EXCLUDE_SET_GORD_NOT_INT}"
    list_of_manglings eti_typedefs
    "TPETRA_INSTANTIATE_GLGN_NO_INT(GO,LO,GO,N)"           TPETRA_ETIMACRO_GLGN_NO_INT)

TRIBITS_ETI_GENERATE_MACROS(
    "${${PACKAGE_NAME}_ETI_FIELDS}" "${${PACKAGE_NAME}_ETI_LIBRARYSET}"
    "${${PACKAGE_NAME}_ETI_EXCLUDE_SET}"
    list_of_manglings eti_typedefs
    "TPETRA_INSTANTIATE_LLGN(LO,LO,GO,N)"           TPETRA_ETIMACRO_LLGN)


# Tpetra::DistObject is templated on "Packet", the type of thing to be
# sent and received.  Tpetra wants to send both Scalar and
# GlobalOrdinal objects directly.  In some cases, it also serializes
# into char (byte).  The TPETRA_DISTOBJECT_INSTANT macro takes care of
# that for us; we don't have to add P=char here.
STRING(REPLACE "S="  "P=" ScalarToPacketSet "${${PACKAGE_NAME}_ETI_LIBRARYSET}")
STRING(REGEX REPLACE "GO=([^{ ]+)" "GO=\\1 P=\\1" GlobalToPacketSet1 "${${PACKAGE_NAME}_ETI_LIBRARYSET}")
STRING(REGEX REPLACE "GO={([^}]+)}" "GO={\\1} P={\\1}" GlobalToPacketSet2 "${${PACKAGE_NAME}_ETI_LIBRARYSET}")
TRIBITS_ETI_GENERATE_MACROS(
    "${${PACKAGE_NAME}_ETI_FIELDS}|P" "${ScalarToPacketSet};${GlobalToPacketSet1};${GlobalToPacketSet2}"
    "${${PACKAGE_NAME}_ETI_EXCLUDE_SET}"
    list_of_manglings eti_typedefs
    "TPETRA_INSTANTIATE_PLGN(P,LO,GO,N)"            TPETRA_ETIMACRO_PLGN)

TRIBITS_ETI_TYPE_EXPANSION(Tpetra_DII   "S=double" "N=${${PACKAGE_NAME}_ETI_NODES}" "LO=int" "GO=int")
TRIBITS_ETI_GENERATE_MACROS(
    "${${PACKAGE_NAME}_ETI_FIELDS}" "${Tpetra_DII}"
    "${${PACKAGE_NAME}_ETI_EXCLUDE_SET}"
    list_of_manglings eti_typedefs
    "TPETRA_INSTANTIATE_DOUBLE_INT_INT_N(S,LO,GO,N)"      TPETRA_ETIMACRO_DII_NODE)

# Generate "mangled" typedefs.  Macros sometimes get grumpy when types
# have spaces, colons, or angle brackets in them.  This includes types
# like "long long" or "std::complex<double>".  Thus, we define
# typedefs that remove the offending characters.  The typedefs also
# get written to the generated header file.
TRIBITS_ETI_GENERATE_TYPEDEF_MACRO(TPETRA_ETI_TYPEDEFS "TPETRA_ETI_MANGLING_TYPEDEFS" "${eti_typedefs}")

# Generate the header file TpetraCore_ETIHelperMacros.h, from the file
# TpetraCore_ETIHelperMacros.h.in (that lives in this directory).  The
# generated header file gets written to the Trilinos build directory,
# in packages/tpetra/core/src/.
CONFIGURE_FILE(
  ${${PACKAGE_NAME}_SOURCE_DIR}/cmake/TpetraCore_ETIHelperMacros.h.in
  ${${PACKAGE_NAME}_ETI_FILE_PATH}
  )
add_custom_target(Tpetra_ETI_generated DEPENDS ${${PACKAGE_NAME}_ETI_FILE_PATH})
