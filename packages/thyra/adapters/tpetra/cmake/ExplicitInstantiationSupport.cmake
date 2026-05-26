include(Join)
MESSAGE(STATUS "${PACKAGE_NAME}: Processing ETI / test support")

# Ifpack2 ETI type fields.  S, LO, GO, N correspond to the four
# template parameters of most Tpetra classes: Scalar, LocalOrdinal,
# GlobalOrdinal, and Node.  ThyraTpetraAdapters shares these with Tpetra, because
# ThyraTpetraAdapters only works with Tpetra linear algebra objects.
SET(ThyraTpetraAdapters_ETI_FIELDS "S|LO|GO|N")

# Set up a pattern that excludes all complex Scalar types.
# TriBITS' ETI system knows how to interpret this pattern.
TRIBITS_ETI_TYPE_EXPANSION(${PACKAGE_NAME}_ETI_EXCLUDE_SET_COMPLEX "S=std::complex<float>|std::complex<double>" "LO=.*" "GO=.*" "N=.*")

# TriBITS' ETI system expects a set of types to be a string, delimited
# by |.  Each template parameter (e.g., Scalar, LocalOrdinal, ...) has
# its own set.  The JOIN commands below set up those lists.  We use
# the following sets that ThyraTpetraAdapters defines:
#
# Scalar: ThyraTpetraAdapters_ETI_SCALARS
# LocalOrdinal: ThyraTpetraAdapters_ETI_LORDS
# GlobalOrdinal: ThyraTpetraAdapters_ETI_GORDS
# Node: ThyraTpetraAdapters_ETI_NODES
#
# Note that the Scalar set from Tpetra includes the Scalar =
# GlobalOrdinal case.  However, ThyraTpetraAdapters's CMake logic excludes this,
# so we don't have to worry about it here.

JOIN(ThyraTpetraAdapters_ETI_SCALARS "|" FALSE ${ThyraTpetraAdapters_ETI_SCALARS})
JOIN(ThyraTpetraAdapters_ETI_LORDS   "|" FALSE ${ThyraTpetraAdapters_ETI_LORDS}  )
JOIN(ThyraTpetraAdapters_ETI_GORDS   "|" FALSE ${ThyraTpetraAdapters_ETI_GORDS}  )
JOIN(ThyraTpetraAdapters_ETI_NODES   "|" FALSE ${ThyraTpetraAdapters_ETI_NODES}  )

MESSAGE(STATUS "Enabled Scalar types:        ${ThyraTpetraAdapters_ETI_SCALARS}")
MESSAGE(STATUS "Enabled LocalOrdinal types:  ${ThyraTpetraAdapters_ETI_LORDS}")
MESSAGE(STATUS "Enabled GlobalOrdinal types: ${ThyraTpetraAdapters_ETI_GORDS}")
MESSAGE(STATUS "Enabled Node types:          ${ThyraTpetraAdapters_ETI_NODES}")

# Construct the "type expansion" string that TriBITS' ETI system
# expects.  Even if ETI is OFF, we will use this to generate macros
# for instantiating tests.
TRIBITS_ETI_TYPE_EXPANSION(SingleScalarInsts
  "S=${ThyraTpetraAdapters_ETI_SCALARS}"
  "N=${ThyraTpetraAdapters_ETI_NODES}"
  "LO=${ThyraTpetraAdapters_ETI_LORDS}"
  "GO=${ThyraTpetraAdapters_ETI_GORDS}")

ASSERT_DEFINED(ThyraTpetraAdapters_ENABLE_EXPLICIT_INSTANTIATION)
IF(ThyraTpetraAdapters_ENABLE_EXPLICIT_INSTANTIATION)
  # mfh 17 Aug 2015: If ETI is ON, it looks like users can set these variables.
  MESSAGE(STATUS "User/Downstream ETI set: ${ThyraTpetraAdapters_ETI_LIBRARYSET}")
  TRIBITS_ADD_ETI_INSTANTIATIONS(ThyraTpetraAdapters ${SingleScalarInsts})
  MESSAGE(STATUS "Excluded type combinations: ${ThyraTpetraAdapters_ETI_EXCLUDE_SET}")
ELSE()
  TRIBITS_ETI_TYPE_EXPANSION(ThyraTpetraAdapters_ETI_LIBRARYSET
    "S=${ThyraTpetraAdapters_ETI_SCALARS}"
    "N=${ThyraTpetraAdapters_ETI_NODES}"
    "LO=${ThyraTpetraAdapters_ETI_LORDS}"
    "GO=${ThyraTpetraAdapters_ETI_GORDS}")
ENDIF()
MESSAGE(STATUS "Set of enabled types, before exclusions: ${${PACKAGE_NAME}_ETI_LIBRARYSET}")

#
# Generate the instantiation macros.  These go into
# ThyraTpetraAdapters_ETIHelperMacros.h, which is generated from
# ThyraTpetraAdapters_ETIHelperMacros.h.in (in this directory).
#
TRIBITS_ETI_GENERATE_MACROS("${ThyraTpetraAdapters_ETI_FIELDS}" "${ThyraTpetraAdapters_ETI_LIBRARYSET}" "${ThyraTpetraAdapters_ETI_EXCLUDE_SET}"
                            list_of_manglings eti_typedefs
                            "THYRATPETRAADAPTERS_INSTANTIATE_L(LO)"         THYRATPETRAADAPTERS_ETIMACRO_L
                            "THYRATPETRAADAPTERS_INSTANTIATE_SL(S,LO)"      THYRATPETRAADAPTERS_ETIMACRO_SL
                            "THYRATPETRAADAPTERS_INSTANTIATE_LG(LO,GO)"         THYRATPETRAADAPTERS_ETIMACRO_LG
                            "THYRATPETRAADAPTERS_INSTANTIATE_SLG(S,LO,GO)"      THYRATPETRAADAPTERS_ETIMACRO_SLG
                            "THYRATPETRAADAPTERS_INSTANTIATE_LGN(S,LO,GO,N)"      THYRATPETRAADAPTERS_ETIMACRO_LGN
                            "THYRATPETRAADAPTERS_INSTANTIATE_SLGN(S,LO,GO,N)"     THYRATPETRAADAPTERS_ETIMACRO_SLGN
                            "THYRATPETRAADAPTERS_INSTANTIATE_N(N)"     THYRATPETRAADAPTERS_ETIMACRO_N
                            )
TRIBITS_ETI_GENERATE_MACROS("${ThyraTpetraAdapters_ETI_FIELDS}" "${ThyraTpetraAdapters_ETI_LIBRARYSET}"
                            "${ThyraTpetraAdapters_ETI_EXCLUDE_SET};${ThyraTpetraAdapters_ETI_EXCLUDE_SET_COMPLEX}"
                            list_of_manglings eti_typedefs
                            "THYRATPETRAADAPTERS_INSTANTIATE_SL_REAL(S,LO,GO)" THYRATPETRAADAPTERS_ETIMACRO_SL_REAL
                            "THYRATPETRAADAPTERS_INSTANTIATE_SLG_REAL(S,LO,GO)" THYRATPETRAADAPTERS_ETIMACRO_SLG_REAL
                            "THYRATPETRAADAPTERS_INSTANTIATE_SLGN_REAL(S,LO,GO,N)" THYRATPETRAADAPTERS_ETIMACRO_SLGN_REAL
                            )

# Generate "mangled" typedefs.  Macros sometimes get grumpy when types
# have spaces, colons, or angle brackets in them.  This includes types
# like "long long" or "std::complex<double>".  Thus, we define
# typedefs that remove the offending characters.  The typedefs also
# get written to the generated header file.
TRIBITS_ETI_GENERATE_TYPEDEF_MACRO(THYRATPETRAADAPTERS_ETI_TYPEDEFS "THYRATPETRAADAPTERS_ETI_MANGLING_TYPEDEFS" "${eti_typedefs}")

# Generate the header file ThyraTpetraAdapters_ETIHelperMacros.h, from the file
# ThyraTpetraAdapters_ETIHelperMacros.h.in (that lives in this directory).  The
# generated header file gets written to the Trilinos build directory,
# in packages/ThyraTpetraAdapters/src/.
CONFIGURE_FILE(${ThyraTpetraAdapters_SOURCE_DIR}/cmake/ThyraTpetraAdapters_ETIHelperMacros.h.in ${ThyraTpetraAdapters_BINARY_DIR}/src/ThyraTpetraAdapters_ETIHelperMacros.h)
