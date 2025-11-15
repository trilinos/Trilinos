include(Join)
MESSAGE(STATUS "${PACKAGE_NAME}: Processing ETI / test support")

# TrilinosCouplings ETI type fields.  S, LO, GO, N correspond to the four
# template parameters of most Tpetra classes: Scalar, LocalOrdinal,
# GlobalOrdinal, and Node.  TrilinosCouplings shares these with Tpetra, because
# TrilinosCouplings only works with Tpetra linear algebra objects.
SET(TrilinosCouplings_ETI_FIELDS "S|LO|GO|N")

# TriBITS' ETI system expects a set of types to be a string, delimited
# by |.  Each template parameter (e.g., Scalar, LocalOrdinal, ...) has
# its own set.  The JOIN commands below set up those lists.  We use
# the following sets that TrilinosCouplings defines:
#
# Scalar: TrilinosCouplings_ETI_SCALARS
# LocalOrdinal: TrilinosCouplings_ETI_LORDS
# GlobalOrdinal: TrilinosCouplings_ETI_GORDS
# Node: TrilinosCouplings_ETI_NODES
#
# Note that the Scalar set from Tpetra includes the Scalar =
# GlobalOrdinal case.  However, TrilinosCouplings's CMake logic excludes this,
# so we don't have to worry about it here.

JOIN(TrilinosCouplings_ETI_SCALARS "|" FALSE ${TrilinosCouplings_ETI_SCALARS})
JOIN(TrilinosCouplings_ETI_LORDS   "|" FALSE ${TrilinosCouplings_ETI_LORDS}  )
JOIN(TrilinosCouplings_ETI_GORDS   "|" FALSE ${TrilinosCouplings_ETI_GORDS}  )
JOIN(TrilinosCouplings_ETI_NODES   "|" FALSE ${TrilinosCouplings_ETI_NODES}  )

MESSAGE(STATUS "Enabled Scalar types:        ${TrilinosCouplings_ETI_SCALARS}")
MESSAGE(STATUS "Enabled LocalOrdinal types:  ${TrilinosCouplings_ETI_LORDS}")
MESSAGE(STATUS "Enabled GlobalOrdinal types: ${TrilinosCouplings_ETI_GORDS}")
MESSAGE(STATUS "Enabled Node types:          ${TrilinosCouplings_ETI_NODES}")

# Construct the "type expansion" string that TriBITS' ETI system
# expects.  Even if ETI is OFF, we will use this to generate macros
# for instantiating tests.
TRIBITS_ETI_TYPE_EXPANSION(SingleScalarInsts
  "S=${TrilinosCouplings_ETI_SCALARS}"
  "N=${TrilinosCouplings_ETI_NODES}"
  "LO=${TrilinosCouplings_ETI_LORDS}"
  "GO=${TrilinosCouplings_ETI_GORDS}")

ASSERT_DEFINED(TrilinosCouplings_ENABLE_EXPLICIT_INSTANTIATION)
IF(TrilinosCouplings_ENABLE_EXPLICIT_INSTANTIATION)
  # mfh 17 Aug 2015: If ETI is ON, it looks like users can set these variables.
  MESSAGE(STATUS "User/Downstream ETI set: ${TrilinosCouplings_ETI_LIBRARYSET}")
  TRIBITS_ADD_ETI_INSTANTIATIONS(TrilinosCouplings ${SingleScalarInsts})
  MESSAGE(STATUS "Excluded type combinations: ${TrilinosCouplings_ETI_EXCLUDE_SET}")
ELSE()
  TRIBITS_ETI_TYPE_EXPANSION(TrilinosCouplings_ETI_LIBRARYSET
    "S=${TrilinosCouplings_ETI_SCALARS}"
    "N=${TrilinosCouplings_ETI_NODES}"
    "LO=${TrilinosCouplings_ETI_LORDS}"
    "GO=${TrilinosCouplings_ETI_GORDS}")
ENDIF()
MESSAGE(STATUS "Set of enabled types, before exclusions: ${${PACKAGE_NAME}_ETI_LIBRARYSET}")

#
# Generate the instantiation macros.  These go into
# TrilinosCouplings_ETIHelperMacros.h, which is generated from
# TrilinosCouplings_ETIHelperMacros.h.in (in this directory).
#
TRIBITS_ETI_GENERATE_MACROS("${TrilinosCouplings_ETI_FIELDS}" "${TrilinosCouplings_ETI_LIBRARYSET}" "${TrilinosCouplings_ETI_EXCLUDE_SET}"
                            list_of_manglings eti_typedefs
                            "TRILINOSCOUPLINGS_INSTANTIATE_SLGN(S,LO,GO,N)"     TRILINOSCOUPLINGS_ETIMACRO_SLGN
                            )

# Generate "mangled" typedefs.  Macros sometimes get grumpy when types
# have spaces, colons, or angle brackets in them.  This includes types
# like "long long" or "std::complex<double>".  Thus, we define
# typedefs that remove the offending characters.  The typedefs also
# get written to the generated header file.
TRIBITS_ETI_GENERATE_TYPEDEF_MACRO(TRILINOSCOUPLINGS_ETI_TYPEDEFS "TRILINOSCOUPLINGS_ETI_MANGLING_TYPEDEFS" "${eti_typedefs}")

# Generate the header file TrilinosCouplings_ETIHelperMacros.h, from the file
# TrilinosCouplings_ETIHelperMacros.h.in (that lives in this directory).  The
# generated header file gets written to the Trilinos build directory,
# in packages/TrilinosCouplings/src/.
CONFIGURE_FILE(${TrilinosCouplings_SOURCE_DIR}/cmake/TrilinosCouplings_ETIHelperMacros.h.in ${TrilinosCouplings_BINARY_DIR}/src/TrilinosCouplings_ETIHelperMacros.h)
