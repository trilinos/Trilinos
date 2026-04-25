include(Join)
MESSAGE(STATUS "${PACKAGE_NAME}: Processing ETI / test support")

# Galeri ETI type fields.  S, LO, GO, N correspond to the four
# template parameters of most Tpetra classes: Scalar, LocalOrdinal,
# GlobalOrdinal, and Node.  Galeri shares these with Tpetra, because
# Galeri only works with Tpetra linear algebra objects.
SET(Galeri_ETI_FIELDS "S|LO|GO|N")

# Set up a pattern that excludes all complex Scalar types.
# TriBITS' ETI system knows how to interpret this pattern.
TRIBITS_ETI_TYPE_EXPANSION(${PACKAGE_NAME}_ETI_EXCLUDE_SET_COMPLEX "S=std::complex<float>|std::complex<double>" "LO=.*" "GO=.*" "N=.*")

# TriBITS' ETI system expects a set of types to be a string, delimited
# by |.  Each template parameter (e.g., Scalar, LocalOrdinal, ...) has
# its own set.  The JOIN commands below set up those lists.  We use
# the following sets that Galeri defines:
#
# Scalar: Galeri_ETI_SCALARS
# LocalOrdinal: Galeri_ETI_LORDS
# GlobalOrdinal: Galeri_ETI_GORDS
# Node: Galeri_ETI_NODES
#
# Note that the Scalar set from Tpetra includes the Scalar =
# GlobalOrdinal case.  However, Galeri's CMake logic excludes this,
# so we don't have to worry about it here.

JOIN(Galeri_ETI_SCALARS "|" FALSE ${Galeri_ETI_SCALARS})
JOIN(Galeri_ETI_LORDS   "|" FALSE ${Galeri_ETI_LORDS}  )
JOIN(Galeri_ETI_GORDS   "|" FALSE ${Galeri_ETI_GORDS}  )
JOIN(Galeri_ETI_NODES   "|" FALSE ${Galeri_ETI_NODES}  )

MESSAGE(STATUS "Enabled Scalar types:        ${Galeri_ETI_SCALARS}")
MESSAGE(STATUS "Enabled LocalOrdinal types:  ${Galeri_ETI_LORDS}")
MESSAGE(STATUS "Enabled GlobalOrdinal types: ${Galeri_ETI_GORDS}")
MESSAGE(STATUS "Enabled Node types:          ${Galeri_ETI_NODES}")

# Construct the "type expansion" string that TriBITS' ETI system
# expects.  Even if ETI is OFF, we will use this to generate macros
# for instantiating tests.
TRIBITS_ETI_TYPE_EXPANSION(SingleScalarInsts 
  "S=${Galeri_ETI_SCALARS}"
  "N=${Galeri_ETI_NODES}"
  "LO=${Galeri_ETI_LORDS}"
  "GO=${Galeri_ETI_GORDS}")

ASSERT_DEFINED(Galeri_ENABLE_EXPLICIT_INSTANTIATION)
IF(Galeri_ENABLE_EXPLICIT_INSTANTIATION)
  # mfh 17 Aug 2015: If ETI is ON, it looks like users can set these variables.
  MESSAGE(STATUS "User/Downstream ETI set: ${Galeri_ETI_LIBRARYSET}")
  TRIBITS_ADD_ETI_INSTANTIATIONS(Galeri ${SingleScalarInsts})
  MESSAGE(STATUS "Excluded type combinations: ${Galeri_ETI_EXCLUDE_SET}")
ELSE()
  TRIBITS_ETI_TYPE_EXPANSION(Galeri_ETI_LIBRARYSET
    "S=${Galeri_ETI_SCALARS}"
    "N=${Galeri_ETI_NODES}"
    "LO=${Galeri_ETI_LORDS}"
    "GO=${Galeri_ETI_GORDS}")
ENDIF()
MESSAGE(STATUS "Set of enabled types, before exclusions: ${${PACKAGE_NAME}_ETI_LIBRARYSET}")

#
# Generate the instantiation macros.  These go into
# Galeri_ETIHelperMacros.h, which is generated from
# Galeri_ETIHelperMacros.h.in (in this directory).
#
TRIBITS_ETI_GENERATE_MACROS("${Galeri_ETI_FIELDS}" "${Galeri_ETI_LIBRARYSET}" "${Galeri_ETI_EXCLUDE_SET}"
                            list_of_manglings eti_typedefs
                            "GALERI_INSTANTIATE_L(LO)"         GALERI_ETIMACRO_L
                            "GALERI_INSTANTIATE_SL(S,LO)"      GALERI_ETIMACRO_SL
                            "GALERI_INSTANTIATE_LG(LO,GO)"         GALERI_ETIMACRO_LG
                            "GALERI_INSTANTIATE_SLG(S,LO,GO)"      GALERI_ETIMACRO_SLG
                            "GALERI_INSTANTIATE_LGN(S,LO,GO,N)"      GALERI_ETIMACRO_LGN
                            "GALERI_INSTANTIATE_SLGN(S,LO,GO,N)"     GALERI_ETIMACRO_SLGN
                            "GALERI_INSTANTIATE_G(GO)"     GALERI_ETIMACRO_G
                            "GALERI_INSTANTIATE_N(N)"     GALERI_ETIMACRO_N
                            )
TRIBITS_ETI_GENERATE_MACROS("${Galeri_ETI_FIELDS}" "${Galeri_ETI_LIBRARYSET}"
                            "${Galeri_ETI_EXCLUDE_SET};${Galeri_ETI_EXCLUDE_SET_COMPLEX}"
                            list_of_manglings eti_typedefs
                            "GALERI_INSTANTIATE_SL_REAL(S,LO,GO)" GALERI_ETIMACRO_SL_REAL
                            "GALERI_INSTANTIATE_SLG_REAL(S,LO,GO)" GALERI_ETIMACRO_SLG_REAL
                            "GALERI_INSTANTIATE_SLGN_REAL(S,LO,GO,N)" GALERI_ETIMACRO_SLGN_REAL
                            )

# Generate "mangled" typedefs.  Macros sometimes get grumpy when types
# have spaces, colons, or angle brackets in them.  This includes types
# like "long long" or "std::complex<double>".  Thus, we define
# typedefs that remove the offending characters.  The typedefs also
# get written to the generated header file.
TRIBITS_ETI_GENERATE_TYPEDEF_MACRO(GALERI_ETI_TYPEDEFS "GALERI_ETI_MANGLING_TYPEDEFS" "${eti_typedefs}")

# Generate the header file Galeri_ETIHelperMacros.h, from the file
# Galeri_ETIHelperMacros.h.in (that lives in this directory).  The
# generated header file gets written to the Trilinos build directory,
# in packages/galeri/src/.
CONFIGURE_FILE(${Galeri_SOURCE_DIR}/cmake/Galeri_ETIHelperMacros.h.in ${Galeri_BINARY_DIR}/src/Galeri_ETIHelperMacros.h)
