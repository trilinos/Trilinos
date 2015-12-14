include(Join)
MESSAGE(STATUS "${PACKAGE_NAME}: Processing ETI / test support")

# Ifpack2 ETI type fields.  S, LO, GO, N correspond to the four
# template parameters of most Tpetra classes: Scalar, LocalOrdinal,
# GlobalOrdinal, and Node.  Ifpack2 shares these with Tpetra, because
# Ifpack2 only works with Tpetra linear algebra objects.
SET(Ifpack2_ETI_FIELDS "S|LO|GO|N")

# Set up a pattern that excludes all complex Scalar types.
# TriBITS' ETI system knows how to interpret this pattern.
TRIBITS_ETI_TYPE_EXPANSION(${PACKAGE_NAME}_ETI_EXCLUDE_SET_COMPLEX "S=std::complex<float>|std::complex<double>" "LO=.*" "GO=.*" "N=.*")

# TriBITS' ETI system expects a set of types to be a string, delimited
# by |.  Each template parameter (e.g., Scalar, LocalOrdinal, ...) has
# its own set.  The JOIN commands below set up those lists.  We use
# the following sets that Ifpack2 defines:
#
# Scalar: Ifpack2_ETI_SCALARS
# LocalOrdinal: Ifpack2_ETI_LORDS
# GlobalOrdinal: Ifpack2_ETI_GORDS
# Node: Ifpack2_ETI_NODES 
#
# Note that the Scalar set from Tpetra includes the Scalar =
# GlobalOrdinal case.  However, Ifpack2's CMake logic excludes this,
# so we don't have to worry about it here.

JOIN(Ifpack2_ETI_SCALARS "|" FALSE ${Ifpack2_ETI_SCALARS})
JOIN(Ifpack2_ETI_LORDS   "|" FALSE ${Ifpack2_ETI_LORDS}  )
JOIN(Ifpack2_ETI_GORDS   "|" FALSE ${Ifpack2_ETI_GORDS}  )
JOIN(Ifpack2_ETI_NODES   "|" FALSE ${Ifpack2_ETI_NODES}  )  

MESSAGE(STATUS "Enabled Scalar types:        ${Ifpack2_ETI_SCALARS}")
MESSAGE(STATUS "Enabled LocalOrdinal types:  ${Ifpack2_ETI_LORDS}")
MESSAGE(STATUS "Enabled GlobalOrdinal types: ${Ifpack2_ETI_GORDS}")
MESSAGE(STATUS "Enabled Node types:          ${Ifpack2_ETI_NODES}")  

# Construct the "type expansion" string that TriBITS' ETI system
# expects.  Even if ETI is OFF, we will use this to generate macros
# for instantiating tests.
TRIBITS_ETI_TYPE_EXPANSION(SingleScalarInsts 
  "S=${Ifpack2_ETI_SCALARS}" 
  "N=${Ifpack2_ETI_NODES}"
  "LO=${Ifpack2_ETI_LORDS}" 
  "GO=${Ifpack2_ETI_GORDS}")

ASSERT_DEFINED(Ifpack2_ENABLE_EXPLICIT_INSTANTIATION)
IF(Ifpack2_ENABLE_EXPLICIT_INSTANTIATION)
  # mfh 17 Aug 2015: If ETI is ON, it looks like users can set these variables.
  MESSAGE(STATUS "User/Downstream ETI set: ${Ifpack2_ETI_LIBRARYSET}")
  TRIBITS_ADD_ETI_INSTANTIATIONS(Ifpack2 ${SingleScalarInsts})
  MESSAGE(STATUS "Excluded type combinations: ${Ifpack2_ETI_EXCLUDE_SET}")
ELSE()
  TRIBITS_ETI_TYPE_EXPANSION(Ifpack2_ETI_LIBRARYSET
    "S=${Ifpack2_ETI_SCALARS}" 
    "N=${Ifpack2_ETI_NODES}"
    "LO=${Ifpack2_ETI_LORDS}" 
    "GO=${Ifpack2_ETI_GORDS}")
ENDIF()
MESSAGE(STATUS "Set of enabled types, before exclusions: ${${PACKAGE_NAME}_ETI_LIBRARYSET}")

#
# Generate the instantiation macros.  These go into
# Ifpack2_ETIHelperMacros.h, which is generated from
# Ifpack2_ETIHelperMacros.h.in (in this directory).
#
TRIBITS_ETI_GENERATE_MACROS("${Ifpack2_ETI_FIELDS}" "${Ifpack2_ETI_LIBRARYSET}" "${Ifpack2_ETI_EXCLUDE_SET}"  
                            list_of_manglings eti_typedefs
                            "IFPACK2_INSTANTIATE_L(LO)"         IFPACK2_ETIMACRO_L
                            "IFPACK2_INSTANTIATE_SL(S,LO)"      IFPACK2_ETIMACRO_SL
                            "IFPACK2_INSTANTIATE_LG(LO,GO)"         IFPACK2_ETIMACRO_LG
                            "IFPACK2_INSTANTIATE_SLG(S,LO,GO)"      IFPACK2_ETIMACRO_SLG
                            "IFPACK2_INSTANTIATE_LGN(S,LO,GO,N)"      IFPACK2_ETIMACRO_LGN
                            "IFPACK2_INSTANTIATE_SLGN(S,LO,GO,N)"     IFPACK2_ETIMACRO_SLGN   
                            "IFPACK2_INSTANTIATE_N(N)"     IFPACK2_ETIMACRO_N                          
                            )
TRIBITS_ETI_GENERATE_MACROS("${Ifpack2_ETI_FIELDS}" "${Ifpack2_ETI_LIBRARYSET}" 
                            "${Ifpack2_ETI_EXCLUDE_SET};${Ifpack2_ETI_EXCLUDE_SET_COMPLEX}"
                            list_of_manglings eti_typedefs
                            "IFPACK2_INSTANTIATE_SL_REAL(S,LO,GO)" IFPACK2_ETIMACRO_SL_REAL
                            "IFPACK2_INSTANTIATE_SLG_REAL(S,LO,GO)" IFPACK2_ETIMACRO_SLG_REAL
                            "IFPACK2_INSTANTIATE_SLGN_REAL(S,LO,GO,N)" IFPACK2_ETIMACRO_SLGN_REAL
                            )

# Generate "mangled" typedefs.  Macros sometimes get grumpy when types
# have spaces, colons, or angle brackets in them.  This includes types
# like "long long" or "std::complex<double>".  Thus, we define
# typedefs that remove the offending characters.  The typedefs also
# get written to the generated header file.
TRIBITS_ETI_GENERATE_TYPEDEF_MACRO(IFPACK2_ETI_TYPEDEFS "IFPACK2_ETI_MANGLING_TYPEDEFS" "${eti_typedefs}")

# Generate the header file Ifpack2_ETIHelperMacros.h, from the file
# Ifpack2_ETIHelperMacros.h.in (that lives in this directory).  The
# generated header file gets written to the Trilinos build directory,
# in packages/ifpack2/src/.
CONFIGURE_FILE(${Ifpack2_SOURCE_DIR}/cmake/Ifpack2_ETIHelperMacros.h.in ${Ifpack2_BINARY_DIR}/src/Ifpack2_ETIHelperMacros.h)
