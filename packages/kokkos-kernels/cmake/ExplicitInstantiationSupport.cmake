message(STATUS "${PACKAGE_NAME}: Processing ETI / test support")

# mfh 11 Oct 2016: Hack to work around #701.  Strip out anything that
# looks like a Kokkos "Node" type (not the same as a Device type!!!)
# from the manglings and typedefs lists.

set(list_of_manglings "")
set(eti_typedefs "")

# This CMake module generates the following header file, which gets
# written to the build directory (like other header files that CMake
# generates).  The file contains macros that do instantiation over a
# finite set of template parameter combinations.  We use the macros
# both for ETI (explicit template instantiation), and for tests.
# Thus, the macros need to be generated even if ETI is OFF.
set(${PACKAGE_NAME}_ETI_FILE ${PACKAGE_NAME}_ETIHelperMacros.h)
set(${PACKAGE_NAME}_ETI_FILE_PATH ${${PACKAGE_NAME}_BINARY_DIR}/src/${${PACKAGE_NAME}_ETI_FILE})

#
# Users have the option to generate the above header file themselves.
# We prefer that users let Trilinos generate the header file.
# However, folks who make intense use of TriBITS sometimes find that
# reusing a previously generated header file shaves a couple minutes
# off their CMake configuration time.  Thus, we give them that option.
#

advanced_set(${PACKAGE_NAME}_USE_STATIC_ETI_MACROS_HEADER_FILE "" CACHE PATH
  "If set, gives the path to a static version of the file ${${PACKAGE_NAME}_ETI_FILE}.  If not set (the default), then the file is generated automatically.")

if(${PACKAGE_NAME}_USE_STATIC_ETI_MACROS_HEADER_FILE)
  # The user wants us to accept their header file and not generate one.
  message("-- NOTE: Skipping generation and using provided static file"
          " '${${PACKAGE_NAME}_USE_STATIC_ETI_MACROS_HEADER_FILE}'")
  kokkos_configure_file(
    ${${PACKAGE_NAME}_USE_STATIC_ETI_MACROS_HEADER_FILE}
    ${${PACKAGE_NAME}_ETI_FILE_PATH}
    COPYONY)
  return()
endif()

#
# The user wants us to generate the header file.  This is the default
# behavior.
#

# Tpetra ETI type fields.  S, LO, D correspond to the template
# parameters Scalar, LocalOrdinal, and DeviceType.  KokkosKernels does
# not need to know about GlobalOrdinal so we omit that.
set(${PACKAGE_NAME}_ETI_FIELDS "S|LO|D")

# Set up a pattern that excludes all Scalar types that are also
# possible GlobalOrdinal types, but includes all other types.
# TriBITS' ETI system knows how to interpret this pattern.
#
# FIXME (mfh 17 Aug 2015, 16 Oct 2015) A better way to do this would
# be to subtract away all enabled GlobalOrdinal types.  Plus, what if
# someone really wants a CrsMatrix<int,...>?
kokkos_eti_type_expansion(
  ${PACKAGE_NAME}_ETI_EXCLUDE_SET_ORDINAL_SCALAR
  "S=short|short int|unsigned short|unsigned short int|int|unsigned|unsigned int|long|long int|unsigned long|unsigned long int|long long|long long int|unsigned long long|unsigned long long int|int16_t|uint16_t|int32_t|uint32_t|int64_t|uint64_t|size_t|ptrdiff_t"
  "LO=.*"
  "D=.*")

# TriBITS' ETI system expects a set of types to be a string, delimited
# by |.  Each template parameter (e.g., Scalar, LocalOrdinal, ...) has
# its own set.  The JOIN commands below set up those lists.  We use
# the following sets that this subpackage defines:
#
# Scalar:       ${PACKAGE_NAME}_ETI_SCALARS
# LocalOrdinal: ${PACKAGE_NAME}_ETI_LORDS
# Device:       ${PACKAGE_NAME}_ETI_DEVICES
#
# Note that the Scalar set from Tpetra includes the Scalar =
# GlobalOrdinal case.  We have to exclude that explicitly in what
# follows.
join(${PACKAGE_NAME}_ETI_SCALARS "|" FALSE ${${PACKAGE_NAME}_ETI_SCALARS})
join(${PACKAGE_NAME}_ETI_LORDS "|" FALSE ${${PACKAGE_NAME}_ETI_LORDS})
join(${PACKAGE_NAME}_ETI_DEVICES "|" FALSE ${${PACKAGE_NAME}_ETI_DEVICES})

message(STATUS "Enabled Scalar types:       ${${PACKAGE_NAME}_ETI_SCALARS}")
message(STATUS "Enabled LocalOrdinal types: ${${PACKAGE_NAME}_ETI_LORDS}")
message(STATUS "Enabled Device types:       ${${PACKAGE_NAME}_ETI_DEVICES}")

# Set up the set of enabled type combinations, in a format that
# TriBITS understands.
#
# mfh 17 Aug 2015, 16 Oct 2015: I don't exactly understand what's
# going on here, but it looks like if ETI is enabled, we let users
# modify ${PACKAGE_NAME}_ETI_LIBRARYSET, and if it's not, we don't.
assert_defined(${PACKAGE_NAME}_ENABLE_EXPLICIT_INSTANTIATION)

# Construct the "type expansion" string that TriBITS' ETI system
# expects.  Even if ETI is OFF, we will use this to generate macros
# for instantiating tests.
kokkos_eti_type_expansion(
  ${PACKAGE_NAME}_ETI_LIBRARYSET
  "S=${${PACKAGE_NAME}_ETI_SCALARS}"
  "LO=${${PACKAGE_NAME}_ETI_LORDS}"
  "D=${${PACKAGE_NAME}_ETI_DEVICES}")

# Construct the "type expansion" string that TriBITS' ETI system
# expects.  Even if ETI is OFF, we will use this to generate macros
# for instantiating tests.
kokkos_eti_type_expansion(
  ${PACKAGE_NAME}_ETI_LIBRARYSET_ORDINAL_SCALAR
  "S=${${PACKAGE_NAME}_ETI_SCALARS_ORDS}"
  "LO=${${PACKAGE_NAME}_ETI_LORDS}"
  "D=${${PACKAGE_NAME}_ETI_DEVICES}")

kokkos_add_eti_instantiations(${PACKAGE_NAME} ${${PACKAGE_NAME}_ETI_LIBRARYSET})

message(STATUS "Set of enabled types, before exclusions: ${${PACKAGE_NAME}_ETI_LIBRARYSET}")

#
# Generate the instantiation macros.  These go into
# ${PACKAGE_NAME}_ETIHelperMacros.h, which is generated from
# ${PACKAGE_NAME}_ETIHelperMacros.h.in (in this directory).
#

# Generate macros that exclude possible ordinal types from the
# list of Scalar types.
kokkos_eti_generate_macros(
  "${${PACKAGE_NAME}_ETI_FIELDS}"
  "${${PACKAGE_NAME}_ETI_LIBRARYSET}"
  "${${PACKAGE_NAME}_ETI_EXCLUDE_SET};${${PACKAGE_NAME}_ETI_EXCLUDE_SET_ORDINAL_SCALAR}"
  list_of_manglings
  eti_typedefs
  "KOKKOSKERNELS_INSTANTIATE_SLD_NO_ORDINAL_SCALAR(S,LO,D)" KOKKOSKERNELS_INSTANTIATE_SLD_NO_ORDINAL_SCALAR
  "KOKKOSKERNELS_INSTANTIATE_SL_NO_ORDINAL_SCALAR(S,LO)"    KOKKOSKERNELS_INSTANTIATE_SL_NO_ORDINAL_SCALAR
  "KOKKOSKERNELS_INSTANTIATE_SD_NO_ORDINAL_SCALAR(S,D)"     KOKKOSKERNELS_INSTANTIATE_SD_NO_ORDINAL_SCALAR
  "KOKKOSKERNELS_INSTANTIATE_S_NO_ORDINAL_SCALAR(S)"        KOKKOSKERNELS_INSTANTIATE_S_NO_ORDINAL_SCALAR)

# Generate macros include ONLY possible ordinal types in the list of Scalar types.
kokkos_eti_generate_macros(
  "${${PACKAGE_NAME}_ETI_FIELDS}"
  "${${PACKAGE_NAME}_ETI_LIBRARYSET_ORDINAL_SCALAR}"
  "${${PACKAGE_NAME}_ETI_EXCLUDE_SET}"
  list_of_manglings
  eti_typedefs
  "KOKKOSKERNELS_INSTANTIATE_SLD_ORDINAL_SCALAR(S,LO,D)" KOKKOSKERNELS_INSTANTIATE_SLD_ORDINAL_SCALAR
  "KOKKOSKERNELS_INSTANTIATE_SL_ORDINAL_SCALAR(S,LO)"    KOKKOSKERNELS_INSTANTIATE_SL_ORDINAL_SCALAR
  "KOKKOSKERNELS_INSTANTIATE_SD_ORDINAL_SCALAR(S,D)"     KOKKOSKERNELS_INSTANTIATE_SD_ORDINAL_SCALAR
  "KOKKOSKERNELS_INSTANTIATE_S_ORDINAL_SCALAR(S)"        KOKKOSKERNELS_INSTANTIATE_S_ORDINAL_SCALAR)

# Generate macros that include all Scalar types (if applicable),
# including possible ordinal types.
kokkos_eti_generate_macros(
  "${${PACKAGE_NAME}_ETI_FIELDS}"
  "${${PACKAGE_NAME}_ETI_LIBRARYSET}"
  "${${PACKAGE_NAME}_ETI_EXCLUDE_SET}"
  list_of_manglings
  eti_typedefs
  "KOKKOSKERNELS_INSTANTIATE_SLD(S,LO,D)" KOKKOSKERNELS_INSTANTIATE_SLD
  "KOKKOSKERNELS_INSTANTIATE_SL(S,LO)"    KOKKOSKERNELS_INSTANTIATE_SL
  "KOKKOSKERNELS_INSTANTIATE_SD(S,D)"     KOKKOSKERNELS_INSTANTIATE_SD
  "KOKKOSKERNELS_INSTANTIATE_S(S)"        KOKKOSKERNELS_INSTANTIATE_S
  "KOKKOSKERNELS_INSTANTIATE_LD(LO,D)"    KOKKOSKERNELS_INSTANTIATE_LD
  "KOKKOSKERNELS_INSTANTIATE_L(LO)"       KOKKOSKERNELS_INSTANTIATE_L
  "KOKKOSKERNELS_INSTANTIATE_D(D)"        KOKKOSKERNELS_INSTANTIATE_D)

# Generate "mangled" typedefs.  Macros sometimes get grumpy when types
# have spaces, colons, or angle brackets in them.  This includes types
# like "long long" or "std::complex<double>".  Thus, we define
# typedefs that remove the offending characters.  The typedefs also
# get written to the generated header file.
kokkos_eti_generate_typedef_macro(KOKKOSKERNELS_ETI_MANGLING_TYPEDEFS "KOKKOSKERNELS_ETI_MANGLING_TYPEDEFS" "${eti_typedefs}")

# Generate the header file ${PACKAGE_NAME}_ETIHelperMacros.h, from the
# file ${PACKAGE_NAME}_ETIHelperMacros.h.in (that lives in this
# directory).  The generated header file gets written to the Trilinos
# build directory, in packages/tpetra/kernels/src/.
kokkos_configure_file(
  ${${PACKAGE_NAME}_SOURCE_DIR}/cmake/${PACKAGE_NAME}_ETIHelperMacros.h.in
  ${${PACKAGE_NAME}_ETI_FILE_PATH})
