
INCLUDE(TribitsPackageMacros)
INCLUDE(TribitsAddOptionAndDefine)

#
# A) Define the package
#

TRIBITS_PACKAGE(Komplex)

#
# B) Set up package-specific options
#

# ...

#
# C) Add the libraries, tests, and examples
#

ADD_SUBDIRECTORY(src)

TRIBITS_ADD_TEST_DIRECTORIES(test)

#TRIBITS_ADD_EXAMPLE_DIRECTORIES(example)

#
# Exclude files for source package.
#

TRIBITS_EXCLUDE_AUTOTOOLS_FILES()

TRIBITS_EXCLUDE_FILES(
  doc/Komplex.*.vsd
  doc/header.tex
  doc/komplex.eps
  doc/komplex.gif
  doc/komplex_user_guide.ps
  example/komplex_hb/README
  example/komplex_hb/blassm.f
  example/komplex_hb/create_vbr.c
  example/komplex_hb/distrib_.*_matrix.c
  example/komplex_hb/formats.f
  example/komplex_hb/iohb.*
  example/komplex_hb/main.c
  example/komplex_hb/prototypes.h
  example/komplex_hb/read_.*
  example/komplex_hb/sc.*
  example/komplex_hb/smsrres.c
  example/komplex_hb/svbrres.c
  example/komplex_hb/unary.f
  example/komplex_hb/write_vec.c
  src/new
  src/stamp-h.in
  test/definition
  )

#
# D) Do standard postprocessing
#

TRIBITS_PACKAGE_POSTPROCESS()
