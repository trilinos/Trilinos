#@HEADER
#@HEADER

#ML_STANDALONE_MODE is set only if ml is being built by itself.
#It controls which linker is used:  in the case that ml is built
#with aztecoo and epetra support enabled, then a c++ linker is invoked
#for ML's interfaces to those two packages.

EXEEXT = .exe

# List the names of all exe's in this directory here
if ML_STANDALONE_MODE

if HAVE_MPI
noinst_PROGRAMS = mlguide mlguide_par
else
noinst_PROGRAMS = mlguide
endif

else

if HAVE_MPI
if ML_BENCHMARKING
noinst_PROGRAMS = mlguide mlguide_par ml_aztec_simple ml_example1d \
                  ml_example1dGS ml_example2d ml_example3d ml_simple_max \
                  ml_aztec_simple_METIS ml_aztec_simple_ParMETIS \
		  mlaz_example \
		  ml_two_level_DD \
                  ml_readfile ml_read_elas ml_recirc \
                  ml_read_salsa ml_read_maxwell \
                  ml_example_epetra_preconditioner
else
noinst_PROGRAMS = mlguide mlguide_par ml_aztec_simple ml_example1d \
                  ml_example1dGS ml_example2d ml_example3d ml_simple_max \
                  ml_aztec_simple_METIS ml_aztec_simple_ParMETIS \
		  mlaz_example \
		  ml_two_level_DD \
		          ml_readfile ml_read_elas ml_recirc \
                  ml_read_salsa \
                  ml_example_epetra_preconditioner
endif
else
noinst_PROGRAMS = mlguide ml_aztec_simple ml_example1d \
                  ml_example1dGS ml_example2d ml_example3d ml_simple_max \
                  ml_aztec_simple_METIS ml_aztec_simple_ParMETIS \
                  mlaz_example \
		  ml_two_level_DD \
		  ml_readfile ml_read_elas ml_recirc \
                  ml_read_salsa \
                  ml_example_epetra_preconditioner
endif

endif
#The following line cleans up bootstrap problems because these programs
#are not built when nothing is defined.
EXTRA_PROGRAMS = ml_read_maxwell mlguide_par

# One SOURCES primary is needed for EACH exe listed above.  (Notice for all 
# of the  following primaries that the pattern is 
# "(Exe_name)_(PRIMARY_NAME) = ").
# List all source files needed to compile the exe.
mlguide_SOURCES = $(srcdir)/mlguide.c
mlguide_par_SOURCES = $(srcdir)/mlguide_par.c
if !ML_STANDALONE_MODE
ml_aztec_simple_SOURCES = $(srcdir)/ml_aztec_simple.c
ml_aztec_simple_METIS_SOURCES = $(srcdir)/ml_aztec_simple_METIS.c
ml_aztec_simple_ParMETIS_SOURCES = $(srcdir)/ml_aztec_simple_ParMETIS.c
mlaz_example_SOURCES = $(srcdir)/mlaz_example.c
ml_two_level_DD_SOURCES = $(srcdir)/ml_two_level_DD.c
ml_example1d_SOURCES = $(srcdir)/ml_example1d.c
ml_example1dGS_SOURCES = $(srcdir)/ml_example1dGS.c
ml_example2d_SOURCES = $(srcdir)/ml_example2d.c
ml_example3d_SOURCES = $(srcdir)/ml_example3d.c
ml_simple_max_SOURCES = $(srcdir)/ml_simple_max.c
ml_readfile_SOURCES = $(srcdir)/ml_readfile.c
ml_read_elas_SOURCES = $(srcdir)/ml_read_elas.c
ml_recirc_SOURCES = $(srcdir)/ml_recirc.c
ml_read_salsa_SOURCES = $(srcdir)/ml_read_salsa.c
ml_example_epetra_preconditioner_SOURCES = $(srcdir)/ml_example_epetra_preconditioner.cpp
if ML_BENCHMARKING
ml_read_maxwell_SOURCES = $(srcdir)/ml_read_maxwell.c
endif
endif

# One DEPENDENCIES primary is needed for each exe listed above.  List any
# libraries that the exe depends on.  (This primary is not absolutely
# essential.)

ml_integrated_dependencies = $(top_builddir)/../aztecoo/src/libaztecoo.a \
                       $(top_builddir)/../epetra/src/libepetra.a

# This is not working.  
# It should be ok to omit it because when _DEPENDENCIES is not
# specified, the value is calculated from the corresponding _LDADD value.
#if ML_STANDALONE_MODE
#mlguide_DEPENDENCIES =
#else
#mlguide_DEPENDENCIES = ${ml_integrated_dependencies}
#endif


ml_common_include_locations = -I$(top_srcdir)/src/Include          \
                   -I$(top_srcdir)/src/Coarsen                  \
                   -I$(top_srcdir)/src/Comm                  \
                   -I$(top_srcdir)/src/FEGrid                  \
                   -I$(top_srcdir)/src/Include                  \
                   -I$(top_srcdir)/src/Krylov                  \
                   -I$(top_srcdir)/src/Main                  \
                   -I$(top_srcdir)/src/Operator                  \
                   -I$(top_srcdir)/src/Smoother                  \
                   -I$(top_srcdir)/src/Utils                  \
                   -I$(top_srcdir)/../aztecoo/src            \
                   -I$(top_builddir)/../aztecoo/src \
                   -I$(top_srcdir)/../epetra/src            \
                   -I$(top_builddir)/../epetra/src \
                   -I$(top_srcdir)/../teuchos/src            \
                   -I$(top_builddir)/../teuchos/src \
                   -I$(top_srcdir)/../triutils/src            \
                   -I$(top_builddir)/../triutils/src

if ML_STANDALONE_MODE
#I'm defining the linker in this way, rather than
#       ml_invididual_example_LINK = ${CC}
#because the above definition will not use the compiler and loader flags
#specified in AM_CFLAGS, CFLAGS, AM_LDFLAGS, and LDFLAGS.
mllink = $(CC) $(AM_CFLAGS) $(CFLAGS) $(AM_LDFLAGS) $(LDFLAGS) -o $@
ml_libs = $(top_builddir)/src/libml.a        \
                $(LIBS)             \
                $(LAPACK_LIBS)                            \
                $(BLAS_LIBS)                            \
                $(FLIBS)
else
mllink = $(CXX) $(AM_CFLAGS) $(CFLAGS) $(AM_LDFLAGS) $(LDFLAGS) -o $@
ml_libs = $(top_builddir)/src/libml.a        \
                $(top_builddir)/../aztecoo/src/libaztecoo.a        \
                $(top_builddir)/../epetra/src/libepetra.a         \
                $(LIBS) \
                $(LAPACK_LIBS)                            \
                $(BLAS_LIBS)                            \
                $(FLIBS)
endif

# One CXXFLAGS primary can be defined for each exe listed above.  Include 
# any -I's necessary for building the exe.

# One LDADD primary can be defined for each exe listed above.  Include any
# LDFLAGS  necessary for building the exe.

mlguide_CFLAGS = ${ml_common_include_locations}
mlguide_LDADD = ${ml_libs}
mlguide_LINK = ${mllink}

mlguide_par_CFLAGS = ${ml_common_include_locations}
mlguide_par_LDADD = ${ml_libs}
mlguide_par_LINK = ${mllink}

if !ML_STANDALONE_MODE
ml_aztec_simple_CFLAGS = ${ml_common_include_locations}
ml_aztec_simple_LDADD = ${ml_libs}
ml_aztec_simple_LINK = ${mllink}

ml_example1d_CFLAGS = ${ml_common_include_locations}
ml_example1d_LDADD = ${ml_libs}
ml_example1d_LINK = ${mllink}

ml_example1dGS_CFLAGS = ${ml_common_include_locations}
ml_example1dGS_LDADD = ${ml_libs}
ml_example1dGS_LINK = ${mllink}

ml_example2d_CFLAGS = ${ml_common_include_locations}
ml_example2d_LDADD = ${ml_libs}
ml_example2d_LINK = ${mllink}

ml_example3d_CFLAGS = ${ml_common_include_locations}
ml_example3d_LDADD = ${ml_libs}
ml_example3d_LINK = ${mllink}

ml_simple_max_CFLAGS = ${ml_common_include_locations}
ml_simple_max_LDADD = ${ml_libs}
ml_simple_max_LINK = ${mllink}

ml_aztec_simple_METIS_CFLAGS = ${ml_common_include_locations}
ml_aztec_simple_METIS_LDADD = ${ml_libs}
ml_aztec_simple_METIS_LINK = ${mllink}

ml_aztec_simple_ParMETIS_CFLAGS = ${ml_common_include_locations}
ml_aztec_simple_ParMETIS_LDADD = ${ml_libs}
ml_aztec_simple_ParMETIS_LINK = ${mllink}

mlaz_example_CFLAGS = ${ml_common_include_locations}
mlaz_example_LDADD = ${ml_libs}
mlaz_example_LINK = ${mllink}

ml_two_level_DD_CFLAGS = ${ml_common_include_locations}
ml_two_level_DD_LDADD = ${ml_libs}
ml_two_level_DD_LINK = ${mllink}

ml_readfile_CFLAGS = ${ml_common_include_locations}
ml_readfile_LDADD = ${ml_libs}
ml_readfile_LINK = ${mllink}

ml_read_elas_CFLAGS = ${ml_common_include_locations}
ml_read_elas_LDADD = ${ml_libs}
ml_read_elas_LINK = ${mllink}

ml_recirc_CFLAGS = ${ml_common_include_locations}
ml_recirc_LDADD = ${ml_libs}
ml_recirc_LINK = ${mllink}

ml_read_salsa_CFLAGS = ${ml_common_include_locations}
ml_read_salsa_LDADD = ${ml_libs}
ml_read_salsa_LINK = ${mllink}

if ML_USING_TEUCHOS

ml_example_epetra_preconditioner_CPPFLAGS = ${ml_common_include_locations}
ml_example_epetra_preconditioner_CXXFLAGS = ${ml_common_include_locations}
ml_example_epetra_preconditioner_LDADD = ${ml_libs} ../../teuchos/src/libteuchos.a ../../triutils/src/libtriutils.a ../../epetra/src/libepetra.a
ml_example_epetra_preconditioner_LINK = ${mllink}
endif

if ML_BENCHMARKING
ml_read_maxwell_CFLAGS = ${ml_common_include_locations}
ml_read_maxwell_LDADD = ${ml_libs}
ml_read_maxwell_LINK = ${mllink}
endif
endif
