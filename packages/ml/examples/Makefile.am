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
endif #HAVE_MPI

else # not standalone

if HAVE_MPI
if ML_BENCHMARKING
noinst_PROGRAMS = mlguide mlguide_par ml_aztec_simple ml_example1d \
                  ml_example1dGS ml_example2d ml_example3d ml_simple_max \
                  ml_aztec_simple_METIS ml_aztec_simple_ParMETIS \
		  mlaz_example \
		  ml_two_level_DD \
                  ml_readfile ml_read_elas ml_recirc \
                  ml_read_salsa ml_read_maxwell \
                  ml_example_MultiLevelPreconditioner \
                  ml_example_MultiLevelPreconditioner_2level \
                  ml_example_MultiLevelOperator \
                  ml_example_MultiLevelPreconditioner_AztecMSR \
                  ml_example_MultiLevelPreconditioner_filtering \
		  ml_example_MultiLevelPreconditioner_vbr \
		  ml_example_MultiLevelPreconditioner_viz \
		  ml_example_MultiLevelPreconditioner_analyze \
		  ml_example_MultiLevelPreconditioner_maxwell 
else
noinst_PROGRAMS = mlguide mlguide_par ml_aztec_simple ml_example1d \
                  ml_example1dGS ml_example2d ml_example3d ml_simple_max \
                  ml_aztec_simple_METIS ml_aztec_simple_ParMETIS \
		  mlaz_example \
		  ml_two_level_DD \
		          ml_readfile ml_read_elas ml_recirc \
                  ml_read_salsa \
                  ml_example_MultiLevelPreconditioner \
                  ml_example_MultiLevelPreconditioner_2level \
                  ml_example_MultiLevelOperator \
                  ml_example_MultiLevelPreconditioner_AztecMSR \
                  ml_example_MultiLevelPreconditioner_filtering \
		  ml_example_MultiLevelPreconditioner_vbr \
		  ml_example_MultiLevelPreconditioner_viz \
		  ml_example_MultiLevelPreconditioner_analyze \
		  ml_example_MultiLevelPreconditioner_maxwell  \
          ml_multiple_rhs ml_aztecoo

endif #ML_BENCHMARKING
else #no mpi
noinst_PROGRAMS = mlguide ml_aztec_simple ml_example1d \
                  ml_example1dGS ml_example2d ml_example3d ml_simple_max \
                  ml_aztec_simple_METIS ml_aztec_simple_ParMETIS \
                  mlaz_example \
		  ml_two_level_DD \
		  ml_readfile ml_read_elas ml_recirc \
                  ml_read_salsa \
                  ml_example_MultiLevelPreconditioner \
                  ml_example_MultiLevelPreconditioner_2level \
                  ml_example_MultiLevelOperator \
                  ml_example_MultiLevelPreconditioner_AztecMSR \
                  ml_example_MultiLevelPreconditioner_filtering \
		  ml_example_MultiLevelPreconditioner_vbr \
		  ml_example_MultiLevelPreconditioner_viz \
		  ml_example_MultiLevelPreconditioner_analyze \
		  ml_example_MultiLevelPreconditioner_maxwell  \
          ml_multiple_rhs ml_aztecoo

endif #HAVE_MPI

endif # ML_STANDALONE_MODE
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
ml_example_MultiLevelPreconditioner_SOURCES = $(srcdir)/ml_example_MultiLevelPreconditioner.cpp
ml_example_MultiLevelPreconditioner_2level_SOURCES = $(srcdir)/ml_example_MultiLevelPreconditioner_2level.cpp
ml_example_MultiLevelOperator_SOURCES = $(srcdir)/ml_example_MultiLevelOperator.cpp
ml_example_MultiLevelPreconditioner_AztecMSR_SOURCES = $(srcdir)/ml_example_MultiLevelPreconditioner_AztecMSR.cpp
ml_example_MultiLevelPreconditioner_filtering_SOURCES = $(srcdir)/ml_example_MultiLevelPreconditioner_filtering.cpp
ml_example_MultiLevelPreconditioner_vbr_SOURCES= $(srcdir)/ml_example_MultiLevelPreconditioner_vbr.cpp
ml_example_MultiLevelPreconditioner_viz_SOURCES= $(srcdir)/ml_example_MultiLevelPreconditioner_viz.cpp
ml_example_MultiLevelPreconditioner_analyze_SOURCES= $(srcdir)/ml_example_MultiLevelPreconditioner_analyze.cpp
ml_example_MultiLevelPreconditioner_maxwell_SOURCES = $(srcdir)/ml_example_MultiLevelPreconditioner_maxwell.cpp
ml_multiple_rhs_SOURCES = $(srcdir)/ml_multiple_rhs.cpp
ml_aztecoo_SOURCES = $(srcdir)/ml_aztecoo.cpp

if ML_BENCHMARKING
ml_read_maxwell_SOURCES = $(srcdir)/ml_read_maxwell.c
endif # ML_BENCHMARKING
endif # !ML_STANDALONE_MODE

# One DEPENDENCIES primary is needed for each exe listed above.  List any
# libraries that the exe depends on.  (This primary is not absolutely
# essential.)

# This is not working.  
# It should be ok to omit it because when _DEPENDENCIES is not
# specified, the value is calculated from the corresponding _LDADD value.
#if ML_STANDALONE_MODE
#mlguide_DEPENDENCIES =
#else
#mlguide_DEPENDENCIES = ${ml_integrated_dependencies}
#endif

if ML_USING_AZTECOO
aztecoo_incdir = -I$(top_srcdir)/../aztecoo/src -I$(top_builddir)/../aztecoo/src
else
aztecoo_incdir = 
endif
if ML_USING_AMESOS
amesos_incdir = -I$(top_srcdir)/../amesos/src -I$(top_builddir)/../amesos/src
else
amesos_incdir = 
endif
if ML_USING_IFPACK
ifpack_incdir = -I$(top_srcdir)/../ifpack/src -I$(top_builddir)/../ifpack/src 
else
ifpack_incdir = 
endif
if ML_USING_TEUCHOS
teuchos_incdir = -I$(top_srcdir)/../teuchos/src -I$(top_builddir)/../teuchos/src
else
teucho_incdir =
endif
if ML_USING_ANASAZI
anasazi_incdir = -I$(top_srcdir)/../anasazi/src -I$(top_builddir)/../anasazi/src
else
anasazi_incdir =
endif
if ML_USING_TRIUTILS
triutils_incdir = -I$(top_srcdir)/../triutils/src -I$(top_builddir)/../triutils/src
else
triutils_incdir =
endif
if ML_USING_EPETRA
epetra_incdir = -I$(top_srcdir)/../epetra/src -I$(top_builddir)/../epetra/src
else
epetra_incdir =
endif

ml_common_include_locations = -I$(top_srcdir)/src/Include \
                   -I$(top_srcdir)/src/Coarsen            \
                   -I$(top_srcdir)/src/Comm               \
                   -I$(top_srcdir)/src/FEGrid             \
                   -I$(top_srcdir)/src/Include            \
                   -I$(top_srcdir)/src/Krylov             \
                   -I$(top_srcdir)/src/Main               \
                   -I$(top_srcdir)/src/Operator           \
                   -I$(top_srcdir)/src/Smoother           \
                   -I$(top_srcdir)/src/Utils              \
	           $(aztecoo_incdir) \
	           $(ifpack_incdir) \
	           $(amesos_incdir) \
	           $(teuchos_incdir) \
	           $(anasazi_incdir) \
	           $(triutils_incdir) \
	           $(epetra_incdir)

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

if ML_USING_AZTECOO
aztecoo_lib = $(top_builddir)/../aztecoo/src/libaztecoo.a
else
aztecoo_lib =
endif

if ML_USING_EPETRA
epetra_lib = $(top_builddir)/../epetra/src/libepetra.a
else
epetra_lib =
endif

if ML_USING_AMESOS
amesos_lib = $(top_builddir)/../amesos/src/libamesos.a
else
amesos_lib =
endif

if ML_USING_ANASAZI
anasazi_lib = $(top_builddir)/../anasazi/src/libanasazi.a
else
anasazi_lib = 
endif

if ML_USING_TEUCHOS 
teuchos_lib = $(top_builddir)/../teuchos/src/libteuchos.a
else
teuchos_lib =
endif

if ML_USING_TRIUTILS
triutils_lib = $(top_builddir)/../triutils/src/libtriutils.a
else
triutils_lib =
endif

if ML_USING_IFPACK
ifpack_lib = $(top_builddir)/../ifpack/src/libifpack.a
else
ifpack_lib =
endif
#
ml_libs = $(top_builddir)/src/libml.a        \
                $(aztecoo_lib) \
                $(ifpack_lib) \
                $(amesos_lib) \
                $(anasazi_lib) \
                $(teuchos_lib) \
                $(triutils_lib) \
                $(epetra_lib) \
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

ml_example_MultiLevelPreconditioner_CPPFLAGS = ${ml_common_include_locations}
ml_example_MultiLevelPreconditioner_CXXFLAGS = ${ml_common_include_locations}
ml_example_MultiLevelPreconditioner_LDADD = ${ml_libs} 
ml_example_MultiLevelPreconditioner_LINK = ${mllink}

ml_example_MultiLevelPreconditioner_2level_CPPFLAGS = ${ml_common_include_locations}
ml_example_MultiLevelPreconditioner_2level_CXXFLAGS = ${ml_common_include_locations}
ml_example_MultiLevelPreconditioner_2level_LDADD = ${ml_libs} 
ml_example_MultiLevelPreconditioner_2level_LINK = ${mllink}

ml_example_MultiLevelPreconditioner_AztecMSR_CPPFLAGS = ${ml_common_include_locations}
ml_example_MultiLevelPreconditioner_AztecMSR_CXXFLAGS = ${ml_common_include_locations}
ml_example_MultiLevelPreconditioner_AztecMSR_LDADD = ${ml_libs} 
ml_example_MultiLevelPreconditioner_AztecMSR_LINK = ${mllink}

ml_example_MultiLevelOperator_CPPFLAGS = ${ml_common_include_locations}
ml_example_MultiLevelOperator_CXXFLAGS = ${ml_common_include_locations}
ml_example_MultiLevelOperator_LDADD = ${ml_libs}
ml_example_MultiLevelOperator_LINK = ${mllink}

ml_example_MultiLevelPreconditioner_filtering_CPPFLAGS = ${ml_common_include_locations}
ml_example_MultiLevelPreconditioner_filtering_CXXFLAGS = ${ml_common_include_locations}
ml_example_MultiLevelPreconditioner_filtering_LDADD = ${ml_libs}
ml_example_MultiLevelPreconditioner_filtering_LINK = ${mllink}

ml_example_MultiLevelPreconditioner_vbr_CPPFLAGS = ${ml_common_include_locations}
ml_example_MultiLevelPreconditioner_vbr_CXXFLAGS = ${ml_common_include_locations}
ml_example_MultiLevelPreconditioner_vbr_LDADD = ${ml_libs}
ml_example_MultiLevelPreconditioner_vbr_LINK = ${mllink}

ml_example_MultiLevelPreconditioner_viz_CPPFLAGS = ${ml_common_include_locations}
ml_example_MultiLevelPreconditioner_viz_CXXFLAGS = ${ml_common_include_locations}
ml_example_MultiLevelPreconditioner_viz_LDADD = ${ml_libs}
ml_example_MultiLevelPreconditioner_viz_LINK = ${mllink}

ml_example_MultiLevelPreconditioner_analyze_CPPFLAGS = ${ml_common_include_locations}
ml_example_MultiLevelPreconditioner_analyze_CXXFLAGS = ${ml_common_include_locations}
ml_example_MultiLevelPreconditioner_analyze_LDADD = ${ml_libs}
ml_example_MultiLevelPreconditioner_analyze_LINK = ${mllink}

ml_example_MultiLevelPreconditioner_maxwell_CPPFLAGS = ${ml_common_include_locations}
ml_example_MultiLevelPreconditioner_maxwell_CXXFLAGS = ${ml_common_include_locations}
ml_example_MultiLevelPreconditioner_maxwell_LDADD = ${ml_libs}
ml_example_MultiLevelPreconditioner_maxwell_LINK = ${mllink}

ml_multiple_rhs_CPPFLAGS = ${ml_common_include_locations}
ml_multiple_rhs_CXXFLAGS = ${ml_common_include_locations}
ml_multiple_rhs_LDADD = ${ml_libs}
ml_multiple_rhs_LINK = ${mllink}

ml_aztecoo_CPPFLAGS = ${ml_common_include_locations}
ml_aztecoo_CXXFLAGS = ${ml_common_include_locations}
ml_aztecoo_LDADD = ${ml_libs}
ml_aztecoo_LINK = ${mllink}

if ML_BENCHMARKING
ml_read_maxwell_CFLAGS = ${ml_common_include_locations}
ml_read_maxwell_LDADD = ${ml_libs}
ml_read_maxwell_LINK = ${mllink}
endif #ML_BENCHMARKING
endif #!ML_STANDALONE_MODE
