

#
# This file contains global-level macros that are specific to Trilinos
#



#
# Macro that defines Trilinos testing support
#

MACRO(TRILINOS_SETUP_TESTING_SUPPORT)

  CONFIGURE_FILE(
    ${Trilinos_SOURCE_DIR}/cmake/ctest/CTestCustom.ctest.in
    ${Trilinos_BINARY_DIR}/CTestCustom.ctest
    )

ENDMACRO()


#
#  Function for helping set up exclude files only for the packages
#  that will not be supporting autotools.
#  Returns a list of the given file name with a path for packages
#  that are not supporting autotools anymore.
#
#  example: PACKAGE_APPLY_TO_NO_AUTOTOOLS_PACKAGES("configure.ac" list)
#    assuming that the packages epetra and teuchos are not supporting 
#    autotools anymore then the return value would be:
#    "epetra/configure.ac;teuchos/configure.ac"
#
#

FUNCTION(APPLY_TO_NO_AUTOTOOLS_PACKAGES FILE_NAME LIST_RETURN)
  SET(NON_AUTOTOOLS_PACKAGES
    /packages/amesos
    /packages/anasazi
    /packages/aztecoo
    /packages/belos
    /packages/didasko
    /packages/epetra
    /packages/epetraext
    /packages/fei
    /packages/galeri
    /packages/ifpack
    /packages/intrepid
    /packages/isorropia
    /packages/kokkos
    /packages/komplex
    /packages/meros
    /packages/ml
    /packages/moertel
    /packages/moocho
    /packages/nox
    /packages/pamgen
    /packages/phalanx
    /packages/phdmesh
    /packages/pliris
    /packages/PyTrilinos
    /packages/rtop
    /packages/rythmos
    /packages/sacado
    /packages/shards
    /packages/stratimikos
    /packages/stk
    /packages/Sundance
    /packages/teuchos
    /packages/ThreadPool
    /packages/thyra
    /packages/tpetra
    /packages/tifpack
    /packages/trilinoscouplings
    /packages/triutils
  )
  
  FOREACH(PACKAGE ${NON_AUTOTOOLS_PACKAGES})
    SET(LIST_RETURN_TMP ${LIST_RETURN_TMP} ${PACKAGE}/${FILE_NAME} ${PACKAGE}/\(.*/\)*${FILE_NAME})
  ENDFOREACH()
  
  SET(${LIST_RETURN} ${LIST_RETURN_TMP} PARENT_SCOPE)
ENDFUNCTION()

#
# Macro that defines Trilinos packaging options:
#

MACRO(TRILINOS_DEFINE_PACKAGING)
  APPLY_TO_NO_AUTOTOOLS_PACKAGES("configure.ac" CONFIGURE_AC_LIST)
  APPLY_TO_NO_AUTOTOOLS_PACKAGES("configure"    CONFIGURE_LIST)
  APPLY_TO_NO_AUTOTOOLS_PACKAGES("Makefile.am"  MAKEFILE_AM_LIST)
  APPLY_TO_NO_AUTOTOOLS_PACKAGES("Makefile.in"  MAKEFILE_AC_LIST)
  APPLY_TO_NO_AUTOTOOLS_PACKAGES(".*.m4"        M4_LIST)
  APPLY_TO_NO_AUTOTOOLS_PACKAGES("bootstrap"    BOOTSTRAP_LIST)
  APPLY_TO_NO_AUTOTOOLS_PACKAGES("config/"      CONFIG_LIST)
 
    
  SET(CPACK_SOURCE_IGNORE_FILES
    /.git/
    ".gitignore"
    classicMakefile
    /Trilinos/cmake/CMakeKitwareBacklog.txt
    /Trilinos/cmake/TODO
    /Trilinos/packages/ITAPS
    /Trilinos/packages/aristos
    /Trilinos/packages/phdmesh
    /Trilinos/packages/claps
    /Trilinos/packages/external
    /Trilinos/packages/jpetra
    /Trilinos/packages/new_package
    /Trilinos/packages/rbgen
    /Trilinos/packages/WebTrilinos
    /Trilinos/packages/cmmlib
    /Trilinos/packages/lyno
    /Trilinos/packages/stalix
    /Trilinos/packages/teko
    /Trilinos/packages/Trios
    /Trilinos/demos/FEApp
    ${CONFIGURE_AC_LIST}
    ${CONFIGURE_LIST}
    ${MAKEFILE_AM_LIST}
    ${MAKEFILE_AC_LIST}
    ${M4_LIST}
    ${BOOTSTRAP_LIST}
    ${CONFIG_LIST}
    /Trilinos/packages/configure.ac
    /Trilinos/packages/configure
    /Trilinos/packages/Makefile.am
    /Trilinos/packages/Makefile.in
    /Trilinos/configure.ac
    /Trilinos/aclocal.m4
    /Trilinos/configure
    /Trilinos/Makefile.am
    /Trilinos/Makefile.in
    /Trilinos/bootstrap
    /Trilinos/config
    /Trilinos/doc/[^b]
    /Trilinos/README_old
    /Trilinos/sampleScripts/old_autotools
    /Trilinos/sampleScripts/git-profiles
    ".*.pyc"
    /Trilinos/SIERRA/
    /Trilinos/commonTools/test/coverage
    /Trilinos/commonTools/test/harness
    /Trilinos/commonTools/test/utilities/README
    /Trilinos/commonTools/test/utilities/dependencies
    /Trilinos/commonTools/test/utilities/packages
    /Trilinos/commonTools/test/utilities/r.*
    /Trilinos/commonTools/scripts
    /Trilinos/commonTools/git
    /Trilinos/packages/PyTrilinos/Notes.txt
    /Trilinos/packages/PyTrilinos/aclocal.m4
    /Trilinos/packages/PyTrilinos/bootstrap
    /Trilinos/packages/PyTrilinos/config
    /Trilinos/packages/PyTrilinos/lib
    /Trilinos/packages/PyTrilinos/macdist
    /Trilinos/packages/PyTrilinos/shared
    /Trilinos/packages/PyTrilinos/src/PyTrilinos_config.h.in
    /Trilinos/packages/PyTrilinos/src/depend
    /Trilinos/packages/PyTrilinos/src/setup.py
    /Trilinos/packages/PyTrilinos/src-boost
    /Trilinos/packages/zoltan/test/ch_brack2_3
    /Trilinos/packages/zoltan/test/ch_bug
    /Trilinos/packages/zoltan/test/ch_degenerate
    /Trilinos/packages/zoltan/test/ch_degenerateAA
    /Trilinos/packages/zoltan/test/ch_drake
    /Trilinos/packages/zoltan/test/ch_ewgt
    /Trilinos/packages/zoltan/test/ch_grid20x19
    /Trilinos/packages/zoltan/test/ch_hammond
    /Trilinos/packages/zoltan/test/ch_hammond2
    /Trilinos/packages/zoltan/test/ch_nograph
    /Trilinos/packages/zoltan/test/ch_onedbug
    /Trilinos/packages/zoltan/test/ch_random
    /Trilinos/packages/zoltan/test/ch_serial
    /Trilinos/packages/zoltan/test/ch_slac
    /Trilinos/packages/zoltan/test/ch_vwgt
    /Trilinos/packages/zoltan/test/ch_vwgt2
    /Trilinos/packages/zoltan/test/hg_cage10
    /Trilinos/packages/zoltan/test/hg_diag500_4
    /Trilinos/packages/zoltan/test/hg_ewgt
    /Trilinos/packages/zoltan/test/hg_felix
    /Trilinos/packages/zoltan/test/hg_ibm03
    /Trilinos/packages/zoltan/test/hg_ml27
    /Trilinos/packages/zoltan/test/hg_nograph
    /Trilinos/packages/zoltan/test/hg_vwgt
    /Trilinos/packages/zoltan/test/nem_ti_20k
    /Trilinos/packages/zoltan/test/nem_ti_4k
    /Trilinos/packages/zoltan/test/misc_siefert
    /Trilinos/packages/zoltan/test/th
    /Trilinos/packages/zoltan/test/bin
    /Trilinos/packages/zoltan/doc/Zoltan_html/tu_html
    /Trilinos/packages/zoltan/src/ZoltanComponent
    /Trilinos/packages/zoltan/src/driver_old
    /Trilinos/packages/zoltan/src/fdriver_old
    /Trilinos/packages/amesos/doc/AmesosOverview
    /Trilinos/packages/amesos/doc/PARA06
    /Trilinos/packages/anasazi/doc/TOMS
    /Trilinos/packages/anasazi/doc/OrthoStudy
    /Trilinos/packages/anasazi/doc/ThyraPerf
    /Trilinos/packages/aztecoo/doc/AZ_capture_matrix_howto.txt
    /Trilinos/packages/aztecoo/doc/Aztec2.0
    /Trilinos/packages/aztecoo/doc/Aztec2.1
    /Trilinos/packages/aztecoo/doc/Managing_conditioning_howto.txt
    /Trilinos/packages/aztecoo/doc/UserGuide
    /Trilinos/packages/aztecoo/doc/azteclogo.gif
    /Trilinos/packages/aztecoo/doc/read_captured_matrix.c
    /Trilinos/packages/aztecoo/example/AztecOO_RecursiveCall
    /Trilinos/packages/aztecoo/example/Epetra_MsrMatrix_AztecOO
    /Trilinos/packages/aztecoo/example/Epetra_MsrMatrix_PowerMethod
    /Trilinos/packages/aztecoo/example/IfpackIctAztecOO
    /Trilinos/packages/aztecoo/example/IfpackAztecOO
    /Trilinos/packages/aztecoo/example/IfpackVbrAztecOO
    /Trilinos/packages/aztecoo/example/MLAztecOO
    /Trilinos/packages/aztecoo/example/azoo_iterate_hb
    /Trilinos/packages/aztecoo/example/aztec_app
    /Trilinos/packages/aztecoo/example/aztec_hb
    /Trilinos/packages/galeri/src-pfem
    /Trilinos/packages/galeri/example-pfem
    /Trilinos/packages/tpetra/doc/CodingGuidelines
    /Trilinos/packages/tpetra/doc/TpetraDesign
    /Trilinos/packages/kokkos/doc
    /Trilinos/packages/aztecoo/example/AztecOO/adapt_main.mk
    /Trilinos/packages/aztecoo/example/AztecOO/vbr_main.mk
    /Trilinos/packages/aztecoo/example/AztecOO_MatlabInput/A.dat
    /Trilinos/packages/aztecoo/example/AztecOO_MatlabInput/Ainv.dat
    /Trilinos/packages/aztecoo/src/AztecOO_string_maps.txt
    /Trilinos/packages/aztecoo/src/AztecOO_string_maps_builder.pl
    /Trilinos/packages/aztecoo/src/az_comm_.*
    /Trilinos/packages/aztecoo/src/md_timer_intel.c
    /Trilinos/packages/aztecoo/src/md_timer_ncube.c
    /Trilinos/packages/aztecoo/src/md_timer_sol.c
    /Trilinos/packages/aztecoo/src/md_timer_sp2.c
    /Trilinos/packages/aztecoo/src/md_timer_sun.c
    /Trilinos/packages/aztecoo/src/md_timer_win2000.c
    /Trilinos/packages/aztecoo/src/md_wrap_intel_c.c
    /Trilinos/packages/aztecoo/src/md_wrap_ncube_c.c
    /Trilinos/packages/aztecoo/src/md_wrap_puma_c.c
    /Trilinos/packages/aztecoo/src/md_wrap_sp2_c.c
    /Trilinos/packages/aztecoo/src/stamp-h.in
    /Trilinos/packages/aztecoo/test/scripts/daily/serial/Ex_AztecOO_UserOpUserMat
    /Trilinos/packages/common/DoxyfilePackageTemplate
    /Trilinos/packages/didasko/examples/teuchos/xml-data
    /Trilinos/packages/didasko/src/.*.eps
    /Trilinos/packages/didasko/src/.*.tex
    /Trilinos/packages/didasko/src/.*.ps
    /Trilinos/packages/didasko/src/DOEbwlogo.pdf
    /Trilinos/packages/didasko/src/Makefile
    /Trilinos/packages/didasko/src/SANDreport.cls
    /Trilinos/packages/didasko/src/Trilinos60Tutorial.pdf
    /Trilinos/packages/didasko/src/Trilinos70Tutorial.pdf
    /Trilinos/packages/didasko/src/Trilinos80Tutorial.pdf
    /Trilinos/packages/didasko/src/TrilinosTutorial_ReviewAndApproval.doc
    /Trilinos/packages/didasko/src/chapterbox.pdf
    /Trilinos/packages/didasko/src/colabarticle.cls
    /Trilinos/packages/didasko/src/snllineblk.pdf
    /Trilinos/packages/didasko/src/tutorial_biblio.bib
    /Trilinos/packages/epetra/doc
    /Trilinos/packages/epetra/example/C_wrappers
    /Trilinos/packages/epetra/example/Fortran
    /Trilinos/packages/epetra/example/ImportExport
    /Trilinos/packages/epetra/example/InverseIteration
    /Trilinos/packages/epetra/example/MapColoring
    /Trilinos/packages/epetra/example/ReducedLinearProblem
    /Trilinos/packages/epetra/example/petra_howle
    /Trilinos/packages/epetra/example/petra_nonlinear
    /Trilinos/packages/epetra/example/petra_transpose
    /Trilinos/packages/epetra/src/Epetra_FastCrsMatrix.cpp
    /Trilinos/packages/epetra/src/Epetra_FastCrsMatrix.h
    /Trilinos/packages/epetra/src/Epetra_InvOperator.cpp
    /Trilinos/packages/epetra/src/Epetra_LinearProblemRedistor.cpp
    /Trilinos/packages/epetra/src/Epetra_LinearProblemRedistor.h
    /Trilinos/packages/epetra/src/Epetra_MpiSmpComm.*
    /Trilinos/packages/epetra/src/stamp-h.in
    /Trilinos/packages/epetra/src/xxemacs
    /Trilinos/packages/epetra/test/BasicPerfTest/runSummary
    /Trilinos/packages/epetra/test/Comm/simple_mpi.cpp
    /Trilinos/packages/epetra/test/Comm/threaded_Makefile
    /Trilinos/packages/epetra/test/Comm/threaded_main.cpp
    /Trilinos/packages/epetra/test/EpetraBenchmarkTest
    /Trilinos/packages/epetra/test/LinearProblemRedistor
    /Trilinos/packages/epetra/test/Makefile.template
    /Trilinos/packages/epetra/test/Map/c_main.c
    /Trilinos/packages/epetra/test/MultiVector/Makefile.purify
    /Trilinos/packages/epetra/test/OSKI
    /Trilinos/packages/epetra/test/VbrMatrix/Suppressions.in
    /Trilinos/packages/epetra/test/Vector/Makefile.purify
    /Trilinos/packages/epetra/test/testAll.*
    /Trilinos/packages/epetraext/doc/UserGuide
    /Trilinos/packages/epetraext/doc/inout
    /Trilinos/packages/epetraext/doc/matlab.README
    /Trilinos/packages/epetraext/example/MapColoring/sample_map
    /Trilinos/packages/epetraext/example/MapColoring/sample_matrix
    /Trilinos/packages/epetraext/example/inout/build
    /Trilinos/packages/epetraext/example/model_evaluator/GLpApp/Parallel2DMeshGeneratorFormat.pdf
    /Trilinos/packages/epetraext/example/model_evaluator/GLpApp/README
    /Trilinos/packages/epetraext/example/model_evaluator/GLpApp/data/square/generate-serial-meshes-1-2.sh
    /Trilinos/packages/epetraext/example/model_evaluator/GLpApp/data/square/generate-serial-meshes-1-2.sh.out
    /Trilinos/packages/epetraext/example/model_evaluator/GLpApp/data/square/square.1.2
    /Trilinos/packages/epetraext/example/model_evaluator/GLpApp/data/square/square.1.2.000
    /Trilinos/packages/epetraext/example/model_evaluator/GLpApp/data/square/square.1.2.001
    /Trilinos/packages/epetraext/example/model_evaluator/GLpApp/data/square/square.1.2.edge
    /Trilinos/packages/epetraext/example/model_evaluator/GLpApp/data/square/square.1.2.ele
    /Trilinos/packages/epetraext/example/model_evaluator/GLpApp/data/square/square.1.2.epart.2
    /Trilinos/packages/epetraext/example/model_evaluator/GLpApp/data/square/square.1.2.node
    /Trilinos/packages/epetraext/example/model_evaluator/GLpApp/data/square/square.1.2.npart.2
    /Trilinos/packages/epetraext/example/model_evaluator/GLpApp/data/square/square.1.2.poly
    /Trilinos/packages/epetraext/example/model_evaluator/GLpApp/data/square/square.1.edge
    /Trilinos/packages/epetraext/example/model_evaluator/GLpApp/data/square/square.1.ele
    /Trilinos/packages/epetraext/example/model_evaluator/GLpApp/data/square/square.1.node
    /Trilinos/packages/epetraext/example/model_evaluator/GLpApp/data/square/square.1.poly
    /Trilinos/packages/epetraext/example/model_evaluator/GLpApp/data/square/square.2.edge
    /Trilinos/packages/epetraext/example/model_evaluator/GLpApp/data/square/square.2.ele
    /Trilinos/packages/epetraext/example/model_evaluator/GLpApp/data/square/square.2.node
    /Trilinos/packages/epetraext/example/model_evaluator/GLpApp/data/square/square.2.poly
    /Trilinos/packages/epetraext/example/model_evaluator/GLpApp/data/square/square.poly
    /Trilinos/packages/epetraext/example/model_evaluator/GLpApp/from-triangle-to-serial-input-mesh.pl
    /Trilinos/packages/epetraext/src/btf/pothen/btf_notes
    /Trilinos/packages/epetraext/src/btf/pothen/main.f
    /Trilinos/packages/epetraext/src/distdir
    /Trilinos/packages/epetraext/src/transform/EpetraExt_Dirichlet_.*
    /Trilinos/packages/epetraext/src/transform/EpetraExt_StaticCondensation_LinearProblem..*
    /Trilinos/packages/epetraext/src/transform/EpetraExt_SubCopy_CrsMatrix..*
    /Trilinos/packages/epetraext/src/zoltan/EpetraExt_ZoltanMpi.*
    /Trilinos/packages/epetraext/test/Copy
    /Trilinos/packages/epetraext/test/Makefile.template
    /Trilinos/packages/epetraext/test/Zoltan/Dummy
    /Trilinos/packages/epetraext/test/inout/build
    /Trilinos/packages/epetraext/test/testAll.*
    /Trilinos/packages/triutils/src/stamp-h.in
    /stamp-h.in
    /Trilinos/packages/galeri/doc/AdvDiffSquare.png
    /Trilinos/packages/galeri/doc/L.*.png
    /Trilinos/packages/galeri/example-fem/TwoSquares.cpp
    /Trilinos/packages/galeri/src-fem/Galeri_FileGrid.h
    /Trilinos/packages/ifpack/doc/UsersGuide
    /Trilinos/packages/ifpack/example/Ifpack_ex_ScalarLaplacian_FEM.cpp
    /Trilinos/packages/ifpack/example/Ifpack_ex_VectorLaplacian_FEM.cpp
    /Trilinos/packages/ifpack/example/ifpack_hb
    /Trilinos/packages/ifpack/example/ifpack_threaded_hb
    /Trilinos/packages/ifpack/src/Ifpack_CrsGraph.h
    /Trilinos/packages/ifpack/src/Ifpack_CrsIlut.cpp
    /Trilinos/packages/ifpack/src/Ifpack_CrsIlut.h
    /Trilinos/packages/ifpack/src/Ifpack_CrsRick.cpp
    /Trilinos/packages/ifpack/src/Ifpack_CrsRick.h
    /Trilinos/packages/ifpack/src/Ifpack_HashTable.cpp
    /Trilinos/packages/ifpack/src/Ifpack_OverlapFactor.*
    /Trilinos/packages/ifpack/src/Ifpack_OverlapSolveObject..*
    /Trilinos/packages/ifpack/src/Ifpack_PerturbedMatrix.h
    /Trilinos/packages/ifpack/src/az_ifpack.*
    /Trilinos/packages/ifpack/src/ifp_Block
    /Trilinos/packages/ifpack/src/ifp_DenseMat.*
    /Trilinos/packages/ifpack/src/ifp_GlobalPrecon.h
    /Trilinos/packages/ifpack/src/ifp_Local.*
    /Trilinos/packages/ifpack/src/ifp_Matrix.h
    /Trilinos/packages/ifpack/src/ifp_Precon..*
    /Trilinos/packages/ifpack/src/ifp_SparseUtil..*
    /Trilinos/packages/ifpack/src/ifp_arch.h
    /Trilinos/packages/ifpack/src/ifp_b.*
    /Trilinos/packages/ifpack/src/ifp_c_wrappers..*
    /Trilinos/packages/ifpack/src/ifp_ifpack.h
    /Trilinos/packages/ifpack/src/ifp_lapackd.h
    /Trilinos/packages/ifpack/src/ifp_sp.*
    /Trilinos/packages/ifpack/src/old.Makefile
    /Trilinos/packages/ifpack/src/stamp-h.in
    /Trilinos/packages/ifpack/src/xxemacs
    /Trilinos/packages/ifpack/test/PointPreconditioner
    /Trilinos/packages/ifpack/test/scripts
    /Trilinos/packages/ifpack/test/scripts/run-tests.sh
    /Trilinos/packages/komplex/doc/Komplex.*.vsd
    /Trilinos/packages/komplex/doc/header.tex
    /Trilinos/packages/komplex/doc/komplex.eps
    /Trilinos/packages/komplex/doc/komplex.gif
    /Trilinos/packages/komplex/doc/komplex_user_guide.ps
    /Trilinos/packages/komplex/example/komplex_hb/README
    /Trilinos/packages/komplex/example/komplex_hb/blassm.f
    /Trilinos/packages/komplex/example/komplex_hb/create_vbr.c
    /Trilinos/packages/komplex/example/komplex_hb/distrib_.*_matrix.c
    /Trilinos/packages/komplex/example/komplex_hb/formats.f
    /Trilinos/packages/komplex/example/komplex_hb/iohb.*
    /Trilinos/packages/komplex/example/komplex_hb/main.c
    /Trilinos/packages/komplex/example/komplex_hb/prototypes.h
    /Trilinos/packages/komplex/example/komplex_hb/read_.*
    /Trilinos/packages/komplex/example/komplex_hb/sc.*
    /Trilinos/packages/komplex/example/komplex_hb/smsrres.c
    /Trilinos/packages/komplex/example/komplex_hb/svbrres.c
    /Trilinos/packages/komplex/example/komplex_hb/unary.f
    /Trilinos/packages/komplex/example/komplex_hb/write_vec.c
    /Trilinos/packages/komplex/src/new
    /Trilinos/packages/komplex/src/stamp-h.in
    /Trilinos/packages/komplex/test/definition
    /Trilinos/packages/amesos/example/RunParaklete.cpp
    /Trilinos/packages/amesos/example/Thyra_AmesosLinearOpWithSolveFactory.cpp
    /Trilinos/packages/amesos/example/pk.h
    /Trilinos/packages/amesos/example/run_pk.c
    /Trilinos/packages/amesos/example/simpleStratimikosSolve.cpp
    /Trilinos/packages/amesos/example/simpleStratimikosSolve.hpp
    /Trilinos/packages/amesos/example/stratimikos_example.cpp
    /Trilinos/packages/amesos/src/Amesos_BTF.h
    /Trilinos/packages/amesos/src/Amesos_Component.h
    /Trilinos/packages/amesos/src/Amesos_Merikos.h
    /Trilinos/packages/amesos/src/Amesos_BTF.h
    /Trilinos/packages/amesos/src/SuiteSparse/CAMD/Source/amesos_camd_1.c
    /Trilinos/packages/amesos/src/SuiteSparse/CAMD/Source/amesos_camd_aat.c
    /Trilinos/packages/amesos/src/SuiteSparse/CAMD/Source/amesos_camd_control.c
    /Trilinos/packages/amesos/src/SuiteSparse/CAMD/Source/amesos_camd_defaults.c
    /Trilinos/packages/amesos/src/SuiteSparse/CAMD/Source/amesos_camd_dump.c
    /Trilinos/packages/amesos/src/SuiteSparse/CAMD/Source/amesos_camd_info.c
    /Trilinos/packages/amesos/src/SuiteSparse/CAMD/Source/amesos_camd_order.c
    /Trilinos/packages/amesos/src/SuiteSparse/CAMD/Source/amesos_camd_preprocess.c
    /Trilinos/packages/amesos/src/SuiteSparse/CAMD/Source/amesos_camd_valid.c
    /Trilinos/packages/amesos/src/SuiteSparse/CCOLAMD/Source/amesos_ccolamd.c
    /Trilinos/packages/amesos/src/SuiteSparse/CHOLMOD/Cholesky/amesos_cholmod_amd.c
    /Trilinos/packages/amesos/src/SuiteSparse/CHOLMOD/Cholesky/amesos_cholmod_analyze.c
    /Trilinos/packages/amesos/src/SuiteSparse/CHOLMOD/Cholesky/amesos_cholmod_colamd.c
    /Trilinos/packages/amesos/src/SuiteSparse/CHOLMOD/Cholesky/amesos_cholmod_etree.c
    /Trilinos/packages/amesos/src/SuiteSparse/CHOLMOD/Cholesky/amesos_cholmod_factorize.c
    /Trilinos/packages/amesos/src/SuiteSparse/CHOLMOD/Cholesky/amesos_cholmod_postorder.c
    /Trilinos/packages/amesos/src/SuiteSparse/CHOLMOD/Cholesky/amesos_cholmod_rcond.c
    /Trilinos/packages/amesos/src/SuiteSparse/CHOLMOD/Cholesky/amesos_cholmod_resymbol.c
    /Trilinos/packages/amesos/src/SuiteSparse/CHOLMOD/Cholesky/amesos_cholmod_rowcolcounts.c
    /Trilinos/packages/amesos/src/SuiteSparse/CHOLMOD/Cholesky/amesos_cholmod_rowfac.c
    /Trilinos/packages/amesos/src/SuiteSparse/CHOLMOD/Cholesky/amesos_cholmod_solve.c
    /Trilinos/packages/amesos/src/SuiteSparse/CHOLMOD/Cholesky/amesos_cholmod_spsolve.c
    /Trilinos/packages/amesos/src/SuiteSparse/CHOLMOD/Cholesky/amesos_t_cholmod_rowfac.c
    /Trilinos/packages/amesos/src/SuiteSparse/CHOLMOD/Core/amesos_cholmod_aat.c
    /Trilinos/packages/amesos/src/SuiteSparse/CHOLMOD/Core/amesos_cholmod_add.c
    /Trilinos/packages/amesos/src/SuiteSparse/CHOLMOD/Core/amesos_cholmod_band.c
    /Trilinos/packages/amesos/src/SuiteSparse/CHOLMOD/Core/amesos_cholmod_change_factor.c
    /Trilinos/packages/amesos/src/SuiteSparse/CHOLMOD/Core/amesos_cholmod_copy.c
    /Trilinos/packages/amesos/src/SuiteSparse/CHOLMOD/Core/amesos_cholmod_dense.c
    /Trilinos/packages/amesos/src/SuiteSparse/CHOLMOD/Core/amesos_cholmod_factor.c
    /Trilinos/packages/amesos/src/SuiteSparse/CHOLMOD/Partition/amesos_cholmod_camd.c
    /Trilinos/packages/amesos/src/SuiteSparse/CHOLMOD/Partition/amesos_cholmod_ccolamd.c
    /Trilinos/packages/amesos/src/SuiteSparse/CHOLMOD/Partition/amesos_cholmod_csymamd.c
    /Trilinos/packages/amesos/src/SuiteSparse/CHOLMOD/Partition/amesos_cholmod_metis.c
    /Trilinos/packages/amesos/src/SuiteSparse/CHOLMOD/Partition/amesos_cholmod_nesdis.c
    /Trilinos/packages/amesos/src/src-repository
    /Trilinos/packages/amesos/src/stamp-h.in
    /Trilinos/packages/amesos/test/TestOptions/Dummy
    /Trilinos/packages/amesos/test/Test_Basic/NotQuiteDense.triU
    /Trilinos/packages/amesos/test/Test_Performance/In_Test_UmfpackPerformance.csh
    /Trilinos/packages/amesos/test/scripts/daily/mpi/TestBasic.sh
    /Trilinos/packages/amesos/test/scripts/daily/serial/TestAmesos.sh
    /Trilinos/packages/pliris/doc/matrix_.*.gif
    /Trilinos/packages/pliris/src/Make..*
    /Trilinos/packages/pliris/src/clean_code.h
    /Trilinos/packages/pliris/src/init..*
    /Trilinos/packages/pliris/src/malloc.c
    /Trilinos/packages/pliris/src/my_srand48.c
    /Trilinos/packages/meros/README-MEROS
    /Trilinos/packages/meros/doc/UsersGuide
    /Trilinos/packages/meros/example/data/mac
    /Trilinos/packages/meros/example/data/mac-vbr
    /Trilinos/packages/meros/example/data/q1p0
    /Trilinos/packages/meros/example/data/salsa
    /Trilinos/packages/meros/example/data/tmac-vbr
    /Trilinos/package/phdmesh/Make.in
    /Trilinos/packages/phdmesh/\(.*/\)*Make.in
    /Trilinos/packages/phdmesh/build_examples
    /Trilinos/packages/teuchos/config.h.in
    /Trilinos/packages/teuchos/doc/images
    /Trilinos/packages/teuchos/example/config.h.in
    /Trilinos/packages/nox/src-loca/python
    /Trilinos/packages/nox/test/lapack/LOCA_python
    /Trilinos/packages/anasazi/src/ModalAnalysisSolvers
    /Trilinos/packages/ml/util
    /Trilinos/packages/ml/etc
    /Trilinos/packages/ml/test/tmp
    /Trilinos/packages/ml/doc/UsersGuide
    /Trilinos/packages/ml/doc/DevelopersGuide
    /Trilinos/packages/ml/doc/MLAPI
    /Trilinos/packages/ml/python
    /Trilinos/packages/ml/doc/DoxyfileWeb
    /Trilinos/packages/ml/doc/build_docs
    /Trilinos/packages/ml/doc/ml-logo.eps
    /Trilinos/packages/ml/doc/ml-logo.jpg
    /Trilinos/packages/ml/doc/sc2000.ps.gz
    /Trilinos/packages/ml/examples/Makefile-common.include
    /Trilinos/packages/ml/examples/Maxwell/ml_periodic_max.c
    /Trilinos/packages/ml/examples/Other/ml_read_complexmaxwell.c
    /Trilinos/packages/ml/examples/Other/ml_read_maxwell.c
    /Trilinos/packages/ml/examples/Other/ml_star2d.c
    /Trilinos/packages/ml/examples/Other/new_readex.c
    /Trilinos/packages/ml/examples/Other/oldml_readex.c
    /Trilinos/packages/ml/examples/Other/seg_readex.c
    /Trilinos/packages/ml/examples/README.AddingExamples
    /Trilinos/packages/ml/examples/RefMaxwell
    /Trilinos/packages/ml/examples/RefMaxwell/rpc.cpp
    /Trilinos/packages/ml/src/Coarsen/README
    /Trilinos/packages/ml/src/Main/ml_v_cycle.c
    /Trilinos/packages/ml/src/Smoother/README
    /Trilinos/packages/ml/src/Utils/jmpilib.c
    /Trilinos/packages/ml/src/Utils/jostle.h
    /Trilinos/packages/ml/src/Utils/ml_vampir.c
    /Trilinos/packages/ml/src/Utils/tumi.c
    /Trilinos/packages/ml/test/README.runtests
    /Trilinos/packages/ml/test/Zoltan/cxx_main_simple.cpp
    /Trilinos/packages/ml/test/scripts
    /Trilinos/packages/ml/test/scripts/run-tests.sh
    /Trilinos/packages/sacado/example/FEApp/experimental
    /Trilinos/packages/tifpack/src/Tifpack_AMDReordering.cpp
    /Trilinos/packages/tifpack/src/Tifpack_AMDReordering.hpp
    /Trilinos/packages/tifpack/src/Tifpack_AdditiveSchwarz.hpp
    /Trilinos/packages/tifpack/src/Tifpack_Amesos.cpp
    /Trilinos/packages/tifpack/src/Tifpack_Amesos.hpp
    /Trilinos/packages/tifpack/src/Tifpack_BlockRelaxation.hpp
    /Trilinos/packages/tifpack/src/Tifpack_ConstructLevelFillGraph.hpp
    /Trilinos/packages/tifpack/src/Tifpack_Container.hpp
    /Trilinos/packages/tifpack/src/Tifpack_CrsGraph.hpp
    /Trilinos/packages/tifpack/src/Tifpack_DenseContainer.cpp
    /Trilinos/packages/tifpack/src/Tifpack_DenseContainer.hpp
    /Trilinos/packages/tifpack/src/Tifpack_DiagPreconditioner.cpp
    /Trilinos/packages/tifpack/src/Tifpack_DiagPreconditioner.hpp
    /Trilinos/packages/tifpack/src/Tifpack_DiagonalFilter.cpp
    /Trilinos/packages/tifpack/src/Tifpack_DiagonalFilter.hpp
    /Trilinos/packages/tifpack/src/Tifpack_DropFilter.cpp
    /Trilinos/packages/tifpack/src/Tifpack_DropFilter.hpp
    /Trilinos/packages/tifpack/src/Tifpack_EquationPartitioner.cpp
    /Trilinos/packages/tifpack/src/Tifpack_EquationPartitioner.hpp
    /Trilinos/packages/tifpack/src/Tifpack_Graph.hpp
    /Trilinos/packages/tifpack/src/Tifpack_Graph_Tpetra_CrsGraph.cpp
    /Trilinos/packages/tifpack/src/Tifpack_Graph_Tpetra_CrsGraph.hpp
    /Trilinos/packages/tifpack/src/Tifpack_Graph_Tpetra_RowMatrix.cpp
    /Trilinos/packages/tifpack/src/Tifpack_Graph_Tpetra_RowMatrix.hpp
    /Trilinos/packages/tifpack/src/Tifpack_GreedyPartitioner.cpp
    /Trilinos/packages/tifpack/src/Tifpack_GreedyPartitioner.hpp
    /Trilinos/packages/tifpack/src/Tifpack_HashTable.cpp
    /Trilinos/packages/tifpack/src/Tifpack_HashTable.hpp
    /Trilinos/packages/tifpack/src/Tifpack_IC.cpp
    /Trilinos/packages/tifpack/src/Tifpack_IC.hpp
    /Trilinos/packages/tifpack/src/Tifpack_ICT.cpp
    /Trilinos/packages/tifpack/src/Tifpack_ICT.hpp
    /Trilinos/packages/tifpack/src/Tifpack_IC_Utils.cpp
    /Trilinos/packages/tifpack/src/Tifpack_IC_Utils.hpp
    /Trilinos/packages/tifpack/src/Tifpack_IKLU.cpp
    /Trilinos/packages/tifpack/src/Tifpack_IKLU.hpp
    /Trilinos/packages/tifpack/src/Tifpack_IKLU_Utils.cpp
    /Trilinos/packages/tifpack/src/Tifpack_IKLU_Utils.hpp
    /Trilinos/packages/tifpack/src/Tifpack_ILU.cpp
    /Trilinos/packages/tifpack/src/Tifpack_ILU.hpp
    /Trilinos/packages/tifpack/src/Tifpack_LinearPartitioner.cpp
    /Trilinos/packages/tifpack/src/Tifpack_LinearPartitioner.hpp
    /Trilinos/packages/tifpack/src/Tifpack_LocalFilter.cpp
    /Trilinos/packages/tifpack/src/Tifpack_LocalFilter.hpp
    /Trilinos/packages/tifpack/src/Tifpack_METISPartitioner.cpp
    /Trilinos/packages/tifpack/src/Tifpack_METISPartitioner.hpp
    /Trilinos/packages/tifpack/src/Tifpack_METISReordering.cpp
    /Trilinos/packages/tifpack/src/Tifpack_METISReordering.hpp
    /Trilinos/packages/tifpack/src/Tifpack_NodeFilter.cpp
    /Trilinos/packages/tifpack/src/Tifpack_NodeFilter.hpp
    /Trilinos/packages/tifpack/src/Tifpack_OverlapFactor.cpp
    /Trilinos/packages/tifpack/src/Tifpack_OverlapFactorObject.hpp
    /Trilinos/packages/tifpack/src/Tifpack_OverlapSolveObject.cpp
    /Trilinos/packages/tifpack/src/Tifpack_OverlapSolveObject.hpp
    /Trilinos/packages/tifpack/src/Tifpack_OverlappingPartitioner.cpp
    /Trilinos/packages/tifpack/src/Tifpack_OverlappingPartitioner.hpp
    /Trilinos/packages/tifpack/src/Tifpack_OverlappingRowMatrix.cpp
    /Trilinos/packages/tifpack/src/Tifpack_OverlappingRowMatrix.hpp
    /Trilinos/packages/tifpack/src/Tifpack_Partitioner.hpp
    /Trilinos/packages/tifpack/src/Tifpack_PerturbedMatrix.hpp
    /Trilinos/packages/tifpack/src/Tifpack_RCMReordering.cpp
    /Trilinos/packages/tifpack/src/Tifpack_RCMReordering.hpp
    /Trilinos/packages/tifpack/src/Tifpack_ReorderFilter.cpp
    /Trilinos/packages/tifpack/src/Tifpack_ReorderFilter.hpp
    /Trilinos/packages/tifpack/src/Tifpack_Reordering.hpp
    /Trilinos/packages/tifpack/src/Tifpack_SPARSKIT.cpp
    /Trilinos/packages/tifpack/src/Tifpack_SPARSKIT.hpp
    /Trilinos/packages/tifpack/src/Tifpack_SingletonFilter.cpp
    /Trilinos/packages/tifpack/src/Tifpack_SingletonFilter.hpp
    /Trilinos/packages/tifpack/src/Tifpack_SparseContainer.hpp
    /Trilinos/packages/tifpack/src/Tifpack_SparsityFilter.cpp
    /Trilinos/packages/tifpack/src/Tifpack_SparsityFilter.hpp
    /Trilinos/packages/tifpack/src/Tifpack_UserPartitioner.cpp
    /Trilinos/packages/tifpack/src/Tifpack_UserPartitioner.hpp
    /Trilinos/packages/tifpack/src/Tifpack_Utils.cpp
    /Trilinos/packages/tifpack/src/Tifpack_Utils.hpp
  )
  
  IF(${PROJECT_NAME}_VERBOSE_CONFIGURE)
    MESSAGE("Exclude files when building source packages")
    FOREACH(item ${CPACK_SOURCE_IGNORE_FILES})
      MESSAGE(${item})
    ENDFOREACH()
  ENDIF()
  

  SET(CPACK_PACKAGE_DESCRIPTION "Trilinos provides algorithms and technologies for the solution of large-scale, complex multi-physics engineering and scientific problems.")
  SET(CPACK_PACKAGE_FILE_NAME "trilinos-setup-${Trilinos_VERSION}")
  SET(CPACK_PACKAGE_INSTALL_DIRECTORY "Trilinos ${Trilinos_VERSION}")
  SET(CPACK_PACKAGE_REGISTRY_KEY "Trilinos ${Trilinos_VERSION}")
  SET(CPACK_PACKAGE_NAME "trilinos")
  SET(CPACK_PACKAGE_VENDOR "Sandia National Laboratories")
  SET(CPACK_PACKAGE_VERSION "${Trilinos_VERSION}")
  SET(CPACK_RESOURCE_FILE_README "${Trilinos_SOURCE_DIR}/README")
  SET(CPACK_RESOURCE_FILE_LICENSE "${Trilinos_SOURCE_DIR}/README")
  SET(CPACK_SOURCE_GENERATOR "TGZ;TBZ2")
  SET(CPACK_SOURCE_FILE_NAME "trilinos-source-${Trilinos_VERSION}")
  SET(CPACK_COMPONENTS_ALL ${Trilinos_PACKAGES})
  
  PACKAGE_ARCH_GET_ENABLED_LIST( Trilinos_PACKAGES Trilinos ON
    FALSE ENABLED_PACKAGES NUM_ENABLED)
  string(REPLACE " " ";" ENABLED_PACKAGES "${ENABLED_PACKAGES}")
  
  #message("ENABLED PACKAGES: ${ENABLED_PACKAGES} ${NUM_ENABLED}")
  FOREACH(PKG ${ENABLED_PACKAGES})
    IF(NOT "${${PKG}_LIB_REQUIRED_DEP_PACKAGES}" STREQUAL "")
        string(TOUPPER ${PKG} UPPER_PKG)
        #message("${UPPER_PKG} depends on : ${${PKG}_LIB_REQUIRED_DEP_PACKAGES}")
        SET(CPACK_COMPONENT_${UPPER_PKG}_DEPENDS ${${PKG}_LIB_REQUIRED_DEP_PACKAGES})
    ENDIF()
    #message("${PKG} depends on : ${${PKG}_LIB_REQUIRED_DEP_PACKAGES}")
  ENDFOREACH()

  
  IF(WIN32)
    SET(CPACK_GENERATOR "NSIS")
    SET(CPACK_NSIS_MODIFY_PATH OFF)
  ENDIF()
  
  INCLUDE(CPack)

ENDMACRO()
