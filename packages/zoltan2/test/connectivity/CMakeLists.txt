#TRIBITS_ADD_EXECUTABLE(
#   iceSheet
#   SOURCES ${HEADERS} iceSheetMain.cpp 
#   COMM serial mpi
#)

TRIBITS_ADD_EXECUTABLE(
  iceSheet_int
  SOURCES ${HEADERS} iceSheetIntegration.cpp
  COMM serial mpi
)

TRIBITS_ADD_TEST(
  iceSheet_int
  NAME iceSheet_Antarctica_16KM
  ARGS " ./ ant16km ant16km-basal-friction.ascii ant16km-borders.quad.msh ant16km-s1-answers.txt"
  NUM_MPI_PROCS 4
  COMM mpi
  PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAIL"
)

TRIBITS_ADD_TEST(
  iceSheet_int
  NAME iceSheet_Antarctica_8KM
  ARGS " ./ ant8km ant8km-basal-friction.ascii ant8km-borders.quad.msh ant8km-s1-answers.txt"
  NUM_MPI_PROCS 4
  COMM mpi
  PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAIL"  
)

TRIBITS_ADD_TEST(
  iceSheet_int
  NAME iceSheet_Antarctica_4KM
  ARGS " ./ ant4km ant4km-basal-friction.ascii ant4km-borders.quad.msh ant4km-s1-answers.txt"
  NUM_MPI_PROCS 4
  COMM mpi
  PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAIL"
)

TRIBITS_ADD_TEST(
  iceSheet_int
  NAME iceSheet_BigRectangle
  ARGS " ./ bigrectangle bigrectangle-basal-friction.ascii bigrectangle-borders.quad.msh bigrectangle-s1-answers.txt"
  NUM_MPI_PROCS 4
  COMM mpi
  PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAIL"  
)

TRIBITS_ADD_TEST(
  iceSheet_int
  NAME iceSheet_ChainedDoubleHinge
  ARGS " ./ chaineddoublehinge chaineddoublehinge-basal-friction.ascii chaineddoublehinge-borders.quad.msh chaineddoublehinge-answer.txt"
  NUM_MPI_PROCS 4
  COMM mpi
  PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAIL"  
)

TRIBITS_ADD_TEST(
  iceSheet_int
  NAME iceSheet_SimpleHinge
  ARGS " ./ test1 test1-basal-friction.ascii test1-borders.quad.msh test1-s1-answers.txt"
  NUM_MPI_PROCS 4
  COMM mpi
  PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAIL"  
)

TRIBITS_ADD_TEST(
  iceSheet_int
  NAME iceSheet_TicTacToe_KeepTwo
  ARGS " ./ tictactoe tictactoe-basal-friction-two.ascii tictactoe-borders.quad.msh tictactoe-s1-two-answers.txt"
  NUM_MPI_PROCS 4
  COMM mpi
  PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAIL"  
)


TRIBITS_ADD_TEST(
  iceSheet_int
  NAME iceSheet_TicTacToe_KeepAll
  ARGS " ./ tictactoe tictactoe-basal-friction-all.ascii tictactoe-borders.quad.msh tictactoe-s1-all-answers.txt"
  NUM_MPI_PROCS 4
  COMM mpi
  PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAIL"  
)

TRIBITS_ADD_TEST(
  iceSheet_int
  NAME icesheet_Three_OneProcessGetsNone
  ARGS "./ three three-basal-friction.ascii three-borders.quad.msh three-answers.txt"
  NUM_MPI_PROCS 4
  COMM mpi
  PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAIL"
)

TRIBITS_ADD_TEST(
  iceSheet_int
  NAME icesheet_FourVtxNoEdges
  ARGS "./ noedges noedges-basal-friction.ascii noedges-borders.quad.msh noedges-answers.txt"
  NUM_MPI_PROCS 4
  COMM mpi
  PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAIL"
)

TRIBITS_ADD_TEST(
  iceSheet_int
  NAME iceSheet_ManyHinges
  ARGS " ./ worstcase worstcase-basal-friction.ascii worstcase-borders.quad.msh worstcase-s1-answers.txt"
  NUM_MPI_PROCS 4
  COMM mpi
  PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAIL" 
)

IF (${PACKAGE_NAME}_ENABLE_ParMETIS)
  
TRIBITS_ADD_TEST(
  iceSheet_int
  NAME iceSheet_Antarctica_16KM_ParMETIS
  ARGS " ./ ant16km ant16km-basal-friction.ascii ant16km-borders.quad.msh ant16km-s1-answers.txt prepartition"
  NUM_MPI_PROCS 4
  COMM mpi
  PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAIL"
)

TRIBITS_ADD_TEST(
  iceSheet_int
  NAME iceSheet_Antarctica_8KM_ParMETIS
  ARGS " ./ ant8km ant8km-basal-friction.ascii ant8km-borders.quad.msh ant8km-s1-answers.txt prepartition"
  NUM_MPI_PROCS 4
  COMM mpi
  PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAIL"  
)

TRIBITS_ADD_TEST(
  iceSheet_int
  NAME iceSheet_Antarctica_4KM_ParMETIS
  ARGS " ./ ant4km ant4km-basal-friction.ascii ant4km-borders.quad.msh ant4km-s1-answers.txt prepartition"
  NUM_MPI_PROCS 4
  COMM mpi
  PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAIL"
)

TRIBITS_ADD_TEST(
  iceSheet_int
  NAME iceSheet_BigRectangle_ParMETIS
  ARGS " ./ bigrectangle bigrectangle-basal-friction.ascii bigrectangle-borders.quad.msh bigrectangle-s1-answers.txt prepartition"
  NUM_MPI_PROCS 4
  COMM mpi
  PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAIL"  
)

TRIBITS_ADD_TEST(
  iceSheet_int
  NAME iceSheet_ChainedDoubleHinge_ParMETIS
  ARGS " ./ chaineddoublehinge chaineddoublehinge-basal-friction.ascii chaineddoublehinge-borders.quad.msh chaineddoublehinge-answer.txt prepartition"
  NUM_MPI_PROCS 4
  COMM mpi
  PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAIL"  
)

TRIBITS_ADD_TEST(
  iceSheet_int
  NAME iceSheet_SimpleHinge_ParMETIS
  ARGS " ./ test1 test1-basal-friction.ascii test1-borders.quad.msh test1-s1-answers.txt prepartition"
  NUM_MPI_PROCS 4
  COMM mpi
  PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAIL"  
)

TRIBITS_ADD_TEST(
  iceSheet_int
  NAME iceSheet_TicTacToe_KeepTwo_ParMETIS
  ARGS " ./ tictactoe tictactoe-basal-friction-two.ascii tictactoe-borders.quad.msh tictactoe-s1-two-answers.txt prepartition"
  NUM_MPI_PROCS 4
  COMM mpi
  PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAIL"  
)


TRIBITS_ADD_TEST(
  iceSheet_int
  NAME iceSheet_TicTacToe_KeepAll_ParMETIS
  ARGS " ./ tictactoe tictactoe-basal-friction-all.ascii tictactoe-borders.quad.msh tictactoe-s1-all-answers.txt prepartition"
  NUM_MPI_PROCS 4
  COMM mpi
  PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAIL"  
)


TRIBITS_ADD_TEST(
  iceSheet_int
  NAME iceSheet_ManyHinges_ParMETIS
  ARGS " ./ worstcase worstcase-basal-friction.ascii worstcase-borders.quad.msh worstcase-s1-answers.txt prepartition"
  NUM_MPI_PROCS 4
  COMM mpi
  PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAIL" 
)
ENDIF()


TRIBITS_COPY_FILES_TO_BINARY_DIR(copy_files_for_connectivity_tests
  SOURCE_FILES
    test1.quad.msh
    test1.mtx
    test1-borders.quad.msh
    test1-basal-friction.ascii
    test1-s4-answers.txt
    test1-s3-answers.txt
    test1-s2-answers.txt
    test1-s1-answers.txt
    ant16km.quad.msh
    ant16km.mtx
    ant16km-borders.quad.msh
    ant16km-basal-friction.ascii
    ant16km-s4-answers.txt
    ant16km-s3-answers.txt
    ant16km-s2-answers.txt
    ant16km-s1-answers.txt
    ant8km.quad.msh
    ant8km.mtx
    ant8km-borders.quad.msh
    ant8km-basal-friction.ascii
    ant8km-s4-answers.txt
    ant8km-s3-answers.txt
    ant8km-s2-answers.txt
    ant8km-s1-answers.txt
    ant4km.quad.msh
    ant4km.mtx
    ant4km-borders.quad.msh
    ant4km-basal-friction.ascii
    ant4km-s4-answers.txt
    ant4km-s3-answers.txt
    ant4km-s2-answers.txt
    ant4km-s1-answers.txt
    #ant2km.quad.msh
    ant2km-borders.quad.msh
    ant2km-basal-friction.ascii
    ant2km-s4-answers.txt
    ant2km-s3-answers.txt
    ant2km-s2-answers.txt
    ant2km-s1-answers.txt
    #ant1km.quad.msh
    ant1km-borders.quad.msh
    ant1km-basal-friction.ascii
    ant1km-s4-answers.txt
    ant1km-s3-answers.txt
    ant1km-s2-answers.txt
    ant1km-s1-answers.txt
    bigrectangle.quad.msh
    bigrectangle.mtx
    bigrectangle-borders.quad.msh
    bigrectangle-basal-friction.ascii
    bigrectangle-s4-answers.txt
    bigrectangle-s3-answers.txt
    bigrectangle-s2-answers.txt
    bigrectangle-s1-answers.txt
    chaineddoublehinge.quad.msh
    chaineddoublehinge.mtx
    chaineddoublehinge-borders.quad.msh
    chaineddoublehinge-basal-friction.ascii
    chaineddoublehinge-answer.txt
    worstcase.quad.msh
    worstcase.mtx
    worstcase-borders.quad.msh
    worstcase-basal-friction.ascii
    worstcase-s4-answers.txt
    worstcase-s3-answers.txt
    worstcase-s2-answers.txt
    worstcase-s1-answers.txt
    tictactoe.quad.msh
    tictactoe.mtx
    tictactoe-borders.quad.msh
    tictactoe-basal-friction-all.ascii
    tictactoe-basal-friction-two.ascii
    tictactoe-s4-two-answers.txt
    tictactoe-s3-two-answers.txt
    tictactoe-s2-two-answers.txt
    tictactoe-s1-two-answers.txt
    tictactoe-s4-all-answers.txt
    tictactoe-s3-all-answers.txt
    tictactoe-s2-all-answers.txt
    tictactoe-s1-all-answers.txt
    noedges.mtx
    noedges-basal-friction.ascii
    noedges-borders.quad.msh
    noedges-answers.txt
    three.mtx
    three-basal-friction.ascii
    three-borders.quad.msh
    three-answers.txt
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/meshData
  DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
)
