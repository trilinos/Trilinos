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
#TRIBITS_ADD_TEST(
#    iceSheet
#    NAME iceSheet_Antarctica_1K
#    ARGS "meshfilename.txt" 
#    PASS_REGULAR_EXPRESSION "PASS"
#    FAIL_REGULAR_EXPRESSION "FAIL"
#)
#TRIBITS_ADD_TEST(
#    iceSheet
#    NAME iceSheet_TicTacToe_KeepAll_sensitivity4
#    ARGS " tictactoe.quad.msh tictactoe-borders.quad.msh tictactoe-basal-friction-all.ascii tictactoe-s4-all-answers.txt 4 i" 
#    NUM_MPI_PROCS 4
#    COMM mpi
#    PASS_REGULAR_EXPRESSION "PASS"
#    FAIL_REGULAR_EXPRESSION "FAIL"
#)
#TRIBITS_ADD_TEST(
#    iceSheet
#    NAME iceSheet_TicTacToe_KeepAll_sensitivity3
#    ARGS " tictactoe.quad.msh tictactoe-borders.quad.msh tictactoe-basal-friction-all.ascii tictactoe-s3-all-answers.txt 3 i" 
#    NUM_MPI_PROCS 4
#    COMM mpi
#    PASS_REGULAR_EXPRESSION "PASS"
#    FAIL_REGULAR_EXPRESSION "FAIL"
#)
#TRIBITS_ADD_TEST(
#    iceSheet
#    NAME iceSheet_TicTacToe_KeepAll_sensitivity2
#    ARGS " tictactoe.quad.msh tictactoe-borders.quad.msh tictactoe-basal-friction-all.ascii tictactoe-s2-all-answers.txt 2 i" 
#    NUM_MPI_PROCS 4
#    COMM mpi
#    PASS_REGULAR_EXPRESSION "PASS"
#    FAIL_REGULAR_EXPRESSION "FAIL"
#)
#TRIBITS_ADD_TEST(
#    iceSheet
#    NAME iceSheet_TicTacToe_KeepAll_sensitivity1
#    ARGS " tictactoe.quad.msh tictactoe-borders.quad.msh tictactoe-basal-friction-all.ascii tictactoe-s1-all-answers.txt 1 i" 
#    NUM_MPI_PROCS 4
#    COMM mpi
#    PASS_REGULAR_EXPRESSION "PASS"
#    FAIL_REGULAR_EXPRESSION "FAIL"
#)
#TRIBITS_ADD_TEST(
#    iceSheet
#    NAME iceSheet_TicTacToe_KeepTwo_sensitivity4
#    ARGS " tictactoe.quad.msh tictactoe-borders.quad.msh tictactoe-basal-friction-two.ascii tictactoe-s4-two-answers.txt 4 i" 
#    NUM_MPI_PROCS 4
#    COMM mpi
#    PASS_REGULAR_EXPRESSION "PASS"
#    FAIL_REGULAR_EXPRESSION "FAIL"
#)
#TRIBITS_ADD_TEST(
#    iceSheet
#    NAME iceSheet_TicTacToe_KeepTwo_sensitivity3
#    ARGS " tictactoe.quad.msh tictactoe-borders.quad.msh tictactoe-basal-friction-two.ascii tictactoe-s3-two-answers.txt 3 i" 
#    NUM_MPI_PROCS 4
#    COMM mpi
#    PASS_REGULAR_EXPRESSION "PASS"
#    FAIL_REGULAR_EXPRESSION "FAIL"
#)
#TRIBITS_ADD_TEST(
#    iceSheet
#    NAME iceSheet_TicTacToe_KeepTwo_sensitivity2
#    ARGS " tictactoe.quad.msh tictactoe-borders.quad.msh tictactoe-basal-friction-two.ascii tictactoe-s2-two-answers.txt 2 i" 
#    NUM_MPI_PROCS 4
#    COMM mpi
#    PASS_REGULAR_EXPRESSION "PASS"
#    FAIL_REGULAR_EXPRESSION "FAIL"
#)
#TRIBITS_ADD_TEST(
#    iceSheet
#    NAME iceSheet_TicTacToe_KeepTwo_sensitivity1
#    ARGS " tictactoe.quad.msh tictactoe-borders.quad.msh tictactoe-basal-friction-two.ascii tictactoe-s1-two-answers.txt 1 i" 
#    NUM_MPI_PROCS 4
#    COMM mpi
#    PASS_REGULAR_EXPRESSION "PASS"
#    FAIL_REGULAR_EXPRESSION "FAIL"
#)
#
#TRIBITS_ADD_TEST(
#    iceSheet
#    NAME iceSheet_ManyDoubleHinges_sensitivity4
#    ARGS " worstcase.quad.msh worstcase-borders.quad.msh worstcase-basal-friction.ascii worstcase-s4-answers.txt 4 i" 
#    NUM_MPI_PROCS 4
#    COMM mpi
#    PASS_REGULAR_EXPRESSION "PASS"
#    FAIL_REGULAR_EXPRESSION "FAIL"
#)
#
#TRIBITS_ADD_TEST(
#    iceSheet
#    NAME iceSheet_ManyDoubleHinges_sensitivity3
#    ARGS " worstcase.quad.msh worstcase-borders.quad.msh worstcase-basal-friction.ascii worstcase-s3-answers.txt 3 i" 
#    NUM_MPI_PROCS 4
#    COMM mpi
#    PASS_REGULAR_EXPRESSION "PASS"
#    FAIL_REGULAR_EXPRESSION "FAIL"
#)
#
#TRIBITS_ADD_TEST(
#    iceSheet
#    NAME iceSheet_ManyDoubleHinges_sensitivity2
#    ARGS " worstcase.quad.msh worstcase-borders.quad.msh worstcase-basal-friction.ascii worstcase-s2-answers.txt 2 i" 
#    NUM_MPI_PROCS 4
#    COMM mpi
#    PASS_REGULAR_EXPRESSION "PASS"
#    FAIL_REGULAR_EXPRESSION "FAIL"
#)
#TRIBITS_ADD_TEST(
#    iceSheet
#    NAME iceSheet_ManyDoubleHinges_sensitivity1
#    ARGS " worstcase.quad.msh worstcase-borders.quad.msh worstcase-basal-friction.ascii worstcase-s1-answers.txt 1 i" 
#    NUM_MPI_PROCS 4
#    COMM mpi
#    PASS_REGULAR_EXPRESSION "PASS"
#    FAIL_REGULAR_EXPRESSION "FAIL"
#)
#
#TRIBITS_ADD_TEST(
#    iceSheet
#    NAME iceSheet_ChainedDoubleHinge_sensitivity4
#    ARGS " chaineddoublehinge.quad.msh chaineddoublehinge-borders.quad.msh chaineddoublehinge-basal-friction.ascii chaineddoublehinge-answer.txt 4 i" 
#    NUM_MPI_PROCS 4
#    COMM mpi
#    PASS_REGULAR_EXPRESSION "PASS"
#    FAIL_REGULAR_EXPRESSION "FAIL"
#)
#TRIBITS_ADD_TEST(
#    iceSheet
#    NAME iceSheet_ChainedDoubleHinge_sensitivity3
#    ARGS " chaineddoublehinge.quad.msh chaineddoublehinge-borders.quad.msh chaineddoublehinge-basal-friction.ascii chaineddoublehinge-answer.txt 3 i" 
#    NUM_MPI_PROCS 4
#    COMM mpi
#    PASS_REGULAR_EXPRESSION "PASS"
#    FAIL_REGULAR_EXPRESSION "FAIL"
#)
#TRIBITS_ADD_TEST(
#    iceSheet
#    NAME iceSheet_ChainedDoubleHinge_sensitivity2
#    ARGS " chaineddoublehinge.quad.msh chaineddoublehinge-borders.quad.msh chaineddoublehinge-basal-friction.ascii chaineddoublehinge-answer.txt 2 i" 
#    NUM_MPI_PROCS 4
#    COMM mpi
#    PASS_REGULAR_EXPRESSION "PASS"
#    FAIL_REGULAR_EXPRESSION "FAIL"
#)
#TRIBITS_ADD_TEST(
#    iceSheet
#    NAME iceSheet_ChainedDoubleHinge_sensitivity1
#    ARGS " chaineddoublehinge.quad.msh chaineddoublehinge-borders.quad.msh chaineddoublehinge-basal-friction.ascii chaineddoublehinge-answer.txt 1 i" 
#    NUM_MPI_PROCS 4
#    COMM mpi
#    PASS_REGULAR_EXPRESSION "PASS"
#    FAIL_REGULAR_EXPRESSION "FAIL"
#)
#TRIBITS_ADD_TEST(
#    iceSheet
#    NAME iceSheet_BigRectangle_sensitivity4
#    ARGS " bigrectangle.quad.msh bigrectangle-borders.quad.msh bigrectangle-basal-friction.ascii bigrectangle-s4-answers.txt 4 i" 
#    NUM_MPI_PROCS 4
#    COMM mpi
#    PASS_REGULAR_EXPRESSION "PASS"
#    FAIL_REGULAR_EXPRESSION "FAIL"
#)
#TRIBITS_ADD_TEST(
#    iceSheet
#    NAME iceSheet_BigRectangle_sensitivity3
#    ARGS " bigrectangle.quad.msh bigrectangle-borders.quad.msh bigrectangle-basal-friction.ascii bigrectangle-s3-answers.txt 3 i" 
#    NUM_MPI_PROCS 4
#    COMM mpi
#    PASS_REGULAR_EXPRESSION "PASS"
#    FAIL_REGULAR_EXPRESSION "FAIL"
#)
#TRIBITS_ADD_TEST(
#    iceSheet
#    NAME iceSheet_BigRectangle_sensitivity2
#    ARGS " bigrectangle.quad.msh bigrectangle-borders.quad.msh bigrectangle-basal-friction.ascii bigrectangle-s2-answers.txt 2 i" 
#    NUM_MPI_PROCS 4
#    COMM mpi
#    PASS_REGULAR_EXPRESSION "PASS"
#    FAIL_REGULAR_EXPRESSION "FAIL"
#)
#TRIBITS_ADD_TEST(
#    iceSheet
#    NAME iceSheet_BigRectangle_sensitivity1
#    ARGS " bigrectangle.quad.msh bigrectangle-borders.quad.msh bigrectangle-basal-friction.ascii bigrectangle-s1-answers.txt 1 i" 
#    NUM_MPI_PROCS 4
#    COMM mpi
#    PASS_REGULAR_EXPRESSION "PASS"
#    FAIL_REGULAR_EXPRESSION "FAIL"
#)
#TRIBITS_ADD_TEST(
#    iceSheet
#    NAME iceSheet_SimpleHinge_sensitivity4
#    ARGS " test1.quad.msh test1-borders.quad.msh test1-basal-friction.ascii test1-s4-answers.txt 4 i"
#    PASS_REGULAR_EXPRESSION "PASS"
#    FAIL_REGULAR_EXPRESSION "FAIL"
#)
#
#TRIBITS_ADD_TEST(
#    iceSheet
#    NAME iceSheet_SimpleHinge_sensitivity3
#    ARGS " test1.quad.msh test1-borders.quad.msh test1-basal-friction.ascii test1-s3-answers.txt 3 i"
#    PASS_REGULAR_EXPRESSION "PASS"
#    FAIL_REGULAR_EXPRESSION "FAIL"
#)
#
#TRIBITS_ADD_TEST(
#    iceSheet
#    NAME iceSheet_SimpleHinge_sensitivity2
#    ARGS " test1.quad.msh test1-borders.quad.msh test1-basal-friction.ascii test1-s2-answers.txt 2 i"
#    PASS_REGULAR_EXPRESSION "PASS"
#    FAIL_REGULAR_EXPRESSION "FAIL"
#)
#
#
#TRIBITS_ADD_TEST(
#    iceSheet
#    NAME iceSheet_SimpleHinge_sensitivity1
#    ARGS " test1.quad.msh test1-borders.quad.msh test1-basal-friction.ascii test1-s1-answers.txt 1 i"
#    PASS_REGULAR_EXPRESSION "PASS"
#    FAIL_REGULAR_EXPRESSION "FAIL"
#)
#
#TRIBITS_ADD_TEST(
#    iceSheet
#    NAME iceSheet_Antarctica_16km_sensitivity4
#    ARGS " ant16km.quad.msh ant16km-borders.quad.msh ant16km-basal-friction.ascii ant16km-s4-answers.txt 4 i" 
#    PASS_REGULAR_EXPRESSION "PASS"
#    FAIL_REGULAR_EXPRESSION "FAIL"
#)
#
#TRIBITS_ADD_TEST(
#    iceSheet
#    NAME iceSheet_Antarctica_16km_sensitivity3
#    ARGS " ant16km.quad.msh ant16km-borders.quad.msh ant16km-basal-friction.ascii ant16km-s3-answers.txt 3 i"
#    PASS_REGULAR_EXPRESSION "PASS"
#    FAIL_REGULAR_EXPRESSION "FAIL"
#)
#
#TRIBITS_ADD_TEST(
#    iceSheet
#    NAME iceSheet_Antarctica_16km_sensitivity2
#    ARGS " ant16km.quad.msh ant16km-borders.quad.msh ant16km-basal-friction.ascii ant16km-s2-answers.txt 2 i" 
#    PASS_REGULAR_EXPRESSION "PASS"
#    FAIL_REGULAR_EXPRESSION "FAIL"
#)
#
#TRIBITS_ADD_TEST(
#    iceSheet
#    NAME iceSheet_Antarctica_16km_sensitivity1
#    ARGS " ant16km.quad.msh ant16km-borders.quad.msh ant16km-basal-friction.ascii ant16km-s1-answers.txt 1 i" 
#    PASS_REGULAR_EXPRESSION "PASS"
#    FAIL_REGULAR_EXPRESSION "FAIL"
#)
#
#TRIBITS_ADD_TEST(
#    iceSheet
#    NAME iceSheet_Antarctica_8km_sensitivity4
#    ARGS " ant8km.quad.msh ant8km-borders.quad.msh ant8km-basal-friction.ascii ant8km-s4-answers.txt 4 i" 
#    PASS_REGULAR_EXPRESSION "PASS"
#    FAIL_REGULAR_EXPRESSION "FAIL"
#)
#
#TRIBITS_ADD_TEST(
#    iceSheet
#    NAME iceSheet_Antarctica_8km_sensitivity3
#    ARGS " ant8km.quad.msh ant8km-borders.quad.msh ant8km-basal-friction.ascii ant8km-s3-answers.txt 3 i" 
#    PASS_REGULAR_EXPRESSION "PASS"
#    FAIL_REGULAR_EXPRESSION "FAIL"
#)
#
#TRIBITS_ADD_TEST(
#    iceSheet
#    NAME iceSheet_Antarctica_8km_sensitivity2
#    ARGS " ant8km.quad.msh ant8km-borders.quad.msh ant8km-basal-friction.ascii ant8km-s2-answers.txt 2 i" 
#    PASS_REGULAR_EXPRESSION "PASS"
#    FAIL_REGULAR_EXPRESSION "FAIL"
#)
#
#TRIBITS_ADD_TEST(
#    iceSheet
#    NAME iceSheet_Antarctica_8km_sensitivity1
#    ARGS " ant8km.quad.msh ant8km-borders.quad.msh ant8km-basal-friction.ascii ant8km-s1-answers.txt 1 i" 
#    PASS_REGULAR_EXPRESSION "PASS"
#    FAIL_REGULAR_EXPRESSION "FAIL"
#)
#
#TRIBITS_ADD_TEST(
#    iceSheet
#    NAME iceSheet_Antarctica_4km_sensitivity4
#    ARGS " ant4km.quad.msh ant4km-borders.quad.msh ant4km-basal-friction.ascii ant4km-s4-answers.txt 4 i" 
#    NUM_MPI_PROCS 8
#    COMM mpi
#    PASS_REGULAR_EXPRESSION "PASS"
#    FAIL_REGULAR_EXPRESSION "FAIL"
#)
#TRIBITS_ADD_TEST(
#    iceSheet
#    NAME iceSheet_Antarctica_4km_sensitivity3
#    ARGS " ant4km.quad.msh ant4km-borders.quad.msh ant4km-basal-friction.ascii ant4km-s3-answers.txt 3 i" 
#    NUM_MPI_PROCS 8
#    COMM mpi
#    PASS_REGULAR_EXPRESSION "PASS"
#    FAIL_REGULAR_EXPRESSION "FAIL"
#)
#
#TRIBITS_ADD_TEST(
#    iceSheet
#    NAME iceSheet_Antarctica_4km_sensitivity2
#    ARGS " ant4km.quad.msh ant4km-borders.quad.msh ant4km-basal-friction.ascii ant4km-s2-answers.txt 2 i" 
#    NUM_MPI_PROCS 8
#    COMM mpi
#    PASS_REGULAR_EXPRESSION "PASS"
#    FAIL_REGULAR_EXPRESSION "FAIL"
#)
#
#TRIBITS_ADD_TEST(
#    iceSheet
#    NAME iceSheet_Antarctica_4km_sensitivity1
#    ARGS " ant4km.quad.msh ant4km-borders.quad.msh ant4km-basal-friction.ascii ant4km-s1-answers.txt 1 i" 
#    NUM_MPI_PROCS 8
#    COMM mpi
#    PASS_REGULAR_EXPRESSION "PASS"
#    FAIL_REGULAR_EXPRESSION "FAIL"
#)

#TRIBITS_ADD_TEST(
    #iceSheet
    #NAME iceSheet_Antarctica_2km_sensitivity4
    #ARGS " ant2km.quad.msh ant2km-borders.quad.msh ant2km-basal-friction.ascii ant2km-s4-answers.txt 4 i" 
    #NUM_MPI_PROCS 32 
    #COMM mpi
    #PASS_REGULAR_EXPRESSION "PASS"
    #FAIL_REGULAR_EXPRESSION "FAIL"
#)
#TRIBITS_ADD_TEST(
    #iceSheet
    #NAME iceSheet_Antarctica_2km_sensitivity3
    #ARGS " ant2km.quad.msh ant2km-borders.quad.msh ant2km-basal-friction.ascii ant2km-s3-answers.txt 3 i" 
    #NUM_MPI_PROCS 32
    #COMM mpi
    #PASS_REGULAR_EXPRESSION "PASS"
    #FAIL_REGULAR_EXPRESSION "FAIL"
#)
#TRIBITS_ADD_TEST(
    #iceSheet
    #NAME iceSheet_Antarctica_2km_sensitivity2
    #ARGS " ant2km.quad.msh ant2km-borders.quad.msh ant2km-basal-friction.ascii ant2km-s2-answers.txt 2 i" 
    #NUM_MPI_PROCS 32
    #COMM mpi
    #PASS_REGULAR_EXPRESSION "PASS"
    #FAIL_REGULAR_EXPRESSION "FAIL"
#)
#TRIBITS_ADD_TEST(
    #iceSheet
    #NAME iceSheet_Antarctica_2km_sensitivity1
    #ARGS " ant2km.quad.msh ant2km-borders.quad.msh ant2km-basal-friction.ascii ant2km-s1-answers.txt 1 i" 
    #NUM_MPI_PROCS 32
    #COMM mpi
    #PASS_REGULAR_EXPRESSION "PASS"
    #FAIL_REGULAR_EXPRESSION "FAIL"
#)

#TRIBITS_ADD_TEST(
    #iceSheet
    #NAME iceSheet_Antarctica_1km_sensitivity4
    #ARGS " ant1km.quad.msh ant1km-borders.quad.msh ant1km-basal-friction.ascii ant1km-s4-answers.txt 4 i" 
    #NUM_MPI_PROCS 32
    #COMM mpi
    #PASS_REGULAR_EXPRESSION "PASS"
    #FAIL_REGULAR_EXPRESSION "FAIL"
#)
#TRIBITS_ADD_TEST(
    #iceSheet
    #NAME iceSheet_Antarctica_1km_sensitivity3
    #ARGS " ant1km.quad.msh ant1km-borders.quad.msh ant1km-basal-friction.ascii ant1km-s3-answers.txt 3 i" 
    #NUM_MPI_PROCS 32
    #COMM mpi
    #PASS_REGULAR_EXPRESSION "PASS"
    #FAIL_REGULAR_EXPRESSION "FAIL"
#)
#TRIBITS_ADD_TEST(
    #iceSheet
    #NAME iceSheet_Antarctica_1km_sensitivity2
    #ARGS " ant1km.quad.msh ant1km-borders.quad.msh ant1km-basal-friction.ascii ant1km-s2-answers.txt 2 i" 
    #NUM_MPI_PROCS 32
    #COMM mpi
    #PASS_REGULAR_EXPRESSION "PASS"
    #FAIL_REGULAR_EXPRESSION "FAIL"
#)
#TRIBITS_ADD_TEST(
    #iceSheet
    #NAME iceSheet_Antarctica_1km_sensitivity1
    #ARGS " ant1km.quad.msh ant1km-borders.quad.msh ant1km-basal-friction.ascii ant1km-s1-answers.txt 1 i" 
    #NUM_MPI_PROCS 32
    #COMM mpi
    #PASS_REGULAR_EXPRESSION "PASS"
    #FAIL_REGULAR_EXPRESSION "FAIL"
#)
#
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
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/meshData
  DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
)
