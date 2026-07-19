TRIBITS_ADD_EXECUTABLE(dirichlet_boundary_detection SOURCES main.cpp)

TRIBITS_COPY_FILES_TO_BINARY_DIR(DirichletBoundaryDetectionCopyFiles
    DEST_FILES Parameters.xml
    EXEDEPS dirichlet_boundary_detection
)

FUNCTION(ADD_DIRICHLET_BOUNDARY_TEST method detection)
    TRIBITS_ADD_TEST(dirichlet_boundary_detection
        NAME DirichletBoundaryDetection_${method}_${detection}
        ARGS "--method=${method} --detection=${detection}"
        COMM mpi
        NUM_MPI_PROCS 4
    )
ENDFUNCTION()

ADD_DIRICHLET_BOUNDARY_TEST(TwoLevelPreconditioner_GDSW      enabled)
ADD_DIRICHLET_BOUNDARY_TEST(TwoLevelPreconditioner_GDSW      disabled)
ADD_DIRICHLET_BOUNDARY_TEST(TwoLevelBlockPreconditioner_GDSW enabled)
ADD_DIRICHLET_BOUNDARY_TEST(TwoLevelBlockPreconditioner_GDSW disabled)
ADD_DIRICHLET_BOUNDARY_TEST(GDSWPreconditioner               enabled)
ADD_DIRICHLET_BOUNDARY_TEST(GDSWPreconditioner               disabled)
ADD_DIRICHLET_BOUNDARY_TEST(RGDSWPreconditioner              enabled)
ADD_DIRICHLET_BOUNDARY_TEST(RGDSWPreconditioner              disabled)
