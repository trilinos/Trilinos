SET(LIB_REQUIRED_DEP_PACKAGES Teuchos Xpetra KokkosCore KokkosContainers KokkosKernels Tpetra Amesos2 Ifpack2)
SET(LIB_OPTIONAL_DEP_PACKAGES AvatarT Belos Teko
                              Intrepid2
                              Zoltan2Core Stratimikos Thyra ThyraTpetraAdapters)
SET(TEST_REQUIRED_DEP_PACKAGES Galeri)
SET(TEST_OPTIONAL_DEP_PACKAGES Pamgen)
SET(LIB_REQUIRED_DEP_TPLS BLAS LAPACK)
SET(LIB_OPTIONAL_DEP_TPLS Boost MATLAB AmgX CGAL ViennaCL MKL Avatar CUSPARSE MAGMASparse mlpack)
SET(TEST_REQUIRED_DEP_TPLS)
SET(TEST_OPTIONAL_DEP_TPLS HYPRE PETSC)
