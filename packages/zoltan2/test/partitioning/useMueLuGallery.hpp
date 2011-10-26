// @HEADER
// ***********************************************************************
//
//         Zoltan2: Sandia Partitioning Ordering & Coloring Library
//
//                Copyright message goes here.   TODO
//
// ***********************************************************************
//
//  Create a large distributed matrix using MueLu's matrix gallery

#include <Teuchos_Comm.hpp>
#include <Teuchos_CommandLineProcessor.hpp>
#include <MueLu_MatrixFactory.hpp>
#include <MueLu_GalleryParameters.hpp>

#ifndef USEMUELUGALLERYHPP_
#define USEMUELUGALLERYHPP_

using Teuchos::RCP;

template <typename Scalar, typename LNO, typename GNO>
RCP<Tpetra::CrsMatrix<Scalar, LNO, GNO> > 
useMueLuGallery(
  int narg,
  char *arg[],
  const RCP<const Teuchos::Comm<int> > &comm
)
{
  int fail = 0;
  ///////////////////////////////////////////////////
  // options
  int xdim=10;  
  int ydim=10;
  int zdim=10;
  std::string matrixType("Laplace2D");

  Teuchos::CommandLineProcessor clp(false, false);
  clp.setOption("x", &xdim, 
                "number of gridpoints in X dimension for "
                "mesh used to generate matrix.");
  clp.setOption("y", &ydim, 
                "number of gridpoints in Y dimension for "
                "mesh used to generate matrix.");
  clp.setOption("z", &zdim, 
                "number of gridpoints in Z dimension for "
                "mesh used to generate matrix.");
  clp.setOption("matrix", &matrixType, 
                "Matrix type: Laplace1D, Laplace2D, or Laplace3D");

  clp.parse(narg, arg);

  ///////////////////////////////////////////////////
  // Use Muelu to create a distributed matrix.

  typedef Tpetra::CrsMatrix<Scalar, LNO, GNO> tcrsMatrix_t;
  typedef Tpetra::Map<LNO, GNO> tMap_t;

  Teuchos::CommandLineProcessor tclp;
  MueLu::Gallery::Parameters<GNO> params(tclp, xdim, ydim, zdim, matrixType);

  RCP<const tMap_t> map = Teuchos::rcp(new tMap_t(
      params.GetNumGlobalElements(), 0, comm));

  RCP<tcrsMatrix_t> matrix;

  try{
    matrix = MueLu::Gallery::CreateCrsMatrix<Scalar, LNO, GNO,
      Tpetra::Map<LNO, GNO>, Tpetra::CrsMatrix<Scalar, LNO, GNO> >(
        params.GetMatrixType(), map, params.GetParameterList());
  }
  catch (std::exception &e) { 
    std::cerr << comm->getRank() << ": " << e.what();
    fail = 1;   // Probably not enough memory
  }

  // test globally for failures
  int gfail;
  Teuchos::reduceAll<int,int>(*comm, Teuchos::REDUCE_MAX, 1, &fail, &gfail); 
  if (gfail) { 
    if (comm->getRank() == 0) { 
      std::cerr << "Error in " << __func__ << " while creating matrix" 
                << std::endl; 
      std::cout << "FAIL" << std::endl; 
    } 
    return RCP<Tpetra::CrsMatrix<Scalar, LNO, GNO> > (Teuchos::null); 
  } 

  return matrix;
}

#endif
