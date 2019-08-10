// @HEADER
//
// ***********************************************************************
//
//        MueLu: A package for multigrid based preconditioning
//                  Copyright 2012 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
#ifndef MUELU_SETUPREGIONSMOOTHERS_DEF_HPP
#define MUELU_SETUPREGIONSMOOTHERS_DEF_HPP

#include <vector>
#include <iostream>
#include <numeric>

#ifdef HAVE_MPI
#include "mpi.h"
#endif

#include <Kokkos_DefaultNode.hpp>

#include <Teuchos_RCP.hpp>

#include <Xpetra_ConfigDefs.hpp>
#include <Xpetra_Export.hpp>
#include <Xpetra_Import.hpp>
#include <Xpetra_Map.hpp>
#include <Xpetra_MultiVector.hpp>
#include <Xpetra_Vector.hpp>
#include <Xpetra_CrsMatrixWrap.hpp>

#include <MueLu_CreateXpetraPreconditioner.hpp>
#include <MueLu_Utilities.hpp>

#include "SetupRegionMatrix_def.hpp"


#if defined(HAVE_MUELU_TPETRA) && defined(HAVE_MUELU_AMESOS2)
#include <Amesos2_config.h>
#include <Amesos2.hpp>
#endif

using Teuchos::RCP;
using Teuchos::ArrayRCP;
using Teuchos::Array;
using Teuchos::ArrayView;
using Teuchos::ParameterList;

/*! \brief performs Jacobi setup
 *
 * Computes the inverse of the diagonal in region format and with interface scaling
 */
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void jacobiSetup(RCP<Teuchos::ParameterList> params,
                 const int maxRegPerProc,
                 const std::vector<RCP<Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > > revisedRowMapPerGrp,
                 const std::vector<RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > > regionGrpMats,
                 const std::vector<RCP<Xpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node> > > regionInterfaceScaling) {
#include "Xpetra_UseShortNames.hpp"
  const Scalar SC_ZERO = Teuchos::ScalarTraits<Scalar>::zero();
  const Scalar SC_ONE = Teuchos::ScalarTraits<Scalar>::one();

  std::vector<RCP<Vector> > regRes(maxRegPerProc);
  createRegionalVector(regRes, maxRegPerProc, revisedRowMapPerGrp);

  // extract diagonal from region matrices, recover true diagonal values, invert diagonal
  Teuchos::Array<RCP<Vector> > diag(maxRegPerProc);
  for (int j = 0; j < maxRegPerProc; j++) {
    // extract inverse of diagonal from matrix
    diag[j] = VectorFactory::Build(regionGrpMats[j]->getRowMap(), true);
    regionGrpMats[j]->getLocalDiagCopy(*diag[j]);
    diag[j]->elementWiseMultiply(SC_ONE, *diag[j], *regionInterfaceScaling[j], SC_ZERO); // ToDo Does it work to pass in diag[j], but also return into the same variable?
    diag[j]->reciprocal(*diag[j]);
  }

  params->set<Teuchos::Array<RCP<Vector> > >("jacobi: inverse diagonal", diag);
}

/*! \brief Do Jacobi smoothing
 *
 *  Perform Jacobi smoothing in the region layout using the true diagonal value
 *  recovered from the splitted matrix.
 */
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void jacobiIterate(RCP<Teuchos::ParameterList> smootherParams,
                   std::vector<RCP<Xpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node> > >& regX, // left-hand side (or solution)
                   const std::vector<RCP<Xpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node> > > regB, // right-hand side (or residual)
                   const std::vector<RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > > regionGrpMats, // matrices in true region layout
                   const std::vector<RCP<Xpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node> > > regionInterfaceScaling, // recreate on coarse grid by import Add on region vector of ones
                   const int maxRegPerProc, ///< max number of regions per proc [in]
                   const RCP<Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > mapComp, ///< composite map
                   const std::vector<RCP<Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > > rowMapPerGrp, ///< row maps in region layout [in] requires the mapping of GIDs on fine mesh to "filter GIDs"
                   const std::vector<RCP<Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > > revisedRowMapPerGrp, ///< revised row maps in region layout [in] (actually extracted from regionGrpMats)
                   const std::vector<RCP<Xpetra::Import<LocalOrdinal, GlobalOrdinal, Node> > > rowImportPerGrp ///< row importer in region layout [in]
    )
{
#include "Xpetra_UseShortNames.hpp"
  // const Scalar SC_ZERO = Teuchos::ScalarTraits<Scalar>::zero();
  const Scalar SC_ONE = Teuchos::ScalarTraits<Scalar>::one();

  const int maxIter    = smootherParams->get<int>   ("smoother: sweeps");
  const double damping = smootherParams->get<double>("smoother: damping");
  Teuchos::Array<RCP<Vector> > diag_inv = smootherParams->get<Teuchos::Array<RCP<Vector> > >("jacobi: inverse diagonal");

  std::vector<RCP<Vector> > regRes(maxRegPerProc);
  createRegionalVector(regRes, maxRegPerProc, revisedRowMapPerGrp);


  for (int iter = 0; iter < maxIter; ++iter) {

    /* Update the residual vector
     * 1. Compute tmp = A * regX in each region
     * 2. Sum interface values in tmp due to duplication (We fake this by scaling to reverse the basic splitting)
     * 3. Compute r = B - tmp
     */
    for (int j = 0; j < maxRegPerProc; j++) { // step 1

//      Teuchos::RCP<Teuchos::FancyOStream> fos = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
//      regRes[j]->getMap()->describe(*fos, Teuchos::VERB_EXTREME);
//      regionGrpMats[j]->getRangeMap()->describe(*fos, Teuchos::VERB_EXTREME);

//      TEUCHOS_ASSERT(regionGrpMats[j]->getDomainMap()->isSameAs(*regX[j]->getMap()));
//      TEUCHOS_ASSERT(regionGrpMats[j]->getRangeMap()->isSameAs(*regRes[j]->getMap()));

      regionGrpMats[j]->apply(*regX[j], *regRes[j]);
    }

    sumInterfaceValues(regRes, mapComp, maxRegPerProc, rowMapPerGrp,
                       revisedRowMapPerGrp, rowImportPerGrp); // step 2

    for (int j = 0; j < maxRegPerProc; j++) { // step 3
      regRes[j]->update(1.0, *regB[j], -1.0);
    }

    // check for convergence
    {
      RCP<Vector> compRes = VectorFactory::Build(mapComp, true);
      regionalToComposite(regRes, compRes, maxRegPerProc, rowMapPerGrp,
                          rowImportPerGrp, Xpetra::ADD);
      typename Teuchos::ScalarTraits<Scalar>::magnitudeType normRes = compRes->norm2();

      if (normRes < 1.0e-12) {return;}
    }

    for (int j = 0; j < maxRegPerProc; j++) {
      // update solution according to Jacobi's method
      regX[j]->elementWiseMultiply(damping, *diag_inv[j], *regRes[j], SC_ONE);
    }
  }

  return;
} // jacobiIterate

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void smootherSetup(RCP<Teuchos::ParameterList> params,
                   const int maxRegPerProc,
                   const std::vector<RCP<Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > > revisedRowMapPerGrp,
                   const std::vector<RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > > regionGrpMats,
                   const std::vector<RCP<Xpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node> > > regionInterfaceScaling) {
  const std::string type = params->get<std::string>("smoother: type");

  std::map<std::string, int> smootherTypes;
  smootherTypes.insert(std::pair<std::string, int>("None",      0));
  smootherTypes.insert(std::pair<std::string, int>("Jacobi",    1));
  smootherTypes.insert(std::pair<std::string, int>("Gauss",     2));
  smootherTypes.insert(std::pair<std::string, int>("Chebyshev", 3));

  switch(smootherTypes[type]) {
  case 0:
    break;
  case 1:
    jacobiSetup(params, maxRegPerProc, revisedRowMapPerGrp, regionGrpMats, regionInterfaceScaling);
    break;
  case 2:
    std::cout << "Gauss-Seidel smoother not implemented yet no smoother is applied" << std::endl;
    break;
  case 3:
    std::cout << "Chebyshev smoother not implemented yet no smoother is applied" << std::endl;
    break;
  default:
    std::cout << "Unknow smoother: " << type << "!" << std::endl;
    throw;
  }
}

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void smootherApply(RCP<Teuchos::ParameterList> params,
                   const int maxRegPerProc,
                   std::vector<RCP<Xpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node> > >& regX,
                   const std::vector<RCP<Xpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node> > > regB,
                   const std::vector<RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > > regionGrpMats,
                   const std::vector<RCP<Xpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node> > > regionInterfaceScaling,
                   const RCP<Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > mapComp,
                   const std::vector<RCP<Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > > rowMapPerGrp,
                   const std::vector<RCP<Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > > revisedRowMapPerGrp,
                   const std::vector<RCP<Xpetra::Import<LocalOrdinal, GlobalOrdinal, Node> > > rowImportPerGrp) {
  const std::string type = params->get<std::string>("smoother: type");

  std::map<std::string, int> smootherTypes;
  smootherTypes.insert(std::pair<std::string, int>("None",         0));
  smootherTypes.insert(std::pair<std::string, int>("Jacobi",       1));
  smootherTypes.insert(std::pair<std::string, int>("Gauss-Seidel", 2));
  smootherTypes.insert(std::pair<std::string, int>("Chebyshev",    3));

  switch(smootherTypes[type]) {
  case 0:
    break;
  case 1:
    jacobiIterate(params, regX, regB, regionGrpMats, regionInterfaceScaling, maxRegPerProc,
                  mapComp, rowMapPerGrp, revisedRowMapPerGrp, rowImportPerGrp);
  case 2:
    break;
  case 3:
    break;
  }

} // smootherApply

#endif // MUELU_SETUPREGIONSMOOTHERS_DEF_HPP
