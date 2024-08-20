// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_SETUPREGIONSMOOTHERS_DEF_HPP
#define MUELU_SETUPREGIONSMOOTHERS_DEF_HPP

#include <vector>

#include <Tpetra_KokkosCompat_DefaultNode.hpp>

#include <Teuchos_RCP.hpp>

#include <Xpetra_ConfigDefs.hpp>
#include <Xpetra_Export.hpp>
#include <Xpetra_Import.hpp>
#include <Xpetra_Map.hpp>
#include <Xpetra_MultiVector.hpp>
#include <Xpetra_Vector.hpp>
#include <Xpetra_CrsMatrixWrap.hpp>

#include "SetupRegionMatrix_def.hpp"
#include "SetupRegionVector_def.hpp"

using Teuchos::Array;
using Teuchos::ArrayRCP;
using Teuchos::ArrayView;
using Teuchos::ParameterList;
using Teuchos::RCP;

/*! \brief Create list of valid smoother types
 *
 * Create this list here in a routine, such that everyone can use exactly the same list.
 * This avoids copy-and-paste errors.
 *
 * ToDo: replace this list by an enum when we migrate to actual code.
 */
std::map<std::string, int> getListOfValidSmootherTypes() {
  std::map<std::string, int> smootherTypes;
  smootherTypes.insert(std::pair<std::string, int>("None", 0));
  smootherTypes.insert(std::pair<std::string, int>("Jacobi", 1));
  smootherTypes.insert(std::pair<std::string, int>("Gauss", 2));
  smootherTypes.insert(std::pair<std::string, int>("SymmetricGauss", 3));
  smootherTypes.insert(std::pair<std::string, int>("Chebyshev", 4));

  return smootherTypes;
}

/*! \brief Compute inverse of diagonal of the operator
 *
 * Computes the inverse of the diagonal in region format and with interface scaling
 */
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void computeInverseDiagonal(RCP<Teuchos::ParameterList> params,
                            const RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > revisedRowMap,
                            const RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > regionMats,
                            const RCP<Xpetra::Import<LocalOrdinal, GlobalOrdinal, Node> > rowImport)  ///< row importer in region layout [in]
{
#include "Xpetra_UseShortNames.hpp"
  using Teuchos::TimeMonitor;
  RCP<TimeMonitor> tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("Region Relaxation Setup")));

  RCP<Vector> diagReg = VectorFactory::Build(revisedRowMap, true);

  // extract inverse of diagonal from matrix
  diagReg = VectorFactory::Build(regionMats->getRowMap(), true);
  regionMats->getLocalDiagCopy(*diagReg);

  sumInterfaceValues(diagReg, revisedRowMap, rowImport);

  diagReg->reciprocal(*diagReg);

  params->set<RCP<Vector> >("smoothers: inverse diagonal", diagReg);
}

/*! \brief Do Jacobi smoothing
 *
 *  Perform Jacobi smoothing in the region layout using the true diagonal value
 *  recovered from the splitted matrix.
 */
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void jacobiIterate(RCP<Teuchos::ParameterList> smootherParams,
                   RCP<Xpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node> >& regX,             // left-hand side (or solution)
                   const RCP<Xpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node> > regB,        // right-hand side (or residual)
                   const RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > regionMats,  // matrices in true region layout
                   const RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > revisedRowMap,    ///< revised row maps in region layout [in] (actually extracted from regionMats)
                   const RCP<Xpetra::Import<LocalOrdinal, GlobalOrdinal, Node> > rowImport,           ///< row importer in region layout [in]
                   bool& zeroInitGuess) {
#include "Xpetra_UseShortNames.hpp"
  using Teuchos::TimeMonitor;
  RCP<TimeMonitor> tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("Region Jacobi Iterate")));

  // const Scalar SC_ZERO = Teuchos::ScalarTraits<Scalar>::zero();
  const Scalar SC_ONE = Teuchos::ScalarTraits<Scalar>::one();

  const int maxIter    = smootherParams->get<int>("smoother: sweeps");
  const double damping = smootherParams->get<double>("smoother: damping");
  RCP<Vector> diag_inv = smootherParams->get<RCP<Vector> >("smoothers: inverse diagonal");

  RCP<Vector> regRes = VectorFactory::Build(revisedRowMap, true);

  for (int iter = 0; iter < maxIter; ++iter) {
    // Update the residual vector
    if (zeroInitGuess) {
      regX->elementWiseMultiply(damping, *diag_inv, *regB, SC_ONE);
    } else {
      computeResidual(regRes, regX, regB, regionMats, *smootherParams);

      // update solution according to Jacobi's method
      regX->elementWiseMultiply(damping, *diag_inv, *regRes, SC_ONE);
    }
    zeroInitGuess = false;
  }

  return;
}  // jacobiIterate

/*! \brief Do Gauss-Seidel smoothing
 *
 *  Perform Gauss-Seidel smoothing in the region layout using the true diagonal value
 *  recovered from the splitted matrix. Off-diagonal values are just taken as they are
 *  in region format.
 */
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void GSIterate(RCP<Teuchos::ParameterList> smootherParams,
               RCP<Xpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node> >& regX,             // left-hand side (or solution)
               const RCP<Xpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node> > regB,        // right-hand side (or residual)
               const RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > regionMats,  // matrices in true region layout
               const RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > revisedRowMap,    ///< revised row maps in region layout [in] (actually extracted from regionMats)
               const RCP<Xpetra::Import<LocalOrdinal, GlobalOrdinal, Node> > rowImport,           ///< row importer in region layout [in]
               bool& zeroInitGuess,
               bool sgs = false) {
#include "Xpetra_UseShortNames.hpp"
  using Teuchos::TimeMonitor;
  RCP<TimeMonitor> tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("Region Gauss-Seidel Iterate")));

  // Extract user-given and pre-computed data from paremter list
  const int maxIter    = smootherParams->get<int>("smoother: sweeps");
  const double damping = smootherParams->get<double>("smoother: damping");
  RCP<Vector> diag_inv = smootherParams->get<RCP<Vector> >("smoothers: inverse diagonal");

  RCP<Vector> regRes = VectorFactory::Build(revisedRowMap, true);

  const size_t numRows = regionMats->getLocalNumRows();

  // GS iteration loop
  for (int iter = 0; iter < maxIter; ++iter) {
    // Update the residual vector
    if (!zeroInitGuess) {
      computeResidual(regRes, regX, regB, regionMats, *smootherParams);
    }

    // update the solution and the residual

    using MT               = typename Teuchos::ScalarTraits<SC>::magnitudeType;
    RCP<Vector> delta      = VectorFactory::Build(regionMats->getRowMap(), true);
    ArrayRCP<SC> ldelta    = delta->getDataNonConst(0);
    ArrayRCP<SC> OneregX   = regX->getDataNonConst(0);
    ArrayRCP<SC> OneregRes = regRes->getDataNonConst(0);
    if (zeroInitGuess) {  // copy regB to regRes
      ArrayRCP<SC> rhs = regB->getDataNonConst(0);
      for (size_t k = 0; k < numRows; ++k) OneregRes[k] = rhs[k];
    }
    ArrayRCP<const SC> Onediag = diag_inv->getData(0);

    // Loop over all rows in the region matrix
    for (size_t k = 0; k < numRows; ++k) {
      // Extract a single row
      ArrayView<const LO> AAcols;
      ArrayView<const SC> AAvals;
      regionMats->getLocalRowView(k, AAcols, AAvals);
      const int* Acols = AAcols.getRawPtr();
      const SC* Avals  = AAvals.getRawPtr();
      const LO RowLeng = AAvals.size();

      // Loop over entries in row k and perform GS iteration
      for (LO kk = 0; kk < RowLeng; kk++) {
        OneregRes[k] -= Avals[kk] * ldelta[Acols[kk]];
      }
      ldelta[k] = damping * Onediag[k] * OneregRes[k];
      OneregX[k] += ldelta[k];
    }
    zeroInitGuess = false;

    if (sgs) {
      for (size_t k = numRows; k--;) {
        // Extract a single row
        ArrayView<const LO> AAcols;
        ArrayView<const SC> AAvals;
        regionMats->getLocalRowView(k, AAcols, AAvals);
        const int* Acols = AAcols.getRawPtr();
        const SC* Avals  = AAvals.getRawPtr();
        const LO RowLeng = AAvals.size();

        // Loop over entries in row k and perform GS iteration
        for (LO kk = 0; kk < RowLeng; kk++) {
          OneregRes[k] -= Avals[kk] * ldelta[Acols[kk]];
        }
        ldelta[k] = damping * Onediag[k] * OneregRes[k];
        OneregX[k] += ldelta[k];
      }
    }
  }

  return;
}  // GS

//! Transfer region vector to composite format and compute its 2-norm
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
typename Teuchos::ScalarTraits<Scalar>::magnitudeType
calcNorm2(RCP<Xpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node> >& regVec,
          const RCP<Xpetra::Import<LocalOrdinal, GlobalOrdinal, Node> > rowImport) {
#include "Xpetra_UseShortNames.hpp"
  const RCP<const Map> mapComp = rowImport->getSourceMap();
  RCP<Vector> compVec          = VectorFactory::Build(mapComp, true);
  regionalToComposite(regVec, compVec, rowImport);
  typename Teuchos::ScalarTraits<Scalar>::magnitudeType norm = compVec->norm2();

  return norm;
}  // calcNorm2

/*! Compute inner product of two region vectors
 *
 * First, we transform the region vectors to the composite layout. Then, we utilize
 * the Xpetra::Vector::dot() capability to compute the inner product.
 *
 * @param[in] regX First region vector
 * @param[in] regY Second region vector
 * @param[in] rowImport Importer to transfer region vectors to composite layout
 *
 * @return Inner product of regX and regY
 */
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Scalar
dotProd(RCP<Xpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node> >& regX,
        RCP<Xpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node> >& regY,
        const RCP<Xpetra::Import<LocalOrdinal, GlobalOrdinal, Node> > rowImport) {
#include "Xpetra_UseShortNames.hpp"
  const RCP<const Map> mapComp = rowImport->getSourceMap();
  RCP<Vector> compX            = VectorFactory::Build(mapComp, true);
  RCP<Vector> compY            = VectorFactory::Build(mapComp, true);
  regionalToComposite(regX, compX, rowImport);
  regionalToComposite(regY, compY, rowImport);
  SC dotVal = compX->dot(*compY);

  return dotVal;
}  // dotProd

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Scalar
powerMethod(RCP<Teuchos::ParameterList> params,
            const RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > regionMats,
            const RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > revisedRowMap,
            const RCP<Xpetra::Import<LocalOrdinal, GlobalOrdinal, Node> > rowImport,
            const int numIters) {
#include "Xpetra_UseShortNames.hpp"

  RCP<Vector> diag_inv = params->get<RCP<Vector> >("smoothers: inverse diagonal");
  const SC SC_ZERO     = Teuchos::ScalarTraits<Scalar>::zero();
  const SC SC_ONE      = Teuchos::ScalarTraits<Scalar>::one();
  SC lambdaMax         = SC_ZERO;
  SC RQ_top, RQ_bottom, norm;

  RCP<Vector> regX = VectorFactory::Build(revisedRowMap, true);
  RCP<Vector> regY = VectorFactory::Build(revisedRowMap, true);

  regX->randomize();

  norm = calcNorm2(regX, rowImport);
  regX->scale(SC_ONE / norm);

  for (int iter = 0; iter < numIters; ++iter) {
    regionMats->apply(*regX, *regY);                     // A.apply (x, y);
    sumInterfaceValues(regY, revisedRowMap, rowImport);  // step 2

    // Scale by inverse of diagonal
    regY->elementWiseMultiply(SC_ONE, *diag_inv, *regY, SC_ZERO);

    RQ_top    = dotProd(regY, regX, rowImport);
    RQ_bottom = dotProd(regX, regX, rowImport);
    lambdaMax = RQ_top / RQ_bottom;

    norm = calcNorm2(regY, rowImport);

    if (norm == SC_ZERO) {  // Return something reasonable.
      return SC_ZERO;
    }
    regX->update(SC_ONE / norm, *regY, SC_ZERO);
  }

  return lambdaMax;
}  // powerMethod

/*! \brief Performs Chebyshev specific setup
 *
 * Use power method to estimate lambdaMx
 */
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void chebyshevSetup(RCP<Teuchos::ParameterList> params,
                    const RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > regionMats,
                    const RCP<Xpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node> > regionInterfaceScaling,
                    const RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > revisedRowMap,
                    const RCP<Xpetra::Import<LocalOrdinal, GlobalOrdinal, Node> > rowImport) {
#include "Xpetra_UseShortNames.hpp"
  using Teuchos::TimeMonitor;
  RCP<TimeMonitor> tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("Region Chebyshev Setup")));

  // Calculate lambdaMax
  Scalar lambdaMax = 1;
  lambdaMax        = powerMethod(params,
                                 regionMats,
                                 revisedRowMap,
                                 rowImport,
                                 10);
  params->set<Scalar>("chebyshev: lambda max", lambdaMax);

}  // chebyshevSetup

/*! \brief The textbook Chebyshev algorithm from Ifpack2 translated into the region format
 */
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void chebyshevIterate(RCP<Teuchos::ParameterList> smootherParams,
                      RCP<Xpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node> >& regX,             ///< left-hand side (or solution)
                      const RCP<Xpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node> > regB,        ///< right-hand side (or residual)
                      const RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > regionMats,  ///< matrices in true region layout
                      const RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > revisedRowMap,    ///< revised row maps in region layout [in] (actually extracted from regionMats)
                      const RCP<Xpetra::Import<LocalOrdinal, GlobalOrdinal, Node> > rowImport,           ///< row importer in region layout [in]
                      bool& zeroInitGuess                                                                ///< Use a zero vector as initial guess?
) {
#include "Xpetra_UseShortNames.hpp"
  using Teuchos::TimeMonitor;
  RCP<TimeMonitor> tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("Region Chebyshev Iterate")));

  // Extract input data from parameter list
  const int maxIter        = smootherParams->get<int>("smoother: sweeps");
  const Scalar eigRatio    = smootherParams->get<double>("smoother: Chebyshev eigRatio");
  const Scalar lambdaMax   = smootherParams->get<Scalar>("chebyshev: lambda max");
  const Scalar boostFactor = smootherParams->get<double>("smoother: Chebyshev boost factor");
  RCP<Vector> diag_inv     = smootherParams->get<RCP<Vector> >("smoothers: inverse diagonal");

  // Define some constants for convenience
  const Scalar SC_ZERO = Teuchos::ScalarTraits<Scalar>::zero();
  const Scalar SC_ONE  = Teuchos::ScalarTraits<Scalar>::one();
  const Scalar SC_TWO  = Teuchos::as<Scalar>(2);

  // Estimation of eigenvalue interval of interest: [alpha, beta]
  const Scalar alpha = lambdaMax / eigRatio;     // lower bound (estimate via given max-to-min ratio)
  const Scalar beta  = boostFactor * lambdaMax;  // upper bound (estimated via boost factor)

  // Algorithmic constants
  const Scalar delta = SC_TWO / (beta - alpha);
  const Scalar theta = (beta + alpha) / SC_TWO;
  const Scalar s1    = theta * delta;

  // Algorithmic parameters
  Scalar dtemp1 = SC_ZERO;
  Scalar dtemp2 = SC_ZERO;
  Scalar rhokp1 = SC_ZERO;
  Scalar rhok   = SC_ONE / s1;

  RCP<Vector> regRes = VectorFactory::Build(revisedRowMap, true);

  RCP<Vector> regP = VectorFactory::Build(revisedRowMap, true);
  RCP<Vector> regZ = VectorFactory::Build(revisedRowMap, true);

  // First Iteration
  if (zeroInitGuess) {
    regZ->elementWiseMultiply(SC_ONE, *diag_inv, *regB, SC_ZERO);  // Z = D_inv * b
    regP->update(SC_ONE / theta, *regZ, SC_ZERO);                  // P = 1/theta Z
    regX->update(SC_ONE, *regP, SC_ZERO);                          // X = 0 + P
  } else {
    // Compute residual vector
    computeResidual(regRes, regX, regB, regionMats, *smootherParams);

    regZ->elementWiseMultiply(SC_ONE, *diag_inv, *regRes, SC_ZERO);  // z = D_inv * R, that is, D \ R.
    regP->update(SC_ONE / theta, *regZ, SC_ZERO);                    // P = 1/theta Z
    regX->update(SC_ONE, *regP, SC_ONE);                             // X = X + P
  }

  // The rest of the iterations
  for (int i = 1; i < maxIter; ++i) {
    // Compute residual vector
    computeResidual(regRes, regX, regB, regionMats, *smootherParams);

    // z = D_inv * R, that is, D \ R.
    regZ->elementWiseMultiply(SC_ONE, *diag_inv, *regRes, SC_ZERO);

    rhokp1 = SC_ONE / (SC_TWO * s1 - rhok);
    dtemp1 = rhokp1 * rhok;
    dtemp2 = SC_TWO * rhokp1 * delta;
    rhok   = rhokp1;
    regP->update(dtemp2, *regZ, dtemp1);  // P = dtemp2*Z + dtemp1*P
    regX->update(SC_ONE, *regP, SC_ONE);  // X = X + P

    // If we compute the residual here, we could either do R = B -
    // A*X, or R = R - alpha*A*P.  Since we choose the former, we
    // can move the computeResidual call to the top of the loop.
  }

  zeroInitGuess = false;
}  // chebyshevIterate

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void smootherSetup(RCP<Teuchos::ParameterList> params,
                   const RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > revisedRowMap,
                   const RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > regionMats,
                   const RCP<Xpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node> > regionInterfaceScaling,
                   const RCP<Xpetra::Import<LocalOrdinal, GlobalOrdinal, Node> > rowImport)  ///< row importer in region layout [in]
{
  using Teuchos::TimeMonitor;
  RCP<TimeMonitor> tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("Region Smoother: 1 - Setup")));

  const std::string type = params->get<std::string>("smoother: type");

  std::map<std::string, int> smootherTypes = getListOfValidSmootherTypes();

  switch (smootherTypes[type]) {
    case 0:  // None
    {
      break;
    }
    case 1:  // Jacobi
    case 2:  // Gauss-Seidel
    case 3:  // Symmetric Gauss-Seidel
    {
      computeInverseDiagonal(params, revisedRowMap, regionMats, rowImport);
      break;
    }
    case 4:  // Chebyshev
    {
      computeInverseDiagonal(params, revisedRowMap, regionMats, rowImport);
      chebyshevSetup(params, regionMats, regionInterfaceScaling, revisedRowMap, rowImport);
      break;
    }
    default: {
      std::cout << "Unknown smoother: " << type << "!" << std::endl;
      throw;
    }
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void smootherApply(RCP<Teuchos::ParameterList> params,
                   RCP<Xpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node> >& regX,
                   const RCP<Xpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node> > regB,
                   const RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > regionMats,
                   const RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > revisedRowMap,
                   const RCP<Xpetra::Import<LocalOrdinal, GlobalOrdinal, Node> > rowImport,
                   bool& zeroInitGuess) {
  using Teuchos::TimeMonitor;
  RCP<TimeMonitor> tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("Region Smoother: 2 - Apply")));

  const std::string type = params->get<std::string>("smoother: type");

  std::map<std::string, int> smootherTypes = getListOfValidSmootherTypes();

  switch (smootherTypes[type]) {
    case 0:  // None
    {
      break;
    }
    case 1:  // Jacobi
    {
      jacobiIterate(params, regX, regB, regionMats, revisedRowMap, rowImport, zeroInitGuess);
      break;
    }
    case 2:  // Gauss-Seidel
    {
      GSIterate(params, regX, regB, regionMats, revisedRowMap, rowImport, zeroInitGuess);
      break;
    }
    case 3:  // Symmetric Gauss-Seidel
    {
      GSIterate(params, regX, regB, regionMats, revisedRowMap, rowImport, zeroInitGuess, true);
      break;
    }
    case 4:  // Chebyshev
      chebyshevIterate(params, regX, regB, regionMats, revisedRowMap, rowImport, zeroInitGuess);
      {
        break;
      }
    default: {
      std::cout << "Unknown smoother: " << type << "!" << std::endl;
      throw;
    }
  }

}  // smootherApply

#endif  // MUELU_SETUPREGIONSMOOTHERS_DEF_HPP
