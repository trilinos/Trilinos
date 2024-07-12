// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_SETUPREGIONVECTOR_DEF_HPP
#define MUELU_SETUPREGIONVECTOR_DEF_HPP

#include <vector>
#include <iostream>
#include <numeric>

#ifdef HAVE_MPI
#include "mpi.h"
#endif

#include <Tpetra_KokkosCompat_DefaultNode.hpp>

#include <Xpetra_ConfigDefs.hpp>
#include <Xpetra_Export.hpp>
#include <Xpetra_Import.hpp>
#include <Xpetra_Map.hpp>
#include <Xpetra_MultiVector.hpp>
#include <Xpetra_Vector.hpp>

using Teuchos::Array;
using Teuchos::ArrayRCP;
using Teuchos::RCP;

/*! \brief Transform composite vector to regional layout
 *
 *  Starting from a vector in composite layout, we
 *  1. import it into an auxiliary vector in the quasiRegional layout
 *  2. replace the quasiRegional map of the auxiliary vector with the regional map
 */
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void compositeToRegional(RCP<Xpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node> > compVec,         ///< Vector in composite layout [in]
                         RCP<Xpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node> >& quasiRegVecs,   ///< Vector in quasiRegional layout [in/out]
                         RCP<Xpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node> >& regVecs,        ///< Vector in regional layout [in/out]
                         const RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > revisedRowMap,  ///< revised row maps in region layout [in]
                         const RCP<Xpetra::Import<LocalOrdinal, GlobalOrdinal, Node> > rowImport          ///< row importer in region layout [in]
) {
#include "Xpetra_UseShortNames.hpp"

  // quasiRegional layout
  // create empty vectors and fill it by extracting data from composite vector
  quasiRegVecs = VectorFactory::Build(rowImport->getTargetMap(), true);
  TEUCHOS_ASSERT(!quasiRegVecs.is_null());
  quasiRegVecs->doImport(*compVec, *(rowImport), Xpetra::INSERT);

  // regional layout
  // create regVecs vector (copy from quasiRegVecs and swap the map)
  regVecs = quasiRegVecs;  // assignment operator= does deep copy in Xpetra
  TEUCHOS_ASSERT(!regVecs.is_null());
  regVecs->replaceMap(revisedRowMap);

  return;
}  // compositeToRegional

/*! \brief Transform composite vector to regional layout
 *
 *  Starting from a vector in composite layout, we
 *  1. import it into an auxiliary vector in the quasiRegional layout
 *  2. replace the quasiRegional map of the auxiliary vector with the regional map
 */
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void compositeToRegional(RCP<Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> > compVec,        ///< Vector in composite layout [in]
                         RCP<Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> >& quasiRegVecs,  ///< Vector in quasiRegional layout [in/out]
                         RCP<Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> >& regVecs,       ///< Vector in regional layout [in/out]
                         const RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > revisedRowMap,      ///< revised row maps in region layout [in]
                         const RCP<Xpetra::Import<LocalOrdinal, GlobalOrdinal, Node> > rowImport              ///< row importer in region layout [in]
) {
#include "Xpetra_UseShortNames.hpp"

  // quasiRegional layout
  // create empty vectors and fill it by extracting data from composite vector
  quasiRegVecs = MultiVectorFactory::Build(rowImport->getTargetMap(), compVec->getNumVectors(), true);
  TEUCHOS_ASSERT(!quasiRegVecs.is_null());
  quasiRegVecs->doImport(*compVec, *(rowImport), Xpetra::INSERT);

  // regional layout
  // create regVecs vector (copy from quasiRegVecs and swap the map)
  regVecs = quasiRegVecs;  // assignment operator= does deep copy in Xpetra
  TEUCHOS_ASSERT(!regVecs.is_null());
  regVecs->replaceMap(revisedRowMap);

  return;
}  // compositeToRegional

/*! \brief Transform regional vector to composite layout
 *
 *  Starting from a \c Epetra_Vector in regional layout, we
 *  1. replace the regional map with the quasiRegional map
 *  2. export it into a vector with composite layout using the \c CombineMode \c Xpetra::ADD.
 *     Note: on-process values also have to be added to account for region interfaces inside a process.
 *
 *  \note We also need the capability to add processor-local values. This is not supported by
 *  available CombineMode options in Xpetra/Tpetra, so we use a manual implementation here.
 */
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void regionalToComposite(const RCP<Xpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node> >& regVec,  ///< Vector in region layout [in]
                         RCP<Xpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node> > compVec,        ///< Vector in composite layout [in/out]
                         const RCP<Xpetra::Import<LocalOrdinal, GlobalOrdinal, Node> > rowImport         ///< row importer in region layout [in]
) {
  /* Let's fake an ADD combine mode that also adds local values by
   * 1. exporting quasiRegional vectors to auxiliary composite vectors (1 per group)
   * 2. add all auxiliary vectors together
   */
#include "Xpetra_UseShortNames.hpp"
  using Teuchos::TimeMonitor;

  RCP<TimeMonitor> tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("regionalToComposite: 1 - compVec setup")));

  const SC SC_ZERO = Teuchos::ScalarTraits<Scalar>::zero();

  const size_t compVecLocalLength = compVec->getLocalLength();
  compVec->putScalar(SC_ZERO);

  tm = Teuchos::null;

  {
    RCP<Vector> quasiRegVec;
    tm          = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("regionalToComposite: 2 - quasiRegVec")));
    quasiRegVec = regVec;
    TEUCHOS_ASSERT(Teuchos::nonnull(quasiRegVec));
    quasiRegVec->replaceMap(rowImport->getTargetMap());

    tm = Teuchos::null;
    tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("regionalToComposite: 3 - partialCompVec")));

    RCP<Vector> partialCompVec = VectorFactory::Build(rowImport->getSourceMap(), true);
    TEUCHOS_ASSERT(Teuchos::nonnull(partialCompVec));
    TEUCHOS_ASSERT(partialCompVec->getLocalLength() == compVecLocalLength);
    partialCompVec->doExport(*quasiRegVec, *(rowImport), Xpetra::ADD);

    tm = Teuchos::null;
    tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("regionalToComposite: 4 - compVec->sumIntoLocalValue")));

    ArrayRCP<const SC> partialCompVecData = partialCompVec->getData(0);
    ArrayRCP<SC> compVecData              = compVec->getDataNonConst(0);
    for (size_t entryIdx = 0; entryIdx < compVecLocalLength; ++entryIdx) {
      compVecData[entryIdx] += partialCompVecData[entryIdx];
    }

    tm = Teuchos::null;
  }

  return;
}  // regionalToComposite

/*! \brief Transform regional vector to composite layout
 *
 *  Starting from a \c Epetra_Vector in regional layout, we
 *  1. replace the regional map with the quasiRegional map
 *  2. export it into a vector with composite layout using the \c CombineMode \c Xpetra::ADD.
 *     Note: on-process values also have to be added to account for region interfaces inside a process.
 *
 *  \note We also need the capability to add processor-local values. This is not supported by
 *  available CombineMode options in Xpetra/Tpetra, so we use a manual implementation here.
 */
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void regionalToComposite(const RCP<Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> >& regVec,  ///< Vector in region layout [in]
                         RCP<Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> > compVec,        ///< Vector in composite layout [in/out]
                         const RCP<Xpetra::Import<LocalOrdinal, GlobalOrdinal, Node> > rowImport              ///< row importer in region layout [in]
) {
  /* Let's fake an ADD combine mode that also adds local values by
   * 1. exporting quasiRegional vectors to auxiliary composite vectors (1 per group)
   * 2. add all auxiliary vectors together
   */
#include "Xpetra_UseShortNames.hpp"
  using Teuchos::TimeMonitor;

  RCP<TimeMonitor> tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("regionalToComposite: 1 - compVec setup")));

  const SC SC_ZERO = Teuchos::ScalarTraits<Scalar>::zero();

  const size_t compVecLocalLength = compVec->getLocalLength();
  compVec->putScalar(SC_ZERO);

  tm = Teuchos::null;

  {
    RCP<MultiVector> quasiRegVec;
    tm          = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("regionalToComposite: 2 - quasiRegVec")));
    quasiRegVec = regVec;
    TEUCHOS_ASSERT(Teuchos::nonnull(quasiRegVec));
    quasiRegVec->replaceMap(rowImport->getTargetMap());

    tm = Teuchos::null;
    tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("regionalToComposite: 3 - partialCompVec")));

    RCP<MultiVector> partialCompVec = MultiVectorFactory::Build(rowImport->getSourceMap(), quasiRegVec->getNumVectors(), true);
    TEUCHOS_ASSERT(Teuchos::nonnull(partialCompVec));
    TEUCHOS_ASSERT(partialCompVec->getLocalLength() == compVecLocalLength);
    partialCompVec->doExport(*quasiRegVec, *(rowImport), Xpetra::ADD);

    tm = Teuchos::null;
    tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("regionalToComposite: 4 - compVec->sumIntoLocalValue")));

    for (LO vecIdx = 0; vecIdx < static_cast<LO>(partialCompVec->getNumVectors()); ++vecIdx) {
      ArrayRCP<const SC> partialCompVecData = partialCompVec->getData(vecIdx);
      ArrayRCP<SC> compVecData              = compVec->getDataNonConst(vecIdx);
      for (size_t entryIdx = 0; entryIdx < compVecLocalLength; ++entryIdx) {
        compVecData[entryIdx] += partialCompVecData[entryIdx];
      }
    }

    tm = Teuchos::null;
  }

  return;
}  // regionalToComposite

/*! \brief Sum region interface values
 *
 *  Sum values of interface GIDs using the underlying Export() routines. Technically, we perform the
 *  exchange/summation of interface data by exporting a regional vector to the composite layout and
 *  then immediately importing it back to the regional layout. The Export() involved when going to the
 *  composite layout takes care of the summation of interface values.
 */
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void sumInterfaceValues(RCP<Xpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node> >& regVec,
                        const RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > revisedRowMap,  ///< revised row maps in region layout [in]
                        const RCP<Xpetra::Import<LocalOrdinal, GlobalOrdinal, Node> > rowImport          ///< row importer in region layout [in])
) {
#include "Xpetra_UseShortNames.hpp"
  using Teuchos::TimeMonitor;

  // Composite map is the same in every group, so just take the first one.
  const RCP<const Map> compMap = rowImport->getSourceMap();

  RCP<Vector> compVec = VectorFactory::Build(compMap, true);
  TEUCHOS_ASSERT(!compVec.is_null());

  RCP<TimeMonitor> tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("sumInterfaceValues: 1 - regionalToComposite")));
  regionalToComposite(regVec, compVec, rowImport);

  tm = Teuchos::null;
  tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("sumInterfaceValues: 2 - compositeToRegional")));

  Teuchos::RCP<Vector> quasiRegVec;
  compositeToRegional(compVec, quasiRegVec, regVec,
                      revisedRowMap, rowImport);

  tm = Teuchos::null;

  return;
}  // sumInterfaceValues

/*! \brief Sum region interface values
 *
 *  Sum values of interface GIDs using the underlying Export() routines. Technically, we perform the
 *  exchange/summation of interface data by exporting a regional vector to the composite layout and
 *  then immediately importing it back to the regional layout. The Export() involved when going to the
 *  composite layout takes care of the summation of interface values.
 */
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void sumInterfaceValues(RCP<Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> >& regVec,
                        const RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > revisedRowMap,  ///< revised row maps in region layout [in]
                        const RCP<Xpetra::Import<LocalOrdinal, GlobalOrdinal, Node> > rowImport          ///< row importer in region layout [in])
) {
#include "Xpetra_UseShortNames.hpp"
  using Teuchos::TimeMonitor;

  // Composite map is the same in every group, so just take the first one.
  const RCP<const Map> compMap = rowImport->getSourceMap();

  RCP<MultiVector> compVec = MultiVectorFactory::Build(compMap, regVec->getNumVectors(), true);
  TEUCHOS_ASSERT(!compVec.is_null());

  RCP<TimeMonitor> tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("sumInterfaceValues: 1 - regionalToComposite")));
  regionalToComposite(regVec, compVec, rowImport);

  tm = Teuchos::null;
  tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("sumInterfaceValues: 2 - compositeToRegional")));

  Teuchos::RCP<MultiVector> quasiRegVec;
  compositeToRegional(compVec, quasiRegVec, regVec,
                      revisedRowMap, rowImport);

  tm = Teuchos::null;

  return;
}  // sumInterfaceValues

/*! \brief Apply scaling to interface DOFs
 *
 * The vector scalingFactors contains the number of adjacent regions for every DOFs.
 * In the interior of a region, this is always 1. On region interfaces, this is >1.
 *
 * We often need to scale interface entries by 1/N with N being the number of adjacent regions.
 * This can be achieved by setting \c inverseScaling to \c true.
 */
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void scaleInterfaceDOFs(RCP<Xpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node> >& regVec,                ///< Vector to be scaled
                        const RCP<Xpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node> >& scalingFactors,  ///< Vector with scaling factors
                        bool inverseScaling                                                                     ///< Divide by scaling factors (yes/no?)
) {
  using Vector        = Xpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
  using VectorFactory = Xpetra::VectorFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>;

  const Scalar zero = Teuchos::ScalarTraits<Scalar>::zero();
  const Scalar one  = Teuchos::ScalarTraits<Scalar>::one();

  if (inverseScaling) {
    RCP<Vector> inverseScalingFactors = VectorFactory::Build(scalingFactors->getMap());
    inverseScalingFactors->reciprocal(*scalingFactors);
    regVec->elementWiseMultiply(one, *regVec, *inverseScalingFactors, zero);
  } else {
    regVec->elementWiseMultiply(one, *regVec, *scalingFactors, zero);
  }
}  // scaleInterfaceDOFs

#endif  // MUELU_SETUPREGIONVECTOR_DEF_HPP
