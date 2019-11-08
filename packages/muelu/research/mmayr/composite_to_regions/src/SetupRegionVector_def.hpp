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
#ifndef MUELU_SETUPREGIONVECTOR_DEF_HPP
#define MUELU_SETUPREGIONVECTOR_DEF_HPP

#include <vector>
#include <iostream>
#include <numeric>

#ifdef HAVE_MPI
#include "mpi.h"
#endif

#include <Kokkos_DefaultNode.hpp>

#include <Xpetra_ConfigDefs.hpp>
#include <Xpetra_Export.hpp>
#include <Xpetra_Import.hpp>
#include <Xpetra_Map.hpp>
#include <Xpetra_MultiVector.hpp>
#include <Xpetra_Vector.hpp>

using Teuchos::RCP;
using Teuchos::ArrayRCP;
using Teuchos::Array;

//! Create an empty vector in the regional layout
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void createRegionalVector(Teuchos::Array<RCP<Xpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node> > >& regVecs, ///< regional vector to be filled
                          const int maxRegPerProc, ///< max number of regions per process
                          const std::vector<RCP<Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > > revisedRowMapPerGrp ///< regional map
                          )
{
#include "Xpetra_UseShortNames.hpp"
  regVecs.resize(maxRegPerProc);
  for (int j = 0; j < maxRegPerProc; j++)
    regVecs[j] = VectorFactory::Build(revisedRowMapPerGrp[j], true);

  return;
} // createRegionalVector

/*! \brief Transform composite vector to regional layout
 *
 *  Starting from a vector in composite layout, we
 *  1. import it into an auxiliary vector in the quasiRegional layout
 *  2. replace the quasiRegional map of the auxiliary vector with the regional map
 */
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void compositeToRegional(RCP<Xpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node> > compVec, ///< Vector in composite layout [in]
                         Array<RCP<Xpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node> > >& quasiRegVecs, ///< Vector in quasiRegional layout [in/out]
                         Array<RCP<Xpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node> > >& regVecs, ///< Vector in regional layout [in/out]
                         const int maxRegPerProc, ///< max number of regions per proc [in]
                         const std::vector<RCP<Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > > rowMapPerGrp, ///< row maps in region layout [in]
                         const std::vector<RCP<Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > > revisedRowMapPerGrp, ///< revised row maps in region layout [in]
                         const std::vector<RCP<Xpetra::Import<LocalOrdinal, GlobalOrdinal, Node> > > rowImportPerGrp ///< row importer in region layout [in]
                         )
{
#include "Xpetra_UseShortNames.hpp"
  // quasiRegional layout
  for (int j = 0; j < maxRegPerProc; j++) {
    // create empty vectors and fill it by extracting data from composite vector
    quasiRegVecs[j] = VectorFactory::Build(rowMapPerGrp[j], true);
    TEUCHOS_ASSERT(!quasiRegVecs[j].is_null());
    quasiRegVecs[j]->doImport(*compVec, *(rowImportPerGrp[j]), Xpetra::INSERT);
  }

  // regional layout
  for (int j = 0; j < maxRegPerProc; j++) {
    // create regVecs vector (copy from quasiRegVecs and swap the map)
    regVecs[j] = quasiRegVecs[j]; // assignment operator= does deep copy in Xpetra
    TEUCHOS_ASSERT(!regVecs[j].is_null());
    regVecs[j]->replaceMap(revisedRowMapPerGrp[j]);
  }

  return;
} // compositeToRegional

/*! \brief Transform regional vector to composite layout
 *
 *  Starting from a \c Epetra_Vector in regional layout, we
 *  1. replace the regional map with the quasiRegional map
 *  2. export it into a vector with composite layout using the \c CombineMode \c Add.
 *     Note: on-process values also have to be added to account for region interfaces inside a process.
 *
 *  \note We also need the capability to add processor-local values. This is not supported by
 *  available CombineMode options in Xpetra/Tpetra, so we use a manual implementation here.
 */
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void regionalToComposite(const Array<RCP<Xpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node> > >& regVec, ///< Vector in region layout [in]
                         RCP<Xpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node> > compVec, ///< Vector in composite layout [in/out]
                         const int maxRegPerProc, ///< max number of regions per proc
                         const std::vector<RCP<Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > > rowMapPerGrp, ///< row maps in quasiRegion layout [in]
                         const std::vector<RCP<Xpetra::Import<LocalOrdinal, GlobalOrdinal, Node> > > rowImportPerGrp, ///< row importer in region layout [in]
                         const Xpetra::CombineMode combineMode ///< Combine mode for import/export [in]
                         )
{
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
    for(int grpIdx = 0; grpIdx < maxRegPerProc; ++grpIdx) {
      tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("regionalToComposite: 2 - quasiRegVec")));
      quasiRegVec = regVec[grpIdx];
      TEUCHOS_ASSERT(Teuchos::nonnull(quasiRegVec));
      quasiRegVec->replaceMap(rowImportPerGrp[grpIdx]->getTargetMap());

      tm = Teuchos::null;
      tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("regionalToComposite: 3 - partialCompVec")));

      RCP<Vector> partialCompVec = VectorFactory::Build(rowImportPerGrp[0]->getSourceMap(), true);
      TEUCHOS_ASSERT(Teuchos::nonnull(partialCompVec));
      TEUCHOS_ASSERT(partialCompVec->getLocalLength() == compVecLocalLength);
      // ToDo (mayr.mt) Use input variable 'combineMode'
      partialCompVec->doExport(*quasiRegVec, *(rowImportPerGrp[grpIdx]), Xpetra::ADD);

      tm = Teuchos::null;
      tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("regionalToComposite: 4 - compVec->sumIntoLocalValue")));

      ArrayRCP<const SC> partialCompVecData = partialCompVec->getData(0);
      ArrayRCP<SC>       compVecData        = compVec->getDataNonConst(0);
      for(size_t entryIdx = 0; entryIdx < compVecLocalLength; ++entryIdx) {
        compVecData[entryIdx] += partialCompVecData[entryIdx];
      }

      tm = Teuchos::null;
    }
  }

  return;
} // regionalToComposite

/*! \brief Sum region interface values
 *
 *  Sum values of interface GIDs using the underlying Export() routines. Technically, we perform the
 *  exchange/summation of interface data by exporting a regional vector to the composite layout and
 *  then immediately importing it back to the regional layout. The Export() involved when going to the
 *  composite layout takes care of the summation of interface values.
 */
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void sumInterfaceValues(Array<RCP<Xpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node> > >& regVec,
                        const RCP<Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > compMap,
                        const int maxRegPerProc, ///< max number of regions per proc [in]
                        const std::vector<RCP<Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > > rowMapPerGrp,///< row maps in region layout [in]
                        const std::vector<RCP<Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > > revisedRowMapPerGrp,///< revised row maps in region layout [in]
                        const std::vector<RCP<Xpetra::Import<LocalOrdinal, GlobalOrdinal, Node> > > rowImportPerGrp ///< row importer in region layout [in])
                        )
{
#include "Xpetra_UseShortNames.hpp"
  using Teuchos::TimeMonitor;

  RCP<Vector> compVec = VectorFactory::Build(compMap, true);
  TEUCHOS_ASSERT(!compVec.is_null());

  Array<Teuchos::RCP<Vector> > quasiRegVec(maxRegPerProc);

  RCP<TimeMonitor> tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("sumInterfaceValues: 1 - regionalToComposite")));
  regionalToComposite(regVec, compVec, maxRegPerProc, rowMapPerGrp,
                      rowImportPerGrp, Xpetra::ADD);

  tm = Teuchos::null;
  tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("sumInterfaceValues: 2 - compositeToRegional")));

  compositeToRegional(compVec, quasiRegVec, regVec, maxRegPerProc,
                      rowMapPerGrp, revisedRowMapPerGrp, rowImportPerGrp);

  tm = Teuchos::null;

  return;
} // sumInterfaceValues


/*! \brief Apply scaling to interface DOFs
 *
 * The vector scalingFactors contains the number of adjacent regions for every DOFs.
 * In the interior of a region, this is always 1. On region interfaces, this is >1.
 *
 * We often need to scale interface entries by 1/N with N being the number of adjacent regions.
 * This can be achieved by setting \c inverseScaling to \c true.
 */
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void scaleInterfaceDOFs(Array<RCP<Xpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node> > >& regVec, ///< Vector to be scaled
                        const Array<RCP<Xpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node> > >& scalingFactors, ///< Vector with scaling factors
                        bool inverseScaling ///< Divide by scaling factors (yes/no?)
                        )
{
  using Vector = Xpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
  using VectorFactory = Xpetra::VectorFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>;

  const Scalar zero = Teuchos::ScalarTraits<Scalar>::zero();
  const Scalar one = Teuchos::ScalarTraits<Scalar>::one();

  for (int j = 0; j < regVec.size(); j++)
  {
    if (inverseScaling)
    {
      RCP<Vector> inverseScalingFactors = VectorFactory::Build(scalingFactors[j]->getMap());
      inverseScalingFactors->reciprocal(*scalingFactors[j]);
      regVec[j]->elementWiseMultiply(one, *regVec[j], *inverseScalingFactors, zero);
    }
    else
    {
      regVec[j]->elementWiseMultiply(one, *regVec[j], *scalingFactors[j], zero);
    }
  }
}

#endif // MUELU_SETUPREGIONVECTOR_DEF_HPP
