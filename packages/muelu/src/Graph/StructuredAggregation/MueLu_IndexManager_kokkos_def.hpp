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
//                    Jonathan Hu        (jhu@sandia.gov)
//                    Ray Tuminaro       (rstumin@sandia.gov)
//                    Luc Berger-Vergiat (lberge@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
#ifndef MUELU_INDEXMANAGER_DEF_KOKKOS_HPP
#define MUELU_INDEXMANAGER_DEF_KOKKOS_HPP

#include <utility>

#include "Teuchos_OrdinalTraits.hpp"

#include <Xpetra_MapFactory.hpp>

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_Types.hpp"
#include "MueLu_BaseClass.hpp"
#include "MueLu_Exceptions.hpp"
#include <MueLu_IndexManager_kokkos_decl.hpp>

/*****************************************************************************

****************************************************************************/

namespace MueLu {

  template<class LocalOrdinal, class GlobalOrdinal, class Node>
  IndexManager_kokkos<LocalOrdinal, GlobalOrdinal, Node>::
  IndexManager_kokkos(const int NumDimensions,
                      const int interpolationOrder,
                      const int MyRank,
                      const ArrayView<const LO> LFineNodesPerDir,
                      const ArrayView<const int> CoarseRate) :
    myRank(MyRank), coarseRate("coarsening rate"), endRate("endRate"),
    lFineNodesPerDir("lFineNodesPerDir"), coarseNodesPerDir("lFineNodesPerDir") {

    RCP<Teuchos::FancyOStream> out;
    if(const char* dbg = std::getenv("MUELU_INDEXMANAGER_DEBUG")) {
      out = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
      out->setShowAllFrontMatter(false).setShowProcRank(true);
    } else {
      out = Teuchos::getFancyOStream(rcp(new Teuchos::oblackholestream()));
    }

    setupIM(NumDimensions, interpolationOrder, CoarseRate, LFineNodesPerDir);

    *out << "Done setting up the IndexManager" << std::endl;

    computeMeshParameters();

    *out << "Computed Mesh Parameters" << std::endl;

  } // IndexManager_kokkos Constructor

  template<class LocalOrdinal, class GlobalOrdinal, class Node>
  void IndexManager_kokkos<LocalOrdinal, GlobalOrdinal, Node>::
  setupIM(const int NumDimensions, const int interpolationOrder,
          const ArrayView<const int> CoarseRate, const ArrayView<const LO> LFineNodesPerDir) {

    numDimensions       = NumDimensions;
    interpolationOrder_ = interpolationOrder;

    TEUCHOS_TEST_FOR_EXCEPTION((LFineNodesPerDir.size() != 3)
                               && (LFineNodesPerDir.size() != numDimensions),
                               Exceptions::RuntimeError,
                               "LFineNodesPerDir has to be of size 3 or of size numDimensions!");

    typename Kokkos::View<LO[3], device_type>::HostMirror lFineNodesPerDir_h = Kokkos::create_mirror_view(lFineNodesPerDir);
    Kokkos::deep_copy(lFineNodesPerDir_h, lFineNodesPerDir);
    typename Kokkos::View<LO[3], device_type>::HostMirror coarseRate_h = Kokkos::create_mirror_view(coarseRate);
    Kokkos::deep_copy(coarseRate_h, coarseRate);

    // Load coarse rate, being careful about formating
    // Also load lFineNodesPerDir
    for(int dim = 0; dim < 3; ++dim) {
      if(dim < getNumDimensions()) {
        lFineNodesPerDir_h(dim) = LFineNodesPerDir[dim];
        if(CoarseRate.size() == 1) {
          coarseRate_h(dim) = CoarseRate[0];
        } else if(CoarseRate.size() == getNumDimensions()) {
          coarseRate_h(dim) = CoarseRate[dim];
        }
      } else {
        lFineNodesPerDir_h(dim) = 1;
        coarseRate_h(dim) = 1;
      }
    }

    Kokkos::deep_copy(lFineNodesPerDir, lFineNodesPerDir_h);
    Kokkos::deep_copy(coarseRate, coarseRate_h);

  } // setupIM

  template<class LocalOrdinal, class GlobalOrdinal, class Node>
  void IndexManager_kokkos<LocalOrdinal, GlobalOrdinal, Node>::computeMeshParameters() {

    RCP<Teuchos::FancyOStream> out;
    if(const char* dbg = std::getenv("MUELU_INDEXMANAGER_DEBUG")) {
      out = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
      out->setShowAllFrontMatter(false).setShowProcRank(true);
    } else {
      out = Teuchos::getFancyOStream(rcp(new Teuchos::oblackholestream()));
    }

    typename Kokkos::View<int[3], device_type>::HostMirror coarseRate_h = Kokkos::create_mirror_view(coarseRate);
    typename Kokkos::View<int[3], device_type>::HostMirror endRate_h = Kokkos::create_mirror_view(endRate);


    typename Kokkos::View<LO[3], device_type>::HostMirror lFineNodesPerDir_h = Kokkos::create_mirror_view(lFineNodesPerDir);
    typename Kokkos::View<LO[3], device_type>::HostMirror coarseNodesPerDir_h = Kokkos::create_mirror_view(coarseNodesPerDir);
    Kokkos::deep_copy(lFineNodesPerDir_h, lFineNodesPerDir);
    Kokkos::deep_copy(coarseRate_h, coarseRate);

    lNumFineNodes10 = lFineNodesPerDir_h(1)*lFineNodesPerDir_h(0);
    lNumFineNodes   = lFineNodesPerDir_h(2)*lNumFineNodes10;
    for(int dim = 0; dim < 3; ++dim) {
      if(dim < numDimensions) {
        endRate_h(dim) = (lFineNodesPerDir_h(dim) - 1) % coarseRate_h(dim);
        if(endRate_h(dim) == 0) {endRate_h(dim) = coarseRate_h(dim);}

      } else { // Default value for dim >= numDimensions
        endRate_h(dim) = 1;
      }
    }

    *out << "lFineNodesPerDir: {" << lFineNodesPerDir_h(0) << ", " << lFineNodesPerDir_h(1) << ", "
         << lFineNodesPerDir_h(2) << "}" << std::endl;
    *out << "endRate: {" << endRate_h(0) << ", " << endRate_h(1) << ", "
         << endRate_h(2) << "}"   << std::endl;

    // Here one element can represent either the degenerate case of one node or the more general
    // case of two nodes, i.e. x---x is a 1D element with two nodes and x is a 1D element with
    // one node. This helps generating a 3D space from tensorial products...
    // A good way to handle this would be to generalize the algorithm to take into account the
    // discretization order used in each direction, at least in the FEM sense, since a 0 degree
    // discretization will have a unique node per element. This way 1D discretization can be
    // viewed as a 3D problem with one 0 degree element in the y direction and one 0 degre
    // element in the z direction.
    // !!! Operations below are aftecting both local and global values that have two         !!!
    // different orientations. Orientations can be interchanged using mapDirG2L and mapDirL2G.
    // coarseRate, endRate and offsets are in the global basis, as well as all the variables
    // starting with a g.
    // !!! while the variables starting with an l are in the local basis.                    !!!
    for(int dim = 0; dim < 3; ++dim) {
      if(dim < numDimensions) {
        // Check whether the partition includes the "end" of the mesh which means that endRate
        // will apply. Also make sure that endRate is not 0 which means that the mesh does not
        // require a particular treatment at the boundaries.
        coarseNodesPerDir_h(dim) = (lFineNodesPerDir_h(dim) - endRate_h(dim) - 1)
          / coarseRate_h(dim) + 2;

      } else { // Default value for dim >= numDimensions
        // endRate[dim] = 1;
        coarseNodesPerDir_h(dim) = 1;
      } // if (dim < numDimensions)

      // This would happen if the rank does not own any nodes but in that case a subcommunicator
      // should be used so this should really not be a concern.
      if(lFineNodesPerDir_h(dim) < 1) {coarseNodesPerDir_h(dim) = 0;}
    } // Loop for dim=0:3

    // Compute cummulative values
    numCoarseNodes10 = coarseNodesPerDir_h(0)*coarseNodesPerDir_h(1);
    numCoarseNodes   = numCoarseNodes10*coarseNodesPerDir_h(2);

    *out << "coarseNodesPerDir: {" << coarseNodesPerDir_h(0) << ", "
         << coarseNodesPerDir_h(1) << ", " << coarseNodesPerDir_h(2) << "}" << std::endl;
    *out << "numCoarseNodes=" << numCoarseNodes << std::endl;

    // Copy Host data to Device.
    Kokkos::deep_copy(coarseRate, coarseRate_h);
    Kokkos::deep_copy(endRate, endRate_h);
    Kokkos::deep_copy(lFineNodesPerDir, lFineNodesPerDir_h);
    Kokkos::deep_copy(coarseNodesPerDir, coarseNodesPerDir_h);
  }

  template<class LocalOrdinal, class GlobalOrdinal, class Node>
  Array<LocalOrdinal> IndexManager_kokkos<LocalOrdinal, GlobalOrdinal, Node>::
  getCoarseNodesPerDirArray() const {
    typename LOTupleView::HostMirror coarseNodesPerDir_h = Kokkos::create_mirror_view(coarseNodesPerDir);
    Kokkos::deep_copy(coarseNodesPerDir_h, coarseNodesPerDir);
    Array<LO> coarseNodesPerDirArray(3);

    for(int dim = 0; dim < 3; ++dim) {
      coarseNodesPerDirArray[dim] = coarseNodesPerDir_h(dim);
    }

    return coarseNodesPerDirArray;
  } // getCoarseNodesData

} //namespace MueLu

#define MUELU_INDEXMANAGER_KOKKOS_SHORT
#endif // MUELU_INDEXMANAGER_DEF_KOKKOS_HPP
