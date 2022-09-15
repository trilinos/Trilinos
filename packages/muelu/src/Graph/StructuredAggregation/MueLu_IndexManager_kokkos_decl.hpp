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
#ifndef MUELU_INDEXMANAGER_KOKKOS_DECL_HPP
#define MUELU_INDEXMANAGER_KOKKOS_DECL_HPP

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_Types.hpp"

#ifdef HAVE_MUELU_KOKKOS_REFACTOR
#include <KokkosCompat_ClassicNodeAPI_Wrapper.hpp>

#include "Teuchos_OrdinalTraits.hpp"

#include <Xpetra_Map_fwd.hpp>
#include <Xpetra_Vector_fwd.hpp>
#include <Xpetra_VectorFactory_fwd.hpp>

#include "MueLu_BaseClass.hpp"
#include "MueLu_IndexManager_kokkos_fwd.hpp"

/*****************************************************************************

****************************************************************************/

namespace MueLu {

/*!
    @class IndexManager_kokkos
    @brief Container class for mesh layout and indices calculation.

    @ingroup Aggregation

    Structure holding mesh parameters for structured mesh. Based on these
    parameters the IndexManager_kokkos computes indices in different index
    spaces and it also provides utilites for coarsening.
*/

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  class IndexManager_kokkos : public BaseClass {
#undef MUELU_INDEXMANAGER_KOKKOS_SHORT
#include "MueLu_UseShortNamesOrdinal.hpp"

  public:
    using execution_space = typename Node::execution_space;
    using memory_space    = typename Node::memory_space;
    using device_type     = Kokkos::Device<execution_space, memory_space>;
    using intTupleView    = typename Kokkos::View<int[3], device_type>;
    using LOTupleView     = typename Kokkos::View<LO[3],  device_type>;

  private:

    const int meshLayout = UNCOUPLED;
    int myRank     = -1;
    int numDimensions;                 ///< Number of spacial dimensions in the problem
    int interpolationOrder_;           ///< Interpolation order used by grid transfer operators using these aggregates.
    intTupleView coarseRate;           ///< coarsening rate in each direction
    intTupleView endRate;              ///< adapted coarsening rate at the edge of the mesh in each direction.

    LO lNumFineNodes;                  ///< local number of nodes.
    LO lNumFineNodes10;                ///< local number of nodes per 0-1 slice.
    LOTupleView lFineNodesPerDir;      ///< local number of nodes per direction.

    LO numCoarseNodes;                 ///< local number of nodes remaining after coarsening.
    LO numCoarseNodes10;               ///< local number of nodes per 0-1 slice remaining after coarsening.
    LOTupleView coarseNodesPerDir;     ///< local number of nodes per direction remaing after coarsening.

  public:

    //! Default constructor, return empty object
    IndexManager_kokkos() = default;

    //! Constructs for uncoupled meshes
    IndexManager_kokkos(const int NumDimensions,
                        const int interpolationOrder,
                        const int MyRank,
                        const ArrayView<const LO> LFineNodesPerDir,
                        const ArrayView<const int> CoarseRate);

    virtual ~IndexManager_kokkos() {}

    //! Common setup pattern used for all the different types of undelying mesh
    void setupIM(const int NumDimensions,
                 const int interpolationOrder,
                 const ArrayView<const int> coarseRate,
                 const ArrayView<const LO> LFineNodesPerDir);

    //! Sets basic parameters used to compute indices on the mesh.
    //! This method requires you to have set this->coarseRate.
    void computeMeshParameters();

    int getNumDimensions()      const {return numDimensions;}

    int getInterpolationOrder() const {return interpolationOrder_;}

    LO getNumLocalFineNodes()   const {return lNumFineNodes;}

    LO getNumCoarseNodes()      const {return numCoarseNodes;}

    KOKKOS_INLINE_FUNCTION
    intTupleView getCoarseningRates()      const {return coarseRate;}

    KOKKOS_INLINE_FUNCTION
    intTupleView getCoarseningEndRates()   const {return endRate;}

    KOKKOS_INLINE_FUNCTION
    LOTupleView  getLocalFineNodesPerDir() const {return lFineNodesPerDir;}

    KOKKOS_INLINE_FUNCTION
    LOTupleView  getCoarseNodesPerDir()    const {return coarseNodesPerDir;}

    Array<LO> getCoarseNodesPerDirArray() const;

    KOKKOS_INLINE_FUNCTION
    void getFineLID2FineTuple(const LO myLID, LO (&tuple)[3]) const {
      LO tmp;
      tuple[2] = myLID / (lFineNodesPerDir(1)*lFineNodesPerDir(0));
      tmp      = myLID % (lFineNodesPerDir(1)*lFineNodesPerDir(0));
      tuple[1] = tmp   / lFineNodesPerDir(0);
      tuple[0] = tmp   % lFineNodesPerDir(0);
    } // getFineNodeLocalTuple

    KOKKOS_INLINE_FUNCTION
    void getFineTuple2FineLID(const LO tuple[3], LO& myLID) const {
      myLID = tuple[2]*lNumFineNodes10 + tuple[1]*lFineNodesPerDir[0] + tuple[0];
    } // getFineNodeLID

    KOKKOS_INLINE_FUNCTION
    void getCoarseLID2CoarseTuple(const LO myLID, LO (&tuple)[3]) const {
      LO tmp;
      tuple[2] = myLID / numCoarseNodes10;
      tmp      = myLID % numCoarseNodes10;
      tuple[1] = tmp   / coarseNodesPerDir[0];
      tuple[0] = tmp   % coarseNodesPerDir[0];
    } // getCoarseNodeLocalTuple

    KOKKOS_INLINE_FUNCTION
    void getCoarseTuple2CoarseLID(const LO i, const LO j, const LO k, LO& myLID) const {
      myLID = k*numCoarseNodes10 + j*coarseNodesPerDir[0] + i;
    } // getCoarseNodeLID

  };

} //namespace MueLu

#define MUELU_INDEXMANAGER_KOKKOS_SHORT
#endif // HAVE_MUELU_KOKKOS_REFACTOR
#endif // MUELU_INDEXMANAGER_KOKKOS_DECL_HPP
