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
//                    Jeremie Gaidamour (jngaida@sandia.gov)
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
#ifndef MUELU_TENTATIVEPFACTORY_DECL_HPP
#define MUELU_TENTATIVEPFACTORY_DECL_HPP

#include <Teuchos_ScalarTraits.hpp> 

#include <Xpetra_Matrix_fwd.hpp>
#include <Xpetra_MultiVector_fwd.hpp>
#include <Xpetra_MapFactory_fwd.hpp>
#include <Xpetra_Map_fwd.hpp>
#include <Xpetra_Matrix_fwd.hpp>
#include <Xpetra_MultiVector_fwd.hpp>
#include <Xpetra_MultiVectorFactory_fwd.hpp>
#include <Xpetra_ImportFactory_fwd.hpp>
#include <Xpetra_CrsMatrixWrap_fwd.hpp>

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_PFactory.hpp"
#include "MueLu_TentativePFactory_fwd.hpp"
#include "MueLu_CoarseMapFactory_fwd.hpp"

#include "MueLu_Level_fwd.hpp"
#include "MueLu_Aggregates_fwd.hpp"
#include "MueLu_AmalgamationFactory_fwd.hpp"
#include "MueLu_AmalgamationInfo_fwd.hpp"

namespace MueLu {

  /*!
    @class TentativePFactory class.
    @brief Factory for building tentative prolongator.

    Factory for creating tentative prolongator.   Nullspace vectors are split across aggregates so that they
    have local support.  The vectors with local support are factored via LAPACK QR.  The Q becomes the
    tentative prolongator, and the R becomes the coarse nullspace. 

    @ingroup MueLuTransferClasses
  */

  template <class Scalar = double, class LocalOrdinal = int, class GlobalOrdinal = LocalOrdinal, class Node = Kokkos::DefaultNode::DefaultNodeType, class LocalMatOps = typename Kokkos::DefaultKernels<void,LocalOrdinal,Node>::SparseOps> //TODO: or BlockSparseOp ?
  class TentativePFactory : public PFactory {
#undef MUELU_TENTATIVEPFACTORY_SHORT
#include "MueLu_UseShortNames.hpp"

  public:
    //! @name Constructors/Destructors.
    //@{
    
    /*! @brief Constructor.
      \param aggregatesFact -- (optional) factory that creates aggregates.
      \param amalgFact -- (optional) factory creates (Un)Amalgamation information from matrix A
      \param nullspaceFact -- (optional) factory that creates (fine level) null space approximation
      \param AFact -- (optional) factory that creates level matrix A
    */
    TentativePFactory(RCP<const FactoryBase> aggregatesFact = Teuchos::null, RCP<const FactoryBase> amalgFact = Teuchos::null, RCP<const FactoryBase> nullspaceFact = Teuchos::null, RCP<const FactoryBase> AFact = Teuchos::null, RCP<const FactoryBase> coarseMapFact = Teuchos::null);
    
    //! Destructor.
    virtual ~TentativePFactory();
    //@}

    //! Input
    //@{

    void DeclareInput(Level & fineLevel, Level & coarseLevel) const;

    //@}

    //! @name Set/Get methods.
    //@{

    void TentativeWithQR(bool value);

    bool TentativeWithQR();

    //@}

    //! @name Build methods.
    //@{

    void Build(Level & fineLevel, Level & coarseLevel) const;

    void BuildP(Level & fineLevel, Level & coarseLevel) const;

    //@}

    //! @name Get/Set functions

    /*! @brief setDomainMapOffset
     sets offset for domain map DOF Gids in tentative prolongation operator.
     offset must not be smaller than zero. Note: Direct solvers (Amesos/Amesos2) are not properly working with offset > 0.
     */
    void setDomainMapOffset(GlobalOrdinal offset);

    /*! @brief getDomainMapOffset
     * returns offset of the domain DOF map (=coarse map).
     */
    GlobalOrdinal getDomainMapOffset() const;

    /*! @brief getStridingData
     * returns vector with size of striding blocks in the domain DOF map (= coarse map). 
     * e.g. for 2 velocity dofs and 1 pressure dof the vector is (2,1)
     */
    std::vector<size_t> getStridingData() { return stridingInfo_; }

    /*! @brief setStridingData
     * set striding vector for the domain DOF map (= coarse map),
     * e.g. (2,1) for 2 velocity dofs and 1 pressure dof
     */
    void setStridingData(std::vector<size_t> stridingInfo) { stridingInfo_ = stridingInfo; }

    /*! @brief getFixedBlockSize
     * returns the full block size (number of DOFs per node) of the domain DOF map (= coarse map). 
     * This is the sum of all entries in the striding vector.
     * e.g. for 2 velocity dofs and 1 pressure dof the return value is 3.
     */ // TODO remove me, code has been moved to CoarseMapFactory
    size_t getFixedBlockSize() const {
      // sum up size of all strided blocks (= number of dofs per node)
      size_t blkSize = 0;
      std::vector<size_t>::const_iterator it;
      for(it = stridingInfo_.begin(); it != stridingInfo_.end(); ++it) {
        blkSize += *it;
      }
      return blkSize;
    }

    /*! @brief getStridedBlockId
     * returns strided block id for the domain DOF map of Ptent (= coarse map)
     * or -1 if full strided map is stored in the domain map of Ptent (= coarse map)
     */
    LocalOrdinal getStridedBlockId() { return stridedBlockId_; }

    /*! @brief setStridedBlockId
     * set strided block id for the domain DOF map of Ptent (= coarse map)
     * or -1 if full strided map is stored in the domain map of Ptent (= coarse map)
     */
    void setStridedBlockId(LocalOrdinal stridedBlockId) {
      stridedBlockId_ = stridedBlockId;
    }

    //@}
  private:
    //! @name Static methods.
    //@{
    
    typedef typename Teuchos::ScalarTraits<SC>::magnitudeType Magnitude;
    
    /*! @brief Make tentative prolongator with QR.

    We note that the implementation would have been *much* easier
    if Ptent were allowed to have a row map based upon aggregates, i.e., a row
    map such that all DoF's in an aggregate were local and consecutive.
    In this case, we could still have used the matrix A's domain map to FillComplete Ptent.
    However, the prolongator smoothing step Ptent - A*Ptent would then be impossible
    because neither Epetra nor Tpetra support adding matrices with differing row maps.

    The following is a high-level view of how this method is implemented.  The result is a tentative
    prolongator such that Ptent has the same row map as A.

    1) The nullspace NS is communicated in such a way that if proc A owns aggregate k, the part of NS
    corresponding to k is replicated on A.  This will happen if aggregate k contains DoFs belonging
    to proc A and other procs.

    2) Each processor A does a local QR on the parts of the NS corresponding to aggregates that A owns.

    3) The rows of Q that A owns (i.e., rows of Q corresponding to DoFs that A owns) are inserted immediately into Ptentative.

    4) Any rows of Q that A doesn't own  are set  to the owning processor as follows:
    (This step is necessary because Epetra does not allow insertion of rows owned by other processors.)

    a) Q is stored as a CSR matrix.  We know each row has at most dim(NS) nonzeros.
    For each row of Q, we must transmit both the values and columns.  We do this by communicating separate vectors,
    each of length dim(NS), for the values and column indices.

    b) We try to make this communication more efficient by first communicating just the global row numbers
    that each processor should expect to receive. A "reduced" map is created from just the expected row numbers.

    c) We then communicate the rows of Q themselves using this "reduced" map, i.e., the target of the Importer is the reduced map.
    Otherwise, the only alternative was to base the Importer on A's rowmap, which is probably much larger than the reduced map.

    5) Once received, the rows are inserted by the owning processes and Ptent is fillCompleted.

    - FIXME There is no attempt to detect if the aggregate is too small to support the NS.
    */
    void MakeTentative(const Matrix& fineA, const Aggregates& aggregates, const AmalgamationInfo& amalgInfo, const MultiVector & fineNullspace, RCP<const Map> coarseMap, //-> INPUT
                       RCP<MultiVector> & coarseNullspace, RCP<Matrix> & Ptentative) const;                  //-> OUTPUT

  private:
    RCP<const FactoryBase> aggregatesFact_; //! Factory that creates aggregates
    RCP<const FactoryBase> amalgFact_;      //! Factory that (Un)Amalgamation info from A
    RCP<const FactoryBase> nullspaceFact_;  //! Factory creating the nullspace
    RCP<const FactoryBase> AFact_;          //! Define which matrix A is used in this factory
    RCP<const FactoryBase> coarseMapFact_;  //! coarseMap factory creates coarse map of prolongation transfer operator

    bool QR_; //! use QR decomposition for improving nullspace information per default

    // TODO to be removed: add a domainGidOffset to strided maps
    GlobalOrdinal domainGidOffset_; //! offset for domain gids (coarse gids) of tentative prolongator  (default = 0). The GIDs for the domain dofs of Ptent start with domainGidOffset, are contiguous and distributed equally over the procs (unless some reordering is done).

    mutable std::vector<size_t> stridingInfo_;   // vector with size of strided blocks (dofs)
    LocalOrdinal stridedBlockId_;        // member variable denoting which dofs are stored in map
                                         // stridedBlock == -1: the full map (with all strided block dofs)
                                         // stridedBlock >  -1: only dofs of strided block with index "stridedBlockId" are stored in this map


  }; //class TentativePFactory

} //namespace MueLu

//TODO: noQR_

#define MUELU_TENTATIVEPFACTORY_SHORT
#endif // MUELU_TENTATIVEPFACTORY_DECL_HPP
