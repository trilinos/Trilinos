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
//                    Andrey Prokopenko (aprokop@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
#ifndef MUELU_BLOCKEDRAPFACTORY_DECL_HPP
#define MUELU_BLOCKEDRAPFACTORY_DECL_HPP

#include <Xpetra_Matrix_fwd.hpp>
#include <Xpetra_CrsMatrix_fwd.hpp>
#include <Xpetra_CrsMatrixWrap_fwd.hpp>
#include <Xpetra_BlockedCrsMatrix_fwd.hpp>
#include <Xpetra_MatrixFactory_fwd.hpp>
#include <Xpetra_Vector_fwd.hpp>
#include <Xpetra_VectorFactory_fwd.hpp>

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_RAPFactory_fwd.hpp"

#include "MueLu_Level_fwd.hpp"
#include "MueLu_FactoryBase_fwd.hpp"
#include "MueLu_PerfUtils_fwd.hpp"
#include "MueLu_TwoLevelFactoryBase.hpp"
#include "MueLu_Utilities_fwd.hpp"

namespace MueLu {
  /*!
    @class BlockedRAPFactory
    @brief Factory for building coarse matrices.
  */
  template <class Scalar = double, class LocalOrdinal = int, class GlobalOrdinal = LocalOrdinal, class Node = KokkosClassic::DefaultNode::DefaultNodeType>
  class BlockedRAPFactory : public TwoLevelFactoryBase {
#undef MUELU_BLOCKEDRAPFACTORY_SHORT
#include "MueLu_UseShortNames.hpp"

  public:
    //! @name Constructors/Destructors.
    //@{

    BlockedRAPFactory();

    virtual ~BlockedRAPFactory() { }
    //@}

    //! @name Input
    //@{

    RCP<const ParameterList> GetValidParameterList() const;

    void DeclareInput(Level &fineLevel, Level &coarseLevel) const;

    //@}

    //! @name Build methods.
    //@{
    void Build(Level &fineLevel, Level &coarseLevel) const;
    //@}

    //! @name Handling of user-defined transfer factories
    //@{

    //! Indicate that zero entries on the diagonal of Ac shall be repaired (i.e. if A(i,i) == 0.0 set A(i,i) = 1.0)
    void SetRepairZeroDiagonal(bool const &repair) {
      repairZeroDiagonals_ = repair;
      if(repair) checkAc_ = true; // make sure that plausibility check is performed. Otherwise SetRepairZeroDiagonal(true) has no effect.
    }

    //! Indicate that a simple plausibility check shall be done for Ac after building RAP
    void SetPlausibilityCheck(bool const &check) {
      checkAc_ = check;
    }

    //@{
    /*! @brief Add transfer factory in the end of list of transfer factories in RepartitionAcFactory.

    Transfer factories are derived from TwoLevelFactoryBase and project some data from the fine level to
    the next coarser level.
    */
    void AddTransferFactory(const RCP<const FactoryBase>& factory);

    // TODO add a function to remove a specific transfer factory?

    //! Returns number of transfer factories.
    size_t NumTransferFactories() const { return transferFacts_.size(); }

    //@}

  private:

    //! @name internal plausibility check methods
    //! checks main diagonal entries of (0,0) block. Does not affect entries in (1,1) block!
    static void CheckMainDiagonal(RCP<BlockedCrsMatrix> & bAc, bool repairZeroDiagonals = false);

    //! If true, perform a basic plausibility check on Ac (default = false)
    //! note, that the repairZeroDiagonals_ flag only is valid for checkAc_ == true
    bool checkAc_;

    //! If true, the CheckMainDiagonal routine automatically repairs zero entries on main diagonal (default = false)
    //! i.e. if A(i,i) == 0.0 set A(i,i) = 1.0
    //! note, that the repairZeroDiagonals_ flag only is valid for checkAc_ == true
    bool repairZeroDiagonals_;

    //@{

    //! list of user-defined transfer Factories
    std::vector<RCP<const FactoryBase> > transferFacts_;

    //@}

  }; //class BlockedRAPFactory

} //namespace MueLu

#define MUELU_BLOCKEDRAPFACTORY_SHORT
#endif // MUELU_BLOCKEDRAPFACTORY_DECL_HPP
