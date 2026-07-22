// @HEADER
// *****************************************************************************
//      Teko: A package for block and physics based preconditioning
//
// Copyright 2010 NTESS and the Teko contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __Teko_HierarchicalGaussSeidelPreconditionerFactory_hpp__
#define __Teko_HierarchicalGaussSeidelPreconditionerFactory_hpp__

#include "Teuchos_RCP.hpp"

#include "Teko_BlockPreconditionerFactory.hpp"
#include "Teko_BlockImplicitLinearOp.hpp"
#include "Teko_Utilities.hpp"
#include <map>
#include <vector>

namespace Teko {

class NestedBlockGS : public BlockImplicitLinearOp {
 public:
  NestedBlockGS(const std::map<int, std::vector<int>>& blockToRow_,
                const std::map<int, LinearOp>& blockToInvOp_, BlockedLinearOp& A_,
                bool useLowerTriangle_ = false);

  VectorSpace range() const override { return productRange_; }
  VectorSpace domain() const override { return productDomain_; }

  void implicitApply(const BlockedMultiVector& r, BlockedMultiVector& z, const double alpha = 1.0,
                     const double beta = 0.0) const override;

 private:
  void upperTriangularImplicitApply(std::vector<BlockedMultiVector>& r,
                                    std::vector<BlockedMultiVector>& z, const double alpha = 1.0,
                                    const double beta = 0.0) const;

  void lowerTriangularImplicitApply(std::vector<BlockedMultiVector>& r,
                                    std::vector<BlockedMultiVector>& z, const double alpha = 1.0,
                                    const double beta = 0.0) const;

  // block operators
  std::map<int, std::vector<int>> blockToRow;
  std::map<int, LinearOp> blockToInvOp;
  std::vector<LinearOp> invOps;
  BlockedLinearOp A;
  std::vector<BlockedLinearOp> Ab;
  bool useLowerTriangle = false;

  Teuchos::RCP<const Thyra::ProductVectorSpaceBase<double>> productRange_;  ///< Range vector space.
  Teuchos::RCP<const Thyra::ProductVectorSpaceBase<double>>
      productDomain_;  ///< Domain vector space.
};

class HierarchicalGaussSeidelPreconditionerFactory : public BlockPreconditionerFactory {
 public:
  ~HierarchicalGaussSeidelPreconditionerFactory() override = default;
  HierarchicalGaussSeidelPreconditionerFactory();

  LinearOp buildPreconditionerOperator(BlockedLinearOp& blo,
                                       BlockPreconditionerState& state) const override;

 protected:
  void initializeFromParameterList(const Teuchos::ParameterList& pl) override;
  using BlockPreconditionerFactory::buildPreconditionerOperator;

 private:
  LinearOp buildBlockInverse(const InverseFactory& invFact,
                             const Teuchos::RCP<InverseFactory>& precFact,
                             const BlockedLinearOp& matrix, BlockPreconditionerState& state,
                             int hierarchicalBlockNum) const;
  template <typename LinearOpType>
  LinearOp buildInverseImpl(const InverseFactory& invFact,
                            const Teuchos::RCP<InverseFactory>& precFact,
                            const LinearOpType& matrix, BlockPreconditionerState& state,
                            int hierarchicalBlockNum) const {
    std::stringstream ss;
    ss << "hierarchical_block_" << hierarchicalBlockNum;

    ModifiableLinearOp& invOp  = state.getModifiableOp(ss.str());
    ModifiableLinearOp& precOp = state.getModifiableOp("prec_" + ss.str());

    if (precFact != Teuchos::null) {
      if (precOp == Teuchos::null) {
        precOp = precFact->buildInverse(matrix);
        state.addModifiableOp("prec_" + ss.str(), precOp);
      } else {
        Teko::rebuildInverse(*precFact, matrix, precOp);
      }
    }

    if (invOp == Teuchos::null)
      if (precOp.is_null())
        invOp = Teko::buildInverse(invFact, matrix);
      else
        invOp = Teko::buildInverse(invFact, matrix, precOp);
    else {
      if (precOp.is_null())
        Teko::rebuildInverse(invFact, matrix, invOp);
      else
        Teko::rebuildInverse(invFact, matrix, precOp, invOp);
    }

    return invOp;
  }

  std::map<int, std::vector<int>> blockToRow;
  std::map<int, Teuchos::RCP<InverseFactory>> blockToInverse;
  std::map<int, Teuchos::RCP<InverseFactory>> blockToPreconditioner;
  mutable std::map<int, LinearOp> blockToInvOp;

  bool useLowerTriangle = false;
};
}  // end namespace Teko

#endif
