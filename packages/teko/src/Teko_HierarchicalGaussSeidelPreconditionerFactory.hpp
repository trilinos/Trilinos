/*
// @HEADER
//
// ***********************************************************************
//
//      Teko: A package for block and physics based preconditioning
//                  Copyright 2010 Sandia Corporation
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
// Questions? Contact Eric C. Cyr (eccyr@sandia.gov)
//
// ***********************************************************************
//
// @HEADER

*/

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

  std::map<int, std::vector<int>> blockToRow;
  std::map<int, Teuchos::RCP<InverseFactory>> blockToInverse;
  std::map<int, Teuchos::RCP<InverseFactory>> blockToPreconditioner;
  mutable std::map<int, LinearOp> blockToInvOp;

  bool useLowerTriangle = false;
};
}  // end namespace Teko

#endif
