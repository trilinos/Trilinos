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

#ifndef __Teko_BlockInvDiagonalStrategy_hpp__
#define __Teko_BlockInvDiagonalStrategy_hpp__

#include <vector>

// Teuchos includes
#include "Teuchos_RCP.hpp"

// Thyra includes
#include "Thyra_LinearOpBase.hpp"

// Teko includes
#include "Teko_Utilities.hpp"
#include "Teko_InverseFactory.hpp"
#include "Teko_BlockPreconditionerFactory.hpp"

namespace Teko {

/** this should be paired with a BlockJacobiPreconditionerFactory
 * or BlockGSPreconditionerFactory.  The idea is that this object
 * provides an approximate inverse operator for each of the diagonal
 * blocks.  Then, the [whatever]PreconditionerFactory can easily
 * construct an approximate inverse. The system under consideration
 * is
 *
 *    A = [ D0  U01 U02 ...]
 *        [ L10  D1 U12 ...]
 *        [ L20 L21  D2 ...]
 *        [  .   .   .  ...]
 *
 * where inverses of D0,D1,D2...DN are needed.
 */
class BlockInvDiagonalStrategy {
 public:
  virtual ~BlockInvDiagonalStrategy() {}

  //! returns an (approximate) inverse of the diagonal blocks of A
  virtual void getInvD(const BlockedLinearOp &A, BlockPreconditionerState &state,
                       std::vector<LinearOp> &invDiag) const = 0;
};

/** This is a simple strategy for a [whatever]PreconditionerFactory
 * it simply returns statically set RCP pointers to the passed in
 * inv(D0) and inv(D1) operators. Not this will _not_ permit
 * efficient implementations when the preconditioner has to be rebuilt
 * or reused often.
 */
class StaticInvDiagStrategy : public BlockInvDiagonalStrategy {
 public:
  StaticInvDiagStrategy(const LinearOp &invD0, const LinearOp &invD1) {
    invDiag_.push_back(invD0);
    invDiag_.push_back(invD1);
  }

  StaticInvDiagStrategy(const std::vector<LinearOp> &invD) : invDiag_(invD) {}

  virtual ~StaticInvDiagStrategy() {}

  /** returns an (approximate) inverse of the diagonal blocks of A
   * where A is closely related to the original source for invD0 and invD1
   */
  virtual void getInvD(const BlockedLinearOp & /* A */, BlockPreconditionerState & /* state */,
                       std::vector<LinearOp> &invDiag) const {
    invDiag.clear();
    invDiag = invDiag_;
  }

 protected:
  // stored inverse operators
  std::vector<Teuchos::RCP<const Thyra::LinearOpBase<double> > > invDiag_;
};

/** A simple class that takes a vector of the <code>InverseFactory</code> objects
 * and pairs each with the diagonal element of the block matrix. This provides
 * the operators needed to use Gauss-Seidel or Jacobi.
 */
class InvFactoryDiagStrategy : public BlockInvDiagonalStrategy {
 public:
  /** Constructor accepting a single inverse factory that will be used
   * to invert all diagonal blocks.
   *
   * \param[in] factory Factory to be used to invert each diagonal block.
   */
  InvFactoryDiagStrategy(const Teuchos::RCP<InverseFactory> &factory);

  /** Constructor that lets the inverse of each block be set individually.
   *
   * \param[in] factories A vector of <code>InverseFactory</code> objects
   *                      which should be the same length as the number of
   *                      diagonal blocks.
   * \param[in] defaultFact The default factory to use if none is specified.
   */
  InvFactoryDiagStrategy(const std::vector<Teuchos::RCP<InverseFactory> > &factories,
                         const Teuchos::RCP<InverseFactory> &defaultFact = Teuchos::null);

  InvFactoryDiagStrategy(
      const std::vector<Teuchos::RCP<InverseFactory> > &inverseFactories,
      const std::vector<Teuchos::RCP<InverseFactory> > &preconditionerFactories,
      const Teuchos::RCP<InverseFactory> &defaultInverseFact        = Teuchos::null,
      const Teuchos::RCP<InverseFactory> &defaultPreconditionerFact = Teuchos::null);

  virtual ~InvFactoryDiagStrategy() {}

  /** returns an (approximate) inverse of the diagonal blocks of A
   * where A is closely related to the original source for invD0 and invD1
   *
   * \param[in]     A       Operator to extract the block diagonals from.
   * \param[in]     state   State object for this operator.
   * \param[in,out] invDiag Vector eventually containing the inverse operators
   */
  virtual void getInvD(const BlockedLinearOp &A, BlockPreconditionerState &state,
                       std::vector<LinearOp> &invDiag) const;

  //! Get factories for testing purposes.
  const std::vector<Teuchos::RCP<InverseFactory> > &getFactories() const { return invDiagFact_; }

 protected:
  // stored inverse operators
  std::vector<Teuchos::RCP<InverseFactory> > invDiagFact_;
  std::vector<Teuchos::RCP<InverseFactory> > precDiagFact_;
  Teuchos::RCP<InverseFactory> defaultInvFact_;
  Teuchos::RCP<InverseFactory> defaultPrecFact_;

  //! Conveinence function for building inverse operators
  LinearOp buildInverse(const InverseFactory &invFact, Teuchos::RCP<InverseFactory> &precFact,
                        const LinearOp &matrix, BlockPreconditionerState &state,
                        const std::string &opPrefix, int i) const;

 private:
  InvFactoryDiagStrategy();
  InvFactoryDiagStrategy(const InvFactoryDiagStrategy &);
};

}  // end namespace Teko

#endif
