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

#include "Teko_BlockPreconditionerFactory.hpp"
#include "Teko_Utilities.hpp"
#include "Teko_InverseFactory.hpp"
#include "Teko_BlockLowerTriInverseOp.hpp"
#include "Teko_BlockUpperTriInverseOp.hpp"

using Teuchos::RCP;

// Declaration of the example preconditioner state
class ExamplePreconditionerState /*@ \label{lne2:begin-decl-state} @*/
    : public Teko::BlockPreconditionerState {
 public:
  // Default constructor
  ExamplePreconditionerState() : Teko::BlockPreconditionerState() {}

  // class members
  Teko::LinearOp P_;  // store P = A_00+\alpha A_01

}; /*@ \label{lne2:end-decl-state} @*/

// Declaration of the preconditioner factory
class ExamplePreconditionerFactory /*@ \label{lne2:begin-decl} @*/
    : public Teko::BlockPreconditionerFactory {
 public:
  // Constructor
  ExamplePreconditionerFactory(const RCP<const Teko::InverseFactory> &inverse, double alpha);

  // Function inherited from Teko::BlockPreconditionerFactory
  Teko::LinearOp buildPreconditionerOperator(Teko::BlockedLinearOp &blo,
                                             Teko::BlockPreconditionerState &state) const;

  // Function that returns the correct type of state object
  virtual RCP<Teko::BlockPreconditionerState> buildPreconditionerState() const;

 protected:
  // class members
  RCP<const Teko::InverseFactory> inverse_;
  double alpha_;
  std::string invP_str_;

}; /*@ \label{lne2:end-decl} @*/

// Constructor definition
ExamplePreconditionerFactory /*@ \label{lne2:begin-constructor} @*/
    ::ExamplePreconditionerFactory(const RCP<const Teko::InverseFactory> &inverse, double alpha)
    : inverse_(inverse), alpha_(alpha) {
  // store the string name for retrieving the inverse
  invP_str_ = "invP";
} /*@ \label{lne2:end-constructor} @*/

// Use the factory to build the preconditioner (this is where the work goes)
Teko::LinearOp ExamplePreconditionerFactory /*@ \label{lne2:begin-bpo} @*/
    ::buildPreconditionerOperator(Teko::BlockedLinearOp &blockOp,
                                  Teko::BlockPreconditionerState &state) const {
  int rows = Teko::blockRowCount(blockOp); /*@ \label{lne2:begin-extraction} @*/
  int cols = Teko::blockColCount(blockOp);

  TEUCHOS_ASSERT(rows == 2);  // sanity checks
  TEUCHOS_ASSERT(cols == 2);

  // get the casted version of the state object
  ExamplePreconditionerState &exampleState /*@ \label{lne2:state-cast} @*/
      = dynamic_cast<ExamplePreconditionerState &>(state);

  // extract subblocks
  const Teko::LinearOp A_00 = Teko::getBlock(0, 0, blockOp);
  const Teko::LinearOp A_01 = Teko::getBlock(0, 1, blockOp);
  const Teko::LinearOp A_10 = Teko::getBlock(1, 0, blockOp);
  const Teko::LinearOp A_11 = Teko::getBlock(1, 1, blockOp); /*@ \label{lne2:end-extraction} @*/

  // get inverse of diag(A11)
  const Teko::LinearOp invH = Teko::getInvDiagonalOp(A_11); /*@ \label{lne2:invH} @*/

  // build or rebuild inverse P /*@ \label{lne2:invP} @*/
  Teko::InverseLinearOp invP = exampleState.getInverse(invP_str_);
  if (invP == Teuchos::null) {
    // build 0,0 block in the preconditioner
    exampleState.P_ = Teko::explicitAdd(A_00, Teko::scale(alpha_, A_01)); /*@ \label{lne2:P} @*/

    invP = Teko::buildInverse(*inverse_, exampleState.P_);  // build inverse P
    exampleState.addInverse(invP_str_, invP);               // add inverse operator to state
  }

  // build lower triangular inverse matrix
  Teko::BlockedLinearOp L = Teko::zeroBlockedOp(blockOp); /*@ \label{lne2:begin-trisolve} @*/
  Teko::setBlock(1, 0, L, A_10);
  Teko::endBlockFill(L);

  std::vector<Teko::LinearOp> invDiag(
      2);  // vector storing inverses /*@ \label{lne2:begin-invdiags} @*/
  invDiag[0] = invP;
  invDiag[1] = invH; /*@ \label{lne2:end-invdiags} @*/

  Teko::LinearOp invTildeA =
      Teko::createBlockLowerTriInverseOp(L, invDiag); /*@ \label{lne2:invLower} @*/

  // tell the state object it has been initialized for this operator
  exampleState.setInitialized(true);

  // return fully constructed preconditioner
  return invTildeA; /*@ \label{lne2:end-trisolve} @*/
} /*@ \label{lne2:end-bpo} @*/

// Function that returns the correct type of state object
RCP<Teko::BlockPreconditionerState> ExamplePreconditionerFactory /*@ \label{lne2:begin-bps} @*/
    ::buildPreconditionerState() const {
  // build the state object
  return Teuchos::rcp(new ExamplePreconditionerState());

} /*@ \label{lne2:end-bps} @*/
