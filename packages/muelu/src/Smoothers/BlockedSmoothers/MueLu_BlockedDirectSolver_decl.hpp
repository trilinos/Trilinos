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
/*
 * MueLu_BlockedDirectSolver_decl.hpp
 *
 *  Created on: 09.02.2014
 *      Author: tobias
 */

#ifndef MUELU_BLOCKEDDIRECTSOLVER_DECL_HPP_
#define MUELU_BLOCKEDDIRECTSOLVER_DECL_HPP_

#include "MueLu_ConfigDefs.hpp"

#include <Teuchos_ParameterList.hpp>

#include <Xpetra_MapExtractor_fwd.hpp>

#include "MueLu_BlockedDirectSolver_fwd.hpp"
#include "MueLu_SmootherPrototype.hpp"
#include "MueLu_DirectSolver_fwd.hpp"
#include "MueLu_MergedBlockedMatrixFactory_fwd.hpp"
#include "MueLu_FactoryBase_fwd.hpp"

namespace MueLu {

/*!
  @class BlockedDirectSolver
  @brief direct solver for nxn blocked matrices

  The nxn block matrix A as input is automatically merged and then solved
  by a direct solver.
    */

template <class Scalar        = SmootherPrototype<>::scalar_type,
          class LocalOrdinal  = typename SmootherPrototype<Scalar>::local_ordinal_type,
          class GlobalOrdinal = typename SmootherPrototype<Scalar, LocalOrdinal>::global_ordinal_type,
          class Node          = typename SmootherPrototype<Scalar, LocalOrdinal, GlobalOrdinal>::node_type>
class BlockedDirectSolver : public SmootherPrototype<Scalar, LocalOrdinal, GlobalOrdinal, Node> {
  typedef Xpetra::MapExtractor<Scalar, LocalOrdinal, GlobalOrdinal, Node> MapExtractorClass;

#undef MUELU_BLOCKEDDIRECTSOLVER_SHORT
#include "MueLu_UseShortNames.hpp"

 public:
  //! @name Constructors / destructors
  //@{

  /*! @brief Constructor
   */
  BlockedDirectSolver(const std::string& type = "", const Teuchos::ParameterList& paramList = Teuchos::ParameterList());

  //! Destructor
  virtual ~BlockedDirectSolver() {}
  //@}

  //! Input
  //@{
  RCP<const ParameterList> GetValidParameterList() const;

  void DeclareInput(Level& currentLevel) const;
  //@}

  //! @name Setup and Apply methods.
  //@{

  /*! @brief Setup routine
   * Call the underlaying Setup routine of the nested direct solver once the input block matrix has been merged
   */
  void Setup(Level& currentLevel);

  /*! @brief Apply the direct solver.
  Solves the linear system <tt>AX=B</tt> using the constructed solver.
  @param X initial guess
  @param B right-hand side
  @param InitialGuessIsZero This option has no effect.
  */
  void Apply(MultiVector& X, const MultiVector& B, bool InitialGuessIsZero = false) const;
  //@}

  RCP<SmootherPrototype> Copy() const;

  //! @name Overridden from Teuchos::Describable
  //@{

  //! Return a simple one-line description of this object.
  std::string description() const;

  //! Print the object with some verbosity level to an FancyOStream object.
  // using MueLu::Describable::describe; // overloading, not hiding
  void print(Teuchos::FancyOStream& out, const VerbLevel verbLevel = Default) const;

  //! Get a rough estimate of cost per iteration
  size_t getNodeSmootherComplexity() const;

  //@}

 private:
  //! smoother type
  std::string type_;

  //! Factory to generate merged block matrix
  RCP<MergedBlockedMatrixFactory> MergedAFact_;

  //! Direct solver
  RCP<DirectSolver> s_;  // solver object

  //! block operator
  RCP<Matrix> A_;  // < ! internal blocked operator "A" generated by AFact_

};  // class BlockedDirectSolver

}  // namespace MueLu

#define MUELU_BLOCKEDDIRECTSOLVER_SHORT

#endif /* MUELU_BLOCKEDDIRECTSOLVER_DECL_HPP_ */
