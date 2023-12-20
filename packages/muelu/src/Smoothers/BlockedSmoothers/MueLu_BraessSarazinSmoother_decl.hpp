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

#ifndef MUELU_BRAESSSARAZINSMOOTHER_DECL_HPP_
#define MUELU_BRAESSSARAZINSMOOTHER_DECL_HPP_
#include "MueLu_ConfigDefs.hpp"

#include <Teuchos_ParameterList.hpp>

// Xpetra
#include <Xpetra_MapExtractor_fwd.hpp>
#include <Xpetra_MultiVectorFactory_fwd.hpp>
#include <Xpetra_Matrix_fwd.hpp>

// MueLu
#include "MueLu_BraessSarazinSmoother_fwd.hpp"
#include "MueLu_SmootherPrototype.hpp"
#include "MueLu_FactoryBase_fwd.hpp"
#include "MueLu_FactoryManagerBase_fwd.hpp"
#include "MueLu_SmootherBase_fwd.hpp"
#include "MueLu_Utilities_fwd.hpp"

#include "MueLu_FactoryManager_fwd.hpp"

namespace MueLu {

/*!
  @class BraessSarazinSmoother
  @brief BraessSarazin smoother for 2x2 block matrices

*/

template <class Scalar        = SmootherPrototype<>::scalar_type,
          class LocalOrdinal  = typename SmootherPrototype<Scalar>::local_ordinal_type,
          class GlobalOrdinal = typename SmootherPrototype<Scalar, LocalOrdinal>::global_ordinal_type,
          class Node          = typename SmootherPrototype<Scalar, LocalOrdinal, GlobalOrdinal>::node_type>
class BraessSarazinSmoother : public SmootherPrototype<Scalar, LocalOrdinal, GlobalOrdinal, Node> {
#undef MUELU_BRAESSSARAZINSMOOTHER_SHORT
#include "MueLu_UseShortNames.hpp"

 public:
  //! Input
  //@{

  RCP<const ParameterList> GetValidParameterList() const;

  void DeclareInput(Level &currentLevel) const;

  //! Add a factory manager for BraessSarazin internal SchurComplement handling
  void AddFactoryManager(RCP<const FactoryManagerBase> FactManager, int pos = 0);

  //@}

  //! @name Setup and Apply methods.
  //@{

  /*! @brief Setup routine

  Setup routine can be summarized in 4 steps:
  - set the map extractors
  - set the blocks
  - create and set the inverse of the diagonal of F
  - set the smoother for the Schur Complement
  */
  void Setup(Level &currentLevel);

  /*! @brief Apply the Braess Sarazin smoother.
  @param X initial guess
  @param B right-hand side
  @param InitialGuessIsZero TODO This option has no effect.
  */
  void Apply(MultiVector &X, const MultiVector &B, bool InitialGuessIsZero = false) const;
  //@}

  RCP<SmootherPrototype> Copy() const;

  //! @name Overridden from Teuchos::Describable
  //@{

  //! Return a simple one-line description of this object.
  std::string description() const;

  /*!
  @brief Print the object with some verbosity level to an FancyOStream object

  Using MueLu::Describable::describe; overloading, not hiding
  */
  void print(Teuchos::FancyOStream &out, const VerbLevel verbLevel = Default) const;

  //! Get a rough estimate of cost per iteration
  size_t getNodeSmootherComplexity() const;

  //@}

 private:
  //! smoother type
  std::string type_ = "Braess Sarazin";

  RCP<const FactoryBase> AFact_;               //!< A Factory
  RCP<const FactoryManagerBase> FactManager_;  //!< Factory manager for creating the Schur Complement

  //! block operator
  RCP<Matrix> A_ = Teuchos::null;  // < ! internal blocked operator "A" generated by AFact_

  RCP<const MapExtractor> rangeMapExtractor_;   //!< range  map extractor (from A_ generated by AFact)
  RCP<const MapExtractor> domainMapExtractor_;  //!< domain map extractor (from A_ generated by AFact)

  //! matrices
  RCP<Matrix> A00_;  //!< Block (0,0) [typically, fluid operator]
  RCP<Matrix> A01_;  //!< Block (0,1) [typically, pressure gradient operator]
  RCP<Matrix> A10_;  //!< Block (1,0) [typically, divergence operator]
  RCP<Matrix> A11_;  //!< Block (1,1) [typically, pressure stabilization term or null block]
  RCP<Matrix> S_;    //!< Schur complement
  RCP<Vector> D_;    //!< Inverse to approximation to block (0,0). Here, D_ = omega*inv(diag(A(0,0)))

  Teuchos::RCP<SmootherBase> smoo_;  //!< Smoother for SchurComplement equation

};  // class Amesos2Smoother

}  // namespace MueLu

#define MUELU_BRAESSSARAZINSMOOTHER_SHORT
#endif /* MUELU_BRAESSSARAZINSMOOTHER_DECL_HPP_ */
