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
 * MueLu_IndefBlockedDiagonalSmoother_decl.hpp
 *
 *  Created on: 13 May 2014
 *      Author: wiesner
 */

#ifndef MUELU_INDEFBLOCKEDDIAGONALSMOOTHER_DECL_HPP_
#define MUELU_INDEFBLOCKEDDIAGONALSMOOTHER_DECL_HPP_

#include "MueLu_ConfigDefs.hpp"

#include <Teuchos_ParameterList.hpp>

// Xpetra
#include <Xpetra_MapExtractor_fwd.hpp>
#include <Xpetra_MultiVectorFactory_fwd.hpp>
#include <Xpetra_Matrix_fwd.hpp>

// MueLu
#include "MueLu_IndefBlockedDiagonalSmoother_fwd.hpp"
#include "MueLu_SmootherPrototype.hpp"
#include "MueLu_FactoryBase_fwd.hpp"
#include "MueLu_FactoryManagerBase_fwd.hpp"
#include "MueLu_SmootherBase_fwd.hpp"

#include "MueLu_SchurComplementFactory_fwd.hpp"
#include "MueLu_FactoryManager_fwd.hpp"

namespace MueLu {

/*!
  @class IndefBlockedDiagonalSmoother
  @brief Cheap Blocked diagonal smoother for indefinite 2x2 block matrices

  Uses the original upper left block and the Schur Complement block on the diagonal blocks.
  Instead of solving the block equations exactly, we apply some sweeps with cheap iterative smoothers
  (e.g. Gauss-Seidel iterations)

*/

template <class Scalar        = SmootherPrototype<>::scalar_type,
          class LocalOrdinal  = typename SmootherPrototype<Scalar>::local_ordinal_type,
          class GlobalOrdinal = typename SmootherPrototype<Scalar, LocalOrdinal>::global_ordinal_type,
          class Node          = typename SmootherPrototype<Scalar, LocalOrdinal, GlobalOrdinal>::node_type>
class IndefBlockedDiagonalSmoother : public SmootherPrototype<Scalar, LocalOrdinal, GlobalOrdinal, Node> {
  typedef Xpetra::MapExtractor<Scalar, LocalOrdinal, GlobalOrdinal, Node> MapExtractorClass;

#undef MUELU_INDEFBLOCKEDDIAGONALSMOOTHER_SHORT
#include "MueLu_UseShortNames.hpp"

 public:
  //! @name Constructors / destructors
  //@{

  /*! @brief Constructor
   */
  IndefBlockedDiagonalSmoother();

  //! Destructor
  virtual ~IndefBlockedDiagonalSmoother();
  //@}

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

  //! Print the object with some verbosity level to an FancyOStream object.
  // using MueLu::Describable::describe; // overloading, not hiding
  void print(Teuchos::FancyOStream &out, const VerbLevel verbLevel = Default) const;

  //! Get a rough estimate of cost per iteration
  size_t getNodeSmootherComplexity() const;

  //@}

 private:
  //! smoother type
  std::string type_;

  RCP<const FactoryBase> AFact_;  //!< A Factory

  //! block operator
  RCP<Matrix> A_;  // < ! internal blocked operator "A" generated by AFact_
  RCP<Matrix> F_;  //!< fluid operator
  RCP<Matrix> Z_;  //!< pressure stabilization term or null block

  RCP<const MapExtractorClass> rangeMapExtractor_;   //!< range  map extractor (from A_ generated by AFact)
  RCP<const MapExtractorClass> domainMapExtractor_;  //!< domain map extractor (from A_ generated by AFact)

  //! Block smoothers
  Teuchos::RCP<SmootherBase> velPredictSmoo_;  //!< smoother for velocity prediction
  Teuchos::RCP<SmootherBase> schurCompSmoo_;   //!< smoother for SchurComplement equation

  //! vector of factory managers
  /*!
   * vector of internal factory managers
   * FactManager_[0] holds the factory manager for the predicting the primary variable
   * FactManager_[1] stores the factory manager used for the SchurComplement correction step.
   */
  std::vector<Teuchos::RCP<const FactoryManagerBase> > FactManager_;
};

}  // namespace MueLu

#define MUELU_INDEFBLOCKEDDIAGONALSMOOTHER_SHORT

#endif /* MUELU_INDEFBLOCKEDDIAGONALSMOOTHER_DECL_HPP_ */
