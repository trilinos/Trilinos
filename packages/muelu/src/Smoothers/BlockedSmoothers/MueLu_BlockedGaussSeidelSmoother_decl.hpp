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

#ifndef MUELU_BLOCKEDGAUSSSEIDELSMOOTHER_DECL_HPP_
#define MUELU_BLOCKEDGAUSSSEIDELSMOOTHER_DECL_HPP_

#include "MueLu_ConfigDefs.hpp"

#include <Teuchos_ParameterList.hpp>

#include <Xpetra_Matrix_fwd.hpp>
#include <Xpetra_BlockedCrsMatrix_fwd.hpp>
#include <Xpetra_ReorderedBlockedCrsMatrix_fwd.hpp>
#include <Xpetra_MultiVectorFactory_fwd.hpp>
#include <Xpetra_MapExtractor_fwd.hpp>

#include "MueLu_BlockedGaussSeidelSmoother_fwd.hpp"

#include "MueLu_SmootherPrototype.hpp"
#include "MueLu_FactoryBase_fwd.hpp"
#include "MueLu_FactoryManagerBase_fwd.hpp"
#include "MueLu_SmootherBase_fwd.hpp"

namespace MueLu {

/*!
  @class BlockedGaussSeidelSmoother
  @brief block Gauss-Seidel method for blocked matrices

  Implementation of a block Gauss-Seidel methods for blocked matrices
  @param LocalOrdinal sweeps = 1: number of BGS sweeps
  @param Scalar omega = 1.0: damping parameter
  @param RCP<FactoryBase> AFact = Teuchos::null: factory for blocked "A" operator

  Use the AddFactoryManager routine to declare the subsmoothers/subsolvers for the block Gauss-Seidel method for
  the block rows. The corresponding factory manager has to provide a variable "A" (pointing to the subblock of the blocked A operator)
  and a smoother object (variable: "PreSmoother").

  Example
  \code

  // prototypes for direct solvers for blocks 1 and 2
  RCP<SmootherPrototype> smoProto11 = rcp(new DirectSolver("", Teuchos::ParameterList(), A11Fact));
  RCP<SmootherPrototype> smoProto22 = rcp(new DirectSolver("", Teuchos::ParameterList(), A22Fact));
  RCP<SmootherFactory> Smoo11Fact = rcp(new SmootherFactory(smoProto11));
  RCP<SmootherFactory> Smoo22Fact = rcp(new SmootherFactory(smoProto22));

  // define factory manager objects for sublocks
  RCP<FactoryManager> M11 = rcp(new FactoryManager());
  M11->SetFactory("A", A11Fact);
  M11->SetFactory("Smoother", Smoo11Fact);

  RCP<FactoryManager> M22 = rcp(new FactoryManager());
  M22->SetFactory("A", A22Fact);
  M22->SetFactory("Smoother", Smoo22Fact);

  // create blocked Gauss-Seidel smoother for 2x2 blocked matrix
  RCP<BlockedGaussSeidelSmoother> smootherPrototype = rcp(new BlockedGaussSeidelSmoother(2,1.0));
  smootherPrototype->AddFactoryManager(M11);
  smootherPrototype->AddFactoryManager(M22);
  RCP<SmootherFactory> smootherFact = rcp( new SmootherFactory(smootherPrototype) );

  // use smootherFact in main-factory manager
  \endcode
*/

template <class Scalar        = SmootherPrototype<>::scalar_type,
          class LocalOrdinal  = typename SmootherPrototype<Scalar>::local_ordinal_type,
          class GlobalOrdinal = typename SmootherPrototype<Scalar, LocalOrdinal>::global_ordinal_type,
          class Node          = typename SmootherPrototype<Scalar, LocalOrdinal, GlobalOrdinal>::node_type>
class BlockedGaussSeidelSmoother : public SmootherPrototype<Scalar, LocalOrdinal, GlobalOrdinal, Node> {
  typedef Xpetra::MapExtractor<Scalar, LocalOrdinal, GlobalOrdinal, Node> MapExtractorClass;

#undef MUELU_BLOCKEDGAUSSSEIDELSMOOTHER_SHORT
#include "MueLu_UseShortNames.hpp"

 public:
  //! @name Constructors / destructors
  //@{

  /*! @brief Constructor
   */
  BlockedGaussSeidelSmoother();

  //! Destructor
  virtual ~BlockedGaussSeidelSmoother();
  //@}

  //! Input
  //@{
  RCP<const ParameterList> GetValidParameterList() const;

  void DeclareInput(Level &currentLevel) const;

  //! Add a factory manager
  // void AddFactoryManager(RCP<const FactoryManagerBase> FactManager);

  //! Add a factory manager at a specific position
  void AddFactoryManager(RCP<const FactoryManagerBase> FactManager, int pos);

  //@}

  //! @name Setup and Apply methods.
  //@{

  /*! @brief Setup routine
   * In the Setup method the Inverse_ vector is filled with the corresponding
   * SmootherBase objects. Without the Inverse_ vector being filled we cannot call
   * BlockedGaussSeidelSmoother::Apply.
   */
  void Setup(Level &currentLevel);

  /*! @brief Apply the direct solver.
  Solves the linear system <tt>AX=B</tt> using the constructed solver.
  @param X initial guess
  @param B right-hand side
  @param InitialGuessIsZero This option has no effect.
  */
  void Apply(MultiVector &X, const MultiVector &B, bool InitialGuessIsZero = false) const;
  //@}

  RCP<SmootherPrototype> Copy() const;

  //! @name Overridden from Teuchos::Describable
  //@{

  //! Return a simple one-line description of this object.
  std::string description() const;

  /*! \brief Print the object with some verbosity level \c verbLevel to an FancyOStream object \c out
   *
   * - use MueLu::Describable::describe;
   * - overloading, not hiding
   */
  void print(Teuchos::FancyOStream &out, const VerbLevel verbLevel = Default) const;

  //! Get a rough estimate of cost per iteration
  size_t getNodeSmootherComplexity() const;

  //@}

 private:
  //! smoother type
  std::string type_;

  //! vector of factory managers
  std::vector<Teuchos::RCP<const FactoryManagerBase> > FactManager_;

  //! vector of smoother/solver factories
  std::vector<Teuchos::RCP<const SmootherBase> > Inverse_;

  //! vector storing whether sub-block is a blocked operator (needed for nested blocked smoothers using Thyra GIDs)
  std::vector<bool> bIsBlockedOperator_;

  //! A Factory
  RCP<FactoryBase> AFact_;

  //! internal blocked operator "A" generated by AFact_
  RCP<Matrix> A_;

  //! range  map extractor (from A_ generated by AFact)
  RCP<const MapExtractorClass> rangeMapExtractor_;

  //! domain map extractor (from A_ generated by AFact)
  RCP<const MapExtractorClass> domainMapExtractor_;
};  // class Amesos2Smoother

}  // namespace MueLu

#define MUELU_BLOCKEDGAUSSSEIDELSMOOTHER_SHORT

#endif /* MUELU_BLOCKEDGAUSSSEIDELSMOOTHER_DECL_HPP_ */
