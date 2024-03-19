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
#ifndef MUELU_DIRECTSOLVER_DECL_HPP
#define MUELU_DIRECTSOLVER_DECL_HPP

#include <Teuchos_ParameterList.hpp>

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_DirectSolver_fwd.hpp"

#include "MueLu_SmootherPrototype.hpp"

#include "MueLu_FactoryBase_fwd.hpp"
#include "MueLu_Amesos2Smoother_fwd.hpp"
#include "MueLu_AmesosSmoother_fwd.hpp"
#include "MueLu_BelosSmoother_fwd.hpp"
#include "MueLu_StratimikosSmoother_fwd.hpp"
#include "MueLu_RefMaxwellSmoother_fwd.hpp"

// Note: DirectSolver is a SmootherPrototype that cannot be turned into a smoother using Setup().
//       When this prototype is cloned using Copy(), the clone is an Amesos or an Amesos2 smoother.
//       The clone can be used as a smoother after calling Setup().

namespace MueLu {

/*!
  @class DirectSolver
  @ingroup MueLuSmootherClasses
  @brief Class that encapsulates direct solvers. Autoselection of AmesosSmoother or Amesos2Smoother according to the compile time configuration of Trilinos
*/

template <class Scalar        = SmootherPrototype<>::scalar_type,
          class LocalOrdinal  = typename SmootherPrototype<Scalar>::local_ordinal_type,
          class GlobalOrdinal = typename SmootherPrototype<Scalar, LocalOrdinal>::global_ordinal_type,
          class Node          = typename SmootherPrototype<Scalar, LocalOrdinal, GlobalOrdinal>::node_type>
class DirectSolver : public SmootherPrototype<Scalar, LocalOrdinal, GlobalOrdinal, Node> {
#undef MUELU_DIRECTSOLVER_SHORT
#include "MueLu_UseShortNames.hpp"

 public:
  //! @name Constructors / destructors
  //@{

  //! @brief Constructor
  //! Note: only parameters shared by Amesos and Amesos2 should be used for type and paramList (example: type= "Klu", "Superlu", paramList = <empty>) .
  DirectSolver(const std::string& type = "", const Teuchos::ParameterList& paramList = Teuchos::ParameterList());

  //! Destructor
  virtual ~DirectSolver() {}

  //@}

  //! Input
  //@{

  void DeclareInput(Level& currentLevel) const;

  //@}

  //! @name Setup and Apply methods.
  //@{

  //! DirectSolver cannot be turned into a smoother using Setup(). Setup() always returns a RuntimeError exception.
  void Setup(Level& currentLevel);

  //! DirectSolver cannot be applied. Apply() always returns a RuntimeError exception.
  void Apply(MultiVector& X, const MultiVector& B, bool InitialGuessIsZero = false) const;

  //@}

  //! Custom SetFactory
  void SetFactory(const std::string& varName, const RCP<const FactoryBase>& factory);

  //! When this prototype is cloned using Copy(), the clone is an Amesos or an Amesos2 smoother.
  RCP<SmootherPrototype> Copy() const;

  //! @name Overridden from Teuchos::Describable
  //@{

  //! Return a simple one-line description of this object.
  std::string description() const;

  void print(Teuchos::FancyOStream& out, const VerbLevel verbLevel = Default) const;

  //! Get a rough estimate of cost per iteration
  size_t getNodeSmootherComplexity() const {
    if (!s_.is_null())
      return s_->getNodeSmootherComplexity();
    else
      return 0.0;
  }

  //@}

 private:
  //
  // Parameters
  //

  //! amesos1/2-specific key phrase that denote smoother type
  std::string type_;

  //
  // Underlying Smoother
  //

  //! Smoother
  RCP<SmootherPrototype> sEpetra_, sTpetra_, sBelos_, sStratimikos_, sRefMaxwell_;
  mutable RCP<SmootherPrototype> s_;

  // Records for the case if something goes wrong
  bool triedEpetra_, triedTpetra_, triedBelos_, triedStratimikos_, triedRefMaxwell_;
  std::string errorEpetra_, errorTpetra_, errorBelos_, errorStratimikos_, errorRefMaxwell_;

};  // class DirectSolver

}  // namespace MueLu

#define MUELU_DIRECTSOLVER_SHORT
#endif  // MUELU_DIRECTSOLVER_DECL_HPP
