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
#ifndef MUELU_REFMAXWELLSMOOTHER_DECL_HPP
#define MUELU_REFMAXWELLSMOOTHER_DECL_HPP

#include <Teuchos_ParameterList.hpp>

#include <Xpetra_Matrix_fwd.hpp>

#include "MueLu_ConfigDefs.hpp"

#include "MueLu_RefMaxwellSmoother_fwd.hpp"

#include "MueLu_Level_fwd.hpp"
#include "MueLu_SmootherPrototype.hpp"
#include "MueLu_RefMaxwell.hpp"

namespace MueLu {

/*!
  @class RefMaxwellSmoother
  @ingroup MueLuSmootherClasses
  @brief Class that encapsulates Operator smoothers.

  */

template <class Scalar        = SmootherPrototype<>::scalar_type,
          class LocalOrdinal  = typename SmootherPrototype<Scalar>::local_ordinal_type,
          class GlobalOrdinal = typename SmootherPrototype<Scalar, LocalOrdinal>::global_ordinal_type,
          class Node          = typename SmootherPrototype<Scalar, LocalOrdinal, GlobalOrdinal>::node_type>
class RefMaxwellSmoother : public SmootherPrototype<Scalar, LocalOrdinal, GlobalOrdinal, Node> {
#undef MUELU_REFMAXWELLSMOOTHER_SHORT
#include "MueLu_UseShortNames.hpp"

 public:
  //! @name Constructors / destructors
  //@{
  /*! @brief Constructor

  @param type smoother type
  @param list options for the particular smoother (e.g., fill factor or damping parameter)

  */

  RefMaxwellSmoother(const std::string type, const Teuchos::ParameterList& paramList);

  //! Destructor
  virtual ~RefMaxwellSmoother() = default;

  //@}

  void SetParameterList(const Teuchos::ParameterList& paramList);

  //! Input
  //@{

  void DeclareInput(Level& currentLevel) const;

  //@}

  //! @name Computational methods.
  //@{

  /*! @brief Set up the smoother.

  This creates the underlying Operator smoother object, copies any parameter list options
  supplied to the constructor to the Operator object, and computes the preconditioner.
  */
  void Setup(Level& currentLevel);

  /*! @brief Apply the preconditioner.

  Solves the linear system <tt>AX=B</tt> using the constructed smoother.

  @param X initial guess
  @param B right-hand side
  @param InitialGuessIsZero (optional) If false, some work can be avoided. Whether this actually saves any work depends on the underlying Operator implementation.
  */
  void Apply(MultiVector& X, const MultiVector& B, bool InitialGuessIsZero = false) const;

  //@}

  //! @name Utilities
  //@{

  RCP<SmootherPrototype> Copy() const;

  //@}

  //! @name Overridden from Teuchos::Describable
  //@{

  //! Return a simple one-line description of this object.
  std::string description() const;

  //! Print the object with some verbosity level to an FancyOStream object.
  // using MueLu::Describable::describe; // overloading, not hiding
  // void describe(Teuchos::FancyOStream &out, const VerbLevel verbLevel = Default) const
  void print(Teuchos::FancyOStream& out, const VerbLevel verbLevel = Default) const;

  //@}

  //! Get a rough estimate of cost per iteration
  size_t getNodeSmootherComplexity() const;

 private:
  void SetupRefMaxwell(Level& currentLevel);

 private:
  std::string type_;
  std::string cachedDescription_;

  //! matrix, used in apply if solving residual equation
  RCP<Operator> op_;

};  // class RefMaxwellSmoother

}  // namespace MueLu

#define MUELU_REFMAXWELLSMOOTHER_SHORT
#endif  // MUELU_REFMAXWELLSMOOTHER_DECL_HPP
