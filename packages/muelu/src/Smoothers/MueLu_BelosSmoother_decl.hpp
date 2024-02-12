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
#ifndef MUELU_BELOSSMOOTHER_DECL_HPP
#define MUELU_BELOSSMOOTHER_DECL_HPP

#include <Teuchos_ParameterList.hpp>

#include <Xpetra_Matrix_fwd.hpp>
#include <Xpetra_MultiVectorFactory_fwd.hpp>

#include "MueLu_ConfigDefs.hpp"

#if defined(HAVE_MUELU_BELOS)

#include "MueLu_BelosSmoother_fwd.hpp"

#ifdef HAVE_XPETRA_TPETRA
#include "BelosTpetraAdapter.hpp"
#include <Tpetra_Operator.hpp>
#include <Tpetra_MultiVector.hpp>
#endif

#include "MueLu_Level_fwd.hpp"
#include "MueLu_SmootherPrototype.hpp"
#include "MueLu_Utilities_fwd.hpp"

#include <BelosSolverManager.hpp>

namespace MueLu {

/*!
  @class BelosSmoother
  @ingroup MueLuSmootherClasses
  @brief Class that encapsulates Belos smoothers.

  This class creates an Belos preconditioner factory. The factory creates a smoother based
  on the type and ParameterList passed into the constructor. See the constructor for more
  information.
  */

template <class Scalar        = SmootherPrototype<>::scalar_type,
          class LocalOrdinal  = typename SmootherPrototype<Scalar>::local_ordinal_type,
          class GlobalOrdinal = typename SmootherPrototype<Scalar, LocalOrdinal>::global_ordinal_type,
          class Node          = typename SmootherPrototype<Scalar, LocalOrdinal, GlobalOrdinal>::node_type>
class BelosSmoother : public SmootherPrototype<Scalar, LocalOrdinal, GlobalOrdinal, Node> {
#undef MUELU_BELOSSMOOTHER_SHORT
#include "MueLu_UseShortNames.hpp"

 public:
  //! @name Constructors / destructors
  //@{
  // TODO: update doc for Belos.
  /*! @brief Constructor

  The options passed into BelosSmoother are those given in the Belos user's manual.

  @param type smoother type
  @param list options for the particular smoother (e.g., fill factor or damping parameter)

  */

#ifndef _MSC_VER
  // Avoid error C3772: invalid friend template declaration
  template <class Scalar2, class LocalOrdinal2, class GlobalOrdinal2, class Node2>
  friend class BelosSmoother;
#endif

  BelosSmoother(const std::string type, const Teuchos::ParameterList& paramList = Teuchos::ParameterList());

  //! Destructor
  virtual ~BelosSmoother() = default;

  //@}

  void SetParameterList(const Teuchos::ParameterList& paramList);

  //! Input
  //@{

  void DeclareInput(Level& currentLevel) const;

  //@}

  //! @name Computational methods.
  //@{

  /*! @brief Set up the smoother.

  This creates the underlying Belos smoother object, copies any parameter list options
  supplied to the constructor to the Belos object, and computes the preconditioner.
  */
  void Setup(Level& currentLevel);

  /*! @brief Apply the preconditioner.

  Solves the linear system <tt>AX=B</tt> using the constructed smoother.

  @param X initial guess
  @param B right-hand side
  @param InitialGuessIsZero (optional) If false, some work can be avoided. Whether this actually saves any work depends on the underlying Belos implementation.
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

  //! For diagnostic purposes
  // RCP<Ifpack2::Preconditioner<Scalar,LocalOrdinal,GlobalOrdinal,Node> > getPreconditioner(){return prec_;}

  //! Get a rough estimate of cost per iteration
  size_t getNodeSmootherComplexity() const;

 private:
  void SetupBelos(Level& currentLevel);

 private:
  std::string type_;

  typedef Tpetra::MultiVector<SC, LO, GO, NO> tMV;
  typedef Tpetra::Operator<SC, LO, GO, NO> tOP;
  RCP<Belos::LinearProblem<Scalar, tMV, tOP> > tBelosProblem_;
  RCP<Belos::SolverManager<Scalar, tMV, tOP> > tSolver_;

  //! matrix, used in apply if solving residual equation
  RCP<Matrix> A_;

};  // class BelosSmoother

#ifdef HAVE_MUELU_EPETRA

#if ((defined(EPETRA_HAVE_OMP) && (!defined(HAVE_TPETRA_INST_OPENMP) || !defined(HAVE_TPETRA_INST_INT_INT))) || \
     (!defined(EPETRA_HAVE_OMP) && (!defined(HAVE_TPETRA_INST_SERIAL) || !defined(HAVE_TPETRA_INST_INT_INT))))
// Stub specialization for missing Epetra template args
template <>
class BelosSmoother<double, int, int, Xpetra::EpetraNode> : public SmootherPrototype<double, int, int, Xpetra::EpetraNode> {
  typedef double Scalar;
  typedef int LocalOrdinal;
  typedef int GlobalOrdinal;
  typedef Xpetra::EpetraNode Node;
#undef MUELU_BELOSSMOOTHER_SHORT
#include "MueLu_UseShortNames.hpp"

 public:
#ifndef _MSC_VER
  // Avoid error C3772: invalid friend template declaration
  template <class Scalar2, class LocalOrdinal2, class GlobalOrdinal2, class Node2>
  friend class BelosSmoother;
#endif

  BelosSmoother(const std::string& type, const Teuchos::ParameterList& paramList = Teuchos::ParameterList(), const LocalOrdinal& overlap = 0) {
    MUELU_TPETRA_ETI_EXCEPTION("BelosSmoother<double,int,int,EpetraNode>", "BelosSmoother<double,int,int,EpetraNode>", "int");
  };

  virtual ~BelosSmoother() {}

  void SetParameterList(const Teuchos::ParameterList& paramList) {}
  void DeclareInput(Level& currentLevel) const {}
  void Setup(Level& currentLevel) {}
  void Apply(MultiVector& X, const MultiVector& B, bool InitialGuessIsZero = false) const {}
  RCP<SmootherPrototype> Copy() const { return Teuchos::null; }

  std::string description() const { return std::string(""); }
  void print(Teuchos::FancyOStream& out, const VerbLevel verbLevel = Default) const {}

  //! Get a rough estimate of cost per iteration
  size_t getNodeSmootherComplexity() const {
    size_t cplx = 0;
    return cplx;
  };
};
#endif

#endif  // HAVE_MUELU_EPETRA

}  // namespace MueLu

#define MUELU_BELOSSMOOTHER_SHORT
#endif  // HAVE_MUELU_BELOS
#endif  // MUELU_BELOSSMOOTHER_DECL_HPP
