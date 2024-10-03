// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_STRATIMIKOSSMOOTHER_DECL_HPP
#define MUELU_STRATIMIKOSSMOOTHER_DECL_HPP

#include <Teuchos_ParameterList.hpp>

#include <Xpetra_CrsMatrixWrap_fwd.hpp>
#include <Xpetra_Matrix_fwd.hpp>

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_StratimikosSmoother_fwd.hpp"

#if defined(HAVE_MUELU_STRATIMIKOS) && defined(HAVE_MUELU_THYRA)

#include <Tpetra_CrsMatrix.hpp>

#include "MueLu_FactoryBase_fwd.hpp"
#include "MueLu_Level_fwd.hpp"
#include "MueLu_SmootherPrototype.hpp"
#include "MueLu_Utilities_fwd.hpp"

#include <Thyra_LinearOpWithSolveBase.hpp>

namespace MueLu {

/*!
  @class StratimikosSmoother
  @ingroup MueLuSmootherClasses
  @brief Class that encapsulates Stratimikos smoothers.

  This class creates an Stratimikos preconditioner factory. The factory creates a smoother based
  on the type and ParameterList passed into the constructor. See the constructor for more
  information.
  */

template <class Scalar        = SmootherPrototype<>::scalar_type,
          class LocalOrdinal  = typename SmootherPrototype<Scalar>::local_ordinal_type,
          class GlobalOrdinal = typename SmootherPrototype<Scalar, LocalOrdinal>::global_ordinal_type,
          class Node          = typename SmootherPrototype<Scalar, LocalOrdinal, GlobalOrdinal>::node_type>
class StratimikosSmoother : public SmootherPrototype<Scalar, LocalOrdinal, GlobalOrdinal, Node> {
#undef MUELU_STRATIMIKOSSMOOTHER_SHORT
#include "MueLu_UseShortNames.hpp"

 public:
  //! @name Constructors / destructors
  //@{
  // TODO: update doc for Stratimikos.
  /*! @brief Constructor

  The options passed into StratimikosSmoother are those given in the Stratimikos user's manual.

  @param type smoother type
  @param list options for the particular smoother (e.g., fill factor or damping parameter)

  */

#ifndef _MSC_VER
  // Avoid error C3772: invalid friend template declaration
  template <class Scalar2, class LocalOrdinal2, class GlobalOrdinal2, class Node2>
  friend class StratimikosSmoother;
#endif

  StratimikosSmoother(const std::string type, const Teuchos::ParameterList& paramList = Teuchos::ParameterList()){
      TEUCHOS_TEST_FOR_EXCEPTION(true, Exceptions::RuntimeError, "Stratimikos only works for Scalar=double.")};

  //! Destructor
  virtual ~StratimikosSmoother() = default;

  //@}

  void SetParameterList(const Teuchos::ParameterList& paramList){};

  //! Input
  //@{

  void DeclareInput(Level& currentLevel) const {};

  //@}

  //! @name Computational methods.
  //@{

  /*! @brief Set up the smoother.

  This creates the underlying Stratimikos smoother object, copies any parameter list options
  supplied to the constructor to the Stratimikos object, and computes the preconditioner.
  */
  void Setup(Level& currentLevel){};

  /*! @brief Apply the preconditioner.

  Solves the linear system <tt>AX=B</tt> using the constructed smoother.

  @param X initial guess
  @param B right-hand side
  @param InitialGuessIsZero (optional) If false, some work can be avoided. Whether this actually saves any work depends on the underlying Stratimikos implementation.
  */
  void Apply(MultiVector& X, const MultiVector& B, bool InitialGuessIsZero = false) const {};

  //@}

  //! @name Utilities
  //@{

  RCP<SmootherPrototype> Copy() const { return Teuchos::null; };

  //@}

  //! @name Overridden from Teuchos::Describable
  //@{

  //! Return a simple one-line description of this object.
  std::string description() const { return std::string(""); };

  //! Print the object with some verbosity level to an FancyOStream object.
  // using MueLu::Describable::describe; // overloading, not hiding
  // void describe(Teuchos::FancyOStream &out, const VerbLevel verbLevel = Default) const
  void print(Teuchos::FancyOStream& out, const VerbLevel verbLevel = Default) const {};

  //@}

  //! For diagnostic purposes
  // RCP<Ifpack2::Preconditioner<Scalar,LocalOrdinal,GlobalOrdinal,Node> > getPreconditioner(){return prec_;}

  //! Get a rough estimate of cost per iteration
  size_t getNodeSmootherComplexity() const { return Teuchos::OrdinalTraits<size_t>::invalid(); };

};  // class StratimikosSmoother

template <class LocalOrdinal, class GlobalOrdinal, class Node>
struct StratimikosSmoother<double, LocalOrdinal, GlobalOrdinal, Node> : public SmootherPrototype<double, LocalOrdinal, GlobalOrdinal, Node> {
  typedef double Scalar;
#undef MUELU_STRATIMIKOSSMOOTHER_SHORT
#include "MueLu_UseShortNames.hpp"

  StratimikosSmoother(const std::string type, const Teuchos::ParameterList& paramList = Teuchos::ParameterList());

  //! Destructor
  virtual ~StratimikosSmoother() = default;

  //@}

  void SetParameterList(const Teuchos::ParameterList& paramList);

  //! Input
  //@{

  void DeclareInput(Level& currentLevel) const;

  //@}

  //! @name Computational methods.
  //@{

  /*! @brief Set up the smoother.

  This creates the underlying Stratimikos smoother object, copies any parameter list options
  supplied to the constructor to the Stratimikos object, and computes the preconditioner.
  */
  void Setup(Level& currentLevel);

  /*! @brief Apply the preconditioner.

  Solves the linear system <tt>AX=B</tt> using the constructed smoother.

  @param X initial guess
  @param B right-hand side
  @param InitialGuessIsZero (optional) If false, some work can be avoided. Whether this actually saves any work depends on the underlying Stratimikos implementation.
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
  void SetupStratimikos(Level& currentLevel);
  /*!
    @brief Filter out any matrix connections that corresponds to coupling along a vertical grid line

    Builds a special filtered matrix for meshes that are structured in at least one direction (the z direction). Specifically
    the new matrix has all the vertical connections removed.
    */
  void ExperimentalDropVertConnections(RCP<Matrix>& filteredA, Level& currentLevel);

  //@}

  std::string type_;

  Teuchos::RCP<Thyra::LinearOpWithSolveBase<Scalar> > solver_;

  //! matrix, used in apply if solving residual equation
  RCP<Matrix> A_;

  bool recurMgOnFilteredA_ = false;
};  // class StratimikosSmoother

}  // namespace MueLu

#define MUELU_STRATIMIKOSSMOOTHER_SHORT
#endif  // HAVE_MUELU_STRATIMIKOS
#endif  // MUELU_STRATIMIKOSSMOOTHER_DECL_HPP
