// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*
 * MueLu_PermutingSmoother_decl.hpp
 *
 *  Created on: Nov 28, 2012
 *      Author: wiesner
 */

#ifndef MUELU_PERMUTINGSMOOTHER_DECL_HPP
#define MUELU_PERMUTINGSMOOTHER_DECL_HPP

#include <Teuchos_ParameterList.hpp>
#include <Xpetra_Matrix_fwd.hpp>
#include <Xpetra_Vector_fwd.hpp>
#include <Xpetra_VectorFactory_fwd.hpp>
#include <Xpetra_MultiVector_fwd.hpp>
#include <Xpetra_MultiVectorFactory_fwd.hpp>

#include "MueLu_ConfigDefs.hpp"

#include "MueLu_SmootherPrototype.hpp"
#include "MueLu_TrilinosSmoother_fwd.hpp"
#include "MueLu_FactoryBase_fwd.hpp"
#include "MueLu_PermutationFactory_fwd.hpp"

namespace MueLu {

/*!
  @class PermutingSmoother
  @ingroup MueLuSmootherClasses
  @brief This class first calculates row- and column permutation operators and applies a smoother to the permuted linear system.

*/

template <class Scalar        = SmootherPrototype<>::scalar_type,
          class LocalOrdinal  = typename SmootherPrototype<Scalar>::local_ordinal_type,
          class GlobalOrdinal = typename SmootherPrototype<Scalar, LocalOrdinal>::global_ordinal_type,
          class Node          = typename SmootherPrototype<Scalar, LocalOrdinal, GlobalOrdinal>::node_type>
class PermutingSmoother : public SmootherPrototype<Scalar, LocalOrdinal, GlobalOrdinal, Node> {
#undef MUELU_PERMUTINGSMOOTHER_SHORT
#include "MueLu_UseShortNames.hpp"

 public:
  //! @name Constructors / destructors
  //@{

  //! @brief Constructor
  //!    @param[in] mapName   Name of map object in level class, which rows/cols can be permuted
  //!    @param[in] mapFact   generating factory of map with name mapName
  //!    @param[in] type      string that contains type of smoother (e.g. "RELAXATION" or "ILU")
  //!    @param[in] paramList parameter list with parameters for smoother (default: empty)
  //!    @param[in] overlap   LocalOrdinal with overlap inforation (default: 0)
  //!    @param[in] permFact  factory, generating permutation and scaling matrices (default: Teuchos::null -> use internal PermutationFactory instance)
  PermutingSmoother(std::string const& mapName, const RCP<const FactoryBase>& mapFact, std::string const& type = "", const Teuchos::ParameterList& paramList = Teuchos::ParameterList(), LO const& overlap = 0, RCP<FactoryBase> permFact = Teuchos::null);

  //! Destructor
  virtual ~PermutingSmoother();
  //@}

  //! Input
  //@{

  void DeclareInput(Level& currentLevel) const;

  //@}

  //! @name Setup and Apply methods.
  //@{

  //! @brief Set up the direct solver.
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
  size_t getNodeSmootherComplexity() const { return s_->getNodeSmootherComplexity(); }

  //@}

 private:
  //! ifpack1/2-specific key phrase that denote smoother type
  std::string type_;

  //! overlap when using the smoother in additive Schwarz mode
  LO overlap_;

  //! Permutation Factory
  RCP<FactoryBase> permFact_;

  //! permQT matrix object
  RCP<Matrix> permQT_;

  //! permP matrix object
  RCP<Matrix> permP_;

  //! scaling matrix object
  Teuchos::RCP<Matrix> diagScalingOp_;

  //
  // Underlying Smoother
  //

  //! Smoother
  RCP<SmootherPrototype> s_;  // TrilinosSmoother object

};  // class PermutingSmoother

}  // namespace MueLu

#define MUELU_PERMUTINGSMOOTHER_SHORT
#endif /* MUELU_PERMUTINGSMOOTHER_DECL_HPP */
