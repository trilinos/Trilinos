// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_PROJECTORSMOOTHER_DECL_HPP
#define MUELU_PROJECTORSMOOTHER_DECL_HPP

#include <Xpetra_Matrix_fwd.hpp>
#include <Xpetra_MultiVector_fwd.hpp>

#include "MueLu_ConfigDefs.hpp"

#include "MueLu_SmootherPrototype.hpp"
#include "MueLu_Utilities_fwd.hpp"

namespace MueLu {

/*!
  @class ProjectorSmoother
  @ingroup MueLuSmootherClasses
  @brief This class enables the elimination of the nullspace component of the solution through the use of projection

  The solution of the coarsest level system may have a significant nullspace component. We can try to eliminate it through the use of nullspaces.
  Due to our construction algorithms, we always have some nullspace vectors on the coarsest level. However, we do not know if they are true
  nullspace vectors or are only approximations to them.

  We will use the Rayleigh quotiont to separate the true nullspace vectors. If Rayleigh quotient is less than specified epsilon, we consider the
  vector to be a true nullspace vector.

  After we separate true nullspace vectors, we project our solution to the space orthogonal to them. To do that in a computationally efficient manner,
  in the Setup stage we orthonormalize the selected nullspace components.
*/

template <class Scalar        = SmootherPrototype<>::scalar_type,
          class LocalOrdinal  = typename SmootherPrototype<Scalar>::local_ordinal_type,
          class GlobalOrdinal = typename SmootherPrototype<Scalar, LocalOrdinal>::global_ordinal_type,
          class Node          = typename SmootherPrototype<Scalar, LocalOrdinal, GlobalOrdinal>::node_type>
class ProjectorSmoother : public SmootherPrototype<Scalar, LocalOrdinal, GlobalOrdinal, Node> {
#undef MUELU_PROJECTORSMOOTHER_SHORT
#include "MueLu_UseShortNames.hpp"

 public:
  //! @name Constructors / destructors
  //@{

  //! @brief Constructor
  ProjectorSmoother(RCP<SmootherPrototype> coarseSolver);

  //! Destructor
  virtual ~ProjectorSmoother();
  //@}

  //! Input
  //@{

  void DeclareInput(Level &currentLevel) const;

  //@}

  //! @name Setup and Apply methods.
  //@{

  //! @brief Set up the direct solver.
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

  //! Print the object with some verbosity level to an FancyOStream object.
  // using MueLu::Describable::describe; // overloading, not hiding
  void print(Teuchos::FancyOStream &out, const VerbLevel verbLevel = Default) const;

  //! Get a rough estimate of cost per iteration
  size_t getNodeSmootherComplexity() const { return coarseSolver_->getNodeSmootherComplexity(); }

  //@}

 private:
  RCP<MultiVector> Borth_;
  RCP<SmootherPrototype> coarseSolver_;

};  // class ProjectorSmoother

}  // namespace MueLu

#define MUELU_PROJECTORSMOOTHER_SHORT
#endif  // MUELU_PROJECTORSMOOTHER_DECL_HPP
