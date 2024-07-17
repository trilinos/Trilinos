// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_IFPACKSMOOTHER_HPP
#define MUELU_IFPACKSMOOTHER_HPP

#include "MueLu_ConfigDefs.hpp"
#if defined(HAVE_MUELU_EPETRA) && defined(HAVE_MUELU_IFPACK)

#include <Teuchos_ParameterList.hpp>
#include <Xpetra_Matrix.hpp>
#include <Xpetra_MultiVectorFactory_fwd.hpp>

class Ifpack_Preconditioner;

#include "MueLu_SmootherPrototype.hpp"
#include "MueLu_Exceptions.hpp"

#include "MueLu_Level_fwd.hpp"
#include "MueLu_FactoryBase_fwd.hpp"
#include "MueLu_Utilities_fwd.hpp"
#include "MueLu_Aggregates_fwd.hpp"

namespace MueLu {

/*!
  @class IfpackSmoother
  @ingroup MueLuSmootherClasses
  @brief Class that encapsulates Ifpack smoothers.

  This class creates an Ifpack preconditioner factory. The factory creates a smoother based on the
  type and ParameterList passed into the constructor. See the constructor for more information.
*/
template <class Node = typename SmootherPrototype<double, int, int>::node_type>
class IfpackSmoother : public MueLu::SmootherPrototype<double, int, int, Node> {
  typedef double Scalar;
  typedef int LocalOrdinal;
  typedef int GlobalOrdinal;
#undef MUELU_IFPACKSMOOTHER_SHORT
#include "MueLu_UseShortNames.hpp"

 public:
  //! @name Constructors / destructors
  //@{

  /*! @brief Constructor

  The options passed into IfpackSmoother are those given in the Ifpack user's manual.

  @param type smoother type
  @param list options for the particular smoother (e.g., fill factor or damping parameter)

  Here is how to select some of the most common smoothers.

  - Gauss-Seidel
  - <tt>type</tt> = <tt>point relaxation stand-alone</tt>
  - parameter list options
  - <tt>relaxation: type</tt> = <tt>Gauss-Seidel</tt>
  - <tt>relaxation: damping factor</tt>
  - symmetric Gauss-Seidel
  - <tt>type</tt> = <tt>point relaxation stand-alone</tt>
  - parameter list options
  - <tt>relaxation: type</tt> = <tt>symmetric Gauss-Seidel</tt>
  - <tt>relaxation: damping factor</tt>
  - Chebyshev
  - <tt>type</tt> = <tt>Chebyshev</tt>
  - parameter list options
  - <tt>chebyshev: ratio eigenvalue</tt>
  - <tt>chebyshev: min eigenvalue</tt>
  - <tt>chebyshev: max eigenvalue</tt>
  - <tt>chebyshev: degree</tt>
  - <tt>chebyshev: zero starting solution</tt> (defaults to <tt>true</tt>)
  - ILU
  - <tt>type</tt> = <tt>ILU</tt>
  - parameter list options
  - <tt>fact: level-of-fill</tt>

  See also Ifpack_PointRelaxation, Ifpack_Chebyshev, Ifpack_ILU.
  */
  IfpackSmoother(std::string const& type, Teuchos::ParameterList const& paramList = Teuchos::ParameterList(), LO const& overlap = 0);  // TODO: empty paramList valid for Ifpack??

  //! Destructor
  virtual ~IfpackSmoother() {}

  //@}

  void SetParameterList(const Teuchos::ParameterList& paramList);

  //! Input
  //@{

  void DeclareInput(Level& currentLevel) const;

  //@}

  //! @name Computational methods.
  //@{

  /*! @brief Set up the smoother.

  This creates the underlying Ifpack smoother object, copies any parameter list options
  supplied to the constructor to the Ifpack object, and computes the preconditioner.
  */
  void Setup(Level& currentLevel);

  /*! @brief Apply the preconditioner.

  Solves the linear system <tt>AX=B</tt> using the constructed smoother.

  @param X initial guess
  @param B right-hand side
  @param InitialGuessIsZero (optional) If false, some work can be avoided.  Whether this actually saves any work depends on the underlying Ifpack implementation.
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

  //! Get a rough estimate of cost per iteration
  size_t getNodeSmootherComplexity() const;

  //@}

 private:
  void SetPrecParameters(const Teuchos::ParameterList& list = Teuchos::ParameterList()) const;

  void SetupAggregate(Level& currentLevel);

 private:
  //! ifpack-specific key phrase that denote smoother type
  std::string type_;

  //! overlap when using the smoother in additive Schwarz mode
  LO overlap_;

  //! Matrix. Not used directly, but held inside of prec_. So we have to keep an RCP pointer to it!
  RCP<Matrix> A_;

  //! pointer to Ifpack solver object
  // Note: prec_ must be destroyed before A_, so declaration of prec_ appears after declaration of A_
  RCP<Ifpack_Preconditioner> prec_;

};  // class IfpackSmoother

//! Non-member templated function GetIfpackSmoother() returns a new IfpackSmoother object when <Scalar, LocalOrdinal, GlobalOrdinal> == <double, int, int>. Otherwise, an exception is thrown.
//! This function simplifies the usage of IfpackSmoother objects inside of templates as templates do not have to be specialized for <double, int, int> (see DirectSolver for an example).
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<MueLu::SmootherPrototype<Scalar, LocalOrdinal, GlobalOrdinal, Node> >
GetIfpackSmoother(const std::string& /* type */                 = "",
                  const Teuchos::ParameterList& /* paramList */ = Teuchos::ParameterList(),
                  const LocalOrdinal& /* overlap */             = 0) {
  TEUCHOS_TEST_FOR_EXCEPTION(true, Exceptions::RuntimeError, "IfpackSmoother cannot be used with Scalar != double, LocalOrdinal != int, GlobalOrdinal != int");
  TEUCHOS_UNREACHABLE_RETURN(Teuchos::null);
}

// Specialization for serial node (used for Epetra)
#if defined(HAVE_MUELU_EPETRA)
template <>
inline RCP<MueLu::SmootherPrototype<double, int, int, Xpetra::EpetraNode> >
GetIfpackSmoother<double, int, int, Xpetra::EpetraNode>(const std::string& type, const Teuchos::ParameterList& paramList, const int& overlap) {
  return rcp(new MueLu::IfpackSmoother<Xpetra::EpetraNode>(type, paramList, overlap));
}
#endif

}  // namespace MueLu

#define MUELU_IFPACKSMOOTHER_SHORT
#endif  // HAVE_MUELU_EPETRA && HAVE_MUELU_IFPACK
#endif  // MUELU_IFPACKSMOOTHER_HPP
