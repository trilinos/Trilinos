// @HEADER
// ***********************************************************************
// 
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//                 Copyright (2004) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#ifndef THYRA_DEFAULT_NONLINEAR_SOLVER_BUILDER_HPP
#define THYRA_DEFAULT_NONLINEAR_SOLVER_BUILDER_HPP

#include "Thyra_NonlinearSolverBuilderBase.hpp"
#include "Teuchos_AbstractFactory.hpp"


namespace Thyra {


/** \brief Concrete subclass of <tt>Thyra::NonlinearSolverBuilderBase</tt> for
 * creating <tt>NonlinearSolverBase</tt> objects and
 * <tt>PreconditionerFactoryBase</tt> object on demand given configured
 * factory objects.
 *
 * ToDo: Finish Documentation!
 *
 * \ingroup Thyra_Nonlin_ME_solvers_grp
 */
class DefaultNonlinearSolverBuilder
  : public Thyra::NonlinearSolverBuilderBase<double>
{
public:

  /** @name Constructors/Initializers/Accessors */
  //@{

  /** \brief . */
  DefaultNonlinearSolverBuilder();

  /** \brief . */
  ~DefaultNonlinearSolverBuilder();

  /** \brief Set a new NonlinearSolverBase factory object. */
  void setNonlinearSolverFactory(
    const RCP<const AbstractFactory<Thyra::NonlinearSolverBase<double> > >
    &nonlinearSolverFactory,
    const std::string &nonlinearSolverTypeName
    );
  
  /** \brief Get the name of the NonlinearSolver type that will be created on
   * the next call to <tt>this->createNonlinearSolver()</tt>.
   */
  std::string getNonlinearSolverName() const;

  //@}

  /** @name Overridden from ParameterListAcceptor */
  //@{

  /** \brief . */
  void setParameterList(RCP<ParameterList> const& paramList);
  /** \brief . */
  RCP<ParameterList> getNonconstParameterList();
  /** \brief . */
  RCP<ParameterList> unsetParameterList();
  /** \brief . */
  RCP<const ParameterList> getParameterList() const;
  /** \brief . */
  RCP<const ParameterList> getValidParameters() const;

  //@}
  
  /** \name Overridden from NonlinearSolverBuilderBase. */
  //@{

  /** \brief . */
  virtual Teuchos::RCP<NonlinearSolverBase<Scalar> >
  createNonlinearSolver(const std::string &nonlinearSolverTypeName) const;

  //@}

private:

  // //////////////////////////////////////
  // Private types

  typedef RCP<const AbstractFactory<Thyra::NonlinearSolverBase<double> > >
  ns_fcty_t;

  // //////////////////////////////////////
  // Private data members
  
  RCP<ParameterList> paramList_;
  mutable RCP<const ParameterList> validParamList_;
  Array<std::string> validNonlinearSolverNames_;
  Array<ns_fcty_t> nonlinearSolverArray_;
  std::string defaultNonlinearSolverName_;

  // //////////////////////////////////////
  // Private member functions

  void initializeDefaults();

};


} // namespace Thyra


/** \brief Inject a new solver type into a DefaultNonlinearSolverBuilder
 * object.
 */
template<class NonlinearSolverType, class Scalar>
void setNonlinearSolverFactory(
  const std::string &nonlinearSolverTypeName,
  const Ptr<DefaultNonlinearSolverBuilder<Scalar> > &defaultNonlinearSolverBuilder
  )
{
  TEST_FOR_EXCEPT(true);
}


#endif // THYRA_DEFAULT_NONLINEAR_SOLVER_BUILDER_HPP
