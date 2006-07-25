
#ifndef THYRA_DEFAULT_REAL_LINEAR_SOLVER_BUILDING_BASE
#define THYRA_DEFAULT_REAL_LINEAR_SOLVER_BUILDING_BASE

#include "Stratimikos_Config.h"
#include "Thyra_LinearSolverBuilderBase.hpp"
#include "Teuchos_AbstractFactory.hpp"

namespace Thyra {

/** \brief Abstract interface for an object that can create
 * <tt>LinearOpWithSolveFactoryBase</tt> objects on demand.
 *
 * ToDo: Finish documentation!
 */
class DefaultRealLinearSolverBuilder : public LinearSolverBuilderBase<double>
{
public:

  /** @name Constructors/Initializers/Accessors */
  //@{

  /** \brief Construct without a parameter list. */
  DefaultRealLinearSolverBuilder();

  /** \brief Construct given a parameter list. */
  DefaultRealLinearSolverBuilder(
    Teuchos::RefCountPtr<Teuchos::ParameterList> const& paramList
    );

  /** \brief Set a new linear solver strategy factory object. */
  void setLinearSolveStrategyFactory(
    const Teuchos::RefCountPtr<const Teuchos::AbstractFactory<LinearOpWithSolveFactoryBase<double> > >  &solveStrategyFactory
    ,const std::string                                                                                  &solveStrategyName
    );

  /** \brief Set a new preconditioner strategy factory object. */
  void setPreconditioningStrategyFactory(
    const Teuchos::RefCountPtr<const Teuchos::AbstractFactory<PreconditionerFactoryBase<double> > >     &precStrategyFactory
    ,const std::string                                                                                  &precStrategyName
    );

  /** \brief Get the name of the linear solver strategy that will be created. */
  std::string getLinearSolveStrategyName() const;

  /** \brief Get the name of the preconditioner strategy that will be created. */
  std::string getPreconditionerStrategyName() const;

  //@}

  /** @name Overridden from ParameterListAcceptor */
  //@{

  /** \brief . */
  void setParameterList(Teuchos::RefCountPtr<Teuchos::ParameterList> const& paramList);
  /** \brief . */
  Teuchos::RefCountPtr<Teuchos::ParameterList> getParameterList();
  /** \brief . */
  Teuchos::RefCountPtr<Teuchos::ParameterList> unsetParameterList();
  /** \brief . */
  Teuchos::RefCountPtr<const Teuchos::ParameterList> getParameterList() const;
  /** \brief . */
  Teuchos::RefCountPtr<const Teuchos::ParameterList> getValidParameters() const;

  //@}
  
  /** \name Overridden from LinearSolverBuilderBase. */
  //@{

  /** \brief . */
  Teuchos::RefCountPtr<LinearOpWithSolveFactoryBase<double> >
  createLinearSolveStrategy(
    const std::string &linearSolveStrategyName
    ) const;
  /** \brief . */
  Teuchos::RefCountPtr<PreconditionerFactoryBase<double> >
  createPreconditioningStrategy(
    const std::string &preconditioningStrategyName
    ) const;

  //@}

private:

  // //////////////////////////////////////
  // Private types

  typedef std::map<std::string,Teuchos::RefCountPtr<const Teuchos::AbstractFactory<LinearOpWithSolveFactoryBase<double> > > >  lowsf_map_t;
  typedef std::map<std::string,Teuchos::RefCountPtr<const Teuchos::AbstractFactory<PreconditionerFactoryBase<double> > > >     pf_map_t;

  // //////////////////////////////////////
  // Private data members
  
  Teuchos::RefCountPtr<Teuchos::ParameterList>                 paramList_;
  mutable Teuchos::RefCountPtr<const Teuchos::ParameterList>   validParamList_;
  lowsf_map_t                                                  lowsf_map_;
  std::string                                                  defaultLOWSF_;
  pf_map_t                                                     pf_map_;
  std::string                                                  defaultPF_;

  // //////////////////////////////////////
  // Private member functions

  void initializeDefaults();

};

} // namespace Thyra

#endif // THYRA_DEFAULT_REAL_LINEAR_SOLVER_BUILDING_BASE
