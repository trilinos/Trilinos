#ifndef TEMPUS_INTEGRATORBASIC_HPP
#define TEMPUS_INTEGRATORBASIC_HPP

// Teuchos
#include "Teuchos_VerboseObject.hpp"
#include "Teuchos_Describable.hpp"
#include "Teuchos_ParameterList.hpp"
// Tempus
#include "Tempus_Integrator.hpp"

#include <string>

namespace tempus {

/** \brief Basic time integrator
 */
template<class Scalar>
class IntegratorBasic : public tempus::Integrator {
public:

  /** \brief Constructor with ParameterList, models, initial conditions
   *  and optional solvers. */
  IntegratorBasic(
    RCP<ParameterList>                     parameterList,
    const RCP<Thyra::VectorBase<Scalar> >& x,
    const RCP<Thyra::VectorBase<Scalar> >& xdot=Teuchos::null,
    const RCP<Thyra::VectorBase<Scalar> >& xdotdot=Teuchos::null );

  /// Destructor
  virtual ~IntegratorBasic() {}

  //! Unique name for this integrator.
  virtual std::string name() const = 0;

  /// \name Basic integrator methods
  //@{
    /// Advance the solution to time, and return true if successful.
    virtual bool advanceTime(const Scalar time) = 0;
  //@}

  /// \name Overridden from Teuchos::ParameterListAcceptor
  //@{
    virtual void setParameterList(RCP<ParameterList> const& pl);
    virtual RCP<ParameterList> getNonconstParameterList();
    virtual RCP<ParameterList> unsetParameterList();
    virtual RCP<const ParameterList> getValidParameters() const;
  //@}

  /// \name Overridden from Teuchos::Describable
  //@{
    virtual std::string description() const;
    virtual void describe(Teuchos::FancyOStream        & out,
                          const Teuchos::EVerbosityLevel verbLevel) const;
  //@}

  /// \name Undo type capabilities
  //@{
    /// Only accept step after meeting time step criteria.
    virtual bool acceptStep() {return false;}
  //@}

};
} // namespace tempus
#endif // TEMPUS_INTEGRATORBASIC_HPP
