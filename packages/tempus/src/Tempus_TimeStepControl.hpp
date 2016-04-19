#ifndef TEMPUS_TIMESTEPCONTROL_HPP
#define TEMPUS_TIMESTEPCONTROL_HPP

// Teuchos
#include "Teuchos_VerboseObject.hpp"
#include "Teuchos_Describable.hpp"
#include "Teuchos_ParameterListAcceptorDefaultBase.hpp"
// Tempus
#include "Tempus_StepType.hpp"
#include "Tempus_SolutionStateMetaData.hpp"

namespace Tempus {

/** \brief TimeStepControl manages the time step size.
 *  There several mechanicisms that effect the time step size and
 *  handled with this class:
 *   - Maximum and minimum time
 *   - Maximum and minimum time index
 *   - Maximum and minimum time step size
 *   - Maximum and minimum error
 *   - Maximum and minimum order
 *   - Startup considerations (e.g., ramping)
 *   - Solution and/or diagnostic output
 *  Additional step control can be added through the step control observer,
 *  or inheriting from this class.
 *   - Stability limits (e.g., CFL number)
 */
template<class Scalar>
class TimeStepControl
  : virtual public Teuchos::Describable,
    virtual public Teuchos::ParameterListAcceptor,
    virtual public Teuchos::VerboseObject<Tempus::TimeStepControl<Scalar> >
{
public:

  /// Destructor
  virtual ~TimeStepControl() {};

  /** \brief Default constructor. */
  TimeStepControl();

  /** \brief. */
  TimeStepControl(Teuchos::RCP<Teuchos::ParameterList> pList_ = Teuchos::null);

  /** \brief. */
  // This is a copy constructor
  TimeStepControl(const TimeStepControl<Scalar>& tsc_);

  /** \brief. .*/
  virtual void getNextTimeStep(Teuchos::RCP<SolutionStateMetaData<Scalar> > metaData,
    bool stepperStatus, bool integratorStatus) const;

  /** \brief. .*/
  virtual bool timeInRange(const Scalar time) const;

  /** \brief. .*/
  virtual bool indexInRange(const int iStep) const;


  Scalar timeMin;         ///< Minimum simulation time
  Scalar timeMax;         ///< Maximum simulation time
  Scalar dtMin;           ///< Minimum time step
  Scalar dtMax;           ///< Maximum time step
  int    iStepMin;        ///< Minimum time step index
  int    iStepMax;        ///< Maximum time step index
  Scalar errorMaxAbs;     ///< Maximum absolute error
  Scalar errorMaxRel;     ///< Maximum relative error
  unsigned int orderMin;  ///< Minimum time integration order
  unsigned int orderMax;  ///< Maximum time integration order

  StepType stepType;      ///< Step type for step control
  Scalar dtConstant;      ///< Constant time step if stepType=CONSTANT_STEP_SIZE

  std::vector<int>    outputIndices;  ///< Vector of output indices.
  std::vector<Scalar> outputTimes;    ///< Vector of output times.

  unsigned int nFailuresMax;            ///< Maximum number of stepper failures
  unsigned int nConsecutiveFailuresMax; ///< Maximum number of consecutive stepper failures

  Teuchos::RCP<Teuchos::ParameterList> pList;

  /// \name Overridden from Teuchos::ParameterListAcceptor
  //@{
    virtual void setParameterList(Teuchos::RCP<Teuchos::ParameterList> const& pl);
    virtual Teuchos::RCP<const Teuchos::ParameterList> getValidParameters() const;
    virtual Teuchos::RCP<const Teuchos::ParameterList> getParameterList() const;
    virtual Teuchos::RCP<Teuchos::ParameterList> getNonconstParameterList();
    virtual Teuchos::RCP<Teuchos::ParameterList> unsetParameterList();
  //@}

  /// \name Overridden from Teuchos::Describable
  //@{
    virtual std::string description() const;
    virtual void describe(Teuchos::FancyOStream          &out,
                          const Teuchos::EVerbosityLevel verbLevel) const;
  //@}

};
} // namespace Tempus
#endif // TEMPUS_TIMESTEPCONTROL_HPP
