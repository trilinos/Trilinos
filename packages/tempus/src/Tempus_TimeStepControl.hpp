#ifndef TEMPUS_TIMESTEPCONTROL_HPP
#define TEMPUS_TIMESTEPCONTROL_HPP

// Teuchos
#include "Teuchos_VerboseObject.hpp"
#include "Teuchos_Describable.hpp"
#include "Teuchos_ParameterListAcceptorDefaultBase.hpp"
// Tempus
#include "Tempus_StepType.hpp"
#include "Tempus_SolutionHistory.hpp"

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

  /** \brief Default constructor. */
  TimeStepControl();

  /** \brief Construct from ParameterList */
  TimeStepControl(Teuchos::RCP<Teuchos::ParameterList> pList_ = Teuchos::null,
                  const Scalar dtConstant_ = 0.0);

  /// This is a copy constructor
  TimeStepControl(const TimeStepControl<Scalar>& tsc_);

  /// Destructor
  virtual ~TimeStepControl() {};

  /** \brief Determine the time step size.*/
  virtual void getNextTimeStep(
    const RCP<SolutionHistory<Scalar> > & solutionHistory,
    const bool stepperStatus, bool & integratorStatus) const;

  /** \brief Check if time is within minimum and maximum time. */
  virtual bool timeInRange(const Scalar time) const;

  /** \brief Check if time step index is within minimum and maximum index. */
  virtual bool indexInRange(const int iStep) const;

  /// \name Overridden from Teuchos::ParameterListAcceptor
  //@{
    void setParameterList(const Teuchos::RCP<Teuchos::ParameterList> & pl);
    Teuchos::RCP<const Teuchos::ParameterList> getValidParameters() const;
    Teuchos::RCP<Teuchos::ParameterList> getNonconstParameterList();
    Teuchos::RCP<Teuchos::ParameterList> unsetParameterList();
  //@}

  /// \name Overridden from Teuchos::Describable
  //@{
    std::string description() const;
    void describe(Teuchos::FancyOStream          &out,
                  const Teuchos::EVerbosityLevel verbLevel) const;
  //@}

  Scalar timeMin;         ///< Minimum simulation time
  Scalar timeMax;         ///< Maximum simulation time
  Scalar dtMin;           ///< Minimum time step
  Scalar dtMax;           ///< Maximum time step
  int    iStepMin;        ///< Minimum time step index
  int    iStepMax;        ///< Maximum time step index
  Scalar errorMaxAbs;     ///< Maximum absolute error
  Scalar errorMaxRel;     ///< Maximum relative error
  int orderMin;  ///< Minimum time integration order
  int orderMax;  ///< Maximum time integration order

  StepType stepType;      ///< Step type for step control
  Scalar dtConstant;      ///< Constant time step if stepType=CONSTANT_STEP_SIZE

  std::vector<int>    outputIndices;  ///< Vector of output indices.
  std::vector<Scalar> outputTimes;    ///< Vector of output times.

  int nFailuresMax;            ///< Maximum number of stepper failures
  int nConsecutiveFailuresMax; ///< Maximum number of consecutive stepper failures

  Teuchos::RCP<Teuchos::ParameterList> pList;
};
} // namespace Tempus

#include "Tempus_TimeStepControl_impl.hpp"

#endif // TEMPUS_TIMESTEPCONTROL_HPP
