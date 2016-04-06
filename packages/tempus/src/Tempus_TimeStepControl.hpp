#ifndef TEMPUS_TIMESTEPCONTROL_HPP
#define TEMPUS_TIMESTEPCONTROL_HPP

#include "Teuchos_VerboseObject.hpp"
#include "Teuchos_Describable.hpp"

#include "Tempus_StepType.hpp"

namespace tempus {

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
    virtual public Teuchos::VerboseObject<tempus::TimeStepControl<Scalar> >
{
public:

  typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType ScalarMag;

  /// Destructor
  virtual ~TimeStepControl() {};

  /** \brief Default constructor. */
  TimeStepControl();

  /** \brief. */
  TimeStepControl( RCP<Teuchos::ParameterList> paramList_ = Teuchos::null );

  /** \brief. */
  // This is a copy constructor
  TimeStepControl(const TimeStepControl<Scalar>& tsc_);


  /** \brief. Basic time-step control.*/
  virtual void TimeStepControl::setNextTimeStep(const Scalar time,
    const int iStep, const Scalar errorAbs, const Scalar errorRel,
    const int order, Scalar dt, bool output) const;

  Scalar timeMin;         ///< Minimum simulation time
  Scalar timeMax;         ///< Maximum simulation time
  Scalar dtMin;           ///< Minimum time step
  Scalar dtMax;           ///< Maximum time step
  int    iStepMin;        ///< Minimum time step index
  int    iStepMax;        ///< Maximum time step index
  Scalar errorMaxAbs;     ///< Maximum absolute error
  Scalar errorMaxRel;     ///< Maximum relative error
  int    orderMin;        ///< Minimum time integration order
  int    orderMax;        ///< Maximum time integration order

  StepType stepType;      ///< Step type for step control
  Scalar dtConstant;      ///< Constant time step if stepType=CONSTANT_STEP_SIZE

  std::vector<int>    outputIndices;  ///< Vector of output indices.
  std::vector<Scalar> outputTimes;    ///< Vector of output times.

  RCP<Teuchos::ParameterList> paramList;

  /// Inherited from Describable:
  /** \brief . */
  virtual std::string description() const;

  /** \brief . */
  virtual void describe( Teuchos::FancyOStream          &out,
                         const Teuchos::EVerbosityLevel verbLevel) const;

  /// Redefined from Teuchos::ParameterListAcceptor
  /** \brief . */
  void setParameterList(RCP<Teuchos::ParameterList> const& paramList);

  /** \brief . */
  RCP<Teuchos::ParameterList> getNonconstParameterList();

  /** \brief . */
  RCP<Teuchos::ParameterList> unsetParameterList();

  RCP<const Teuchos::ParameterList> getValidParameters() const;

};
} // namespace tempus
#endif // TEMPUS_TIMESTEPCONTROL_HPP
