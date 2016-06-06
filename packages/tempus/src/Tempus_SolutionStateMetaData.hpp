#ifndef TEMPUS_SOLUTIONSTATEMETADATA_HPP
#define TEMPUS_SOLUTIONSTATEMETADATA_HPP

// Teuchos
#include "Teuchos_VerboseObject.hpp"
#include "Teuchos_Describable.hpp"
// Tempus
#include "Tempus_Types.hpp"


namespace Tempus {

/** \brief Solution state meta data.
 */
template<class Scalar>
class SolutionStateMetaData :
  public Teuchos::Describable,
  public Teuchos::VerboseObject<Tempus::SolutionStateMetaData<Scalar> >
{
public:

  /// Default constructor.
  SolutionStateMetaData();

  /// Constructor
  SolutionStateMetaData(
    const Scalar time_,
    const int    iStep_,
    const Scalar dt_,
    const Scalar errorAbs_,
    const Scalar errorRel_,
    const int    order_,
    const int    nFailures_,
    const int    nConsecutiveFailures_,
    const Status solutionStatus_,
    const bool   output_,
    const bool   isRestartable_,
    const bool   isInterpolated_,
    const Scalar accuracy_);

  /// Copy constructor
  SolutionStateMetaData(const SolutionStateMetaData<Scalar>& ssmd_);

  /// Clone
  Teuchos::RCP<SolutionStateMetaData<Scalar> > clone();

  /// Destructor
  virtual ~SolutionStateMetaData() {};

  Scalar time;              ///< Time of solution
  int    iStep;             ///< Time step index for this solution
  Scalar dt;                ///< Time step for this solution
  Scalar errorAbs;          ///< Absolute local truncation error
  Scalar errorRel;          ///< Relative local truncation error
  int order;                ///< Order of this solution
  int nFailures;            ///< Total number of stepper failures
  int nConsecutiveFailures; ///< Consecutive number of stepper failures

  /** The solutionStatus is used to indicate
      - if the solution is still being worked on; WORKING
      - if the solution is accepted and completed (e.g., past solutions
        in SolutionHistory); PASSED.
      - if the time step has FAILED.  This may be caused by the Stepper
        failing, or Integrator not accepting the time step.
  */
  Status solutionStatus;
  bool   output;            ///< SolutionState should be or has been outputted
  bool   isRestartable;     ///< T - soln can be used as a restart
  bool   isInterpolated;    ///< F - soln is time integrated; T - soln is interpolated
  Scalar accuracy;          ///< Interpolation accuracy of solution

  /// \name Overridden from Teuchos::Describable
  //@{
    virtual std::string description() const;
    virtual void describe(Teuchos::FancyOStream          &out,
                          const Teuchos::EVerbosityLevel verbLevel) const;
  //@}

};
} // namespace Tempus

#include "Tempus_SolutionStateMetaData_impl.hpp"

#endif // TEMPUS_SOLUTIONSTATEMETADATA_HPP
