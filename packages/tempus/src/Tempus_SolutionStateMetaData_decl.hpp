#ifndef TEMPUS_SOLUTIONSTATEMETADATA_DECL_HPP
#define TEMPUS_SOLUTIONSTATEMETADATA_DECL_HPP

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
    const Scalar time,
    const int    iStep,
    const Scalar dt,
    const Scalar errorAbs,
    const Scalar errorRel,
    const int    order,
    const int    nFailures,
    const int    nConsecutiveFailures,
    const Status solutionStatus,
    const bool   output,
    const bool   outputScreen,
    const bool   isRestartable,
    const bool   isInterpolated,
    const Scalar accuracy);

  /// Copy constructor
  SolutionStateMetaData(const SolutionStateMetaData<Scalar>& ssmd);

  /// Clone constructor
  Teuchos::RCP<SolutionStateMetaData<Scalar> > clone();

  /// This is a deep copy
  void copy(Teuchos::RCP<SolutionStateMetaData<Scalar> > ssmd);

  /// Destructor
  virtual ~SolutionStateMetaData() {};

  /// \name Accessor methods
  //@{
    Scalar getTime()                 const {return time_;}
    int    getIStep()                const {return iStep_;}
    Scalar getDt()                   const {return dt_;}
    Scalar getErrorAbs()             const {return errorAbs_;}
    Scalar getErrorRel()             const {return errorRel_;}
    int    getOrder()                const {return order_;}
    int    getNFailures()            const {return nFailures_;}
    int    getNConsecutiveFailures() const {return nConsecutiveFailures_;}
    Status getSolutionStatus()       const {return solutionStatus_;}
    bool   getOutput()               const {return output_;}
    bool   getOutputScreen()         const {return outputScreen_;}
    bool   getIsRestartable()        const {return isRestartable_;}
    bool   getIsInterpolated()       const {return isInterpolated_;}
    Scalar getAccuracy()             const {return accuracy_;}

    void setTime(Scalar time) {time_ = time;}
    void setIStep(int iStep) {iStep_ = iStep;}
    void setDt(Scalar dt) {dt_ = dt;}
    void setErrorAbs (Scalar errorAbs){errorAbs_ = errorAbs;}
    void setErrorRel (Scalar errorRel){errorRel_ = errorRel;}
    void setOrder(int order) {order_ = order;}
    void setNFailures(int nFailures) {nFailures_ = nFailures;}
    void setNConsecutiveFailures(int nConsecutiveFailures)
      {nConsecutiveFailures_ = nConsecutiveFailures;}
    void setSolutionStatus(Status solutionStatus)
      {solutionStatus_ = solutionStatus;}
    void setOutput(bool output) {output_ = output;}
    void setOutputScreen(bool outputScreen) {outputScreen_ = outputScreen;}
    void setIsRestartable(bool isRestartable) {isRestartable_=isRestartable;}
    void setIsInterpolated(bool isInterpolated)
      {isInterpolated_ = isInterpolated;}
    void setAccuracy(Scalar accuracy) {accuracy_ = accuracy;}
  //@}

  /// \name Overridden from Teuchos::Describable
  //@{
    virtual std::string description() const;
    virtual void describe(Teuchos::FancyOStream          &out,
                          const Teuchos::EVerbosityLevel verbLevel) const;
  //@}

protected:
  Scalar time_;              ///< Time of solution
  int    iStep_;             ///< Time step index for this solution
  Scalar dt_;                ///< Time step for this solution
  Scalar errorAbs_;          ///< Absolute local truncation error
  Scalar errorRel_;          ///< Relative local truncation error
  int order_;                ///< Order of this solution
  int nFailures_;            ///< Total number of stepper failures
  int nConsecutiveFailures_; ///< Consecutive number of stepper failures

  /** The solutionStatus is used to indicate
      - if the solution is still being worked on; WORKING
      - if the solution is accepted and completed (e.g., past solutions
        in SolutionHistory); PASSED.
      - if the time step has FAILED.  This may be caused by the Stepper
        failing, or Integrator not accepting the time step.
  */
  Status solutionStatus_;
  bool   output_;            ///< SolutionState should be or has been outputted
  bool   outputScreen_;      ///< Output screen dump
  bool   isRestartable_;     ///< T - soln can be used as a restart
  bool   isInterpolated_;    ///< F - soln is time integrated; T - soln is interpolated
  Scalar accuracy_;          ///< Interpolation accuracy of solution

};
} // namespace Tempus

#endif // TEMPUS_SOLUTIONSTATEMETADATA_HPP
