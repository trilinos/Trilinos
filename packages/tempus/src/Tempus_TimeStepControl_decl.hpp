// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef Tempus_TimeStepControl_decl_hpp
#define Tempus_TimeStepControl_decl_hpp

// Teuchos
#include "Teuchos_VerboseObject.hpp"
#include "Teuchos_Describable.hpp"
#include "Teuchos_ParameterListAcceptorDefaultBase.hpp"
// Tempus
#include "Tempus_config.hpp"
#include "Tempus_SolutionHistory.hpp"
#include "Tempus_TimeStepControlStrategyComposite.hpp"

#include <iostream>
#include <iterator>
#include <sstream>


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

  /// Constructor
  TimeStepControl(Teuchos::RCP<Teuchos::ParameterList> pList = Teuchos::null);

  /// This is a copy constructor
  TimeStepControl(const TimeStepControl<Scalar>& tsc);

  /// Destructor
  virtual ~TimeStepControl() {}

  virtual void initialize(Teuchos::RCP<Teuchos::ParameterList> pList =
    Teuchos::null) { this->setParameterList(pList); }

  /** \brief Determine the time step size.*/
  virtual void getNextTimeStep(
    const Teuchos::RCP<SolutionHistory<Scalar> > & solutionHistory,
    Status & integratorStatus);

  /** \brief Check if time is within minimum and maximum time. */
  virtual bool timeInRange(const Scalar time) const;

  /** \brief Check if time step index is within minimum and maximum index. */
  virtual bool indexInRange(const int iStep) const;

  /** \brief Set the TimeStepControlStrategy. */
  virtual void setTimeStepControlStrategy(
        Teuchos::RCP<TimeStepControlStrategy<Scalar> > tscs = Teuchos::null);

  /// \name Overridden from Teuchos::ParameterListAccepto{}
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

  /// \name Get ParameterList values
  //@{
    virtual Scalar getInitTime() const
      { return tscPL_->get<double>("Initial Time"); }
    virtual Scalar getFinalTime() const
      { return tscPL_->get<double>("Final Time"); }
    virtual Scalar getMinTimeStep() const
      { return tscPL_->get<double>("Minimum Time Step"); }
    virtual Scalar getInitTimeStep() const
      { return tscPL_->get<double>("Initial Time Step"); }
    virtual Scalar getMaxTimeStep() const
      { return tscPL_->get<double>("Maximum Time Step"); }
    virtual int getInitIndex() const
      { return tscPL_->get<int>   ("Initial Time Index"); }
    virtual int getFinalIndex() const
      { return tscPL_->get<int>   ("Final Time Index"); }
    virtual Scalar getMaxAbsError() const
      { return tscPL_->get<double>("Maximum Absolute Error"); }
    virtual Scalar getMaxRelError() const
      { return tscPL_->get<double>("Maximum Relative Error"); }
    virtual int getMinOrder() const
      { return tscPL_->get<int>   ("Minimum Order"); }
    virtual int getInitOrder() const
      { return tscPL_->get<int>   ("Initial Order"); }
    virtual int getMaxOrder() const
      { return tscPL_->get<int>   ("Maximum Order"); }
    virtual std::string getStepType() const
      { return tscPL_->get<std::string>("Integrator Step Type"); }
    virtual std::vector<int> getOutputIndices() const
      { return outputIndices_; }
    virtual std::vector<Scalar> getOutputTimes() const
      { return outputTimes_; }
    virtual int getMaxFailures() const
      { return tscPL_->get<int>("Maximum Number of Stepper Failures"); }
    virtual int getMaxConsecFailures() const
      { return tscPL_->
               get<int>("Maximum Number of Consecutive Stepper Failures"); }
    virtual int getNumTimeSteps() const
      { return tscPL_->get<int>("Number of Time Steps"); }
    virtual Teuchos::RCP<TimeStepControlStrategyComposite<Scalar>>
       getTimeStepControlStrategy() const { return stepControlStrategy_;}
  //@}

  /// \name Set ParameterList values
  //@{
    virtual void setInitTime(Scalar InitTime)
      { tscPL_->set<double>("Initial Time"             , InitTime    ); }
    virtual void setFinalTime(Scalar FinalTime)
      { tscPL_->set<double>("Final Time"               , FinalTime   ); }
    virtual void setMinTimeStep(Scalar MinTimeStep)
      { tscPL_->set<double>("Minimum Time Step"        , MinTimeStep ); }
    virtual void setInitTimeStep(Scalar InitTimeStep)
      { tscPL_->set<double>("Initial Time Step"        , InitTimeStep); }
    virtual void setMaxTimeStep(Scalar MaxTimeStep)
      { tscPL_->set<double>("Maximum Time Step"        , MaxTimeStep ); }
    virtual void setInitIndex(int InitIndex)
      { tscPL_->set<int>   ("Initial Time Index"       , InitIndex   ); }
    virtual void setFinalIndex(int FinalIndex)
      { tscPL_->set<int>   ("Final Time Index"         , FinalIndex  ); }
    virtual void setMaxAbsError(Scalar MaxAbsError)
      { tscPL_->set<double>("Maximum Absolute Error"   , MaxAbsError ); }
    virtual void setMaxRelError(Scalar MaxRelError)
      { tscPL_->set<double>("Maximum Relative Error"   , MaxRelError ); }
    virtual void setMinOrder(int MinOrder)
     { tscPL_->set<int>   ("Minimum Order"             , MinOrder    ); }
    virtual void setInitOrder(int InitOrder)
      { tscPL_->set<int>   ("Initial Order"            , InitOrder   ); }
    virtual void setMaxOrder(int MaxOrder)
      { tscPL_->set<int>   ("Maximum Order"            , MaxOrder    ); }
    virtual void setStepType(std::string StepType)
      { tscPL_->set<std::string>("Integrator Step Type", StepType    ); }
    virtual void setOutputIndices(std::vector<int> OutputIndices)
      { outputIndices_ = OutputIndices;
        std::ostringstream ss;
        std::copy(OutputIndices.begin(), OutputIndices.end()-1,
                  std::ostream_iterator<int>(ss, ","));
        ss << OutputIndices.back();
        tscPL_->set<std::string>("Output Index List", ss.str());
      }
    virtual void setOutputTimes(std::vector<Scalar> OutputTimes)
      { outputTimes_ = OutputTimes;
        std::ostringstream ss;
        std::copy(OutputTimes.begin(), OutputTimes.end()-1,
                  std::ostream_iterator<Scalar>(ss, ","));
        ss << OutputTimes.back();
        tscPL_->set<std::string>("Output Time List", ss.str());
      }
    virtual void setMaxFailures(int MaxFailures)
      { tscPL_->set<int>("Maximum Number of Stepper Failures", MaxFailures); }
    virtual void setMaxConsecFailures(int MaxConsecFailures)
      { tscPL_->set<int>
        ("Maximum Number of Consecutive Stepper Failures", MaxConsecFailures); }
    virtual void setNumTimeSteps(int numTimeSteps);
  //@}

protected:

  Teuchos::RCP<Teuchos::ParameterList> tscPL_;

  std::vector<int>    outputIndices_;  ///< Vector of output indices.
  std::vector<Scalar> outputTimes_;    ///< Vector of output times.

  bool outputAdjustedDt_; ///< Flag indicating that dt was adjusted for output.
  Scalar dtAfterOutput_;  ///< dt to reinstate after output step.

  Teuchos::RCP<TimeStepControlStrategyComposite<Scalar>> stepControlStrategy_;

};
} // namespace Tempus

#endif // Tempus_TimeStepControl_decl_hpp
