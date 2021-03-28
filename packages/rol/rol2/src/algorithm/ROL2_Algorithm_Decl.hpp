#pragma once
#ifndef ROL2_ALGORITHM_DECL_HPP
#define ROL2_ALGORITHM_DECL_HPP

namespace ROL2 {

/** \class ROL2::Algorithm
    \brief Abstract interface for optimization algorithms
*/
template<class Real>
class Algorithm {
public:

  //-------------------------------------------------------------------
  /** \enum ROL2::Algorithm::ExitStatus
      \brief Encodes the reason why an algorithm has terminated
  */
  enum class ExitStatus : std::int16_t {
    Converged = 0,
    MaxIter,
    StepTolMet,
    EncounteredNaN,
    UserDefined,
    Last
  };

  static EnumMap<ExitStatus> exitstatus_dict;

  //-------------------------------------------------------------------
  /** \class ROL2::Algorithm::State
      \brief Basic interface for the container types used by algorithms
  */
  class State {
  public:

    State() = default;
    virtual ~State() = default;

    virtual void reset() = 0;

    /** \brief Returns true if vector memory has been allocated
    */
    virtual bool is_initialized() const;

    Ptr<Vector<Real>> iterateVec_;
    Ptr<Vector<Real>> minIterVec_;
    Real value_     = ROL_INF<Real>;
    Real minValue_  = ROL_INF<Real>;
    Real gnorm_     = ROL_INF<Real>;
    Real snorm_     = ROL_INF<Real>;
    int iter_;
    int minIter_;
    int nfval_;
    int ngrad_;
    bool flag_      = false;
    ExitStatus statusFlag_ = ExitStatus::Last;

  }; // ROL2::Algorithm<Real>::State

  //-------------------------------------------------------------------
  // Public Interface
 
  Algorithm() = default;

  virtual ~Algorithm() = default;

  virtual void setState( const Ptr<State>& state ) {}

  virtual void setStatusTest( const Ptr<StatusTest<Real>>& status, 
                                    bool                   combineStatus = false ) = 0;

  /** \brief Apply algorithm to an optimization problem
  */ 
//  virtual void run( OptimizationProblem<Real>&, std::ostream& ) = 0;             

  virtual void writeHeader( std::ostream& ) const = 0;

  virtual void writeName( std::ostream& os ) const { os << "ROL >> "; } 

  virtual void writeOutput( std::ostream&, bool ) const = 0;

  virtual const State&            getState()  const = 0;
  virtual const StatusTest<Real>& getStatus() const = 0;

  virtual void reset();

protected:

  void initialize( const Vector<Real>& );

  virtual State&            getState()  = 0;
  virtual StatusTest<Real>& getStatus() = 0;
  
}; // Algorithm

template<class Real>
EnumMap<typename Algorithm<Real>::ExitStatus>
Algorithm<Real>::exitstatus_dict = { "Converged",
                                     "Iteration Limit Exceeded", 
                                     "Step Tolerance Met", 
                                     "Step and/or Gradient Returned NaN", 
                                     "User Defined" };
                                     


} // namesapce ROL2

#endif // ROL2_ALGORITHM_HPP

