#pragma once
#ifndef ROL2_ALGORITHM_DECL_HPP
#define ROL2_ALGORITHM_DECL_HPP

namespace ROL2 {

/** \class ROL2::Algorithm
    \brief Abstract interface for optimization algorithms
*/
template<class Real>
class Algorithm {

  //-------------------------------------------------------------------
  /** \enum ROL2::Algorithm::ExitStatus
      \brief Encodes the reason why an algorithm has terminated
  */
  enum class ExitStatus : std::int16_t {
    Converged = 0,
    MaxIter,
    StepTolMet,
    EncounteredNan,
    UserDefined,
    Last
  };

  //-------------------------------------------------------------------
  /** \class ROL2::Algorithm::State
      \brief Basic interface for the container types used by algorithms
  */
  class State {
  public:

    State() = default;

    virtual ~State() = default;

    virtual void reset() = 0;

    virtual void initialize( const Vector<Real>& ) = 0;

    /** \brief Returns true if vector memory has been allocated
    */
    virtual bool is_initialized() const = 0;

  }; 

  //-------------------------------------------------------------------
  // Public Interface
 
  Algorithm() = default;

  virtual ~Algorithm() = default;

  virtual void setStatusTest( const std::shared_ptr<StatusTest<Real>>&, bool  ) = 0;

  /** \brief Apply algorithm to an optimization problem
  */ 
  virtual void run( OptimizationProblem<Real>&, std::ostream& ) = 0;             

  virtual void writeHeader( std::ostream& ) const = 0;

  virtual void writeName( std::ostream& os ) const { os << "ROL >> "; } 

  virtual void writeOutput( std::ostream& ) const = 0;

  virtual const State&            getState()  const = 0;
  virtual const StatusTest<Real>& getStatus() const = 0;

  virtual void reset() = 0;

protected:

  virtual State&            getState()  = 0;
  virtual StatusTest<Real>& getStatus() = 0;
  
}; // Algorithm

} // namesapce ROL2

#endif // ROL2_ALGORITHM_HPP

