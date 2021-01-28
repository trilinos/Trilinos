#pragma once
#ifndef ROL2_TYPEU_ALGORITHM_DECL_HPP
#define ROL2_TYPEU_ALGORITHM_DECL_HPP

#include "ROL2_Algorithm_Decl.hpp"

namespace ROL2 {
namespace TypeU {

/** \class ROL2::TypeU::Algorithm
    \brief Common interface for unconstrained optimization problems.
*/

template<class Real>
class Algorithm : public ROL2::Algorithm<Real> {
public:

  using ROL2::Algorithm<Real>::ExitStatus;

  /** \class ROL2::Algorithm_U::State 
      \brief Common container for quantities used by unconstrained algorithms.
  */
  struct State : public ROL2::Algorithm<Real>::State {

    State();

    virtual State() = default;

    virtual void reset() override;

    virtual bool is_initialized() const override;

    virtual void initialize( const Vector<Real>& );

    Ptr<Vector<Real>> iterateVec_;     ///< Optimization vector
    Ptr<Vector<Real>> minIterVec_;     ///< Optimization vector that has produced the minimal objective value
    Ptr<Vector<Real>> stepVec_;        ///< Step Vector
    Ptr<Vector<Real>> gradientVec_;    ///< Gradient vector

    Real value_      = ROL_INF<Real>;  ///< Current objective value
    Real minValue_   = ROL_INF<Real>;  ///< Current minimum objective value computed
    Real gnorm_      = ROL_INF<Real>;  ///< Norm of the gradient vector
    Real snorm_      = ROL_INF<Real>;  ///< Norm of the step vector
    Real searchSize_ = ROL_ONE<Real>;  /**< Characteristic size used in optimization algorithms. For example,
                                            the step length of a line search or the radius of a trust region. */

    int iter    = 0; ///< Current iterate number
    int minIter = 0; ///< Iterate number at which minimal objective value was found
    int nfval   = 0; ///< Number of objective function evaluations
    int ngrad   = 0; ///< Number of gradient evaluations
    int nhvec   = 0; ///< Number of Hessian-vector product evaluations
    
    ExitStatus exitStatus_ = ExitStatus::Last;

    bool is_initialized_ = false;

  }; // State

  Algorithm() = default;
  virtual ~Algorithm() = default;

  virtual void setStatusTest( const Ptr<StatusTest<Real>& status,
                              bool combineStatus = false ) override;

  virtual void run( OptimizationProblem<Real>& problem,
                    std::ostream&              outStream = std::cout ) override;

  /** \brief Run algorithm on unconstrained problems (Type-U).
             This is the primary Type-U interface.
  */
  virtual void run( Vector<Real>&    x,
                    Objective<Real>& obj,
                    std::ostream&    outStream = std::cout );

  /** \brief Run algorithm on unconstrained problems with explicit
             linear equality constraints (Type-U).
             This is the primary Type-U with explicit linear equality
             constraints interface.
  */
  virtual void run( Vector<Real>&     x,
                    Objective<Real>&  obj,
                    Constraint<Real>& linear_con,
                    Vector<Real>&     linear_mul,
                    std::ostream&     outStream = std::cout );

  /** \brief Run algorithm on unconstrained problems with explicit
             linear equality constraints (Type-U).
             This general interface supports the use of dual optimization vector spaces,
             where the user does not define the dual() method.
  */
  virtual voi run(       Vector<Real>&     x,
                   const Vector<Real>&     g,
                         Objective<Real>&  obj,
                         Constraint<Real>& linear_con,
                         Vector<Real>&     linear_mul,
                   const Vector<Real>&     linear_c,
                         std::ostream&     outStream = std::cout );

  /** \brief Run algorithm on unconstrained problems (Type-U).
             This general interface supports the use of dual optimization vector spaces,
             where the user does not define the dual() method.
  */
  virtual void run( Vector<Real>&       x,
                    const Vector<Real>& g,
                    Objective<Real>&    obj,
                    std::ostream&       outStream = std::cout ) = 0;

  virtual void writeHeader( std::ostream& ) const override;

  virtual void writeName( std::ostream& ) const override;

  virtual void writeOutput( std::ostream&, bool ) const override;

  virtual const State&            getState()  const override;
  virtual const StatusTest<Real>& getStatus() const override;

  virtual void reset() override;

protected:

  void initialize( const Vector<Real>& x );

  virtual State&            getState()  override;
  virtual StatusTest<Real>& getStatus() override;

private:

  const Ptr<State>                    state_;
  const Ptr<CombinedStatusTest<Real>> status_;

}; 

} // namespace TypeU
} // namespace ROL2

#endif // ROL2_TYPEU_ALGORITHM_DECL_HPP

