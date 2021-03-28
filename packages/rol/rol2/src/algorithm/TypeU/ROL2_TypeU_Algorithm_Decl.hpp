#pragma once
#ifndef ROL2_TYPEU_ALGORITHM_DECL_HPP
#define ROL2_TYPEU_ALGORITHM_DECL_HPP

namespace ROL2 {
namespace TypeU {

/** \class ROL2::TypeU::Algorithm
    \brief Common interface for unconstrained optimization problems.
*/

template<class Real>
class Algorithm : public ROL2::Algorithm<Real> {
public:

  using typename ROL2::Algorithm<Real>::ExitStatus;

  enum class Type : std::int16_t {
    Bundle = 0,
    LineSearch,
    TrustRegion,
    Last
  };

  static EnumMap<Type> type_dict;

  /** \class ROL2::Algorithm_U::State 
      \brief Common container for quantities used by unconstrained algorithms.
  */
  struct State : public ROL2::Algorithm<Real>::State {

    State() = default;

    virtual ~State() = default;

    virtual void reset() override;

    virtual bool is_initialized() const override;

    Ptr<Vector<Real>> stepVec_;        ///< Step Vector
    Ptr<Vector<Real>> gradientVec_;    ///< Gradient vector

    Real searchSize_ = ROL_ONE<Real>;  /**< Characteristic size used in optimization algorithms. For example,
                                            the step length of a line search or the radius of a trust region. */

    ExitStatus exitStatus_ = ExitStatus::Last;

  }; // State

  Algorithm();

  virtual ~Algorithm() = default;

  virtual void setState( const Ptr<State>& state ) {}// { state_ = state; }

  virtual void setStatusTest( const Ptr<StatusTest<Real>>& status,
                                    bool                   combineStatus = false ) override;
//  virtual void run( OptimizationProblem<Real>& problem,
//                    std::ostream&              outStream = std::cout ) override;

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
//  virtual void run( Vector<Real>&     x,
//                    Objective<Real>&  obj,
//                    Constraint<Real>& linear_con,
//                    Vector<Real>&     linear_mul,
//                    std::ostream&     outStream = std::cout );

  /** \brief Run algorithm on unconstrained problems with explicit
             linear equality constraints (Type-U).
             This general interface supports the use of dual optimization vector spaces,
             where the user does not define the dual() method.
  */
//  virtual void run(       Vector<Real>&     x,
//                    const Vector<Real>&     g,
//                          Objective<Real>&  obj,
//                          Constraint<Real>& linear_con,
//                          Vector<Real>&     linear_mul,
//                    const Vector<Real>&     linear_c,
//                          std::ostream&     outStream = std::cout );

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

  virtual const State&                    getState()  const override { return *state_; }
  virtual const CombinedStatusTest<Real>& getStatus() const override { return *status_; }

  virtual void reset() override;

protected:

  void initialize( const Vector<Real>& x,
                   const Vector<Real>& g );

  virtual State&                    getState()  override { return *state_;  }
  virtual CombinedStatusTest<Real>& getStatus() override { return *status_; }

private:

  Ptr<State>                    state_;
  Ptr<CombinedStatusTest<Real>> status_;

}; // class Algorithm

template<class Real> 
EnumMap<typename Algorithm<Real>::Type>
Algorithm<Real>::type_dict = { "Bundle",
                               "Line Search",
                               "Trust Region" };



} // namespace TypeU
} // namespace ROL2

#endif // ROL2_TYPEU_ALGORITHM_DECL_HPP

