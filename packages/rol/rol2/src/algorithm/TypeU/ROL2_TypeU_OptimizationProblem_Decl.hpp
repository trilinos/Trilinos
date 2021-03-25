#pragma once
#ifndef ROL2_TYPEU_OPTIMIZATIONPROBLEM_DECL_HPP
#define ROL2_TYPEU_OPTIMIZATIONPROBLEM_DECL_HPP

namespace ROL2 {
namespace TypeU {

template<class Real>
class OptimizationProblem {
public:

  struct CheckData { 
    std::vector<Real>              solution;
    std::vector<std::vector<Real>> gradient;
    std::vector<std::vector<Real>> hessVec;
    std::vector<Real>              hessSym;
  }; // CheckData

  virtual ~OptimizationProblem() = default;

  OptimizationProblem() = default;

  OptimizationProblem( const Ptr<Objective<Real>>& obj,
                       const Ptr<Vector<Real>>&    x  );

  const Ptr<Objective<Real>> getObjective();

  const Ptr<Vector<Real>> getSolutionVector();

  void checkObjective( CheckData&    data,
                       Vector<Real>& x,
                       Vector<Real>& u,
                       Vector<Real>& v, 
                       std::ostream& os = std::cout,
                       int           numSteps = UFD::NUM_CHECKDERIV_STEPS,
                       int           fd_order = 1 );

  void checkObjective( Vector<Real>& x,
                       Vector<Real>& u,
                       Vector<Real>& v, 
                       std::ostream& os = std::cout,
                       int           numSteps = UFD::NUM_CHECKDERIV_STEPS,
                       int           fd_order = 1 );

  void checkSolutionVector( CheckData&    data,
                            Vector<Real>& x,
                            Vector<Real>& y,
                            Vector<Real>& v, 
                            std::ostream& os = std::cout );  

  void checkSolutionVector( Vector<Real>& x,
                            Vector<Real>& y,
                            Vector<Real>& v, 
                            std::ostream& os = std::cout );

  void check( CheckData&    data, 
              std::ostream& os = std::cout,
              int           numSteps = UFD::NUM_CHECKDERIV_STEPS,
              int           fd_order = 1 );

  void check( std::ostream& os = std::cout,
              int           numSteps = UFD::NUM_CHECKDERIV_STEPS,
              int           fd_order = 1 );

protected:

  Ptr<Objective<Real>> obj_;
  Ptr<Vector<Real>>    sol_;
 
}; // ROL2::TypeU::OptimizationProblem

} // namespace TypeU
} // namespace ROL2

#endif // ROL2_TYPEU_OPTIMIZATIONPROBLEM_DECL_HPP

