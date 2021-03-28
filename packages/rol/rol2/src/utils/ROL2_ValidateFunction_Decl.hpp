#pragma once
#ifndef ROL2_VALIDATEFUNCTION_DECL_HPP
#define ROL2_VALIDATEFUNCTION_DECL_HPP

namespace ROL2 {

/** \file  ROL2_ValidateFunction_Decl.hpp
    \brief Provides a set of tools for validating the behavior of several
           function types that are commonly used in ROL.

           - Finite difference check of derivatives
           - Symmetry check of linear operators
           - Adjoint consistency check for linear operators
           - Inverse identity check for linear operators
*/

template<typename Real>
class ValidateFunction {
public:

  ValidateFunction( int           order = 1,
                    int           numSteps = 13,
                    int           precision = 11,
                    bool          printToStream = true,
                    std::ostream& os = std::cout );

  ValidateFunction( const std::vector<Real>& steps,
                          int                order = 1,
                          int                precision = 11,
                          bool               printToStream = true,
                          std::ostream&      os = std::cout );

  ValidateFunction( ROL2::ParameterList& parlist, 
                    std::ostream&        os = std::cout );

  virtual ~ValidateFunction() = default;

  /* Check the directional derivative of a function that maps vectors to scalars */
  template<typename ValueFunction,
           typename DerivativeFunction,
           typename UpdateFunction,
           typename PrimalVector,
           typename DualVector>
  std::vector<std::vector<Real>>
  scalar_derivative_check( ValueFunction       f_value,
                           DerivativeFunction  f_derivative,
                           UpdateFunction      f_update,
                           const DualVector&   g,
                           const PrimalVector& v,
                           const PrimalVector& x,
                           const std::string&  label ) const;


 /* Check the directional derivative of a function that maps vectors to scalars */
  template<typename ValueFunction,
           typename DerivativeFunction,
           typename UpdateFunction,
           typename DomainVectorType,
           typename RangeVectorType>
  std::vector<std::vector<Real>>
  vector_derivative_check( ValueFunction           f_value,
                           DerivativeFunction      f_derivative,
                           UpdateFunction          f_update,
                           const RangeVectorType&  c,
                           const DomainVectorType& x,
                           const DomainVectorType& v,
                           const std::string&      label ) const;

  std::ostream& getStream() const { return os_; }

private:

  std::vector<Real> steps_;         // Set of step sizes of FD approximation
  std::ostream&     os_;            // Output stream
  int               order_;         // Finite difference order
  int               numSteps_;      // Number of evalutions of different step sizes
  int               width_;         // For print formatting
  int               precision_;     // Number of digits to display
  bool              printToStream_; // False will suppress output

}; // class ValidateFunction

} // namespace ROL2

#endif // ROL2_VALIDATEFUNCTION_DECL_HPP

