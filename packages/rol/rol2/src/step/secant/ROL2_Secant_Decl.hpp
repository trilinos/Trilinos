// @HEADER
// ************************************************************************
//
//               Rapid Optimization Library (ROL2) Package
//                 Copyright (2014) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact lead developers:
//              Drew Kouri   (dpkouri@sandia.gov) and
//              Denis Ridzal (dridzal@sandia.gov)
//
// ************************************************************************
// @HEADER

#pragma once 
#ifndef ROL2_SECANT_DECL_HPP
#define ROL2_SECANT_DECL_HPP

/** \class ROL2::Secant
    \brief Provides interface for and implements limited-memory secant operators.
*/

namespace ROL2 {

template<class Real>
class Secant : public LinearOperator<Real> {
public:

  enum class Mode : std::int16_t {
    Forward = 0,
    Inverse,
    Both
  };

  enum class Type : std::int16_t {
    LBFGS = 0,
    LDFP,
    LSR1,
    BarzilaiBorwein,
    UserDefined,
    Last
  };

  //----------------------------------------------------------
  struct State {
    Ptr<Vector<Real>>              iterate_;
    std::vector<Ptr<Vector<Real>>> iterDiff_; // Step Storage
    std::vector<Ptr<Vector<Real>>> gradDiff_; // Gradient Storage
    std::vector<Real>              product_;  // Step-Gradient Inner Product Storage
    std::vector<Real>              product2_; // Step-Gradient Inner Product Storage
    int storage_;                             // Storage Size
    int current_ = -1;                        // Current Storage Size
    int iter_     = 0;                        // Current Optimization Iteration
    Mode mode_;                               // Intended application mode

    State(int M, Mode sm) : storage(M), current(-1), iter(0), mode(sm) {}
  };  // struct Secant<Real>::State
  //----------------------------------------------------------


  virtual ~Secant() = default;

  // Constructor
  Secant( int  M = 10, 
          bool useDefaultScaling = true, 
          Real Bscaling = Real(1), 
          Mode mode = Mode::Both );

  // Update Secant Approximation
  virtual void updateStorage( const Vector<Real>& x,  
                              const Vector<Real>& grad,
                              const Vector<Real>& gp, 
                              const Vector<Real>& s,
                                    Real          snorm,      
                                    int           iter );

  // Apply Secant Approximate Inverse Hessian
  virtual void applyH(       Vector<Real>& Hv, 
                       const Vector<Real>& v ) const = 0;

  // Apply Initial Secant Approximate Inverse Hessian
  virtual void applyH0(       Vector<Real>& Hv, 
                        const Vector<Real>& v ) const;

  // Apply Secant Approximate Hessian
  virtual void applyB(       Vector<Real>& Bv, 
                       const Vector<Real>& v ) const = 0;

  // Apply Initial Secant Approximate Hessian 
  virtual void applyB0(       Vector<Real>& Bv, 
                        const Vector<Real>& v ) const;

  // Test Secant Approximations 
  void test( std::ostream &stream = std::cout ) const;

  void apply(       Vector<Real>& Hv, 
              const Vector<Real>& v, 
                    Real&         tol) const {
    applyB(Hv,v);
  }

  void applyInverse(       Vector<Real>& Hv, 
                     const Vector<Real>& v, 
                           Real&         tol) const {
    applyH(Hv,v);
  }

  const State& getState() const { return *state_; }

  static Ptr<Secant> create( ParameterList& parlist, Mode mode = Mode::Both );
  static Ptr<Secant> create( Type type, int L = 10, int BBtype = 1 );

  static EnumMap<Type> type_dict;
  static EnumMap<Mode> mode_dict;

protected:

  State& getState() { return *state_; }

  Ptr<Vector<Real>> y_;
  Real              Bscaling_;
  Mode              mode_;
  bool              useDefaultScaling_;
  bool              isInitialized_ = false;

private:
  Ptr<State>        state_; 

}; // Secant

template<typename Real>
inline std::string enumToStr( Secant<Real>::Type e ) {
  return Secant<Real>::type_dict[e];
}

template<class Real> 
EnumMap<Secant<Real>>::Mode>
Secant<Real>::mode_dict = { "Forward",
                            "Inverse",
                            "Both" };

template<class Real> 
EnumMap<Secant<Real>>::Type>
Secant<Real>::type_dict = { "Limited-Memory BFGS",
                            "Limited-Memory DFP",
                            "Limited-Memory SR1", 
                            "Barzilai-Borwein",
                            "User Defined" };


} // namespace ROL2

#endif // ROL2_SECANT_DECL_HPP
