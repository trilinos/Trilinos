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

#ifndef ROL_OBJECTIVE_FROM_BOUND_CONSTRAINT_H
#define ROL_OBJECTIVE_FROM_BOUND_CONSTRAINT_H

#include "ROL_Objective.hpp"
#include "ROL_BoundConstraint.hpp"

#include "Teuchos_ParameterList.hpp"

namespace ROL {


/** @ingroup func_group
    \class ROL::ObjectiveFromBoundConstraint 
    \brief Create a penalty objective from upper and lower bound vectors
 */

template <class Real> 
class ObjectiveFromBoundConstraint : public Objective<Real> {

  typedef Vector<Real> V;

  typedef Elementwise::Axpy<Real>              Axpy;
  typedef Elementwise::Aypx<Real>              Aypx;
  typedef Elementwise::Scale<Real>             Scale;
  typedef Elementwise::Reciprocal<Real>        Reciprocal;
  typedef Elementwise::Power<Real>             Power;
  typedef Elementwise::Logarithm<Real>         Logarithm;
  typedef Elementwise::Multiply<Real>          Multiply; 
  typedef Elementwise::Heaviside<Real>         Heaviside;
  typedef Elementwise::ThresholdUpper<Real>    ThresholdUpper;
  typedef Elementwise::ThresholdLower<Real>    ThresholdLower;
  typedef Elementwise::ReductionSum<Real>      Sum;
  typedef Elementwise::UnaryFunction<Real>     UnaryFunction;
 


  enum EBarrierType {
    BARRIER_LOGARITHM = 0,
    BARRIER_QUADRATIC,     
    BARRIER_DOUBLEWELL,  
    BARRIER_LAST
  } eBarrierType_; 

  inline std::string EBarrierToString( EBarrierType type ) {
    std::string retString;
    switch(type) {
      case BARRIER_LOGARITHM:
        retString = "Logarithmic"; 
        break;
      case BARRIER_QUADRATIC:
        retString = "Quadratic";
        break;
      case BARRIER_DOUBLEWELL:
        retString = "Double Well";
        break;   
      case BARRIER_LAST:
        retString = "Last Type (Dummy)";
        break;
      default:
        retString = "Invalid EBarrierType";
        break;
    }
    return retString;
  }

  inline EBarrierType StringToEBarrierType( std::string s ) {
    s = removeStringFormat(s);
    EBarrierType type = BARRIER_LOGARITHM;
    for( int to = BARRIER_LOGARITHM; to != BARRIER_LAST; ++to ) {
      type = static_cast<EBarrierType>(to);
      if( !s.compare(removeStringFormat(EBarrierToString(type))) ) {
        return type;
      }
    }
    return type;
  }   




private:
  Teuchos::RCP<V> lo_;
  Teuchos::RCP<V> up_;
  Teuchos::RCP<V> a_;     // scratch vector
  Teuchos::RCP<V> b_;     // scratch vector
  EBarrierType    btype_;

public:
 
  ObjectiveFromBoundConstraint( const BoundConstraint<Real> &bc,
                                Teuchos::ParameterList &parlist ) :
    lo_( bc.getLowerVectorRCP() ),
    up_( bc.getUpperVectorRCP() ) {

    a_ = lo_->clone();
    b_ = up_->clone();

    std::string bfstring = parlist.sublist("Barrier Function").get("Type","Logarithmic");
    btype_ = StringToEBarrierType(bfstring); 
  }

  ObjectiveFromBoundConstraint( const BoundConstraint<Real> &bc ) :
    lo_( bc.getLowerVectorRCP() ),
    up_( bc.getUpperVectorRCP() )
    { 
      a_ = lo_->clone();
      b_ = up_->clone();
    }


  Real value( const Vector<Real> &x, Real &tol ) {

    Teuchos::RCP<UnaryFunction> func;

    switch(btype_) {
      case BARRIER_LOGARITHM:

        a_->set(*lo_);                         // a = l
        a_->applyBinary(Aypx(-1.0),x);         // a = x-l
        a_->applyUnary(Logarithm());           // a = log(x-l)
        a_->applyUnary(Scale(-1.0));           // a = -log(x-l)  
               
        b_->set(x);                            // b = x
        b_->applyBinary(Aypx(-1.0),*up_);      // b = u-x
        b_->applyUnary(Logarithm());           // b = log(u-x)
        b_->applyUnary(Scale(-1.0));           // b = -log(u-x)

        b_->plus(*a_);                         // b = -log(x-l)-log(u-x)

        break;

      case BARRIER_QUADRATIC:
        
        a_->set(*lo_);                         // a = l
        a_->applyBinary(Aypx(-1.0),x);         // a = x-l
        a_->applyUnary(ThresholdLower(0.0));   // a = min(x-l,0)
        a_->applyUnary(Power(2.0));            // a = min(x-l,0)^2

        b_->set(*up_);                         // b = u
        b_->applyBinary(Aypx(-1.0),x);         // b = x-u
        b_->applyUnary(ThresholdUpper(0.0));   // b = max(x-u,0)
        b_->applyUnary(Power(2.0));            // b = max(x-u,0)^2

        b_->plus(*a_);                         // b = min(x-l,0)^2 + max(x-u,0)^2  
  
        break;

      case BARRIER_DOUBLEWELL:

        a_->set(*lo_);                        // a = l
        a_->applyBinary(Aypx(-1.0),x);        // a = x-l
        a_->applyUnary(Power(2.0));           // a = (x-l)^2

        b_->set(x);                           // b = x
        b_->applyBinary(Aypx(-1.0),*up_);     // b = u-x
        b_->applyUnary(Power(2.0));           // b = (u-x)^2
       
        b_->applyBinary(Multiply(),*a_);      // b = (x-l)^2*(u-x)^2 

        break;

      default:
        TEUCHOS_TEST_FOR_EXCEPTION(true,std::invalid_argument,
          ">>>(ObjectiveFromBoundConstraint::value): Undefined barrier function type!"); 

        break;
    }

    Real result = b_->reduce(Sum()); 
    return result;

  }

  void gradient( Vector<Real> &g, const Vector<Real> &x, Real &tol ) {

    switch(btype_) {
      case BARRIER_LOGARITHM:

        a_->set(*lo_);                     // a = l
        a_->applyBinary(Aypx(-1.0),x);     // a = x-l
        a_->applyUnary(Reciprocal());      // a = 1/(x-l)
        a_->applyUnary(Scale(-1.0));       // a = -1/(x-l)  
               
        b_->set(x);                        // b = x
        b_->applyBinary(Aypx(-1.0),*up_);  // b = u-x
        b_->applyUnary(Reciprocal());      // b = 1/(u-x)

        b_->plus(*a_);                     // b = -1/(x-l)+1/(u-x)

        break;

      case BARRIER_QUADRATIC:

        a_->set(*lo_);                         // a = l
        a_->applyBinary(Aypx(-1.0),x);         // a = x-l
        a_->applyUnary(ThresholdLower(0.0));   // a = min(x-l,0)

        b_->set(*up_);                         // b = u
        b_->applyBinary(Aypx(-1.0),x);         // b = x-u
        b_->applyUnary(ThresholdUpper(0.0));   // b = max(x-u,0)

        b_->plus(*a_);                         // b = max(x-u,0) + min(x-l,0)
        b_->scale(2.0);                        // b = 2*max(x-u,0) + 2*min(x-l,0)
        break;

      case BARRIER_DOUBLEWELL:

        a_->set(*lo_);                     // a = l
        a_->applyBinary(Aypx(-1.0),x);     // a = x-l

        b_->set(x);                        // b = x
        b_->applyBinary(Aypx(-1.0),*up_);  // b = u-x
       
        a_->applyBinary(Aypx(-1.0),*b_);   // a = (u-x)-(x-l)
        a_->applyUnary(Scale(2.0));        // a = 2*[(u-x)-(x-l)]          
 
        b_->set(x);                        // b = x
        b_->applyBinary(Aypx(-1.0),*up_);  // b = u-x

        a_->applyBinary(Multiply(),*b_);   // a = [(u-x)-(x-l)]*(u-x);
         
        
        b_->set(*lo_);                     // b = l
        b_->applyBinary(Aypx(-1.0),x);     // b = x-l
        b_->applyBinary(Multiply(),*a_);   // b = [(u-x)-(x-l)]*(u-x)*(x-l) 

        break;

      default:
        TEUCHOS_TEST_FOR_EXCEPTION(true,std::invalid_argument,
          ">>>(ObjectiveFromBoundConstraint::gradient): Undefined barrier function type!"); 

        break;
    }

    g.set(*b_);

  }

  void hessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &x, Real &tol ) {

    switch(btype_) {
      case BARRIER_LOGARITHM:

        a_->set(*lo_);                     // a = l
        a_->applyBinary(Axpy(-1.0),x);     // a = x-l
        a_->applyUnary(Reciprocal());      // a = 1/(x-l)
        a_->applyUnary(Power(2.0));        // a = 1/(x-l)^2  
               
        b_->set(x);                        // b = x
        b_->applyBinary(Axpy(-1.0),*up_);  // b = u-x
        b_->applyUnary(Reciprocal());      // b = 1/(u-x)
        b_->applyUnary(Power(2.0));        // b = 1/(u-x)^2

        b_->plus(*a_);                     // b = 1/(x-l)^2 + 1/(u-x)^2

        break;
  
      case BARRIER_QUADRATIC:
 
        a_->set(*lo_);                     // a = l
        a_->applyBinary(Axpy(-1.0),x);     // a = x-l
        a_->scale(-1.0);                   // a = l-x
        a_->applyUnary(Heaviside());       // a = theta(l-x)

        b_->set(*up_);                     // b = u
        b_->applyBinary(Axpy(-1.0),x);     // b = x-u
        b_->applyUnary(Heaviside());       // b = theta(x-u)
        b_->plus(*a_);                     // b = theta(l-x) + theta(x-u)
        b_->scale(2.0);                    // b = 2*theta(l-x) + 2*theta(x-u)

        break;

      case BARRIER_DOUBLEWELL:

        a_->set(*lo_);                      // a = l
        a_->applyBinary(Axpy(-1.0),x);      // a = x-l

        b_->set(x);                         // b = x
        b_->applyBinary(Axpy(-1.0),*up_);   // b = u-x

        b_->applyBinary(Multiply(),*a_);    // b = (u-x)*(x-l)
        b_->applyUnary(Scale(-8.0));        // b = -8*(u-x)*(x-l)
 
        a_->applyUnary(Power(2.0));         // a = (x-l)^2
        a_->applyUnary(Scale(2.0));         // a = 2*(x-l)^2
        
        b_->applyBinary(Axpy(1.0),*a_);     // b = 2*(x-l)^2-8*(u-x)*(x-l)

        a_->set(x);                         // a = x
        a_->applyBinary(Axpy(-1.0),*up_);   // a = u-x
        a_->applyUnary(Power(2.0));         // a = (u-x)^2
        a_->applyUnary(Scale(2.0));         // a = 2*(u-x)^2

        b_->plus(*a_);                      // b = 2*(u-x)^2-8*(u-x)*(x-l)+2*(x-l)         
         
        
        break;

      default:
        TEUCHOS_TEST_FOR_EXCEPTION(true,std::invalid_argument,
          ">>>(ObjectiveFromBoundConstraint::hessVec): Undefined barrier function type!"); 

        break;
    }
   
    hv.set(v);
    hv.applyBinary(Multiply(),*b_); 
    
  }
  
  // For testing purposes
  Teuchos::RCP<Vector<Real> > getBarrierVector(void) {
    return b_;
  } 


}; // class ObjectiveFromBoundConstraint

}


#endif // ROL_OBJECTIVE_FROM_BOUND_CONSTRAINT_H

