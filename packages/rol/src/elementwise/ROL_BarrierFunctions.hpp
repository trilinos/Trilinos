#ifndef ROL_BARRIER_FUNCTIONS_H
#define ROL_BARRIER_FUNCTIONS_H

#include "ROL_UnaryFunctions.hpp"

#include "Teuchos_ParameterList.hpp"


namespace ROL {
namespace Elementwise {

/**
  
 */
template<class Real> 
class BarrierFunction {

private:

  const Teuchos::RCP<const UnaryFunction<Real> > value_;
  const Teuchos::RCP<const UnaryFunction<Real> > deriv_;
  const Teuchos::RCP<const UnaryFunction<Real> > deriv2_;

public:
 
  BarrierFunction( const Teuchos::RCP<const UnaryFunction<Real> > &value,
                   const Teuchos::RCP<const UnaryFunction<Real> > &deriv,
                   const Teuchos::RCP<const UnaryFunction<Real> > &deriv2 ) :
                   value_(value), deriv_(deriv), deriv2_(deriv2) { 
                   }

  Teuchos::RCP<UnaryFunction<Real> > getValue(void) const {
    return value_;
  }

  Teuchos::RCP<UnaryFunction<Real> > getDerivative(void) const {
    return deriv_;
  }
  
  Teuchos::RCP<UnaryFunction<Real> > getSecondDerivative(void) const {
    return deriv2_;
  }

};



template<class Real>
class BarrierFunctionFactory {

private:

  enum EBarrierType {
    BARRIER_LOGARITHM = 0,
    BARRIER_QUADRATIC,
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
    for( EBarrierType to = BARRIER_LOGARITHM; to < BARRIER_LAST; ++to ) {
      if( !s.compare(removeStringFormat(EBarrierToString(to))) ) {
        return to;
      }
    }
    return BARRIER_LOGARITHM;
  }

  Teuchos::RCP<BarrierFunction<Real> > barrier_;

public:
  BarrierFunctionFactory( Teuchos::ParameterList &parlist ) {

    Teuchos::ParameterList &bflist = parlist.sublist("Barrier Function");
    std::string bfstring = bflist.get("Type","Logarithmic");

    eBarrierType_ = StringToEBarrierType(bfstring); 
     
  }

  Teuchos::RCP<BarrierFunction<Real> > getBarrierFunction( void ) const {

    Teuchos::RCP<UnaryFunction<Real> > value;
    Teuchos::RCP<UnaryFunction<Real> > deriv;
    Teuchos::RCP<UnaryFunction<Real> > deriv2;

    Teuchos::RCP<UnaryFunction<Real> > sqr = Teuchos::rcp(new Power<Real>(2.0));

    switch( eBarrierType_ ) {

      case BARRIER_LOGARITHM:  
      {
        Teuchos::RCP<UnaryFunction<Real> > log = Teuchos::rcp(new Logarithm<Real>);
        Teuchos::RCP<UnaryFunction<Real> > inv = Teuchos::rcp(new Reciprocal<Real>);  
        Teuchos::RCP<UnaryFunction<Real> > neg = Teuchos::rcp(new Scale<Real>(-1.0));

        value  = Teuchos::rcp(new Composition<Real>(neg,log));
        deriv  = Teuchos::rcp(new Composition<Real>(neg,inv));
        deriv2 = Teuchos::rcp(new Composition<Real>(sqr,inv)); 
        
        break;
      }  
      case BARRIER_QUADRATIC:
      {
        Teuchos::RCP<UnaryFunction<Real> > thl = Teuchos::rcp(new ThresholdLower<Real>(0.0) );
        Teuchos::RCP<UnaryFunction<Real> > two = Teuchos::rcp(new Scale<Real>(2.0) );
 
        class Indicator : public UnaryFunction<Real> {
        public:
          Real apply( const Real &x ) const {
            Real value;
            if( x<0.0 ) {
              value = 2.0;
            }
            else if( x==0.0 ) {
              value = 1.0;
            }
            else {
              value = 0.0;
            } 
          } 
        }; // class Indicator

        value  = Teuchos::rcp(new Composition<Real>(sqr,thl) );
        deriv  = Teuchos::rcp(new Composition<Real>(two,thl) );
        deriv2 = Teuchos::rcp(new Indicator);
   
        break;
      }
      default:
      {
        TEUCHOS_TEST_FOR_EXCEPTION(true,std::invalid_argument,
          ">>> (BarrierFunctionFactory::getBarrierFunction): Undefined barrier function type!");
      }
    } 
    return Teuchos::rcp( new BarrierFunction<Real>(value, deriv, deriv2) );
  } 


}; // class BarrierFunctionFactory







} // namespace Elementwise
} // namespace ROL











#endif // ROL_BARRIER_FUNCTIONS_H
