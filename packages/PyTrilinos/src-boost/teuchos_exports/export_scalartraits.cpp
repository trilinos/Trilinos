
#include "Teuchos_ScalarTraits.hpp" 

// Boost initialization
#include <boost/python.hpp>
using namespace boost::python;

template <class T>
class STraits
{
public:
	typedef Teuchos::ScalarTraits<T> traits;
	
	static bool my_iscmplx()
	{
		return traits::isComplex;
	}
	static bool my_iscmp()
	{
		return traits::isComparable;
	}
	static bool my_hmp()
	{
		return traits::hasMachineParameters;
	}
	
	static class_< traits > wrap(char* name)
	{

		
		class_< traits > 
		PyTraits = class_< traits >( name , init<>() )
			
		//	static const bool isComplex;
		.add_static_property("isComplex",&my_iscmplx)//	,
//				"Determines if scalar type is complex")

		//	//! Determines if scalar type supports relational operators such as <, >, <=, >=.
		//	static const bool isComparable;
		.add_static_property("isComparable",&my_iscmp)
		
		//    //! Determines if scalar type have machine-specific parameters (i.e. eps(), sfmin(), base(), prec(), t(), rnd(), emin(), rmin(), emax(), rmax() are supported)
		//	static const bool hasMachineParameters;
		.add_static_property("hasMachineParameters",&my_hmp)
					
		//    static inline magnitudeType magnitude(T a)
			.def("magnitude",traits::magnitude,
					"Returns the magnitudeType of the scalar type")
			.staticmethod("magnitude")

		//    static inline T zero()
			.def("zero",traits::zero,
					"Returns representation of zero for this scalar type.")
			.staticmethod("zero")

		//    static inline T one()
			.def("one",traits::one,
					"Returns representation of one for this scalar type.")
			.staticmethod("one")

		//    static inline magnitudeType real(T a) 
			.def("real",traits::real,
					"Returns the real part of the scalar type a.")
			.staticmethod("real")

		//    static inline magnitudeType imag(T a)
			.def("imag",traits::imag,
					"Returns the imaginary part of the scalar type a.")
			.staticmethod("imag")

		//    static inline T conjugate(T a) 
			.def("conjugate",traits::conjugate,
					"Returns the conjugate of the scalar type")
			.staticmethod("conjugate")

		//    static inline void seedrandom(unsigned int s)
			.def("seedrandom",traits::seedrandom,
					"Seed the random number generator returned by random().")
			.staticmethod("seedrandom")
			
		//    static inline T random()                   
			.def("random",traits::random,
					"Returns a random number (between -one() and +one()) of this scalar type.")
			.staticmethod("random")

		//    static inline std::string name()           
			.def("name",traits::name,
					"Returns the name of this scalar type.")
			.staticmethod("name")
		
		//    static inline T squareroot(T x)
			.def("squareroot",traits::squareroot,
					"ReturnsReturns a number of magnitudeType that is \n"
					"the square root of this scalar type x. " )
		 	.staticmethod("squareroot")
		 
		//    //! Returns the result of raising one scalar \c x to the power \c y.
		//    static inline T pow(T x, T y) 
			.def("pow",traits::pow,
					"Returns the result of raising one scalar x to the power y.")
			.staticmethod("pow")
		;
			
		return PyTraits;
			
	}
};
template<class T>
class_< Teuchos::ScalarTraits<T> > addNaninf( class_< Teuchos::ScalarTraits<T> > PyTraits )
{
	typedef Teuchos::ScalarTraits<T> traits;
	
	PyTraits
		//    static inline T nan()
		.def("nan",traits::nan,
				"Returns a number that represents NaN.")
		.staticmethod("nan")
		
		//    static inline bool isnaninf(const T& x)
		.def("isnaninf",traits::isnaninf,
				" Returns True if x is NaN or Inf.")
		.staticmethod("isnaninf")
	;
	return PyTraits;
}

template<class T>
class_< Teuchos::ScalarTraits<T> > addMachineParams( class_< Teuchos::ScalarTraits<T> > PyTraits )
{
	typedef Teuchos::ScalarTraits<T> traits;
	
	PyTraits
		//    static inline magnitudeType eps()   
		.def("eps",traits::eps,
				"Returns relative machine precision.")
		.staticmethod("eps")

		//    static inline magnitudeType sfmin() 
		.def("sfmin",traits::sfmin,
				"Returns safe minimum (sfmin), such that 1/sfmin does not overflow.")
		.staticmethod("sfmin")
				
		//    static inline magnitudeType base()  
		.def("base",traits::base,
				"Returns the base of the machine.")
		.staticmethod("base")

		//    static inline magnitudeType prec()  
		.def("prec",traits::prec,
				"Returns eps*base.")
		.staticmethod("prec")

		//    static inline magnitudeType t()     
		.def("t",traits::t,
				"Returns the number of (base) digits in the mantissa.")
		.staticmethod("t")

		//    static inline magnitudeType rnd()   
		.def("rnd",traits::rnd,
				"Returns 1.0 when rounding occurs in addition, 0.0 otherwise")
		.staticmethod("rnd")
		
		//    static inline magnitudeType emin()  
		.def("emin",traits::emin,
				"Returns the minimum exponent before (gradual) underflow.")
		.staticmethod("emin")
		
		//    static inline magnitudeType rmin()  
		.def("rmin",traits::rmin,
				"Returns the underflow threshold - base^(emin-1).")
		.staticmethod("rmin")
		
		//    static inline magnitudeType emax()  
		.def("emax",traits::emax,
				"Returns the largest exponent before overflow.")
		.staticmethod("emax")
		
		//    static inline magnitudeType rmax()  
		.def("rmax",traits::rmax,
				"Overflow theshold -  (base^emax)*(1-eps)")
		.staticmethod("rmax")
	;
	return PyTraits;

}

void expose_scalartraits()
{
	
	addMachineParams ( 
		addNaninf(
			STraits<float>::wrap("ScalarTraitsFloat") 
		));
	
	addMachineParams (
		addNaninf( 
			STraits<double>::wrap("ScalarTraitsDouble") 
		));
		
	STraits<char>::wrap("ScalarTraitsChar");
	STraits<int>::wrap("ScalarTraitsInt");
}


