#include <ROL_ScalarTraits.hpp>
#include <sstream> // __str__

#include <functional>
#include <pybind11/pybind11.h>
#include <string>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>


#ifndef BINDER_PYBIND11_TYPE_CASTER
	#define BINDER_PYBIND11_TYPE_CASTER
	PYBIND11_DECLARE_HOLDER_TYPE(T, std::shared_ptr<T>)
	PYBIND11_DECLARE_HOLDER_TYPE(T, T*)
	PYBIND11_MAKE_OPAQUE(std::shared_ptr<void>)
#endif

void bind_ROL_ScalarTraits(std::function< pybind11::module &(std::string const &namespace_) > &M)
{
	{ // ROL::ScalarTraits_Magnitude file:ROL_ScalarTraits.hpp line:55
		pybind11::class_<ROL::ScalarTraits_Magnitude<double>, std::shared_ptr<ROL::ScalarTraits_Magnitude<double>>> cl(M("ROL"), "ScalarTraits_Magnitude_double_t", "", pybind11::module_local());
		cl.def( pybind11::init( [](){ return new ROL::ScalarTraits_Magnitude<double>(); } ) );
	}
	{ // ROL::ScalarTraits file:ROL_ScalarTraits.hpp line:66
		pybind11::class_<ROL::ScalarTraits<double>, std::shared_ptr<ROL::ScalarTraits<double>>> cl(M("ROL"), "ScalarTraits_double_t", "", pybind11::module_local());
		cl.def( pybind11::init( [](){ return new ROL::ScalarTraits<double>(); } ) );
		cl.def_static("zero", (double (*)()) &ROL::ScalarTraits<double>::zero, "C++: ROL::ScalarTraits<double>::zero() --> double");
		cl.def_static("half", (double (*)()) &ROL::ScalarTraits<double>::half, "C++: ROL::ScalarTraits<double>::half() --> double");
		cl.def_static("one", (double (*)()) &ROL::ScalarTraits<double>::one, "C++: ROL::ScalarTraits<double>::one() --> double");
		cl.def_static("two", (double (*)()) &ROL::ScalarTraits<double>::two, "C++: ROL::ScalarTraits<double>::two() --> double");
		cl.def_static("eps", (double (*)()) &ROL::ScalarTraits<double>::eps, "C++: ROL::ScalarTraits<double>::eps() --> double");
		cl.def_static("rmin", (double (*)()) &ROL::ScalarTraits<double>::rmin, "C++: ROL::ScalarTraits<double>::rmin() --> double");
		cl.def_static("rmax", (double (*)()) &ROL::ScalarTraits<double>::rmax, "C++: ROL::ScalarTraits<double>::rmax() --> double");
		cl.def_static("two_pi", (double (*)()) &ROL::ScalarTraits<double>::two_pi, "C++: ROL::ScalarTraits<double>::two_pi() --> double");
		cl.def_static("pi", (double (*)()) &ROL::ScalarTraits<double>::pi, "C++: ROL::ScalarTraits<double>::pi() --> double");
		cl.def_static("half_pi", (double (*)()) &ROL::ScalarTraits<double>::half_pi, "C++: ROL::ScalarTraits<double>::half_pi() --> double");
		cl.def_static("quarter_pi", (double (*)()) &ROL::ScalarTraits<double>::quarter_pi, "C++: ROL::ScalarTraits<double>::quarter_pi() --> double");
		cl.def_static("sqrt_two_pi", (double (*)()) &ROL::ScalarTraits<double>::sqrt_two_pi, "C++: ROL::ScalarTraits<double>::sqrt_two_pi() --> double");
		cl.def_static("sqrt_pi", (double (*)()) &ROL::ScalarTraits<double>::sqrt_pi, "C++: ROL::ScalarTraits<double>::sqrt_pi() --> double");
		cl.def_static("sqrt_half_pi", (double (*)()) &ROL::ScalarTraits<double>::sqrt_half_pi, "C++: ROL::ScalarTraits<double>::sqrt_half_pi() --> double");
		cl.def_static("sqrt_two", (double (*)()) &ROL::ScalarTraits<double>::sqrt_two, "C++: ROL::ScalarTraits<double>::sqrt_two() --> double");
	}
}
