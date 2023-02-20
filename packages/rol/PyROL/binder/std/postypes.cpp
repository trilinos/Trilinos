#include <cwchar>
#include <exception>
#include <ios>
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

// std::exception file:bits/exception.h line:60
struct PyCallBack_std_exception : public std::exception {
	using std::exception::exception;

	const char * what() const throw() override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const std::exception *>(this), "what");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<const char *>::value) {
				static pybind11::detail::override_caster_t<const char *> caster;
				return pybind11::detail::cast_ref<const char *>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<const char *>(std::move(o));
		}
		return exception::what();
	}
};

void bind_std_postypes(std::function< pybind11::module &(std::string const &namespace_) > &M)
{
	{ // std::fpos file:bits/postypes.h line:112
		pybind11::class_<std::fpos<__mbstate_t>, std::shared_ptr<std::fpos<__mbstate_t>>> cl(M("std"), "fpos___mbstate_t_t", "", pybind11::module_local());
		cl.def( pybind11::init( [](){ return new std::fpos<__mbstate_t>(); } ) );
		cl.def( pybind11::init<long>(), pybind11::arg("__off") );

		cl.def( pybind11::init( [](std::fpos<__mbstate_t> const &o){ return new std::fpos<__mbstate_t>(o); } ) );
		cl.def("assign", (class std::fpos<__mbstate_t> & (std::fpos<__mbstate_t>::*)(const class std::fpos<__mbstate_t> &)) &std::fpos<__mbstate_t>::operator=, "C++: std::fpos<__mbstate_t>::operator=(const class std::fpos<__mbstate_t> &) --> class std::fpos<__mbstate_t> &", pybind11::return_value_policy::automatic, pybind11::arg(""));
		cl.def("__iadd__", (class std::fpos<__mbstate_t> & (std::fpos<__mbstate_t>::*)(long)) &std::fpos<__mbstate_t>::operator+=, "C++: std::fpos<__mbstate_t>::operator+=(long) --> class std::fpos<__mbstate_t> &", pybind11::return_value_policy::automatic, pybind11::arg("__off"));
		cl.def("__isub__", (class std::fpos<__mbstate_t> & (std::fpos<__mbstate_t>::*)(long)) &std::fpos<__mbstate_t>::operator-=, "C++: std::fpos<__mbstate_t>::operator-=(long) --> class std::fpos<__mbstate_t> &", pybind11::return_value_policy::automatic, pybind11::arg("__off"));
		cl.def("__add__", (class std::fpos<__mbstate_t> (std::fpos<__mbstate_t>::*)(long) const) &std::fpos<__mbstate_t>::operator+, "C++: std::fpos<__mbstate_t>::operator+(long) const --> class std::fpos<__mbstate_t>", pybind11::arg("__off"));
		cl.def("__sub__", (class std::fpos<__mbstate_t> (std::fpos<__mbstate_t>::*)(long) const) &std::fpos<__mbstate_t>::operator-, "C++: std::fpos<__mbstate_t>::operator-(long) const --> class std::fpos<__mbstate_t>", pybind11::arg("__off"));
		cl.def("__sub__", (long (std::fpos<__mbstate_t>::*)(const class std::fpos<__mbstate_t> &) const) &std::fpos<__mbstate_t>::operator-, "C++: std::fpos<__mbstate_t>::operator-(const class std::fpos<__mbstate_t> &) const --> long", pybind11::arg("__other"));
	}
	{ // std::exception file:bits/exception.h line:60
		pybind11::class_<std::exception, std::shared_ptr<std::exception>, PyCallBack_std_exception> cl(M("std"), "exception", "", pybind11::module_local());
		cl.def( pybind11::init( [](){ return new std::exception(); }, [](){ return new PyCallBack_std_exception(); } ) );
		cl.def( pybind11::init( [](PyCallBack_std_exception const &o){ return new PyCallBack_std_exception(o); } ) );
		cl.def( pybind11::init( [](std::exception const &o){ return new std::exception(o); } ) );
		cl.def("assign", (class std::exception & (std::exception::*)(const class std::exception &)) &std::exception::operator=, "C++: std::exception::operator=(const class std::exception &) --> class std::exception &", pybind11::return_value_policy::automatic, pybind11::arg(""));
		cl.def("what", (const char * (std::exception::*)() const) &std::exception::what, "C++: std::exception::what() const --> const char *", pybind11::return_value_policy::automatic);
	}
}
