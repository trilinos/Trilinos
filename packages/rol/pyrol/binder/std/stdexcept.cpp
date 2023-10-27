#include <iterator>
#include <memory>
#include <sstream> // __str__
#include <stdexcept>
#include <string>

#include <functional>
#include <pybind11/pybind11.h>
#include <string>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include <Teuchos_RCP.hpp>


#ifndef BINDER_PYBIND11_TYPE_CASTER
	#define BINDER_PYBIND11_TYPE_CASTER
	PYBIND11_DECLARE_HOLDER_TYPE(T, Teuchos::RCP<T>)
	PYBIND11_DECLARE_HOLDER_TYPE(T, T*)
	PYBIND11_MAKE_OPAQUE(Teuchos::RCP<void>)
#endif

// std::logic_error file:stdexcept line:113
struct PyCallBack_std_logic_error : public std::logic_error {
	using std::logic_error::logic_error;

	const char * what() const noexcept override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const std::logic_error *>(this), "what");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<const char *>::value) {
				static pybind11::detail::override_caster_t<const char *> caster;
				return pybind11::detail::cast_ref<const char *>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<const char *>(std::move(o));
		}
		return logic_error::what();
	}
};

// std::invalid_argument file:stdexcept line:158
struct PyCallBack_std_invalid_argument : public std::invalid_argument {
	using std::invalid_argument::invalid_argument;

	const char * what() const noexcept override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const std::invalid_argument *>(this), "what");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<const char *>::value) {
				static pybind11::detail::override_caster_t<const char *> caster;
				return pybind11::detail::cast_ref<const char *>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<const char *>(std::move(o));
		}
		return logic_error::what();
	}
};

// std::runtime_error file:stdexcept line:197
struct PyCallBack_std_runtime_error : public std::runtime_error {
	using std::runtime_error::runtime_error;

	const char * what() const noexcept override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const std::runtime_error *>(this), "what");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<const char *>::value) {
				static pybind11::detail::override_caster_t<const char *> caster;
				return pybind11::detail::cast_ref<const char *>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<const char *>(std::move(o));
		}
		return runtime_error::what();
	}
};

void bind_std_stdexcept(std::function< pybind11::module &(std::string const &namespace_) > &M)
{
	{ // std::logic_error file:stdexcept line:113
		pybind11::class_<std::logic_error, Teuchos::RCP<std::logic_error>, PyCallBack_std_logic_error, std::exception> cl(M("std"), "logic_error", "", pybind11::module_local());
		cl.def( pybind11::init<const std::string &>(), pybind11::arg("__arg") );

		cl.def( pybind11::init<const char *>(), pybind11::arg("") );

		cl.def( pybind11::init( [](PyCallBack_std_logic_error const &o){ return new PyCallBack_std_logic_error(o); } ) );
		cl.def( pybind11::init( [](std::logic_error const &o){ return new std::logic_error(o); } ) );
		cl.def("what", (const char * (std::logic_error::*)() const) &std::logic_error::what, "C++: std::logic_error::what() const --> const char *", pybind11::return_value_policy::automatic);
		cl.def("assign", (class std::logic_error & (std::logic_error::*)(const class std::logic_error &)) &std::logic_error::operator=, "C++: std::logic_error::operator=(const class std::logic_error &) --> class std::logic_error &", pybind11::return_value_policy::automatic, pybind11::arg(""));
	}
	{ // std::invalid_argument file:stdexcept line:158
		pybind11::class_<std::invalid_argument, Teuchos::RCP<std::invalid_argument>, PyCallBack_std_invalid_argument, std::logic_error> cl(M("std"), "invalid_argument", "", pybind11::module_local());
		cl.def( pybind11::init<const std::string &>(), pybind11::arg("__arg") );

		cl.def( pybind11::init<const char *>(), pybind11::arg("") );

		cl.def( pybind11::init( [](PyCallBack_std_invalid_argument const &o){ return new PyCallBack_std_invalid_argument(o); } ) );
		cl.def( pybind11::init( [](std::invalid_argument const &o){ return new std::invalid_argument(o); } ) );
		cl.def("assign", (class std::invalid_argument & (std::invalid_argument::*)(const class std::invalid_argument &)) &std::invalid_argument::operator=, "C++: std::invalid_argument::operator=(const class std::invalid_argument &) --> class std::invalid_argument &", pybind11::return_value_policy::automatic, pybind11::arg(""));
	}
	{ // std::runtime_error file:stdexcept line:197
		pybind11::class_<std::runtime_error, Teuchos::RCP<std::runtime_error>, PyCallBack_std_runtime_error, std::exception> cl(M("std"), "runtime_error", "", pybind11::module_local());
		cl.def( pybind11::init<const std::string &>(), pybind11::arg("__arg") );

		cl.def( pybind11::init<const char *>(), pybind11::arg("") );

		cl.def( pybind11::init( [](PyCallBack_std_runtime_error const &o){ return new PyCallBack_std_runtime_error(o); } ) );
		cl.def( pybind11::init( [](std::runtime_error const &o){ return new std::runtime_error(o); } ) );
		cl.def("what", (const char * (std::runtime_error::*)() const) &std::runtime_error::what, "C++: std::runtime_error::what() const --> const char *", pybind11::return_value_policy::automatic);
		cl.def("assign", (class std::runtime_error & (std::runtime_error::*)(const class std::runtime_error &)) &std::runtime_error::operator=, "C++: std::runtime_error::operator=(const class std::runtime_error &) --> class std::runtime_error &", pybind11::return_value_policy::automatic, pybind11::arg(""));
	}
}
