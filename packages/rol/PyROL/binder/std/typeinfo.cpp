#include <sstream> // __str__
#include <typeinfo>

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

// std::type_info file:typeinfo line:88
struct PyCallBack_std_type_info : public std::type_info {
	using std::type_info::type_info;

	bool __is_pointer_p() const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const std::type_info *>(this), "__is_pointer_p");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return type_info::__is_pointer_p();
	}
	bool __is_function_p() const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const std::type_info *>(this), "__is_function_p");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return type_info::__is_function_p();
	}
};

// std::bad_cast file:typeinfo line:187
struct PyCallBack_std_bad_cast : public std::bad_cast {
	using std::bad_cast::bad_cast;

	const char * what() const noexcept override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const std::bad_cast *>(this), "what");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<const char *>::value) {
				static pybind11::detail::override_caster_t<const char *> caster;
				return pybind11::detail::cast_ref<const char *>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<const char *>(std::move(o));
		}
		return bad_cast::what();
	}
};

void bind_std_typeinfo(std::function< pybind11::module &(std::string const &namespace_) > &M)
{
	{ // std::type_info file:typeinfo line:88
		pybind11::class_<std::type_info, Teuchos::RCP<std::type_info>, PyCallBack_std_type_info> cl(M("std"), "type_info", "", pybind11::module_local());
		cl.def("name", (const char * (std::type_info::*)() const) &std::type_info::name, "C++: std::type_info::name() const --> const char *", pybind11::return_value_policy::automatic);
		cl.def("before", (bool (std::type_info::*)(const class std::type_info &) const) &std::type_info::before, "C++: std::type_info::before(const class std::type_info &) const --> bool", pybind11::arg("__arg"));
		cl.def("__eq__", (bool (std::type_info::*)(const class std::type_info &) const) &std::type_info::operator==, "C++: std::type_info::operator==(const class std::type_info &) const --> bool", pybind11::arg("__arg"));
		cl.def("__ne__", (bool (std::type_info::*)(const class std::type_info &) const) &std::type_info::operator!=, "C++: std::type_info::operator!=(const class std::type_info &) const --> bool", pybind11::arg("__arg"));
		cl.def("hash_code", (unsigned long (std::type_info::*)() const) &std::type_info::hash_code, "C++: std::type_info::hash_code() const --> unsigned long");
		cl.def("__is_pointer_p", (bool (std::type_info::*)() const) &std::type_info::__is_pointer_p, "C++: std::type_info::__is_pointer_p() const --> bool");
		cl.def("__is_function_p", (bool (std::type_info::*)() const) &std::type_info::__is_function_p, "C++: std::type_info::__is_function_p() const --> bool");
	}
	{ // std::bad_cast file:typeinfo line:187
		pybind11::class_<std::bad_cast, Teuchos::RCP<std::bad_cast>, PyCallBack_std_bad_cast, std::exception> cl(M("std"), "bad_cast", "", pybind11::module_local());
		cl.def( pybind11::init( [](){ return new std::bad_cast(); }, [](){ return new PyCallBack_std_bad_cast(); } ) );
		cl.def( pybind11::init( [](PyCallBack_std_bad_cast const &o){ return new PyCallBack_std_bad_cast(o); } ) );
		cl.def( pybind11::init( [](std::bad_cast const &o){ return new std::bad_cast(o); } ) );
		cl.def("what", (const char * (std::bad_cast::*)() const) &std::bad_cast::what, "C++: std::bad_cast::what() const --> const char *", pybind11::return_value_policy::automatic);
		cl.def("assign", (class std::bad_cast & (std::bad_cast::*)(const class std::bad_cast &)) &std::bad_cast::operator=, "C++: std::bad_cast::operator=(const class std::bad_cast &) --> class std::bad_cast &", pybind11::return_value_policy::automatic, pybind11::arg(""));
	}
}
