#include <ROL_Types.hpp>
#include <iterator>
#include <memory>
#include <sstream> // __str__
#include <stdexcept>
#include <string>

#include <functional>
#include <pybind11/smart_holder.h>
#include <string>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>


#ifndef BINDER_PYBIND11_TYPE_CASTER
	#define BINDER_PYBIND11_TYPE_CASTER
	PYBIND11_DECLARE_HOLDER_TYPE(T, std::shared_ptr<T>)
	PYBIND11_DECLARE_HOLDER_TYPE(T, T*)
	PYBIND11_MAKE_OPAQUE(std::shared_ptr<void>)
#endif

// ROL::Exception::NotImplemented file:ROL_Types.hpp line:914
struct PyCallBack_ROL_Exception_NotImplemented : public ROL::Exception::NotImplemented {
	using ROL::Exception::NotImplemented::NotImplemented;

	const char * what() const throw() override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::Exception::NotImplemented *>(this), "what");
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

void bind_ROL_Types(std::function< pybind11::module &(std::string const &namespace_) > &M)
{
	{ // ROL::Exception::NotImplemented file:ROL_Types.hpp line:914
		pybind11::class_<ROL::Exception::NotImplemented, std::shared_ptr<ROL::Exception::NotImplemented>, PyCallBack_ROL_Exception_NotImplemented, std::logic_error> cl(M("ROL::Exception"), "NotImplemented", "", pybind11::module_local());
		cl.def( pybind11::init<const std::string &>(), pybind11::arg("what_arg") );

		cl.def( pybind11::init( [](PyCallBack_ROL_Exception_NotImplemented const &o){ return new PyCallBack_ROL_Exception_NotImplemented(o); } ) );
		cl.def( pybind11::init( [](ROL::Exception::NotImplemented const &o){ return new ROL::Exception::NotImplemented(o); } ) );
		cl.def("assign", (class ROL::Exception::NotImplemented & (ROL::Exception::NotImplemented::*)(const class ROL::Exception::NotImplemented &)) &ROL::Exception::NotImplemented::operator=, "C++: ROL::Exception::NotImplemented::operator=(const class ROL::Exception::NotImplemented &) --> class ROL::Exception::NotImplemented &", pybind11::return_value_policy::automatic, pybind11::arg(""));
	}
}
