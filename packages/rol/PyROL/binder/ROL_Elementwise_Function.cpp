#include <ROL_Elementwise_Function.hpp>
#include <sstream> // __str__

#include <functional>
#include "PyROL_Smart_Holder.hpp"
#include <string>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>


#ifndef BINDER_PYBIND11_TYPE_CASTER
	#define BINDER_PYBIND11_TYPE_CASTER
	PYBIND11_DECLARE_HOLDER_TYPE(T, std::shared_ptr<T>)
	PYBIND11_DECLARE_HOLDER_TYPE(T, T*)
	PYBIND11_MAKE_OPAQUE(std::shared_ptr<void>)
#endif

// ROL::Elementwise::UnaryFunction file:ROL_Elementwise_Function.hpp line:60
struct PyCallBack_ROL_Elementwise_UnaryFunction_double_t : public ROL::Elementwise::UnaryFunction<double> {
	using ROL::Elementwise::UnaryFunction<double>::UnaryFunction;

	double apply(const double & a0) const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::Elementwise::UnaryFunction<double> *>(this), "apply");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		pybind11::pybind11_fail("Tried to call pure virtual function \"UnaryFunction::apply\"");
	}
};

// ROL::Elementwise::BinaryFunction file:ROL_Elementwise_Function.hpp line:96
struct PyCallBack_ROL_Elementwise_BinaryFunction_double_t : public ROL::Elementwise::BinaryFunction<double> {
	using ROL::Elementwise::BinaryFunction<double>::BinaryFunction;

	double apply(const double & a0, const double & a1) const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::Elementwise::BinaryFunction<double> *>(this), "apply");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		pybind11::pybind11_fail("Tried to call pure virtual function \"BinaryFunction::apply\"");
	}
};

void bind_ROL_Elementwise_Function(std::function< pybind11::module &(std::string const &namespace_) > &M)
{
	{ // ROL::Elementwise::UnaryFunction file:ROL_Elementwise_Function.hpp line:60
		pybind11::class_<ROL::Elementwise::UnaryFunction<double>, std::shared_ptr<ROL::Elementwise::UnaryFunction<double>>, PyCallBack_ROL_Elementwise_UnaryFunction_double_t> cl(M("ROL::Elementwise"), "UnaryFunction_double_t", "", pybind11::module_local());
		cl.def(pybind11::init<PyCallBack_ROL_Elementwise_UnaryFunction_double_t const &>());
		cl.def( pybind11::init( [](){ return new PyCallBack_ROL_Elementwise_UnaryFunction_double_t(); } ) );
		cl.def("apply", (double (ROL::Elementwise::UnaryFunction<double>::*)(const double &) const) &ROL::Elementwise::UnaryFunction<double>::apply, "C++: ROL::Elementwise::UnaryFunction<double>::apply(const double &) const --> double", pybind11::arg("x"));
		cl.def("assign", (class ROL::Elementwise::UnaryFunction<double> & (ROL::Elementwise::UnaryFunction<double>::*)(const class ROL::Elementwise::UnaryFunction<double> &)) &ROL::Elementwise::UnaryFunction<double>::operator=, "C++: ROL::Elementwise::UnaryFunction<double>::operator=(const class ROL::Elementwise::UnaryFunction<double> &) --> class ROL::Elementwise::UnaryFunction<double> &", pybind11::return_value_policy::automatic, pybind11::arg(""));
	}
	{ // ROL::Elementwise::BinaryFunction file:ROL_Elementwise_Function.hpp line:96
		pybind11::class_<ROL::Elementwise::BinaryFunction<double>, std::shared_ptr<ROL::Elementwise::BinaryFunction<double>>, PyCallBack_ROL_Elementwise_BinaryFunction_double_t> cl(M("ROL::Elementwise"), "BinaryFunction_double_t", "", pybind11::module_local());
		cl.def(pybind11::init<PyCallBack_ROL_Elementwise_BinaryFunction_double_t const &>());
		cl.def( pybind11::init( [](){ return new PyCallBack_ROL_Elementwise_BinaryFunction_double_t(); } ) );
		cl.def("apply", (double (ROL::Elementwise::BinaryFunction<double>::*)(const double &, const double &) const) &ROL::Elementwise::BinaryFunction<double>::apply, "C++: ROL::Elementwise::BinaryFunction<double>::apply(const double &, const double &) const --> double", pybind11::arg("x"), pybind11::arg("y"));
		cl.def("assign", (class ROL::Elementwise::BinaryFunction<double> & (ROL::Elementwise::BinaryFunction<double>::*)(const class ROL::Elementwise::BinaryFunction<double> &)) &ROL::Elementwise::BinaryFunction<double>::operator=, "C++: ROL::Elementwise::BinaryFunction<double>::operator=(const class ROL::Elementwise::BinaryFunction<double> &) --> class ROL::Elementwise::BinaryFunction<double> &", pybind11::return_value_policy::automatic, pybind11::arg(""));
	}
}
