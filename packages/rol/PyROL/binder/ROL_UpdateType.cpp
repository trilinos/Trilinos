#include <ROL_Elementwise_Function.hpp>
#include <ROL_Elementwise_Reduce.hpp>
#include <ROL_LinearOperator.hpp>
#include <ROL_PartitionedVector.hpp>
#include <ROL_UpdateType.hpp>
#include <ROL_Vector.hpp>
#include <cwchar>
#include <ios>
#include <iterator>
#include <memory>
#include <ostream>
#include <sstream> // __str__
#include <streambuf>
#include <string>
#include <vector>

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

// ROL::LinearOperator file:ROL_LinearOperator.hpp line:71
struct PyCallBack_ROL_LinearOperator_double_t : public ROL::LinearOperator<double> {
	using ROL::LinearOperator<double>::LinearOperator;

	void update(const class ROL::Vector<double> & a0, bool a1, int a2) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::LinearOperator<double> *>(this), "update");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return LinearOperator::update(a0, a1, a2);
	}
	void apply(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, double & a2) const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::LinearOperator<double> *>(this), "apply");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		pybind11::pybind11_fail("Tried to call pure virtual function \"LinearOperator::apply\"");
	}
	void applyInverse(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, double & a2) const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::LinearOperator<double> *>(this), "applyInverse");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return LinearOperator::applyInverse(a0, a1, a2);
	}
	void applyAdjoint(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, double & a2) const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::LinearOperator<double> *>(this), "applyAdjoint");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return LinearOperator::applyAdjoint(a0, a1, a2);
	}
	void applyAdjointInverse(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, double & a2) const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::LinearOperator<double> *>(this), "applyAdjointInverse");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return LinearOperator::applyAdjointInverse(a0, a1, a2);
	}
};

void bind_ROL_UpdateType(std::function< pybind11::module &(std::string const &namespace_) > &M)
{
	// ROL::UpdateType file:ROL_UpdateType.hpp line:52
	pybind11::enum_<ROL::UpdateType>(M("ROL"), "UpdateType", "", pybind11::module_local())
		.value("Initial", ROL::UpdateType::Initial)
		.value("Accept", ROL::UpdateType::Accept)
		.value("Revert", ROL::UpdateType::Revert)
		.value("Trial", ROL::UpdateType::Trial)
		.value("Temp", ROL::UpdateType::Temp);

;

	// ROL::UpdateTypeToString(const enum ROL::UpdateType &) file:ROL_UpdateType.hpp line:61
	M("ROL").def("UpdateTypeToString", (std::string (*)(const enum ROL::UpdateType &)) &ROL::UpdateTypeToString, "C++: ROL::UpdateTypeToString(const enum ROL::UpdateType &) --> std::string", pybind11::arg("type"));

	{ // ROL::LinearOperator file:ROL_LinearOperator.hpp line:71
		pybind11::class_<ROL::LinearOperator<double>, std::shared_ptr<ROL::LinearOperator<double>>, PyCallBack_ROL_LinearOperator_double_t> cl(M("ROL"), "LinearOperator_double_t", "", pybind11::module_local());
		cl.def(pybind11::init<PyCallBack_ROL_LinearOperator_double_t const &>());
		cl.def( pybind11::init( [](){ return new PyCallBack_ROL_LinearOperator_double_t(); } ) );
		cl.def("update", [](ROL::LinearOperator<double> &o, const class ROL::Vector<double> & a0) -> void { return o.update(a0); }, "", pybind11::arg("x"));
		cl.def("update", [](ROL::LinearOperator<double> &o, const class ROL::Vector<double> & a0, bool const & a1) -> void { return o.update(a0, a1); }, "", pybind11::arg("x"), pybind11::arg("flag"));
		cl.def("update", (void (ROL::LinearOperator<double>::*)(const class ROL::Vector<double> &, bool, int)) &ROL::LinearOperator<double>::update, "C++: ROL::LinearOperator<double>::update(const class ROL::Vector<double> &, bool, int) --> void", pybind11::arg("x"), pybind11::arg("flag"), pybind11::arg("iter"));
		cl.def("apply", (void (ROL::LinearOperator<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, double &) const) &ROL::LinearOperator<double>::apply, "C++: ROL::LinearOperator<double>::apply(class ROL::Vector<double> &, const class ROL::Vector<double> &, double &) const --> void", pybind11::arg("Hv"), pybind11::arg("v"), pybind11::arg("tol"));
		cl.def("applyInverse", (void (ROL::LinearOperator<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, double &) const) &ROL::LinearOperator<double>::applyInverse, "C++: ROL::LinearOperator<double>::applyInverse(class ROL::Vector<double> &, const class ROL::Vector<double> &, double &) const --> void", pybind11::arg("Hv"), pybind11::arg("v"), pybind11::arg("tol"));
		cl.def("applyAdjoint", (void (ROL::LinearOperator<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, double &) const) &ROL::LinearOperator<double>::applyAdjoint, "C++: ROL::LinearOperator<double>::applyAdjoint(class ROL::Vector<double> &, const class ROL::Vector<double> &, double &) const --> void", pybind11::arg("Hv"), pybind11::arg("v"), pybind11::arg("tol"));
		cl.def("applyAdjointInverse", (void (ROL::LinearOperator<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, double &) const) &ROL::LinearOperator<double>::applyAdjointInverse, "C++: ROL::LinearOperator<double>::applyAdjointInverse(class ROL::Vector<double> &, const class ROL::Vector<double> &, double &) const --> void", pybind11::arg("Hv"), pybind11::arg("v"), pybind11::arg("tol"));
		cl.def("assign", (class ROL::LinearOperator<double> & (ROL::LinearOperator<double>::*)(const class ROL::LinearOperator<double> &)) &ROL::LinearOperator<double>::operator=, "C++: ROL::LinearOperator<double>::operator=(const class ROL::LinearOperator<double> &) --> class ROL::LinearOperator<double> &", pybind11::return_value_policy::automatic, pybind11::arg(""));
	}
}
