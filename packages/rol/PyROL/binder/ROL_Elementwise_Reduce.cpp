#include <ROL_Elementwise_Reduce.hpp>
#include <sstream> // __str__

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

// ROL::Elementwise::ReductionOp file:ROL_Elementwise_Reduce.hpp line:68
struct PyCallBack_ROL_Elementwise_ReductionOp_double_t : public ROL::Elementwise::ReductionOp<double> {
	using ROL::Elementwise::ReductionOp<double>::ReductionOp;

	void reduce(const double & a0, double & a1) const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::Elementwise::ReductionOp<double> *>(this), "reduce");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		pybind11::pybind11_fail("Tried to call pure virtual function \"ReductionOp::reduce\"");
	}
	void reduce(const volatile double & a0, volatile double & a1) const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::Elementwise::ReductionOp<double> *>(this), "reduce");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		pybind11::pybind11_fail("Tried to call pure virtual function \"ReductionOp::reduce\"");
	}
	double initialValue() const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::Elementwise::ReductionOp<double> *>(this), "initialValue");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		pybind11::pybind11_fail("Tried to call pure virtual function \"ReductionOp::initialValue\"");
	}
	enum ROL::Elementwise::EReductionType reductionType() const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::Elementwise::ReductionOp<double> *>(this), "reductionType");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<enum ROL::Elementwise::EReductionType>::value) {
				static pybind11::detail::override_caster_t<enum ROL::Elementwise::EReductionType> caster;
				return pybind11::detail::cast_ref<enum ROL::Elementwise::EReductionType>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<enum ROL::Elementwise::EReductionType>(std::move(o));
		}
		pybind11::pybind11_fail("Tried to call pure virtual function \"ReductionOp::reductionType\"");
	}
};

void bind_ROL_Elementwise_Reduce(std::function< pybind11::module &(std::string const &namespace_) > &M)
{
	// ROL::Elementwise::EReductionType file:ROL_Elementwise_Reduce.hpp line:58
	pybind11::enum_<ROL::Elementwise::EReductionType>(M("ROL::Elementwise"), "EReductionType", pybind11::arithmetic(), "")
		.value("REDUCE_SUM", ROL::Elementwise::REDUCE_SUM)
		.value("REDUCE_MIN", ROL::Elementwise::REDUCE_MIN)
		.value("REDUCE_MAX", ROL::Elementwise::REDUCE_MAX)
		.value("REDUCE_AND", ROL::Elementwise::REDUCE_AND)
		.value("REDUCE_BOR", ROL::Elementwise::REDUCE_BOR)
		.export_values();

;

	{ // ROL::Elementwise::ReductionOp file:ROL_Elementwise_Reduce.hpp line:68
		pybind11::class_<ROL::Elementwise::ReductionOp<double>, Teuchos::RCP<ROL::Elementwise::ReductionOp<double>>, PyCallBack_ROL_Elementwise_ReductionOp_double_t> cl(M("ROL::Elementwise"), "ReductionOp_double_t", "");
		cl.def(pybind11::init<PyCallBack_ROL_Elementwise_ReductionOp_double_t const &>());
		cl.def( pybind11::init( [](){ return new PyCallBack_ROL_Elementwise_ReductionOp_double_t(); } ) );
		cl.def("reduce", (void (ROL::Elementwise::ReductionOp<double>::*)(const double &, double &) const) &ROL::Elementwise::ReductionOp<double>::reduce, "C++: ROL::Elementwise::ReductionOp<double>::reduce(const double &, double &) const --> void", pybind11::arg("input"), pybind11::arg("output"));
		cl.def("reduce", (void (ROL::Elementwise::ReductionOp<double>::*)(const volatile double &, volatile double &) const) &ROL::Elementwise::ReductionOp<double>::reduce, "C++: ROL::Elementwise::ReductionOp<double>::reduce(const volatile double &, volatile double &) const --> void", pybind11::arg("input"), pybind11::arg("output"));
		cl.def("initialValue", (double (ROL::Elementwise::ReductionOp<double>::*)() const) &ROL::Elementwise::ReductionOp<double>::initialValue, "C++: ROL::Elementwise::ReductionOp<double>::initialValue() const --> double");
		cl.def("reductionType", (enum ROL::Elementwise::EReductionType (ROL::Elementwise::ReductionOp<double>::*)() const) &ROL::Elementwise::ReductionOp<double>::reductionType, "C++: ROL::Elementwise::ReductionOp<double>::reductionType() const --> enum ROL::Elementwise::EReductionType");
		cl.def("assign", (class ROL::Elementwise::ReductionOp<double> & (ROL::Elementwise::ReductionOp<double>::*)(const class ROL::Elementwise::ReductionOp<double> &)) &ROL::Elementwise::ReductionOp<double>::operator=, "C++: ROL::Elementwise::ReductionOp<double>::operator=(const class ROL::Elementwise::ReductionOp<double> &) --> class ROL::Elementwise::ReductionOp<double> &", pybind11::return_value_policy::automatic, pybind11::arg(""));
	}
}
