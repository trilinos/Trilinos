#include <ROL_Elementwise_Reduce.hpp>
#include <sstream> // __str__

#include <functional>
#include "PyROL_Smart_Holder.hpp"
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

// ROL::Elementwise::ReductionSum file:ROL_Elementwise_Reduce.hpp line:80
struct PyCallBack_ROL_Elementwise_ReductionSum_double_t : public ROL::Elementwise::ReductionSum<double> {
	using ROL::Elementwise::ReductionSum<double>::ReductionSum;

	void reduce(const double & a0, double & a1) const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::Elementwise::ReductionSum<double> *>(this), "reduce");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return ReductionSum::reduce(a0, a1);
	}
	void reduce(const volatile double & a0, volatile double & a1) const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::Elementwise::ReductionSum<double> *>(this), "reduce");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return ReductionSum::reduce(a0, a1);
	}
	double initialValue() const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::Elementwise::ReductionSum<double> *>(this), "initialValue");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return ReductionSum::initialValue();
	}
	enum ROL::Elementwise::EReductionType reductionType() const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::Elementwise::ReductionSum<double> *>(this), "reductionType");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<enum ROL::Elementwise::EReductionType>::value) {
				static pybind11::detail::override_caster_t<enum ROL::Elementwise::EReductionType> caster;
				return pybind11::detail::cast_ref<enum ROL::Elementwise::EReductionType>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<enum ROL::Elementwise::EReductionType>(std::move(o));
		}
		return ReductionSum::reductionType();
	}
};

// ROL::Elementwise::ReductionMin file:ROL_Elementwise_Reduce.hpp line:121
struct PyCallBack_ROL_Elementwise_ReductionMin_double_t : public ROL::Elementwise::ReductionMin<double> {
	using ROL::Elementwise::ReductionMin<double>::ReductionMin;

	void reduce(const double & a0, double & a1) const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::Elementwise::ReductionMin<double> *>(this), "reduce");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return ReductionMin::reduce(a0, a1);
	}
	void reduce(const volatile double & a0, volatile double & a1) const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::Elementwise::ReductionMin<double> *>(this), "reduce");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return ReductionMin::reduce(a0, a1);
	}
	double initialValue() const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::Elementwise::ReductionMin<double> *>(this), "initialValue");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return ReductionMin::initialValue();
	}
	enum ROL::Elementwise::EReductionType reductionType() const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::Elementwise::ReductionMin<double> *>(this), "reductionType");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<enum ROL::Elementwise::EReductionType>::value) {
				static pybind11::detail::override_caster_t<enum ROL::Elementwise::EReductionType> caster;
				return pybind11::detail::cast_ref<enum ROL::Elementwise::EReductionType>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<enum ROL::Elementwise::EReductionType>(std::move(o));
		}
		return ReductionMin::reductionType();
	}
};

// ROL::Elementwise::ReductionMax file:ROL_Elementwise_Reduce.hpp line:148
struct PyCallBack_ROL_Elementwise_ReductionMax_double_t : public ROL::Elementwise::ReductionMax<double> {
	using ROL::Elementwise::ReductionMax<double>::ReductionMax;

	void reduce(const double & a0, double & a1) const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::Elementwise::ReductionMax<double> *>(this), "reduce");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return ReductionMax::reduce(a0, a1);
	}
	void reduce(const volatile double & a0, volatile double & a1) const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::Elementwise::ReductionMax<double> *>(this), "reduce");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return ReductionMax::reduce(a0, a1);
	}
	double initialValue() const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::Elementwise::ReductionMax<double> *>(this), "initialValue");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return ReductionMax::initialValue();
	}
	enum ROL::Elementwise::EReductionType reductionType() const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::Elementwise::ReductionMax<double> *>(this), "reductionType");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<enum ROL::Elementwise::EReductionType>::value) {
				static pybind11::detail::override_caster_t<enum ROL::Elementwise::EReductionType> caster;
				return pybind11::detail::cast_ref<enum ROL::Elementwise::EReductionType>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<enum ROL::Elementwise::EReductionType>(std::move(o));
		}
		return ReductionMax::reductionType();
	}
};

void bind_ROL_Elementwise_Reduce(std::function< pybind11::module &(std::string const &namespace_) > &M)
{
	// ROL::Elementwise::EReductionType file:ROL_Elementwise_Reduce.hpp line:58
	pybind11::enum_<ROL::Elementwise::EReductionType>(M("ROL::Elementwise"), "EReductionType", pybind11::arithmetic(), "", pybind11::module_local())
		.value("REDUCE_SUM", ROL::Elementwise::REDUCE_SUM)
		.value("REDUCE_MIN", ROL::Elementwise::REDUCE_MIN)
		.value("REDUCE_MAX", ROL::Elementwise::REDUCE_MAX)
		.value("REDUCE_AND", ROL::Elementwise::REDUCE_AND)
		.value("REDUCE_BOR", ROL::Elementwise::REDUCE_BOR)
		.export_values();

;

	{ // ROL::Elementwise::ReductionOp file:ROL_Elementwise_Reduce.hpp line:68
		pybind11::class_<ROL::Elementwise::ReductionOp<double>, Teuchos::RCP<ROL::Elementwise::ReductionOp<double>>, PyCallBack_ROL_Elementwise_ReductionOp_double_t> cl(M("ROL::Elementwise"), "ReductionOp_double_t", "", pybind11::module_local());
		cl.def(pybind11::init<PyCallBack_ROL_Elementwise_ReductionOp_double_t const &>());
		cl.def( pybind11::init( [](){ return new PyCallBack_ROL_Elementwise_ReductionOp_double_t(); } ) );
		cl.def("reduce", (void (ROL::Elementwise::ReductionOp<double>::*)(const double &, double &) const) &ROL::Elementwise::ReductionOp<double>::reduce, "C++: ROL::Elementwise::ReductionOp<double>::reduce(const double &, double &) const --> void", pybind11::arg("input"), pybind11::arg("output"));
		cl.def("reduce", (void (ROL::Elementwise::ReductionOp<double>::*)(const volatile double &, volatile double &) const) &ROL::Elementwise::ReductionOp<double>::reduce, "C++: ROL::Elementwise::ReductionOp<double>::reduce(const volatile double &, volatile double &) const --> void", pybind11::arg("input"), pybind11::arg("output"));
		cl.def("initialValue", (double (ROL::Elementwise::ReductionOp<double>::*)() const) &ROL::Elementwise::ReductionOp<double>::initialValue, "C++: ROL::Elementwise::ReductionOp<double>::initialValue() const --> double");
		cl.def("reductionType", (enum ROL::Elementwise::EReductionType (ROL::Elementwise::ReductionOp<double>::*)() const) &ROL::Elementwise::ReductionOp<double>::reductionType, "C++: ROL::Elementwise::ReductionOp<double>::reductionType() const --> enum ROL::Elementwise::EReductionType");
		cl.def("returnReduce", (double (ROL::Elementwise::ReductionOp<double>::*)(const double &, double &) const) &ROL::Elementwise::ReductionOp<double>::returnReduce, "C++: ROL::Elementwise::ReductionOp<double>::returnReduce(const double &, double &) const --> double", pybind11::arg("input"), pybind11::arg("output"));
		cl.def("returnReduce", (double (ROL::Elementwise::ReductionOp<double>::*)(const volatile double &, volatile double &) const) &ROL::Elementwise::ReductionOp<double>::returnReduce, "C++: ROL::Elementwise::ReductionOp<double>::returnReduce(const volatile double &, volatile double &) const --> double", pybind11::arg("input"), pybind11::arg("output"));
		cl.def("assign", (class ROL::Elementwise::ReductionOp<double> & (ROL::Elementwise::ReductionOp<double>::*)(const class ROL::Elementwise::ReductionOp<double> &)) &ROL::Elementwise::ReductionOp<double>::operator=, "C++: ROL::Elementwise::ReductionOp<double>::operator=(const class ROL::Elementwise::ReductionOp<double> &) --> class ROL::Elementwise::ReductionOp<double> &", pybind11::return_value_policy::automatic, pybind11::arg(""));
	}
	{ // ROL::Elementwise::ReductionSum file:ROL_Elementwise_Reduce.hpp line:80
		pybind11::class_<ROL::Elementwise::ReductionSum<double>, Teuchos::RCP<ROL::Elementwise::ReductionSum<double>>, PyCallBack_ROL_Elementwise_ReductionSum_double_t, ROL::Elementwise::ReductionOp<double>> cl(M("ROL::Elementwise"), "ReductionSum_double_t", "", pybind11::module_local());
		cl.def( pybind11::init( [](PyCallBack_ROL_Elementwise_ReductionSum_double_t const &o){ return new PyCallBack_ROL_Elementwise_ReductionSum_double_t(o); } ) );
		cl.def( pybind11::init( [](ROL::Elementwise::ReductionSum<double> const &o){ return new ROL::Elementwise::ReductionSum<double>(o); } ) );
		cl.def( pybind11::init( [](){ return new ROL::Elementwise::ReductionSum<double>(); }, [](){ return new PyCallBack_ROL_Elementwise_ReductionSum_double_t(); } ) );
		cl.def("reduce", (void (ROL::Elementwise::ReductionSum<double>::*)(const double &, double &) const) &ROL::Elementwise::ReductionSum<double>::reduce, "C++: ROL::Elementwise::ReductionSum<double>::reduce(const double &, double &) const --> void", pybind11::arg("input"), pybind11::arg("output"));
		cl.def("reduce", (void (ROL::Elementwise::ReductionSum<double>::*)(const volatile double &, volatile double &) const) &ROL::Elementwise::ReductionSum<double>::reduce, "C++: ROL::Elementwise::ReductionSum<double>::reduce(const volatile double &, volatile double &) const --> void", pybind11::arg("input"), pybind11::arg("output"));
		cl.def("initialValue", (double (ROL::Elementwise::ReductionSum<double>::*)() const) &ROL::Elementwise::ReductionSum<double>::initialValue, "C++: ROL::Elementwise::ReductionSum<double>::initialValue() const --> double");
		cl.def("reductionType", (enum ROL::Elementwise::EReductionType (ROL::Elementwise::ReductionSum<double>::*)() const) &ROL::Elementwise::ReductionSum<double>::reductionType, "C++: ROL::Elementwise::ReductionSum<double>::reductionType() const --> enum ROL::Elementwise::EReductionType");
		cl.def("assign", (class ROL::Elementwise::ReductionSum<double> & (ROL::Elementwise::ReductionSum<double>::*)(const class ROL::Elementwise::ReductionSum<double> &)) &ROL::Elementwise::ReductionSum<double>::operator=, "C++: ROL::Elementwise::ReductionSum<double>::operator=(const class ROL::Elementwise::ReductionSum<double> &) --> class ROL::Elementwise::ReductionSum<double> &", pybind11::return_value_policy::automatic, pybind11::arg(""));
		cl.def("reduce", (void (ROL::Elementwise::ReductionOp<double>::*)(const double &, double &) const) &ROL::Elementwise::ReductionOp<double>::reduce, "C++: ROL::Elementwise::ReductionOp<double>::reduce(const double &, double &) const --> void", pybind11::arg("input"), pybind11::arg("output"));
		cl.def("reduce", (void (ROL::Elementwise::ReductionOp<double>::*)(const volatile double &, volatile double &) const) &ROL::Elementwise::ReductionOp<double>::reduce, "C++: ROL::Elementwise::ReductionOp<double>::reduce(const volatile double &, volatile double &) const --> void", pybind11::arg("input"), pybind11::arg("output"));
		cl.def("initialValue", (double (ROL::Elementwise::ReductionOp<double>::*)() const) &ROL::Elementwise::ReductionOp<double>::initialValue, "C++: ROL::Elementwise::ReductionOp<double>::initialValue() const --> double");
		cl.def("reductionType", (enum ROL::Elementwise::EReductionType (ROL::Elementwise::ReductionOp<double>::*)() const) &ROL::Elementwise::ReductionOp<double>::reductionType, "C++: ROL::Elementwise::ReductionOp<double>::reductionType() const --> enum ROL::Elementwise::EReductionType");
		cl.def("returnReduce", (double (ROL::Elementwise::ReductionOp<double>::*)(const double &, double &) const) &ROL::Elementwise::ReductionOp<double>::returnReduce, "C++: ROL::Elementwise::ReductionOp<double>::returnReduce(const double &, double &) const --> double", pybind11::arg("input"), pybind11::arg("output"));
		cl.def("returnReduce", (double (ROL::Elementwise::ReductionOp<double>::*)(const volatile double &, volatile double &) const) &ROL::Elementwise::ReductionOp<double>::returnReduce, "C++: ROL::Elementwise::ReductionOp<double>::returnReduce(const volatile double &, volatile double &) const --> double", pybind11::arg("input"), pybind11::arg("output"));
		cl.def("assign", (class ROL::Elementwise::ReductionOp<double> & (ROL::Elementwise::ReductionOp<double>::*)(const class ROL::Elementwise::ReductionOp<double> &)) &ROL::Elementwise::ReductionOp<double>::operator=, "C++: ROL::Elementwise::ReductionOp<double>::operator=(const class ROL::Elementwise::ReductionOp<double> &) --> class ROL::Elementwise::ReductionOp<double> &", pybind11::return_value_policy::automatic, pybind11::arg(""));
	}
	{ // ROL::Elementwise::ReductionMin file:ROL_Elementwise_Reduce.hpp line:121
		pybind11::class_<ROL::Elementwise::ReductionMin<double>, Teuchos::RCP<ROL::Elementwise::ReductionMin<double>>, PyCallBack_ROL_Elementwise_ReductionMin_double_t, ROL::Elementwise::ReductionOp<double>> cl(M("ROL::Elementwise"), "ReductionMin_double_t", "", pybind11::module_local());
		cl.def( pybind11::init( [](){ return new ROL::Elementwise::ReductionMin<double>(); }, [](){ return new PyCallBack_ROL_Elementwise_ReductionMin_double_t(); } ) );
		cl.def( pybind11::init( [](PyCallBack_ROL_Elementwise_ReductionMin_double_t const &o){ return new PyCallBack_ROL_Elementwise_ReductionMin_double_t(o); } ) );
		cl.def( pybind11::init( [](ROL::Elementwise::ReductionMin<double> const &o){ return new ROL::Elementwise::ReductionMin<double>(o); } ) );
		cl.def("reduce", (void (ROL::Elementwise::ReductionMin<double>::*)(const double &, double &) const) &ROL::Elementwise::ReductionMin<double>::reduce, "C++: ROL::Elementwise::ReductionMin<double>::reduce(const double &, double &) const --> void", pybind11::arg("input"), pybind11::arg("output"));
		cl.def("reduce", (void (ROL::Elementwise::ReductionMin<double>::*)(const volatile double &, volatile double &) const) &ROL::Elementwise::ReductionMin<double>::reduce, "C++: ROL::Elementwise::ReductionMin<double>::reduce(const volatile double &, volatile double &) const --> void", pybind11::arg("input"), pybind11::arg("output"));
		cl.def("initialValue", (double (ROL::Elementwise::ReductionMin<double>::*)() const) &ROL::Elementwise::ReductionMin<double>::initialValue, "C++: ROL::Elementwise::ReductionMin<double>::initialValue() const --> double");
		cl.def("reductionType", (enum ROL::Elementwise::EReductionType (ROL::Elementwise::ReductionMin<double>::*)() const) &ROL::Elementwise::ReductionMin<double>::reductionType, "C++: ROL::Elementwise::ReductionMin<double>::reductionType() const --> enum ROL::Elementwise::EReductionType");
		cl.def("assign", (class ROL::Elementwise::ReductionMin<double> & (ROL::Elementwise::ReductionMin<double>::*)(const class ROL::Elementwise::ReductionMin<double> &)) &ROL::Elementwise::ReductionMin<double>::operator=, "C++: ROL::Elementwise::ReductionMin<double>::operator=(const class ROL::Elementwise::ReductionMin<double> &) --> class ROL::Elementwise::ReductionMin<double> &", pybind11::return_value_policy::automatic, pybind11::arg(""));
		cl.def("reduce", (void (ROL::Elementwise::ReductionOp<double>::*)(const double &, double &) const) &ROL::Elementwise::ReductionOp<double>::reduce, "C++: ROL::Elementwise::ReductionOp<double>::reduce(const double &, double &) const --> void", pybind11::arg("input"), pybind11::arg("output"));
		cl.def("reduce", (void (ROL::Elementwise::ReductionOp<double>::*)(const volatile double &, volatile double &) const) &ROL::Elementwise::ReductionOp<double>::reduce, "C++: ROL::Elementwise::ReductionOp<double>::reduce(const volatile double &, volatile double &) const --> void", pybind11::arg("input"), pybind11::arg("output"));
		cl.def("initialValue", (double (ROL::Elementwise::ReductionOp<double>::*)() const) &ROL::Elementwise::ReductionOp<double>::initialValue, "C++: ROL::Elementwise::ReductionOp<double>::initialValue() const --> double");
		cl.def("reductionType", (enum ROL::Elementwise::EReductionType (ROL::Elementwise::ReductionOp<double>::*)() const) &ROL::Elementwise::ReductionOp<double>::reductionType, "C++: ROL::Elementwise::ReductionOp<double>::reductionType() const --> enum ROL::Elementwise::EReductionType");
		cl.def("returnReduce", (double (ROL::Elementwise::ReductionOp<double>::*)(const double &, double &) const) &ROL::Elementwise::ReductionOp<double>::returnReduce, "C++: ROL::Elementwise::ReductionOp<double>::returnReduce(const double &, double &) const --> double", pybind11::arg("input"), pybind11::arg("output"));
		cl.def("returnReduce", (double (ROL::Elementwise::ReductionOp<double>::*)(const volatile double &, volatile double &) const) &ROL::Elementwise::ReductionOp<double>::returnReduce, "C++: ROL::Elementwise::ReductionOp<double>::returnReduce(const volatile double &, volatile double &) const --> double", pybind11::arg("input"), pybind11::arg("output"));
		cl.def("assign", (class ROL::Elementwise::ReductionOp<double> & (ROL::Elementwise::ReductionOp<double>::*)(const class ROL::Elementwise::ReductionOp<double> &)) &ROL::Elementwise::ReductionOp<double>::operator=, "C++: ROL::Elementwise::ReductionOp<double>::operator=(const class ROL::Elementwise::ReductionOp<double> &) --> class ROL::Elementwise::ReductionOp<double> &", pybind11::return_value_policy::automatic, pybind11::arg(""));
	}
	{ // ROL::Elementwise::ReductionMax file:ROL_Elementwise_Reduce.hpp line:148
		pybind11::class_<ROL::Elementwise::ReductionMax<double>, Teuchos::RCP<ROL::Elementwise::ReductionMax<double>>, PyCallBack_ROL_Elementwise_ReductionMax_double_t, ROL::Elementwise::ReductionOp<double>> cl(M("ROL::Elementwise"), "ReductionMax_double_t", "", pybind11::module_local());
		cl.def( pybind11::init( [](){ return new ROL::Elementwise::ReductionMax<double>(); }, [](){ return new PyCallBack_ROL_Elementwise_ReductionMax_double_t(); } ) );
		cl.def( pybind11::init( [](PyCallBack_ROL_Elementwise_ReductionMax_double_t const &o){ return new PyCallBack_ROL_Elementwise_ReductionMax_double_t(o); } ) );
		cl.def( pybind11::init( [](ROL::Elementwise::ReductionMax<double> const &o){ return new ROL::Elementwise::ReductionMax<double>(o); } ) );
		cl.def("reduce", (void (ROL::Elementwise::ReductionMax<double>::*)(const double &, double &) const) &ROL::Elementwise::ReductionMax<double>::reduce, "C++: ROL::Elementwise::ReductionMax<double>::reduce(const double &, double &) const --> void", pybind11::arg("input"), pybind11::arg("output"));
		cl.def("reduce", (void (ROL::Elementwise::ReductionMax<double>::*)(const volatile double &, volatile double &) const) &ROL::Elementwise::ReductionMax<double>::reduce, "C++: ROL::Elementwise::ReductionMax<double>::reduce(const volatile double &, volatile double &) const --> void", pybind11::arg("input"), pybind11::arg("output"));
		cl.def("initialValue", (double (ROL::Elementwise::ReductionMax<double>::*)() const) &ROL::Elementwise::ReductionMax<double>::initialValue, "C++: ROL::Elementwise::ReductionMax<double>::initialValue() const --> double");
		cl.def("reductionType", (enum ROL::Elementwise::EReductionType (ROL::Elementwise::ReductionMax<double>::*)() const) &ROL::Elementwise::ReductionMax<double>::reductionType, "C++: ROL::Elementwise::ReductionMax<double>::reductionType() const --> enum ROL::Elementwise::EReductionType");
		cl.def("assign", (class ROL::Elementwise::ReductionMax<double> & (ROL::Elementwise::ReductionMax<double>::*)(const class ROL::Elementwise::ReductionMax<double> &)) &ROL::Elementwise::ReductionMax<double>::operator=, "C++: ROL::Elementwise::ReductionMax<double>::operator=(const class ROL::Elementwise::ReductionMax<double> &) --> class ROL::Elementwise::ReductionMax<double> &", pybind11::return_value_policy::automatic, pybind11::arg(""));
		cl.def("reduce", (void (ROL::Elementwise::ReductionOp<double>::*)(const double &, double &) const) &ROL::Elementwise::ReductionOp<double>::reduce, "C++: ROL::Elementwise::ReductionOp<double>::reduce(const double &, double &) const --> void", pybind11::arg("input"), pybind11::arg("output"));
		cl.def("reduce", (void (ROL::Elementwise::ReductionOp<double>::*)(const volatile double &, volatile double &) const) &ROL::Elementwise::ReductionOp<double>::reduce, "C++: ROL::Elementwise::ReductionOp<double>::reduce(const volatile double &, volatile double &) const --> void", pybind11::arg("input"), pybind11::arg("output"));
		cl.def("initialValue", (double (ROL::Elementwise::ReductionOp<double>::*)() const) &ROL::Elementwise::ReductionOp<double>::initialValue, "C++: ROL::Elementwise::ReductionOp<double>::initialValue() const --> double");
		cl.def("reductionType", (enum ROL::Elementwise::EReductionType (ROL::Elementwise::ReductionOp<double>::*)() const) &ROL::Elementwise::ReductionOp<double>::reductionType, "C++: ROL::Elementwise::ReductionOp<double>::reductionType() const --> enum ROL::Elementwise::EReductionType");
		cl.def("returnReduce", (double (ROL::Elementwise::ReductionOp<double>::*)(const double &, double &) const) &ROL::Elementwise::ReductionOp<double>::returnReduce, "C++: ROL::Elementwise::ReductionOp<double>::returnReduce(const double &, double &) const --> double", pybind11::arg("input"), pybind11::arg("output"));
		cl.def("returnReduce", (double (ROL::Elementwise::ReductionOp<double>::*)(const volatile double &, volatile double &) const) &ROL::Elementwise::ReductionOp<double>::returnReduce, "C++: ROL::Elementwise::ReductionOp<double>::returnReduce(const volatile double &, volatile double &) const --> double", pybind11::arg("input"), pybind11::arg("output"));
		cl.def("assign", (class ROL::Elementwise::ReductionOp<double> & (ROL::Elementwise::ReductionOp<double>::*)(const class ROL::Elementwise::ReductionOp<double> &)) &ROL::Elementwise::ReductionOp<double>::operator=, "C++: ROL::Elementwise::ReductionOp<double>::operator=(const class ROL::Elementwise::ReductionOp<double> &) --> class ROL::Elementwise::ReductionOp<double> &", pybind11::return_value_policy::automatic, pybind11::arg(""));
	}
}
