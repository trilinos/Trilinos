#include <ROL_BinaryFunctions.hpp>
#include <ROL_UnaryFunctions.hpp>
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

// ROL::Elementwise::Fill file:ROL_UnaryFunctions.hpp line:61
struct PyCallBack_ROL_Elementwise_Fill_double_t : public ROL::Elementwise::Fill<double> {
	using ROL::Elementwise::Fill<double>::Fill;

	double apply(const double & a0) const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::Elementwise::Fill<double> *>(this), "apply");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return Fill::apply(a0);
	}
};

// ROL::Elementwise::Shift file:ROL_UnaryFunctions.hpp line:73
struct PyCallBack_ROL_Elementwise_Shift_double_t : public ROL::Elementwise::Shift<double> {
	using ROL::Elementwise::Shift<double>::Shift;

	double apply(const double & a0) const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::Elementwise::Shift<double> *>(this), "apply");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return Shift::apply(a0);
	}
};

// ROL::Elementwise::AbsoluteValue file:ROL_UnaryFunctions.hpp line:95
struct PyCallBack_ROL_Elementwise_AbsoluteValue_double_t : public ROL::Elementwise::AbsoluteValue<double> {
	using ROL::Elementwise::AbsoluteValue<double>::AbsoluteValue;

	double apply(const double & a0) const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::Elementwise::AbsoluteValue<double> *>(this), "apply");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return AbsoluteValue::apply(a0);
	}
};

// ROL::Elementwise::Sign file:ROL_UnaryFunctions.hpp line:104
struct PyCallBack_ROL_Elementwise_Sign_double_t : public ROL::Elementwise::Sign<double> {
	using ROL::Elementwise::Sign<double>::Sign;

	double apply(const double & a0) const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::Elementwise::Sign<double> *>(this), "apply");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return Sign::apply(a0);
	}
};

// ROL::Elementwise::Power file:ROL_UnaryFunctions.hpp line:123
struct PyCallBack_ROL_Elementwise_Power_double_t : public ROL::Elementwise::Power<double> {
	using ROL::Elementwise::Power<double>::Power;

	double apply(const double & a0) const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::Elementwise::Power<double> *>(this), "apply");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return Power::apply(a0);
	}
};

// ROL::Elementwise::SquareRoot file:ROL_UnaryFunctions.hpp line:137
struct PyCallBack_ROL_Elementwise_SquareRoot_double_t : public ROL::Elementwise::SquareRoot<double> {
	using ROL::Elementwise::SquareRoot<double>::SquareRoot;

	double apply(const double & a0) const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::Elementwise::SquareRoot<double> *>(this), "apply");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return SquareRoot::apply(a0);
	}
};

// ROL::Elementwise::UniformlyRandom file:ROL_UnaryFunctions.hpp line:175
struct PyCallBack_ROL_Elementwise_UniformlyRandom_double_t : public ROL::Elementwise::UniformlyRandom<double> {
	using ROL::Elementwise::UniformlyRandom<double>::UniformlyRandom;

	double apply(const double & a0) const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::Elementwise::UniformlyRandom<double> *>(this), "apply");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return UniformlyRandom::apply(a0);
	}
};

// ROL::Elementwise::Multiply file:ROL_BinaryFunctions.hpp line:54
struct PyCallBack_ROL_Elementwise_Multiply_double_t : public ROL::Elementwise::Multiply<double> {
	using ROL::Elementwise::Multiply<double>::Multiply;

	double apply(const double & a0, const double & a1) const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::Elementwise::Multiply<double> *>(this), "apply");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return Multiply::apply(a0, a1);
	}
};

// ROL::Elementwise::Divide file:ROL_BinaryFunctions.hpp line:74
struct PyCallBack_ROL_Elementwise_Divide_double_t : public ROL::Elementwise::Divide<double> {
	using ROL::Elementwise::Divide<double>::Divide;

	double apply(const double & a0, const double & a1) const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::Elementwise::Divide<double> *>(this), "apply");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return Divide::apply(a0, a1);
	}
};

// ROL::Elementwise::DivideAndInvert file:ROL_BinaryFunctions.hpp line:84
struct PyCallBack_ROL_Elementwise_DivideAndInvert_double_t : public ROL::Elementwise::DivideAndInvert<double> {
	using ROL::Elementwise::DivideAndInvert<double>::DivideAndInvert;

	double apply(const double & a0, const double & a1) const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::Elementwise::DivideAndInvert<double> *>(this), "apply");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return DivideAndInvert::apply(a0, a1);
	}
};

// ROL::Elementwise::Min file:ROL_BinaryFunctions.hpp line:120
struct PyCallBack_ROL_Elementwise_Min_double_t : public ROL::Elementwise::Min<double> {
	using ROL::Elementwise::Min<double>::Min;

	double apply(const double & a0, const double & a1) const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::Elementwise::Min<double> *>(this), "apply");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return Min::apply(a0, a1);
	}
};

// ROL::Elementwise::Max file:ROL_BinaryFunctions.hpp line:130
struct PyCallBack_ROL_Elementwise_Max_double_t : public ROL::Elementwise::Max<double> {
	using ROL::Elementwise::Max<double>::Max;

	double apply(const double & a0, const double & a1) const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::Elementwise::Max<double> *>(this), "apply");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return Max::apply(a0, a1);
	}
};

// ROL::Elementwise::ValueSet file:ROL_BinaryFunctions.hpp line:166
struct PyCallBack_ROL_Elementwise_ValueSet_double_t : public ROL::Elementwise::ValueSet<double> {
	using ROL::Elementwise::ValueSet<double>::ValueSet;

	double apply(const double & a0, const double & a1) const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::Elementwise::ValueSet<double> *>(this), "apply");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return ValueSet::apply(a0, a1);
	}
};

void bind_ROL_UnaryFunctions(std::function< pybind11::module &(std::string const &namespace_) > &M)
{
	{ // ROL::Elementwise::Fill file:ROL_UnaryFunctions.hpp line:61
		pybind11::class_<ROL::Elementwise::Fill<double>, Teuchos::RCP<ROL::Elementwise::Fill<double>>, PyCallBack_ROL_Elementwise_Fill_double_t, ROL::Elementwise::UnaryFunction<double>> cl(M("ROL::Elementwise"), "Fill_double_t", "", pybind11::module_local());
		cl.def( pybind11::init<const double &>(), pybind11::arg("value") );

		cl.def( pybind11::init( [](PyCallBack_ROL_Elementwise_Fill_double_t const &o){ return new PyCallBack_ROL_Elementwise_Fill_double_t(o); } ) );
		cl.def( pybind11::init( [](ROL::Elementwise::Fill<double> const &o){ return new ROL::Elementwise::Fill<double>(o); } ) );
		cl.def("apply", (double (ROL::Elementwise::Fill<double>::*)(const double &) const) &ROL::Elementwise::Fill<double>::apply, "C++: ROL::Elementwise::Fill<double>::apply(const double &) const --> double", pybind11::arg("x"));
		cl.def("assign", (class ROL::Elementwise::Fill<double> & (ROL::Elementwise::Fill<double>::*)(const class ROL::Elementwise::Fill<double> &)) &ROL::Elementwise::Fill<double>::operator=, "C++: ROL::Elementwise::Fill<double>::operator=(const class ROL::Elementwise::Fill<double> &) --> class ROL::Elementwise::Fill<double> &", pybind11::return_value_policy::automatic, pybind11::arg(""));
		cl.def("apply", (double (ROL::Elementwise::UnaryFunction<double>::*)(const double &) const) &ROL::Elementwise::UnaryFunction<double>::apply, "C++: ROL::Elementwise::UnaryFunction<double>::apply(const double &) const --> double", pybind11::arg("x"));
		cl.def("assign", (class ROL::Elementwise::UnaryFunction<double> & (ROL::Elementwise::UnaryFunction<double>::*)(const class ROL::Elementwise::UnaryFunction<double> &)) &ROL::Elementwise::UnaryFunction<double>::operator=, "C++: ROL::Elementwise::UnaryFunction<double>::operator=(const class ROL::Elementwise::UnaryFunction<double> &) --> class ROL::Elementwise::UnaryFunction<double> &", pybind11::return_value_policy::automatic, pybind11::arg(""));
	}
	{ // ROL::Elementwise::Shift file:ROL_UnaryFunctions.hpp line:73
		pybind11::class_<ROL::Elementwise::Shift<double>, Teuchos::RCP<ROL::Elementwise::Shift<double>>, PyCallBack_ROL_Elementwise_Shift_double_t, ROL::Elementwise::UnaryFunction<double>> cl(M("ROL::Elementwise"), "Shift_double_t", "", pybind11::module_local());
		cl.def( pybind11::init<const double &>(), pybind11::arg("value") );

		cl.def( pybind11::init( [](PyCallBack_ROL_Elementwise_Shift_double_t const &o){ return new PyCallBack_ROL_Elementwise_Shift_double_t(o); } ) );
		cl.def( pybind11::init( [](ROL::Elementwise::Shift<double> const &o){ return new ROL::Elementwise::Shift<double>(o); } ) );
		cl.def("apply", (double (ROL::Elementwise::Shift<double>::*)(const double &) const) &ROL::Elementwise::Shift<double>::apply, "C++: ROL::Elementwise::Shift<double>::apply(const double &) const --> double", pybind11::arg("x"));
		cl.def("assign", (class ROL::Elementwise::Shift<double> & (ROL::Elementwise::Shift<double>::*)(const class ROL::Elementwise::Shift<double> &)) &ROL::Elementwise::Shift<double>::operator=, "C++: ROL::Elementwise::Shift<double>::operator=(const class ROL::Elementwise::Shift<double> &) --> class ROL::Elementwise::Shift<double> &", pybind11::return_value_policy::automatic, pybind11::arg(""));
		cl.def("apply", (double (ROL::Elementwise::UnaryFunction<double>::*)(const double &) const) &ROL::Elementwise::UnaryFunction<double>::apply, "C++: ROL::Elementwise::UnaryFunction<double>::apply(const double &) const --> double", pybind11::arg("x"));
		cl.def("assign", (class ROL::Elementwise::UnaryFunction<double> & (ROL::Elementwise::UnaryFunction<double>::*)(const class ROL::Elementwise::UnaryFunction<double> &)) &ROL::Elementwise::UnaryFunction<double>::operator=, "C++: ROL::Elementwise::UnaryFunction<double>::operator=(const class ROL::Elementwise::UnaryFunction<double> &) --> class ROL::Elementwise::UnaryFunction<double> &", pybind11::return_value_policy::automatic, pybind11::arg(""));
	}
	{ // ROL::Elementwise::AbsoluteValue file:ROL_UnaryFunctions.hpp line:95
		pybind11::class_<ROL::Elementwise::AbsoluteValue<double>, Teuchos::RCP<ROL::Elementwise::AbsoluteValue<double>>, PyCallBack_ROL_Elementwise_AbsoluteValue_double_t, ROL::Elementwise::UnaryFunction<double>> cl(M("ROL::Elementwise"), "AbsoluteValue_double_t", "", pybind11::module_local());
		cl.def( pybind11::init( [](PyCallBack_ROL_Elementwise_AbsoluteValue_double_t const &o){ return new PyCallBack_ROL_Elementwise_AbsoluteValue_double_t(o); } ) );
		cl.def( pybind11::init( [](ROL::Elementwise::AbsoluteValue<double> const &o){ return new ROL::Elementwise::AbsoluteValue<double>(o); } ) );
		cl.def( pybind11::init( [](){ return new ROL::Elementwise::AbsoluteValue<double>(); }, [](){ return new PyCallBack_ROL_Elementwise_AbsoluteValue_double_t(); } ) );
		cl.def("apply", (double (ROL::Elementwise::AbsoluteValue<double>::*)(const double &) const) &ROL::Elementwise::AbsoluteValue<double>::apply, "C++: ROL::Elementwise::AbsoluteValue<double>::apply(const double &) const --> double", pybind11::arg("x"));
		cl.def("assign", (class ROL::Elementwise::AbsoluteValue<double> & (ROL::Elementwise::AbsoluteValue<double>::*)(const class ROL::Elementwise::AbsoluteValue<double> &)) &ROL::Elementwise::AbsoluteValue<double>::operator=, "C++: ROL::Elementwise::AbsoluteValue<double>::operator=(const class ROL::Elementwise::AbsoluteValue<double> &) --> class ROL::Elementwise::AbsoluteValue<double> &", pybind11::return_value_policy::automatic, pybind11::arg(""));
		cl.def("apply", (double (ROL::Elementwise::UnaryFunction<double>::*)(const double &) const) &ROL::Elementwise::UnaryFunction<double>::apply, "C++: ROL::Elementwise::UnaryFunction<double>::apply(const double &) const --> double", pybind11::arg("x"));
		cl.def("assign", (class ROL::Elementwise::UnaryFunction<double> & (ROL::Elementwise::UnaryFunction<double>::*)(const class ROL::Elementwise::UnaryFunction<double> &)) &ROL::Elementwise::UnaryFunction<double>::operator=, "C++: ROL::Elementwise::UnaryFunction<double>::operator=(const class ROL::Elementwise::UnaryFunction<double> &) --> class ROL::Elementwise::UnaryFunction<double> &", pybind11::return_value_policy::automatic, pybind11::arg(""));
	}
	{ // ROL::Elementwise::Sign file:ROL_UnaryFunctions.hpp line:104
		pybind11::class_<ROL::Elementwise::Sign<double>, Teuchos::RCP<ROL::Elementwise::Sign<double>>, PyCallBack_ROL_Elementwise_Sign_double_t, ROL::Elementwise::UnaryFunction<double>> cl(M("ROL::Elementwise"), "Sign_double_t", "", pybind11::module_local());
		cl.def( pybind11::init( [](){ return new ROL::Elementwise::Sign<double>(); }, [](){ return new PyCallBack_ROL_Elementwise_Sign_double_t(); } ) );
		cl.def( pybind11::init( [](PyCallBack_ROL_Elementwise_Sign_double_t const &o){ return new PyCallBack_ROL_Elementwise_Sign_double_t(o); } ) );
		cl.def( pybind11::init( [](ROL::Elementwise::Sign<double> const &o){ return new ROL::Elementwise::Sign<double>(o); } ) );
		cl.def("apply", (double (ROL::Elementwise::Sign<double>::*)(const double &) const) &ROL::Elementwise::Sign<double>::apply, "C++: ROL::Elementwise::Sign<double>::apply(const double &) const --> double", pybind11::arg("x"));
		cl.def("assign", (class ROL::Elementwise::Sign<double> & (ROL::Elementwise::Sign<double>::*)(const class ROL::Elementwise::Sign<double> &)) &ROL::Elementwise::Sign<double>::operator=, "C++: ROL::Elementwise::Sign<double>::operator=(const class ROL::Elementwise::Sign<double> &) --> class ROL::Elementwise::Sign<double> &", pybind11::return_value_policy::automatic, pybind11::arg(""));
		cl.def("apply", (double (ROL::Elementwise::UnaryFunction<double>::*)(const double &) const) &ROL::Elementwise::UnaryFunction<double>::apply, "C++: ROL::Elementwise::UnaryFunction<double>::apply(const double &) const --> double", pybind11::arg("x"));
		cl.def("assign", (class ROL::Elementwise::UnaryFunction<double> & (ROL::Elementwise::UnaryFunction<double>::*)(const class ROL::Elementwise::UnaryFunction<double> &)) &ROL::Elementwise::UnaryFunction<double>::operator=, "C++: ROL::Elementwise::UnaryFunction<double>::operator=(const class ROL::Elementwise::UnaryFunction<double> &) --> class ROL::Elementwise::UnaryFunction<double> &", pybind11::return_value_policy::automatic, pybind11::arg(""));
	}
	{ // ROL::Elementwise::Power file:ROL_UnaryFunctions.hpp line:123
		pybind11::class_<ROL::Elementwise::Power<double>, Teuchos::RCP<ROL::Elementwise::Power<double>>, PyCallBack_ROL_Elementwise_Power_double_t, ROL::Elementwise::UnaryFunction<double>> cl(M("ROL::Elementwise"), "Power_double_t", "", pybind11::module_local());
		cl.def( pybind11::init<const double &>(), pybind11::arg("exponent") );

		cl.def( pybind11::init( [](PyCallBack_ROL_Elementwise_Power_double_t const &o){ return new PyCallBack_ROL_Elementwise_Power_double_t(o); } ) );
		cl.def( pybind11::init( [](ROL::Elementwise::Power<double> const &o){ return new ROL::Elementwise::Power<double>(o); } ) );
		cl.def("apply", (double (ROL::Elementwise::Power<double>::*)(const double &) const) &ROL::Elementwise::Power<double>::apply, "C++: ROL::Elementwise::Power<double>::apply(const double &) const --> double", pybind11::arg("x"));
		cl.def("assign", (class ROL::Elementwise::Power<double> & (ROL::Elementwise::Power<double>::*)(const class ROL::Elementwise::Power<double> &)) &ROL::Elementwise::Power<double>::operator=, "C++: ROL::Elementwise::Power<double>::operator=(const class ROL::Elementwise::Power<double> &) --> class ROL::Elementwise::Power<double> &", pybind11::return_value_policy::automatic, pybind11::arg(""));
		cl.def("apply", (double (ROL::Elementwise::UnaryFunction<double>::*)(const double &) const) &ROL::Elementwise::UnaryFunction<double>::apply, "C++: ROL::Elementwise::UnaryFunction<double>::apply(const double &) const --> double", pybind11::arg("x"));
		cl.def("assign", (class ROL::Elementwise::UnaryFunction<double> & (ROL::Elementwise::UnaryFunction<double>::*)(const class ROL::Elementwise::UnaryFunction<double> &)) &ROL::Elementwise::UnaryFunction<double>::operator=, "C++: ROL::Elementwise::UnaryFunction<double>::operator=(const class ROL::Elementwise::UnaryFunction<double> &) --> class ROL::Elementwise::UnaryFunction<double> &", pybind11::return_value_policy::automatic, pybind11::arg(""));
	}
	{ // ROL::Elementwise::SquareRoot file:ROL_UnaryFunctions.hpp line:137
		pybind11::class_<ROL::Elementwise::SquareRoot<double>, Teuchos::RCP<ROL::Elementwise::SquareRoot<double>>, PyCallBack_ROL_Elementwise_SquareRoot_double_t, ROL::Elementwise::UnaryFunction<double>> cl(M("ROL::Elementwise"), "SquareRoot_double_t", "", pybind11::module_local());
		cl.def( pybind11::init( [](){ return new ROL::Elementwise::SquareRoot<double>(); }, [](){ return new PyCallBack_ROL_Elementwise_SquareRoot_double_t(); } ) );
		cl.def( pybind11::init( [](PyCallBack_ROL_Elementwise_SquareRoot_double_t const &o){ return new PyCallBack_ROL_Elementwise_SquareRoot_double_t(o); } ) );
		cl.def( pybind11::init( [](ROL::Elementwise::SquareRoot<double> const &o){ return new ROL::Elementwise::SquareRoot<double>(o); } ) );
		cl.def("apply", (double (ROL::Elementwise::SquareRoot<double>::*)(const double &) const) &ROL::Elementwise::SquareRoot<double>::apply, "C++: ROL::Elementwise::SquareRoot<double>::apply(const double &) const --> double", pybind11::arg("x"));
		cl.def("assign", (class ROL::Elementwise::SquareRoot<double> & (ROL::Elementwise::SquareRoot<double>::*)(const class ROL::Elementwise::SquareRoot<double> &)) &ROL::Elementwise::SquareRoot<double>::operator=, "C++: ROL::Elementwise::SquareRoot<double>::operator=(const class ROL::Elementwise::SquareRoot<double> &) --> class ROL::Elementwise::SquareRoot<double> &", pybind11::return_value_policy::automatic, pybind11::arg(""));
		cl.def("apply", (double (ROL::Elementwise::UnaryFunction<double>::*)(const double &) const) &ROL::Elementwise::UnaryFunction<double>::apply, "C++: ROL::Elementwise::UnaryFunction<double>::apply(const double &) const --> double", pybind11::arg("x"));
		cl.def("assign", (class ROL::Elementwise::UnaryFunction<double> & (ROL::Elementwise::UnaryFunction<double>::*)(const class ROL::Elementwise::UnaryFunction<double> &)) &ROL::Elementwise::UnaryFunction<double>::operator=, "C++: ROL::Elementwise::UnaryFunction<double>::operator=(const class ROL::Elementwise::UnaryFunction<double> &) --> class ROL::Elementwise::UnaryFunction<double> &", pybind11::return_value_policy::automatic, pybind11::arg(""));
	}
	{ // ROL::Elementwise::UniformlyRandom file:ROL_UnaryFunctions.hpp line:175
		pybind11::class_<ROL::Elementwise::UniformlyRandom<double>, Teuchos::RCP<ROL::Elementwise::UniformlyRandom<double>>, PyCallBack_ROL_Elementwise_UniformlyRandom_double_t, ROL::Elementwise::UnaryFunction<double>> cl(M("ROL::Elementwise"), "UniformlyRandom_double_t", "", pybind11::module_local());
		cl.def( pybind11::init( [](){ return new ROL::Elementwise::UniformlyRandom<double>(); }, [](){ return new PyCallBack_ROL_Elementwise_UniformlyRandom_double_t(); } ), "doc");
		cl.def( pybind11::init( [](const double & a0){ return new ROL::Elementwise::UniformlyRandom<double>(a0); }, [](const double & a0){ return new PyCallBack_ROL_Elementwise_UniformlyRandom_double_t(a0); } ), "doc");
		cl.def( pybind11::init<const double &, const double &>(), pybind11::arg("lower"), pybind11::arg("upper") );

		cl.def( pybind11::init( [](PyCallBack_ROL_Elementwise_UniformlyRandom_double_t const &o){ return new PyCallBack_ROL_Elementwise_UniformlyRandom_double_t(o); } ) );
		cl.def( pybind11::init( [](ROL::Elementwise::UniformlyRandom<double> const &o){ return new ROL::Elementwise::UniformlyRandom<double>(o); } ) );
		cl.def("apply", (double (ROL::Elementwise::UniformlyRandom<double>::*)(const double &) const) &ROL::Elementwise::UniformlyRandom<double>::apply, "C++: ROL::Elementwise::UniformlyRandom<double>::apply(const double &) const --> double", pybind11::arg("x"));
		cl.def("apply", (double (ROL::Elementwise::UnaryFunction<double>::*)(const double &) const) &ROL::Elementwise::UnaryFunction<double>::apply, "C++: ROL::Elementwise::UnaryFunction<double>::apply(const double &) const --> double", pybind11::arg("x"));
		cl.def("assign", (class ROL::Elementwise::UnaryFunction<double> & (ROL::Elementwise::UnaryFunction<double>::*)(const class ROL::Elementwise::UnaryFunction<double> &)) &ROL::Elementwise::UnaryFunction<double>::operator=, "C++: ROL::Elementwise::UnaryFunction<double>::operator=(const class ROL::Elementwise::UnaryFunction<double> &) --> class ROL::Elementwise::UnaryFunction<double> &", pybind11::return_value_policy::automatic, pybind11::arg(""));
	}
	{ // ROL::Elementwise::Multiply file:ROL_BinaryFunctions.hpp line:54
		pybind11::class_<ROL::Elementwise::Multiply<double>, Teuchos::RCP<ROL::Elementwise::Multiply<double>>, PyCallBack_ROL_Elementwise_Multiply_double_t, ROL::Elementwise::BinaryFunction<double>> cl(M("ROL::Elementwise"), "Multiply_double_t", "", pybind11::module_local());
		cl.def( pybind11::init( [](){ return new ROL::Elementwise::Multiply<double>(); }, [](){ return new PyCallBack_ROL_Elementwise_Multiply_double_t(); } ) );
		cl.def( pybind11::init( [](PyCallBack_ROL_Elementwise_Multiply_double_t const &o){ return new PyCallBack_ROL_Elementwise_Multiply_double_t(o); } ) );
		cl.def( pybind11::init( [](ROL::Elementwise::Multiply<double> const &o){ return new ROL::Elementwise::Multiply<double>(o); } ) );
		cl.def("apply", (double (ROL::Elementwise::Multiply<double>::*)(const double &, const double &) const) &ROL::Elementwise::Multiply<double>::apply, "C++: ROL::Elementwise::Multiply<double>::apply(const double &, const double &) const --> double", pybind11::arg("x"), pybind11::arg("y"));
		cl.def("assign", (class ROL::Elementwise::Multiply<double> & (ROL::Elementwise::Multiply<double>::*)(const class ROL::Elementwise::Multiply<double> &)) &ROL::Elementwise::Multiply<double>::operator=, "C++: ROL::Elementwise::Multiply<double>::operator=(const class ROL::Elementwise::Multiply<double> &) --> class ROL::Elementwise::Multiply<double> &", pybind11::return_value_policy::automatic, pybind11::arg(""));
		cl.def("apply", (double (ROL::Elementwise::BinaryFunction<double>::*)(const double &, const double &) const) &ROL::Elementwise::BinaryFunction<double>::apply, "C++: ROL::Elementwise::BinaryFunction<double>::apply(const double &, const double &) const --> double", pybind11::arg("x"), pybind11::arg("y"));
		cl.def("assign", (class ROL::Elementwise::BinaryFunction<double> & (ROL::Elementwise::BinaryFunction<double>::*)(const class ROL::Elementwise::BinaryFunction<double> &)) &ROL::Elementwise::BinaryFunction<double>::operator=, "C++: ROL::Elementwise::BinaryFunction<double>::operator=(const class ROL::Elementwise::BinaryFunction<double> &) --> class ROL::Elementwise::BinaryFunction<double> &", pybind11::return_value_policy::automatic, pybind11::arg(""));
	}
	{ // ROL::Elementwise::Divide file:ROL_BinaryFunctions.hpp line:74
		pybind11::class_<ROL::Elementwise::Divide<double>, Teuchos::RCP<ROL::Elementwise::Divide<double>>, PyCallBack_ROL_Elementwise_Divide_double_t, ROL::Elementwise::BinaryFunction<double>> cl(M("ROL::Elementwise"), "Divide_double_t", "", pybind11::module_local());
		cl.def( pybind11::init( [](){ return new ROL::Elementwise::Divide<double>(); }, [](){ return new PyCallBack_ROL_Elementwise_Divide_double_t(); } ) );
		cl.def( pybind11::init( [](PyCallBack_ROL_Elementwise_Divide_double_t const &o){ return new PyCallBack_ROL_Elementwise_Divide_double_t(o); } ) );
		cl.def( pybind11::init( [](ROL::Elementwise::Divide<double> const &o){ return new ROL::Elementwise::Divide<double>(o); } ) );
		cl.def("apply", (double (ROL::Elementwise::Divide<double>::*)(const double &, const double &) const) &ROL::Elementwise::Divide<double>::apply, "C++: ROL::Elementwise::Divide<double>::apply(const double &, const double &) const --> double", pybind11::arg("x"), pybind11::arg("y"));
		cl.def("assign", (class ROL::Elementwise::Divide<double> & (ROL::Elementwise::Divide<double>::*)(const class ROL::Elementwise::Divide<double> &)) &ROL::Elementwise::Divide<double>::operator=, "C++: ROL::Elementwise::Divide<double>::operator=(const class ROL::Elementwise::Divide<double> &) --> class ROL::Elementwise::Divide<double> &", pybind11::return_value_policy::automatic, pybind11::arg(""));
		cl.def("apply", (double (ROL::Elementwise::BinaryFunction<double>::*)(const double &, const double &) const) &ROL::Elementwise::BinaryFunction<double>::apply, "C++: ROL::Elementwise::BinaryFunction<double>::apply(const double &, const double &) const --> double", pybind11::arg("x"), pybind11::arg("y"));
		cl.def("assign", (class ROL::Elementwise::BinaryFunction<double> & (ROL::Elementwise::BinaryFunction<double>::*)(const class ROL::Elementwise::BinaryFunction<double> &)) &ROL::Elementwise::BinaryFunction<double>::operator=, "C++: ROL::Elementwise::BinaryFunction<double>::operator=(const class ROL::Elementwise::BinaryFunction<double> &) --> class ROL::Elementwise::BinaryFunction<double> &", pybind11::return_value_policy::automatic, pybind11::arg(""));
	}
	{ // ROL::Elementwise::DivideAndInvert file:ROL_BinaryFunctions.hpp line:84
		pybind11::class_<ROL::Elementwise::DivideAndInvert<double>, Teuchos::RCP<ROL::Elementwise::DivideAndInvert<double>>, PyCallBack_ROL_Elementwise_DivideAndInvert_double_t, ROL::Elementwise::BinaryFunction<double>> cl(M("ROL::Elementwise"), "DivideAndInvert_double_t", "", pybind11::module_local());
		cl.def( pybind11::init( [](){ return new ROL::Elementwise::DivideAndInvert<double>(); }, [](){ return new PyCallBack_ROL_Elementwise_DivideAndInvert_double_t(); } ) );
		cl.def( pybind11::init( [](PyCallBack_ROL_Elementwise_DivideAndInvert_double_t const &o){ return new PyCallBack_ROL_Elementwise_DivideAndInvert_double_t(o); } ) );
		cl.def( pybind11::init( [](ROL::Elementwise::DivideAndInvert<double> const &o){ return new ROL::Elementwise::DivideAndInvert<double>(o); } ) );
		cl.def("apply", (double (ROL::Elementwise::DivideAndInvert<double>::*)(const double &, const double &) const) &ROL::Elementwise::DivideAndInvert<double>::apply, "C++: ROL::Elementwise::DivideAndInvert<double>::apply(const double &, const double &) const --> double", pybind11::arg("x"), pybind11::arg("y"));
		cl.def("assign", (class ROL::Elementwise::DivideAndInvert<double> & (ROL::Elementwise::DivideAndInvert<double>::*)(const class ROL::Elementwise::DivideAndInvert<double> &)) &ROL::Elementwise::DivideAndInvert<double>::operator=, "C++: ROL::Elementwise::DivideAndInvert<double>::operator=(const class ROL::Elementwise::DivideAndInvert<double> &) --> class ROL::Elementwise::DivideAndInvert<double> &", pybind11::return_value_policy::automatic, pybind11::arg(""));
		cl.def("apply", (double (ROL::Elementwise::BinaryFunction<double>::*)(const double &, const double &) const) &ROL::Elementwise::BinaryFunction<double>::apply, "C++: ROL::Elementwise::BinaryFunction<double>::apply(const double &, const double &) const --> double", pybind11::arg("x"), pybind11::arg("y"));
		cl.def("assign", (class ROL::Elementwise::BinaryFunction<double> & (ROL::Elementwise::BinaryFunction<double>::*)(const class ROL::Elementwise::BinaryFunction<double> &)) &ROL::Elementwise::BinaryFunction<double>::operator=, "C++: ROL::Elementwise::BinaryFunction<double>::operator=(const class ROL::Elementwise::BinaryFunction<double> &) --> class ROL::Elementwise::BinaryFunction<double> &", pybind11::return_value_policy::automatic, pybind11::arg(""));
	}
	{ // ROL::Elementwise::Min file:ROL_BinaryFunctions.hpp line:120
		pybind11::class_<ROL::Elementwise::Min<double>, Teuchos::RCP<ROL::Elementwise::Min<double>>, PyCallBack_ROL_Elementwise_Min_double_t, ROL::Elementwise::BinaryFunction<double>> cl(M("ROL::Elementwise"), "Min_double_t", "", pybind11::module_local());
		cl.def( pybind11::init( [](){ return new ROL::Elementwise::Min<double>(); }, [](){ return new PyCallBack_ROL_Elementwise_Min_double_t(); } ) );
		cl.def( pybind11::init( [](PyCallBack_ROL_Elementwise_Min_double_t const &o){ return new PyCallBack_ROL_Elementwise_Min_double_t(o); } ) );
		cl.def( pybind11::init( [](ROL::Elementwise::Min<double> const &o){ return new ROL::Elementwise::Min<double>(o); } ) );
		cl.def("apply", (double (ROL::Elementwise::Min<double>::*)(const double &, const double &) const) &ROL::Elementwise::Min<double>::apply, "C++: ROL::Elementwise::Min<double>::apply(const double &, const double &) const --> double", pybind11::arg("x"), pybind11::arg("y"));
		cl.def("assign", (class ROL::Elementwise::Min<double> & (ROL::Elementwise::Min<double>::*)(const class ROL::Elementwise::Min<double> &)) &ROL::Elementwise::Min<double>::operator=, "C++: ROL::Elementwise::Min<double>::operator=(const class ROL::Elementwise::Min<double> &) --> class ROL::Elementwise::Min<double> &", pybind11::return_value_policy::automatic, pybind11::arg(""));
		cl.def("apply", (double (ROL::Elementwise::BinaryFunction<double>::*)(const double &, const double &) const) &ROL::Elementwise::BinaryFunction<double>::apply, "C++: ROL::Elementwise::BinaryFunction<double>::apply(const double &, const double &) const --> double", pybind11::arg("x"), pybind11::arg("y"));
		cl.def("assign", (class ROL::Elementwise::BinaryFunction<double> & (ROL::Elementwise::BinaryFunction<double>::*)(const class ROL::Elementwise::BinaryFunction<double> &)) &ROL::Elementwise::BinaryFunction<double>::operator=, "C++: ROL::Elementwise::BinaryFunction<double>::operator=(const class ROL::Elementwise::BinaryFunction<double> &) --> class ROL::Elementwise::BinaryFunction<double> &", pybind11::return_value_policy::automatic, pybind11::arg(""));
	}
	{ // ROL::Elementwise::Max file:ROL_BinaryFunctions.hpp line:130
		pybind11::class_<ROL::Elementwise::Max<double>, Teuchos::RCP<ROL::Elementwise::Max<double>>, PyCallBack_ROL_Elementwise_Max_double_t, ROL::Elementwise::BinaryFunction<double>> cl(M("ROL::Elementwise"), "Max_double_t", "", pybind11::module_local());
		cl.def( pybind11::init( [](){ return new ROL::Elementwise::Max<double>(); }, [](){ return new PyCallBack_ROL_Elementwise_Max_double_t(); } ) );
		cl.def( pybind11::init( [](PyCallBack_ROL_Elementwise_Max_double_t const &o){ return new PyCallBack_ROL_Elementwise_Max_double_t(o); } ) );
		cl.def( pybind11::init( [](ROL::Elementwise::Max<double> const &o){ return new ROL::Elementwise::Max<double>(o); } ) );
		cl.def("apply", (double (ROL::Elementwise::Max<double>::*)(const double &, const double &) const) &ROL::Elementwise::Max<double>::apply, "C++: ROL::Elementwise::Max<double>::apply(const double &, const double &) const --> double", pybind11::arg("x"), pybind11::arg("y"));
		cl.def("assign", (class ROL::Elementwise::Max<double> & (ROL::Elementwise::Max<double>::*)(const class ROL::Elementwise::Max<double> &)) &ROL::Elementwise::Max<double>::operator=, "C++: ROL::Elementwise::Max<double>::operator=(const class ROL::Elementwise::Max<double> &) --> class ROL::Elementwise::Max<double> &", pybind11::return_value_policy::automatic, pybind11::arg(""));
		cl.def("apply", (double (ROL::Elementwise::BinaryFunction<double>::*)(const double &, const double &) const) &ROL::Elementwise::BinaryFunction<double>::apply, "C++: ROL::Elementwise::BinaryFunction<double>::apply(const double &, const double &) const --> double", pybind11::arg("x"), pybind11::arg("y"));
		cl.def("assign", (class ROL::Elementwise::BinaryFunction<double> & (ROL::Elementwise::BinaryFunction<double>::*)(const class ROL::Elementwise::BinaryFunction<double> &)) &ROL::Elementwise::BinaryFunction<double>::operator=, "C++: ROL::Elementwise::BinaryFunction<double>::operator=(const class ROL::Elementwise::BinaryFunction<double> &) --> class ROL::Elementwise::BinaryFunction<double> &", pybind11::return_value_policy::automatic, pybind11::arg(""));
	}
	{ // ROL::Elementwise::ValueSet file:ROL_BinaryFunctions.hpp line:166
		pybind11::class_<ROL::Elementwise::ValueSet<double>, Teuchos::RCP<ROL::Elementwise::ValueSet<double>>, PyCallBack_ROL_Elementwise_ValueSet_double_t, ROL::Elementwise::BinaryFunction<double>> cl(M("ROL::Elementwise"), "ValueSet_double_t", "", pybind11::module_local());
		cl.def( pybind11::init( [](const double & a0, const int & a1){ return new ROL::Elementwise::ValueSet<double>(a0, a1); }, [](const double & a0, const int & a1){ return new PyCallBack_ROL_Elementwise_ValueSet_double_t(a0, a1); } ), "doc");
		cl.def( pybind11::init( [](const double & a0, const int & a1, const double & a2){ return new ROL::Elementwise::ValueSet<double>(a0, a1, a2); }, [](const double & a0, const int & a1, const double & a2){ return new PyCallBack_ROL_Elementwise_ValueSet_double_t(a0, a1, a2); } ), "doc");
		cl.def( pybind11::init<const double &, const int, const double &, const double &>(), pybind11::arg("threshold"), pybind11::arg("option"), pybind11::arg("c1"), pybind11::arg("c2") );

		cl.def( pybind11::init( [](PyCallBack_ROL_Elementwise_ValueSet_double_t const &o){ return new PyCallBack_ROL_Elementwise_ValueSet_double_t(o); } ) );
		cl.def( pybind11::init( [](ROL::Elementwise::ValueSet<double> const &o){ return new ROL::Elementwise::ValueSet<double>(o); } ) );
		cl.def("apply", (double (ROL::Elementwise::ValueSet<double>::*)(const double &, const double &) const) &ROL::Elementwise::ValueSet<double>::apply, "C++: ROL::Elementwise::ValueSet<double>::apply(const double &, const double &) const --> double", pybind11::arg("x"), pybind11::arg("y"));
		cl.def("apply", (double (ROL::Elementwise::BinaryFunction<double>::*)(const double &, const double &) const) &ROL::Elementwise::BinaryFunction<double>::apply, "C++: ROL::Elementwise::BinaryFunction<double>::apply(const double &, const double &) const --> double", pybind11::arg("x"), pybind11::arg("y"));
		cl.def("assign", (class ROL::Elementwise::BinaryFunction<double> & (ROL::Elementwise::BinaryFunction<double>::*)(const class ROL::Elementwise::BinaryFunction<double> &)) &ROL::Elementwise::BinaryFunction<double>::operator=, "C++: ROL::Elementwise::BinaryFunction<double>::operator=(const class ROL::Elementwise::BinaryFunction<double> &) --> class ROL::Elementwise::BinaryFunction<double> &", pybind11::return_value_policy::automatic, pybind11::arg(""));
	}
}
