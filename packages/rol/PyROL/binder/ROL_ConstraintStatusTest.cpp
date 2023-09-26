#include <ROL_AugmentedLagrangianObjective.hpp>
#include <ROL_ConstraintStatusTest.hpp>
#include <ROL_Elementwise_Function.hpp>
#include <ROL_Elementwise_Reduce.hpp>
#include <ROL_Objective.hpp>
#include <ROL_Types.hpp>
#include <ROL_UpdateType.hpp>
#include <ROL_Vector.hpp>
#include <Teuchos_ENull.hpp>
#include <Teuchos_FilteredIterator.hpp>
#include <Teuchos_ParameterEntry.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_ParameterListModifier.hpp>
#include <Teuchos_RCPDecl.hpp>
#include <Teuchos_RCPNode.hpp>
#include <Teuchos_StringIndexedOrderedValueObjectContainer.hpp>
#include <cwchar>
#include <deque>
#include <ios>
#include <iterator>
#include <memory>
#include <ostream>
#include <sstream> // __str__
#include <streambuf>
#include <string>
#include <vector>

#include <functional>
#include <pybind11/smart_holder.h>
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

// ROL::ConstraintStatusTest file:ROL_ConstraintStatusTest.hpp line:58
struct PyCallBack_ROL_ConstraintStatusTest_double_t : public ROL::ConstraintStatusTest<double> {
	using ROL::ConstraintStatusTest<double>::ConstraintStatusTest;

	bool check(struct ROL::AlgorithmState<double> & a0) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::ConstraintStatusTest<double> *>(this), "check");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return ConstraintStatusTest::check(a0);
	}
};

// ROL::AugmentedLagrangianObjective file:ROL_AugmentedLagrangianObjective.hpp line:88
struct PyCallBack_ROL_AugmentedLagrangianObjective_double_t : public ROL::AugmentedLagrangianObjective<double> {
	using ROL::AugmentedLagrangianObjective<double>::AugmentedLagrangianObjective;

	void update(const class ROL::Vector<double> & a0, enum ROL::UpdateType a1, int a2) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::AugmentedLagrangianObjective<double> *>(this), "update");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return AugmentedLagrangianObjective::update(a0, a1, a2);
	}
	double value(const class ROL::Vector<double> & a0, double & a1) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::AugmentedLagrangianObjective<double> *>(this), "value");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return AugmentedLagrangianObjective::value(a0, a1);
	}
	void gradient(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, double & a2) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::AugmentedLagrangianObjective<double> *>(this), "gradient");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return AugmentedLagrangianObjective::gradient(a0, a1, a2);
	}
	void hessVec(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2, double & a3) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::AugmentedLagrangianObjective<double> *>(this), "hessVec");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return AugmentedLagrangianObjective::hessVec(a0, a1, a2, a3);
	}
	void update(const class ROL::Vector<double> & a0, bool a1, int a2) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::AugmentedLagrangianObjective<double> *>(this), "update");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Objective::update(a0, a1, a2);
	}
	double dirDeriv(const class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, double & a2) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::AugmentedLagrangianObjective<double> *>(this), "dirDeriv");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return Objective::dirDeriv(a0, a1, a2);
	}
	void invHessVec(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2, double & a3) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::AugmentedLagrangianObjective<double> *>(this), "invHessVec");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Objective::invHessVec(a0, a1, a2, a3);
	}
	void precond(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2, double & a3) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::AugmentedLagrangianObjective<double> *>(this), "precond");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Objective::precond(a0, a1, a2, a3);
	}
};

void bind_ROL_ConstraintStatusTest(std::function< pybind11::module &(std::string const &namespace_) > &M)
{
	{ // ROL::ConstraintStatusTest file:ROL_ConstraintStatusTest.hpp line:58
		pybind11::class_<ROL::ConstraintStatusTest<double>, Teuchos::RCP<ROL::ConstraintStatusTest<double>>, PyCallBack_ROL_ConstraintStatusTest_double_t, ROL::StatusTest<double>> cl(M("ROL"), "ConstraintStatusTest_double_t", "", pybind11::module_local());
		cl.def( pybind11::init<class Teuchos::ParameterList &>(), pybind11::arg("parlist") );

		cl.def( pybind11::init( [](){ return new ROL::ConstraintStatusTest<double>(); }, [](){ return new PyCallBack_ROL_ConstraintStatusTest_double_t(); } ), "doc");
		cl.def( pybind11::init( [](double const & a0){ return new ROL::ConstraintStatusTest<double>(a0); }, [](double const & a0){ return new PyCallBack_ROL_ConstraintStatusTest_double_t(a0); } ), "doc");
		cl.def( pybind11::init( [](double const & a0, double const & a1){ return new ROL::ConstraintStatusTest<double>(a0, a1); }, [](double const & a0, double const & a1){ return new PyCallBack_ROL_ConstraintStatusTest_double_t(a0, a1); } ), "doc");
		cl.def( pybind11::init( [](double const & a0, double const & a1, double const & a2){ return new ROL::ConstraintStatusTest<double>(a0, a1, a2); }, [](double const & a0, double const & a1, double const & a2){ return new PyCallBack_ROL_ConstraintStatusTest_double_t(a0, a1, a2); } ), "doc");
		cl.def( pybind11::init<double, double, double, int>(), pybind11::arg("gtol"), pybind11::arg("ctol"), pybind11::arg("stol"), pybind11::arg("max_iter") );

		cl.def( pybind11::init( [](PyCallBack_ROL_ConstraintStatusTest_double_t const &o){ return new PyCallBack_ROL_ConstraintStatusTest_double_t(o); } ) );
		cl.def( pybind11::init( [](ROL::ConstraintStatusTest<double> const &o){ return new ROL::ConstraintStatusTest<double>(o); } ) );
		cl.def("check", (bool (ROL::ConstraintStatusTest<double>::*)(struct ROL::AlgorithmState<double> &)) &ROL::ConstraintStatusTest<double>::check, "C++: ROL::ConstraintStatusTest<double>::check(struct ROL::AlgorithmState<double> &) --> bool", pybind11::arg("state"));
		cl.def("assign", (class ROL::ConstraintStatusTest<double> & (ROL::ConstraintStatusTest<double>::*)(const class ROL::ConstraintStatusTest<double> &)) &ROL::ConstraintStatusTest<double>::operator=, "C++: ROL::ConstraintStatusTest<double>::operator=(const class ROL::ConstraintStatusTest<double> &) --> class ROL::ConstraintStatusTest<double> &", pybind11::return_value_policy::automatic, pybind11::arg(""));
		cl.def("check", (bool (ROL::StatusTest<double>::*)(struct ROL::AlgorithmState<double> &)) &ROL::StatusTest<double>::check, "C++: ROL::StatusTest<double>::check(struct ROL::AlgorithmState<double> &) --> bool", pybind11::arg("state"));
		cl.def("assign", (class ROL::StatusTest<double> & (ROL::StatusTest<double>::*)(const class ROL::StatusTest<double> &)) &ROL::StatusTest<double>::operator=, "C++: ROL::StatusTest<double>::operator=(const class ROL::StatusTest<double> &) --> class ROL::StatusTest<double> &", pybind11::return_value_policy::automatic, pybind11::arg(""));
	}
	{ // ROL::AugmentedLagrangianObjective file:ROL_AugmentedLagrangianObjective.hpp line:88
		pybind11::class_<ROL::AugmentedLagrangianObjective<double>, Teuchos::RCP<ROL::AugmentedLagrangianObjective<double>>, PyCallBack_ROL_AugmentedLagrangianObjective_double_t, ROL::Objective<double>> cl(M("ROL"), "AugmentedLagrangianObjective_double_t", "", pybind11::module_local());
		cl.def( pybind11::init<const class Teuchos::RCP<class ROL::Objective<double> > &, const class Teuchos::RCP<class ROL::Constraint<double> > &, const double, const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, class Teuchos::ParameterList &>(), pybind11::arg("obj"), pybind11::arg("con"), pybind11::arg("penaltyParameter"), pybind11::arg("dualOptVec"), pybind11::arg("primConVec"), pybind11::arg("dualConVec"), pybind11::arg("parlist") );

		cl.def( pybind11::init<const class Teuchos::RCP<class ROL::Objective<double> > &, const class Teuchos::RCP<class ROL::Constraint<double> > &, const double, const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, const bool, const int>(), pybind11::arg("obj"), pybind11::arg("con"), pybind11::arg("penaltyParameter"), pybind11::arg("dualOptVec"), pybind11::arg("primConVec"), pybind11::arg("dualConVec"), pybind11::arg("scaleLagrangian"), pybind11::arg("HessianApprox") );

		cl.def( pybind11::init( [](PyCallBack_ROL_AugmentedLagrangianObjective_double_t const &o){ return new PyCallBack_ROL_AugmentedLagrangianObjective_double_t(o); } ) );
		cl.def( pybind11::init( [](ROL::AugmentedLagrangianObjective<double> const &o){ return new ROL::AugmentedLagrangianObjective<double>(o); } ) );
		cl.def("update", [](ROL::AugmentedLagrangianObjective<double> &o, const class ROL::Vector<double> & a0, enum ROL::UpdateType const & a1) -> void { return o.update(a0, a1); }, "", pybind11::arg("x"), pybind11::arg("type"));
		cl.def("update", (void (ROL::AugmentedLagrangianObjective<double>::*)(const class ROL::Vector<double> &, enum ROL::UpdateType, int)) &ROL::AugmentedLagrangianObjective<double>::update, "C++: ROL::AugmentedLagrangianObjective<double>::update(const class ROL::Vector<double> &, enum ROL::UpdateType, int) --> void", pybind11::arg("x"), pybind11::arg("type"), pybind11::arg("iter"));
		cl.def("setScaling", [](ROL::AugmentedLagrangianObjective<double> &o) -> void { return o.setScaling(); }, "");
		cl.def("setScaling", [](ROL::AugmentedLagrangianObjective<double> &o, const double & a0) -> void { return o.setScaling(a0); }, "", pybind11::arg("fscale"));
		cl.def("setScaling", (void (ROL::AugmentedLagrangianObjective<double>::*)(const double, const double)) &ROL::AugmentedLagrangianObjective<double>::setScaling, "C++: ROL::AugmentedLagrangianObjective<double>::setScaling(const double, const double) --> void", pybind11::arg("fscale"), pybind11::arg("cscale"));
		cl.def("value", (double (ROL::AugmentedLagrangianObjective<double>::*)(const class ROL::Vector<double> &, double &)) &ROL::AugmentedLagrangianObjective<double>::value, "C++: ROL::AugmentedLagrangianObjective<double>::value(const class ROL::Vector<double> &, double &) --> double", pybind11::arg("x"), pybind11::arg("tol"));
		cl.def("gradient", (void (ROL::AugmentedLagrangianObjective<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, double &)) &ROL::AugmentedLagrangianObjective<double>::gradient, "C++: ROL::AugmentedLagrangianObjective<double>::gradient(class ROL::Vector<double> &, const class ROL::Vector<double> &, double &) --> void", pybind11::arg("g"), pybind11::arg("x"), pybind11::arg("tol"));
		cl.def("hessVec", (void (ROL::AugmentedLagrangianObjective<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, double &)) &ROL::AugmentedLagrangianObjective<double>::hessVec, "C++: ROL::AugmentedLagrangianObjective<double>::hessVec(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, double &) --> void", pybind11::arg("hv"), pybind11::arg("v"), pybind11::arg("x"), pybind11::arg("tol"));
		cl.def("getObjectiveValue", (double (ROL::AugmentedLagrangianObjective<double>::*)(const class ROL::Vector<double> &, double &)) &ROL::AugmentedLagrangianObjective<double>::getObjectiveValue, "C++: ROL::AugmentedLagrangianObjective<double>::getObjectiveValue(const class ROL::Vector<double> &, double &) --> double", pybind11::arg("x"), pybind11::arg("tol"));
		cl.def("getObjectiveGradient", (const class Teuchos::RCP<const class ROL::Vector<double> > (ROL::AugmentedLagrangianObjective<double>::*)(const class ROL::Vector<double> &, double &)) &ROL::AugmentedLagrangianObjective<double>::getObjectiveGradient, "C++: ROL::AugmentedLagrangianObjective<double>::getObjectiveGradient(const class ROL::Vector<double> &, double &) --> const class Teuchos::RCP<const class ROL::Vector<double> >", pybind11::arg("x"), pybind11::arg("tol"));
		cl.def("getConstraintVec", (const class Teuchos::RCP<const class ROL::Vector<double> > (ROL::AugmentedLagrangianObjective<double>::*)(const class ROL::Vector<double> &, double &)) &ROL::AugmentedLagrangianObjective<double>::getConstraintVec, "C++: ROL::AugmentedLagrangianObjective<double>::getConstraintVec(const class ROL::Vector<double> &, double &) --> const class Teuchos::RCP<const class ROL::Vector<double> >", pybind11::arg("x"), pybind11::arg("tol"));
		cl.def("getNumberConstraintEvaluations", (int (ROL::AugmentedLagrangianObjective<double>::*)() const) &ROL::AugmentedLagrangianObjective<double>::getNumberConstraintEvaluations, "C++: ROL::AugmentedLagrangianObjective<double>::getNumberConstraintEvaluations() const --> int");
		cl.def("getNumberFunctionEvaluations", (int (ROL::AugmentedLagrangianObjective<double>::*)() const) &ROL::AugmentedLagrangianObjective<double>::getNumberFunctionEvaluations, "C++: ROL::AugmentedLagrangianObjective<double>::getNumberFunctionEvaluations() const --> int");
		cl.def("getNumberGradientEvaluations", (int (ROL::AugmentedLagrangianObjective<double>::*)() const) &ROL::AugmentedLagrangianObjective<double>::getNumberGradientEvaluations, "C++: ROL::AugmentedLagrangianObjective<double>::getNumberGradientEvaluations() const --> int");
		cl.def("reset", (void (ROL::AugmentedLagrangianObjective<double>::*)(const class ROL::Vector<double> &, const double)) &ROL::AugmentedLagrangianObjective<double>::reset, "C++: ROL::AugmentedLagrangianObjective<double>::reset(const class ROL::Vector<double> &, const double) --> void", pybind11::arg("multiplier"), pybind11::arg("penaltyParameter"));
		cl.def("update", [](ROL::Objective<double> &o, const class ROL::Vector<double> & a0, enum ROL::UpdateType const & a1) -> void { return o.update(a0, a1); }, "", pybind11::arg("x"), pybind11::arg("type"));
		cl.def("update", (void (ROL::Objective<double>::*)(const class ROL::Vector<double> &, enum ROL::UpdateType, int)) &ROL::Objective<double>::update, "Update objective function. \n\n      This function updates the objective function at new iterations. \n      \n\n      is the new iterate. \n      \n\n   is the type of update requested.\n      \n\n   is the outer algorithm iterations count.\n\nC++: ROL::Objective<double>::update(const class ROL::Vector<double> &, enum ROL::UpdateType, int) --> void", pybind11::arg("x"), pybind11::arg("type"), pybind11::arg("iter"));
		cl.def("update", [](ROL::Objective<double> &o, const class ROL::Vector<double> & a0) -> void { return o.update(a0); }, "", pybind11::arg("x"));
		cl.def("update", [](ROL::Objective<double> &o, const class ROL::Vector<double> & a0, bool const & a1) -> void { return o.update(a0, a1); }, "", pybind11::arg("x"), pybind11::arg("flag"));
		cl.def("update", (void (ROL::Objective<double>::*)(const class ROL::Vector<double> &, bool, int)) &ROL::Objective<double>::update, "Update objective function. \n\n      This function updates the objective function at new iterations. \n      \n\n      is the new iterate. \n      \n\n   is true if the iterate has changed.\n      \n\n   is the outer algorithm iterations count.\n\nC++: ROL::Objective<double>::update(const class ROL::Vector<double> &, bool, int) --> void", pybind11::arg("x"), pybind11::arg("flag"), pybind11::arg("iter"));
		cl.def("value", (double (ROL::Objective<double>::*)(const class ROL::Vector<double> &, double &)) &ROL::Objective<double>::value, "C++: ROL::Objective<double>::value(const class ROL::Vector<double> &, double &) --> double", pybind11::arg("x"), pybind11::arg("tol"));
		cl.def("gradient", (void (ROL::Objective<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, double &)) &ROL::Objective<double>::gradient, "Compute gradient.\n\n      This function returns the objective function gradient.\n      \n\n   is the gradient.\n      \n\n   is the current iterate.\n      \n\n is a tolerance for inexact objective function computation.\n\n      The default implementation is a finite-difference approximation based on the function value.\n      This requires the definition of a basis \n\n for the optimization vectors x and\n      the definition of a basis \n\n for the dual optimization vectors (gradient vectors g).\n      The bases must be related through the Riesz map, i.e., \n\n,\n      and this must be reflected in the implementation of the ROL::Vector::dual() method.\n\nC++: ROL::Objective<double>::gradient(class ROL::Vector<double> &, const class ROL::Vector<double> &, double &) --> void", pybind11::arg("g"), pybind11::arg("x"), pybind11::arg("tol"));
		cl.def("dirDeriv", (double (ROL::Objective<double>::*)(const class ROL::Vector<double> &, const class ROL::Vector<double> &, double &)) &ROL::Objective<double>::dirDeriv, "Compute directional derivative.\n\n      This function returns the directional derivative of the objective function in the \n direction.\n      \n\n   is the current iterate.\n      \n\n   is the direction.\n      \n\n is a tolerance for inexact objective function computation.\n\nC++: ROL::Objective<double>::dirDeriv(const class ROL::Vector<double> &, const class ROL::Vector<double> &, double &) --> double", pybind11::arg("x"), pybind11::arg("d"), pybind11::arg("tol"));
		cl.def("hessVec", (void (ROL::Objective<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, double &)) &ROL::Objective<double>::hessVec, "Apply Hessian approximation to vector.\n\n      This function applies the Hessian of the objective function to the vector \n.\n      \n\n  is the the action of the Hessian on \n.\n      \n\n   is the direction vector.\n      \n\n   is the current iterate.\n      \n\n is a tolerance for inexact objective function computation.\n\nC++: ROL::Objective<double>::hessVec(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, double &) --> void", pybind11::arg("hv"), pybind11::arg("v"), pybind11::arg("x"), pybind11::arg("tol"));
		cl.def("invHessVec", (void (ROL::Objective<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, double &)) &ROL::Objective<double>::invHessVec, "Apply inverse Hessian approximation to vector.\n\n      This function applies the inverse Hessian of the objective function to the vector \n.\n      \n\n  is the action of the inverse Hessian on \n.\n      \n\n   is the direction vector.\n      \n\n   is the current iterate.\n      \n\n is a tolerance for inexact objective function computation.\n\nC++: ROL::Objective<double>::invHessVec(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, double &) --> void", pybind11::arg("hv"), pybind11::arg("v"), pybind11::arg("x"), pybind11::arg("tol"));
		cl.def("precond", (void (ROL::Objective<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, double &)) &ROL::Objective<double>::precond, "Apply preconditioner to vector.\n\n      This function applies a preconditioner for the Hessian of the objective function to the vector \n.\n      \n\n  is the action of the Hessian preconditioner on \n.\n      \n\n   is the direction vector.\n      \n\n   is the current iterate.\n      \n\n is a tolerance for inexact objective function computation.\n\nC++: ROL::Objective<double>::precond(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, double &) --> void", pybind11::arg("Pv"), pybind11::arg("v"), pybind11::arg("x"), pybind11::arg("tol"));
		cl.def("assign", (class ROL::Objective<double> & (ROL::Objective<double>::*)(const class ROL::Objective<double> &)) &ROL::Objective<double>::operator=, "C++: ROL::Objective<double>::operator=(const class ROL::Objective<double> &) --> class ROL::Objective<double> &", pybind11::return_value_policy::automatic, pybind11::arg(""));
	}
}
