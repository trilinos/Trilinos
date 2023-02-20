#include <ROL_Constraint_Partitioned.hpp>
#include <ROL_Elementwise_Function.hpp>
#include <ROL_Elementwise_Reduce.hpp>
#include <ROL_FletcherObjectiveBase.hpp>
#include <ROL_FletcherObjectiveE.hpp>
#include <ROL_Objective.hpp>
#include <ROL_PartitionedVector.hpp>
#include <ROL_SlacklessObjective.hpp>
#include <ROL_UpdateType.hpp>
#include <ROL_Vector.hpp>
#include <Teuchos_FilteredIterator.hpp>
#include <Teuchos_ParameterEntry.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_ParameterListModifier.hpp>
#include <Teuchos_RCPDecl.hpp>
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


#ifndef BINDER_PYBIND11_TYPE_CASTER
	#define BINDER_PYBIND11_TYPE_CASTER
	PYBIND11_DECLARE_HOLDER_TYPE(T, std::shared_ptr<T>)
	PYBIND11_DECLARE_HOLDER_TYPE(T, T*)
	PYBIND11_MAKE_OPAQUE(std::shared_ptr<void>)
#endif

// ROL::FletcherObjectiveBase file:ROL_FletcherObjectiveBase.hpp line:60
struct PyCallBack_ROL_FletcherObjectiveBase_double_t : public ROL::FletcherObjectiveBase<double> {
	using ROL::FletcherObjectiveBase<double>::FletcherObjectiveBase;

	void update(const class ROL::Vector<double> & a0, enum ROL::UpdateType a1, int a2) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::FletcherObjectiveBase<double> *>(this), "update");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return FletcherObjectiveBase::update(a0, a1, a2);
	}
	void solveAugmentedSystem(class ROL::Vector<double> & a0, class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2, const class ROL::Vector<double> & a3, const class ROL::Vector<double> & a4, double & a5, bool a6) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::FletcherObjectiveBase<double> *>(this), "solveAugmentedSystem");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4, a5, a6);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		pybind11::pybind11_fail("Tried to call pure virtual function \"FletcherObjectiveBase::solveAugmentedSystem\"");
	}
	void update(const class ROL::Vector<double> & a0, bool a1, int a2) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::FletcherObjectiveBase<double> *>(this), "update");
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
	double value(const class ROL::Vector<double> & a0, double & a1) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::FletcherObjectiveBase<double> *>(this), "value");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		pybind11::pybind11_fail("Tried to call pure virtual function \"Objective::value\"");
	}
	void gradient(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, double & a2) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::FletcherObjectiveBase<double> *>(this), "gradient");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Objective::gradient(a0, a1, a2);
	}
	double dirDeriv(const class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, double & a2) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::FletcherObjectiveBase<double> *>(this), "dirDeriv");
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
	void hessVec(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2, double & a3) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::FletcherObjectiveBase<double> *>(this), "hessVec");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Objective::hessVec(a0, a1, a2, a3);
	}
	void invHessVec(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2, double & a3) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::FletcherObjectiveBase<double> *>(this), "invHessVec");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::FletcherObjectiveBase<double> *>(this), "precond");
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

// ROL::FletcherObjectiveE file:ROL_FletcherObjectiveE.hpp line:53
struct PyCallBack_ROL_FletcherObjectiveE_double_t : public ROL::FletcherObjectiveE<double> {
	using ROL::FletcherObjectiveE<double>::FletcherObjectiveE;

	double value(const class ROL::Vector<double> & a0, double & a1) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::FletcherObjectiveE<double> *>(this), "value");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return FletcherObjectiveE::value(a0, a1);
	}
	void gradient(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, double & a2) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::FletcherObjectiveE<double> *>(this), "gradient");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return FletcherObjectiveE::gradient(a0, a1, a2);
	}
	void hessVec(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2, double & a3) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::FletcherObjectiveE<double> *>(this), "hessVec");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return FletcherObjectiveE::hessVec(a0, a1, a2, a3);
	}
	void solveAugmentedSystem(class ROL::Vector<double> & a0, class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2, const class ROL::Vector<double> & a3, const class ROL::Vector<double> & a4, double & a5, bool a6) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::FletcherObjectiveE<double> *>(this), "solveAugmentedSystem");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4, a5, a6);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return FletcherObjectiveE::solveAugmentedSystem(a0, a1, a2, a3, a4, a5, a6);
	}
	void update(const class ROL::Vector<double> & a0, enum ROL::UpdateType a1, int a2) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::FletcherObjectiveE<double> *>(this), "update");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return FletcherObjectiveBase::update(a0, a1, a2);
	}
	void update(const class ROL::Vector<double> & a0, bool a1, int a2) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::FletcherObjectiveE<double> *>(this), "update");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::FletcherObjectiveE<double> *>(this), "dirDeriv");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::FletcherObjectiveE<double> *>(this), "invHessVec");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::FletcherObjectiveE<double> *>(this), "precond");
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

void bind_ROL_FletcherObjectiveBase(std::function< pybind11::module &(std::string const &namespace_) > &M)
{
	{ // ROL::FletcherObjectiveBase file:ROL_FletcherObjectiveBase.hpp line:60
		pybind11::class_<ROL::FletcherObjectiveBase<double>, std::shared_ptr<ROL::FletcherObjectiveBase<double>>, PyCallBack_ROL_FletcherObjectiveBase_double_t, ROL::Objective<double>> cl(M("ROL"), "FletcherObjectiveBase_double_t", "", pybind11::module_local());
		cl.def( pybind11::init<const class std::shared_ptr<class ROL::Objective<double> > &, const class std::shared_ptr<class ROL::Constraint<double> > &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, class Teuchos::ParameterList &>(), pybind11::arg("obj"), pybind11::arg("con"), pybind11::arg("xprim"), pybind11::arg("xdual"), pybind11::arg("cprim"), pybind11::arg("cdual"), pybind11::arg("parlist") );

		cl.def(pybind11::init<PyCallBack_ROL_FletcherObjectiveBase_double_t const &>());
		cl.def("update", [](ROL::FletcherObjectiveBase<double> &o, const class ROL::Vector<double> & a0, enum ROL::UpdateType const & a1) -> void { return o.update(a0, a1); }, "", pybind11::arg("x"), pybind11::arg("type"));
		cl.def("update", (void (ROL::FletcherObjectiveBase<double>::*)(const class ROL::Vector<double> &, enum ROL::UpdateType, int)) &ROL::FletcherObjectiveBase<double>::update, "C++: ROL::FletcherObjectiveBase<double>::update(const class ROL::Vector<double> &, enum ROL::UpdateType, int) --> void", pybind11::arg("x"), pybind11::arg("type"), pybind11::arg("iter"));
		cl.def("getLagrangianGradient", (class std::shared_ptr<const class ROL::Vector<double> > (ROL::FletcherObjectiveBase<double>::*)(const class ROL::Vector<double> &)) &ROL::FletcherObjectiveBase<double>::getLagrangianGradient, "C++: ROL::FletcherObjectiveBase<double>::getLagrangianGradient(const class ROL::Vector<double> &) --> class std::shared_ptr<const class ROL::Vector<double> >", pybind11::arg("x"));
		cl.def("getConstraintVec", (class std::shared_ptr<const class ROL::Vector<double> > (ROL::FletcherObjectiveBase<double>::*)(const class ROL::Vector<double> &)) &ROL::FletcherObjectiveBase<double>::getConstraintVec, "C++: ROL::FletcherObjectiveBase<double>::getConstraintVec(const class ROL::Vector<double> &) --> class std::shared_ptr<const class ROL::Vector<double> >", pybind11::arg("x"));
		cl.def("getMultiplierVec", (class std::shared_ptr<const class ROL::Vector<double> > (ROL::FletcherObjectiveBase<double>::*)(const class ROL::Vector<double> &)) &ROL::FletcherObjectiveBase<double>::getMultiplierVec, "C++: ROL::FletcherObjectiveBase<double>::getMultiplierVec(const class ROL::Vector<double> &) --> class std::shared_ptr<const class ROL::Vector<double> >", pybind11::arg("x"));
		cl.def("getGradient", (class std::shared_ptr<const class ROL::Vector<double> > (ROL::FletcherObjectiveBase<double>::*)(const class ROL::Vector<double> &)) &ROL::FletcherObjectiveBase<double>::getGradient, "C++: ROL::FletcherObjectiveBase<double>::getGradient(const class ROL::Vector<double> &) --> class std::shared_ptr<const class ROL::Vector<double> >", pybind11::arg("x"));
		cl.def("getObjectiveValue", (double (ROL::FletcherObjectiveBase<double>::*)(const class ROL::Vector<double> &)) &ROL::FletcherObjectiveBase<double>::getObjectiveValue, "C++: ROL::FletcherObjectiveBase<double>::getObjectiveValue(const class ROL::Vector<double> &) --> double", pybind11::arg("x"));
		cl.def("getNumberFunctionEvaluations", (int (ROL::FletcherObjectiveBase<double>::*)() const) &ROL::FletcherObjectiveBase<double>::getNumberFunctionEvaluations, "C++: ROL::FletcherObjectiveBase<double>::getNumberFunctionEvaluations() const --> int");
		cl.def("getNumberGradientEvaluations", (int (ROL::FletcherObjectiveBase<double>::*)() const) &ROL::FletcherObjectiveBase<double>::getNumberGradientEvaluations, "C++: ROL::FletcherObjectiveBase<double>::getNumberGradientEvaluations() const --> int");
		cl.def("getNumberConstraintEvaluations", (int (ROL::FletcherObjectiveBase<double>::*)() const) &ROL::FletcherObjectiveBase<double>::getNumberConstraintEvaluations, "C++: ROL::FletcherObjectiveBase<double>::getNumberConstraintEvaluations() const --> int");
		cl.def("reset", (void (ROL::FletcherObjectiveBase<double>::*)(double, double)) &ROL::FletcherObjectiveBase<double>::reset, "C++: ROL::FletcherObjectiveBase<double>::reset(double, double) --> void", pybind11::arg("sigma"), pybind11::arg("delta"));
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
	{ // ROL::FletcherObjectiveE file:ROL_FletcherObjectiveE.hpp line:53
		pybind11::class_<ROL::FletcherObjectiveE<double>, std::shared_ptr<ROL::FletcherObjectiveE<double>>, PyCallBack_ROL_FletcherObjectiveE_double_t, ROL::FletcherObjectiveBase<double>> cl(M("ROL"), "FletcherObjectiveE_double_t", "", pybind11::module_local());
		cl.def( pybind11::init<const class std::shared_ptr<class ROL::Objective<double> > &, const class std::shared_ptr<class ROL::Constraint<double> > &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, class Teuchos::ParameterList &>(), pybind11::arg("obj"), pybind11::arg("con"), pybind11::arg("xprim"), pybind11::arg("xdual"), pybind11::arg("cprim"), pybind11::arg("cdual"), pybind11::arg("parlist") );

		cl.def( pybind11::init( [](PyCallBack_ROL_FletcherObjectiveE_double_t const &o){ return new PyCallBack_ROL_FletcherObjectiveE_double_t(o); } ) );
		cl.def( pybind11::init( [](ROL::FletcherObjectiveE<double> const &o){ return new ROL::FletcherObjectiveE<double>(o); } ) );
		cl.def("value", (double (ROL::FletcherObjectiveE<double>::*)(const class ROL::Vector<double> &, double &)) &ROL::FletcherObjectiveE<double>::value, "C++: ROL::FletcherObjectiveE<double>::value(const class ROL::Vector<double> &, double &) --> double", pybind11::arg("x"), pybind11::arg("tol"));
		cl.def("gradient", (void (ROL::FletcherObjectiveE<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, double &)) &ROL::FletcherObjectiveE<double>::gradient, "C++: ROL::FletcherObjectiveE<double>::gradient(class ROL::Vector<double> &, const class ROL::Vector<double> &, double &) --> void", pybind11::arg("g"), pybind11::arg("x"), pybind11::arg("tol"));
		cl.def("hessVec", (void (ROL::FletcherObjectiveE<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, double &)) &ROL::FletcherObjectiveE<double>::hessVec, "C++: ROL::FletcherObjectiveE<double>::hessVec(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, double &) --> void", pybind11::arg("hv"), pybind11::arg("v"), pybind11::arg("x"), pybind11::arg("tol"));
		cl.def("update", [](ROL::FletcherObjectiveBase<double> &o, const class ROL::Vector<double> & a0, enum ROL::UpdateType const & a1) -> void { return o.update(a0, a1); }, "", pybind11::arg("x"), pybind11::arg("type"));
		cl.def("update", (void (ROL::FletcherObjectiveBase<double>::*)(const class ROL::Vector<double> &, enum ROL::UpdateType, int)) &ROL::FletcherObjectiveBase<double>::update, "C++: ROL::FletcherObjectiveBase<double>::update(const class ROL::Vector<double> &, enum ROL::UpdateType, int) --> void", pybind11::arg("x"), pybind11::arg("type"), pybind11::arg("iter"));
		cl.def("getLagrangianGradient", (class std::shared_ptr<const class ROL::Vector<double> > (ROL::FletcherObjectiveBase<double>::*)(const class ROL::Vector<double> &)) &ROL::FletcherObjectiveBase<double>::getLagrangianGradient, "C++: ROL::FletcherObjectiveBase<double>::getLagrangianGradient(const class ROL::Vector<double> &) --> class std::shared_ptr<const class ROL::Vector<double> >", pybind11::arg("x"));
		cl.def("getConstraintVec", (class std::shared_ptr<const class ROL::Vector<double> > (ROL::FletcherObjectiveBase<double>::*)(const class ROL::Vector<double> &)) &ROL::FletcherObjectiveBase<double>::getConstraintVec, "C++: ROL::FletcherObjectiveBase<double>::getConstraintVec(const class ROL::Vector<double> &) --> class std::shared_ptr<const class ROL::Vector<double> >", pybind11::arg("x"));
		cl.def("getMultiplierVec", (class std::shared_ptr<const class ROL::Vector<double> > (ROL::FletcherObjectiveBase<double>::*)(const class ROL::Vector<double> &)) &ROL::FletcherObjectiveBase<double>::getMultiplierVec, "C++: ROL::FletcherObjectiveBase<double>::getMultiplierVec(const class ROL::Vector<double> &) --> class std::shared_ptr<const class ROL::Vector<double> >", pybind11::arg("x"));
		cl.def("getGradient", (class std::shared_ptr<const class ROL::Vector<double> > (ROL::FletcherObjectiveBase<double>::*)(const class ROL::Vector<double> &)) &ROL::FletcherObjectiveBase<double>::getGradient, "C++: ROL::FletcherObjectiveBase<double>::getGradient(const class ROL::Vector<double> &) --> class std::shared_ptr<const class ROL::Vector<double> >", pybind11::arg("x"));
		cl.def("getObjectiveValue", (double (ROL::FletcherObjectiveBase<double>::*)(const class ROL::Vector<double> &)) &ROL::FletcherObjectiveBase<double>::getObjectiveValue, "C++: ROL::FletcherObjectiveBase<double>::getObjectiveValue(const class ROL::Vector<double> &) --> double", pybind11::arg("x"));
		cl.def("getNumberFunctionEvaluations", (int (ROL::FletcherObjectiveBase<double>::*)() const) &ROL::FletcherObjectiveBase<double>::getNumberFunctionEvaluations, "C++: ROL::FletcherObjectiveBase<double>::getNumberFunctionEvaluations() const --> int");
		cl.def("getNumberGradientEvaluations", (int (ROL::FletcherObjectiveBase<double>::*)() const) &ROL::FletcherObjectiveBase<double>::getNumberGradientEvaluations, "C++: ROL::FletcherObjectiveBase<double>::getNumberGradientEvaluations() const --> int");
		cl.def("getNumberConstraintEvaluations", (int (ROL::FletcherObjectiveBase<double>::*)() const) &ROL::FletcherObjectiveBase<double>::getNumberConstraintEvaluations, "C++: ROL::FletcherObjectiveBase<double>::getNumberConstraintEvaluations() const --> int");
		cl.def("reset", (void (ROL::FletcherObjectiveBase<double>::*)(double, double)) &ROL::FletcherObjectiveBase<double>::reset, "C++: ROL::FletcherObjectiveBase<double>::reset(double, double) --> void", pybind11::arg("sigma"), pybind11::arg("delta"));
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
