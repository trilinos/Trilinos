#include <ROL_BoundConstraint.hpp>
#include <ROL_Constraint.hpp>
#include <ROL_Elementwise_Function.hpp>
#include <ROL_Elementwise_Reduce.hpp>
#include <ROL_Objective.hpp>
#include <ROL_PartitionedVector.hpp>
#include <ROL_Secant.hpp>
#include <ROL_Step.hpp>
#include <ROL_TrustRegionStep.hpp>
#include <ROL_Types.hpp>
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

// ROL::TrustRegionStep file: line:46
struct PyCallBack_ROL_TrustRegionStep_double_t : public ROL::TrustRegionStep<double> {
	using ROL::TrustRegionStep<double>::TrustRegionStep;

	void initialize(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2, class ROL::Objective<double> & a3, class ROL::BoundConstraint<double> & a4, struct ROL::AlgorithmState<double> & a5) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TrustRegionStep<double> *>(this), "initialize");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4, a5);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return TrustRegionStep::initialize(a0, a1, a2, a3, a4, a5);
	}
	void compute(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::BoundConstraint<double> & a3, struct ROL::AlgorithmState<double> & a4) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TrustRegionStep<double> *>(this), "compute");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return TrustRegionStep::compute(a0, a1, a2, a3, a4);
	}
	void update(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::BoundConstraint<double> & a3, struct ROL::AlgorithmState<double> & a4) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TrustRegionStep<double> *>(this), "update");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return TrustRegionStep::update(a0, a1, a2, a3, a4);
	}
	std::string printHeader() const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TrustRegionStep<double> *>(this), "printHeader");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<std::string>::value) {
				static pybind11::detail::override_caster_t<std::string> caster;
				return pybind11::detail::cast_ref<std::string>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<std::string>(std::move(o));
		}
		return TrustRegionStep::printHeader();
	}
	std::string printName() const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TrustRegionStep<double> *>(this), "printName");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<std::string>::value) {
				static pybind11::detail::override_caster_t<std::string> caster;
				return pybind11::detail::cast_ref<std::string>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<std::string>(std::move(o));
		}
		return TrustRegionStep::printName();
	}
	std::string print(struct ROL::AlgorithmState<double> & a0, bool a1) const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TrustRegionStep<double> *>(this), "print");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1);
			if (pybind11::detail::cast_is_temporary_value_reference<std::string>::value) {
				static pybind11::detail::override_caster_t<std::string> caster;
				return pybind11::detail::cast_ref<std::string>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<std::string>(std::move(o));
		}
		return TrustRegionStep::print(a0, a1);
	}
	void initialize(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::BoundConstraint<double> & a3, struct ROL::AlgorithmState<double> & a4) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TrustRegionStep<double> *>(this), "initialize");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Step::initialize(a0, a1, a2, a3, a4);
	}
	void initialize(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Vector<double> & a2, const class ROL::Vector<double> & a3, class ROL::Objective<double> & a4, class ROL::Constraint<double> & a5, struct ROL::AlgorithmState<double> & a6) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TrustRegionStep<double> *>(this), "initialize");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4, a5, a6);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Step::initialize(a0, a1, a2, a3, a4, a5, a6);
	}
	void initialize(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Vector<double> & a2, const class ROL::Vector<double> & a3, class ROL::Objective<double> & a4, class ROL::Constraint<double> & a5, class ROL::BoundConstraint<double> & a6, struct ROL::AlgorithmState<double> & a7) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TrustRegionStep<double> *>(this), "initialize");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4, a5, a6, a7);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Step::initialize(a0, a1, a2, a3, a4, a5, a6, a7);
	}
	void compute(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2, class ROL::Objective<double> & a3, class ROL::Constraint<double> & a4, struct ROL::AlgorithmState<double> & a5) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TrustRegionStep<double> *>(this), "compute");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4, a5);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Step::compute(a0, a1, a2, a3, a4, a5);
	}
	void update(class ROL::Vector<double> & a0, class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2, class ROL::Objective<double> & a3, class ROL::Constraint<double> & a4, struct ROL::AlgorithmState<double> & a5) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TrustRegionStep<double> *>(this), "update");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4, a5);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Step::update(a0, a1, a2, a3, a4, a5);
	}
	void compute(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2, class ROL::Objective<double> & a3, class ROL::Constraint<double> & a4, class ROL::BoundConstraint<double> & a5, struct ROL::AlgorithmState<double> & a6) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TrustRegionStep<double> *>(this), "compute");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4, a5, a6);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Step::compute(a0, a1, a2, a3, a4, a5, a6);
	}
	void update(class ROL::Vector<double> & a0, class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2, class ROL::Objective<double> & a3, class ROL::Constraint<double> & a4, class ROL::BoundConstraint<double> & a5, struct ROL::AlgorithmState<double> & a6) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TrustRegionStep<double> *>(this), "update");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4, a5, a6);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Step::update(a0, a1, a2, a3, a4, a5, a6);
	}
};

void bind_unknown_unknown_4(std::function< pybind11::module &(std::string const &namespace_) > &M)
{
	{ // ROL::TrustRegionStep file: line:46
		pybind11::class_<ROL::TrustRegionStep<double>, std::shared_ptr<ROL::TrustRegionStep<double>>, PyCallBack_ROL_TrustRegionStep_double_t, ROL::Step<double>> cl(M("ROL"), "TrustRegionStep_double_t", "", pybind11::module_local());
		cl.def( pybind11::init<class Teuchos::ParameterList &>(), pybind11::arg("parlist") );

		cl.def( pybind11::init<class std::shared_ptr<class ROL::Secant<double> > &, class Teuchos::ParameterList &>(), pybind11::arg("secant"), pybind11::arg("parlist") );

		cl.def( pybind11::init( [](PyCallBack_ROL_TrustRegionStep_double_t const &o){ return new PyCallBack_ROL_TrustRegionStep_double_t(o); } ) );
		cl.def( pybind11::init( [](ROL::TrustRegionStep<double> const &o){ return new ROL::TrustRegionStep<double>(o); } ) );
		cl.def("initialize", [](ROL::TrustRegionStep<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Vector<double> & a2, const class ROL::Vector<double> & a3, class ROL::Objective<double> & a4, class ROL::Constraint<double> & a5, class ROL::BoundConstraint<double> & a6, struct ROL::AlgorithmState<double> & a7) -> void { return o.initialize(a0, a1, a2, a3, a4, a5, a6, a7); }, "", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("l"), pybind11::arg("c"), pybind11::arg("obj"), pybind11::arg("con"), pybind11::arg("bnd"), pybind11::arg("algo_state"));
		cl.def("initialize", [](ROL::TrustRegionStep<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Vector<double> & a2, const class ROL::Vector<double> & a3, class ROL::Objective<double> & a4, class ROL::Constraint<double> & a5, struct ROL::AlgorithmState<double> & a6) -> void { return o.initialize(a0, a1, a2, a3, a4, a5, a6); }, "", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("l"), pybind11::arg("c"), pybind11::arg("obj"), pybind11::arg("con"), pybind11::arg("algo_state"));
		cl.def("initialize", [](ROL::TrustRegionStep<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::BoundConstraint<double> & a3, struct ROL::AlgorithmState<double> & a4) -> void { return o.initialize(a0, a1, a2, a3, a4); }, "", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("con"), pybind11::arg("algo_state"));
		cl.def("compute", [](ROL::TrustRegionStep<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2, class ROL::Objective<double> & a3, class ROL::Constraint<double> & a4, class ROL::BoundConstraint<double> & a5, struct ROL::AlgorithmState<double> & a6) -> void { return o.compute(a0, a1, a2, a3, a4, a5, a6); }, "", pybind11::arg("s"), pybind11::arg("x"), pybind11::arg("l"), pybind11::arg("obj"), pybind11::arg("con"), pybind11::arg("bnd"), pybind11::arg("algo_state"));
		cl.def("compute", [](ROL::TrustRegionStep<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2, class ROL::Objective<double> & a3, class ROL::Constraint<double> & a4, struct ROL::AlgorithmState<double> & a5) -> void { return o.compute(a0, a1, a2, a3, a4, a5); }, "", pybind11::arg("s"), pybind11::arg("x"), pybind11::arg("l"), pybind11::arg("obj"), pybind11::arg("con"), pybind11::arg("algo_state"));
		cl.def("update", [](ROL::TrustRegionStep<double> &o, class ROL::Vector<double> & a0, class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2, class ROL::Objective<double> & a3, class ROL::Constraint<double> & a4, class ROL::BoundConstraint<double> & a5, struct ROL::AlgorithmState<double> & a6) -> void { return o.update(a0, a1, a2, a3, a4, a5, a6); }, "", pybind11::arg("x"), pybind11::arg("l"), pybind11::arg("s"), pybind11::arg("obj"), pybind11::arg("con"), pybind11::arg("bnd"), pybind11::arg("algo_state"));
		cl.def("update", [](ROL::TrustRegionStep<double> &o, class ROL::Vector<double> & a0, class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2, class ROL::Objective<double> & a3, class ROL::Constraint<double> & a4, struct ROL::AlgorithmState<double> & a5) -> void { return o.update(a0, a1, a2, a3, a4, a5); }, "", pybind11::arg("x"), pybind11::arg("l"), pybind11::arg("s"), pybind11::arg("obj"), pybind11::arg("con"), pybind11::arg("algo_state"));
		cl.def("initialize", (void (ROL::TrustRegionStep<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, struct ROL::AlgorithmState<double> &)) &ROL::TrustRegionStep<double>::initialize, "Initialize step.\n\n      This function initializes the information necessary to run the trust-region algorithm.\n      \n\n           is the initial guess for the optimization vector.\n      \n\n         is the objective function.\n      \n\n         is the bound constraint.\n      \n\n  is the algorithm state.\n\nC++: ROL::TrustRegionStep<double>::initialize(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, struct ROL::AlgorithmState<double> &) --> void", pybind11::arg("x"), pybind11::arg("s"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("algo_state"));
		cl.def("compute", (void (ROL::TrustRegionStep<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, struct ROL::AlgorithmState<double> &)) &ROL::TrustRegionStep<double>::compute, "Compute step.\n\n      Computes a trial step, \n by solving the trust-region subproblem.  \n      The trust-region subproblem solver is defined by the enum ETrustRegion.  \n      \n\n          is the computed trial step\n      \n\n          is the current iterate\n      \n\n        is the objective function\n      \n\n        are the bound constraints\n      \n\n contains the current state of the algorithm\n\nC++: ROL::TrustRegionStep<double>::compute(class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, struct ROL::AlgorithmState<double> &) --> void", pybind11::arg("s"), pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("algo_state"));
		cl.def("update", (void (ROL::TrustRegionStep<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, struct ROL::AlgorithmState<double> &)) &ROL::TrustRegionStep<double>::update, "Update step, if successful.\n\n      Given a trial step, \n, this function updates \n. \n      This function also updates the secant approximation.\n\n      \n          is the updated iterate\n      \n\n          is the computed trial step\n      \n\n        is the objective function\n      \n\n        are the bound constraints\n      \n\n contains the current state of the algorithm\n\nC++: ROL::TrustRegionStep<double>::update(class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, struct ROL::AlgorithmState<double> &) --> void", pybind11::arg("x"), pybind11::arg("s"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("algo_state"));
		cl.def("printHeader", (std::string (ROL::TrustRegionStep<double>::*)() const) &ROL::TrustRegionStep<double>::printHeader, "Print iterate header.\n\n      This function produces a string containing header information.\n\nC++: ROL::TrustRegionStep<double>::printHeader() const --> std::string");
		cl.def("printName", (std::string (ROL::TrustRegionStep<double>::*)() const) &ROL::TrustRegionStep<double>::printName, "Print step name.\n\n      This function produces a string containing the algorithmic step information.\n\nC++: ROL::TrustRegionStep<double>::printName() const --> std::string");
		cl.def("print", [](ROL::TrustRegionStep<double> const &o, struct ROL::AlgorithmState<double> & a0) -> std::string { return o.print(a0); }, "", pybind11::arg("algo_state"));
		cl.def("print", (std::string (ROL::TrustRegionStep<double>::*)(struct ROL::AlgorithmState<double> &, bool) const) &ROL::TrustRegionStep<double>::print, "Print iterate status.\n\n      This function prints the iteration status.\n\n      \n    is the current state of the algorithm\n      \n\n   if ste to true will print the header at each iteration\n\nC++: ROL::TrustRegionStep<double>::print(struct ROL::AlgorithmState<double> &, bool) const --> std::string", pybind11::arg("algo_state"), pybind11::arg("print_header"));
		cl.def("assign", (class ROL::TrustRegionStep<double> & (ROL::TrustRegionStep<double>::*)(const class ROL::TrustRegionStep<double> &)) &ROL::TrustRegionStep<double>::operator=, "C++: ROL::TrustRegionStep<double>::operator=(const class ROL::TrustRegionStep<double> &) --> class ROL::TrustRegionStep<double> &", pybind11::return_value_policy::automatic, pybind11::arg(""));
		cl.def("initialize", (void (ROL::Step<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, struct ROL::AlgorithmState<double> &)) &ROL::Step<double>::initialize, "C++: ROL::Step<double>::initialize(class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, struct ROL::AlgorithmState<double> &) --> void", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("con"), pybind11::arg("algo_state"));
		cl.def("initialize", (void (ROL::Step<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, struct ROL::AlgorithmState<double> &)) &ROL::Step<double>::initialize, "C++: ROL::Step<double>::initialize(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, struct ROL::AlgorithmState<double> &) --> void", pybind11::arg("x"), pybind11::arg("s"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("con"), pybind11::arg("algo_state"));
		cl.def("initialize", (void (ROL::Step<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::Constraint<double> &, struct ROL::AlgorithmState<double> &)) &ROL::Step<double>::initialize, "C++: ROL::Step<double>::initialize(class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::Constraint<double> &, struct ROL::AlgorithmState<double> &) --> void", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("l"), pybind11::arg("c"), pybind11::arg("obj"), pybind11::arg("con"), pybind11::arg("algo_state"));
		cl.def("initialize", (void (ROL::Step<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::Constraint<double> &, class ROL::BoundConstraint<double> &, struct ROL::AlgorithmState<double> &)) &ROL::Step<double>::initialize, "C++: ROL::Step<double>::initialize(class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::Constraint<double> &, class ROL::BoundConstraint<double> &, struct ROL::AlgorithmState<double> &) --> void", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("l"), pybind11::arg("c"), pybind11::arg("obj"), pybind11::arg("con"), pybind11::arg("bnd"), pybind11::arg("algo_state"));
		cl.def("compute", (void (ROL::Step<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, struct ROL::AlgorithmState<double> &)) &ROL::Step<double>::compute, "C++: ROL::Step<double>::compute(class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, struct ROL::AlgorithmState<double> &) --> void", pybind11::arg("s"), pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("algo_state"));
		cl.def("update", (void (ROL::Step<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, struct ROL::AlgorithmState<double> &)) &ROL::Step<double>::update, "C++: ROL::Step<double>::update(class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, struct ROL::AlgorithmState<double> &) --> void", pybind11::arg("x"), pybind11::arg("s"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("algo_state"));
		cl.def("compute", (void (ROL::Step<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::Constraint<double> &, struct ROL::AlgorithmState<double> &)) &ROL::Step<double>::compute, "C++: ROL::Step<double>::compute(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::Constraint<double> &, struct ROL::AlgorithmState<double> &) --> void", pybind11::arg("s"), pybind11::arg("x"), pybind11::arg("l"), pybind11::arg("obj"), pybind11::arg("con"), pybind11::arg("algo_state"));
		cl.def("update", (void (ROL::Step<double>::*)(class ROL::Vector<double> &, class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::Constraint<double> &, struct ROL::AlgorithmState<double> &)) &ROL::Step<double>::update, "C++: ROL::Step<double>::update(class ROL::Vector<double> &, class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::Constraint<double> &, struct ROL::AlgorithmState<double> &) --> void", pybind11::arg("x"), pybind11::arg("l"), pybind11::arg("s"), pybind11::arg("obj"), pybind11::arg("con"), pybind11::arg("algo_state"));
		cl.def("compute", (void (ROL::Step<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::Constraint<double> &, class ROL::BoundConstraint<double> &, struct ROL::AlgorithmState<double> &)) &ROL::Step<double>::compute, "C++: ROL::Step<double>::compute(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::Constraint<double> &, class ROL::BoundConstraint<double> &, struct ROL::AlgorithmState<double> &) --> void", pybind11::arg("s"), pybind11::arg("x"), pybind11::arg("l"), pybind11::arg("obj"), pybind11::arg("con"), pybind11::arg("bnd"), pybind11::arg("algo_state"));
		cl.def("update", (void (ROL::Step<double>::*)(class ROL::Vector<double> &, class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::Constraint<double> &, class ROL::BoundConstraint<double> &, struct ROL::AlgorithmState<double> &)) &ROL::Step<double>::update, "C++: ROL::Step<double>::update(class ROL::Vector<double> &, class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::Constraint<double> &, class ROL::BoundConstraint<double> &, struct ROL::AlgorithmState<double> &) --> void", pybind11::arg("x"), pybind11::arg("l"), pybind11::arg("s"), pybind11::arg("obj"), pybind11::arg("con"), pybind11::arg("bnd"), pybind11::arg("algo_state"));
		cl.def("printHeader", (std::string (ROL::Step<double>::*)() const) &ROL::Step<double>::printHeader, "C++: ROL::Step<double>::printHeader() const --> std::string");
		cl.def("printName", (std::string (ROL::Step<double>::*)() const) &ROL::Step<double>::printName, "C++: ROL::Step<double>::printName() const --> std::string");
		cl.def("print", [](ROL::Step<double> const &o, struct ROL::AlgorithmState<double> & a0) -> std::string { return o.print(a0); }, "", pybind11::arg("algo_state"));
		cl.def("print", (std::string (ROL::Step<double>::*)(struct ROL::AlgorithmState<double> &, bool) const) &ROL::Step<double>::print, "C++: ROL::Step<double>::print(struct ROL::AlgorithmState<double> &, bool) const --> std::string", pybind11::arg("algo_state"), pybind11::arg("printHeader"));
		cl.def("getStepState", (const class std::shared_ptr<const struct ROL::StepState<double> > (ROL::Step<double>::*)() const) &ROL::Step<double>::getStepState, "C++: ROL::Step<double>::getStepState() const --> const class std::shared_ptr<const struct ROL::StepState<double> >");
		cl.def("reset", [](ROL::Step<double> &o) -> void { return o.reset(); }, "");
		cl.def("reset", (void (ROL::Step<double>::*)(const double)) &ROL::Step<double>::reset, "C++: ROL::Step<double>::reset(const double) --> void", pybind11::arg("searchSize"));
		cl.def("assign", (class ROL::Step<double> & (ROL::Step<double>::*)(const class ROL::Step<double> &)) &ROL::Step<double>::operator=, "C++: ROL::Step<double>::operator=(const class ROL::Step<double> &) --> class ROL::Step<double> &", pybind11::return_value_policy::automatic, pybind11::arg(""));
	}
}
