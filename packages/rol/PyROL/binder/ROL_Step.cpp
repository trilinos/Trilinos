#include <ROL_BoundConstraint.hpp>
#include <ROL_Constraint.hpp>
#include <ROL_DescentDirection_U.hpp>
#include <ROL_Elementwise_Function.hpp>
#include <ROL_Elementwise_Reduce.hpp>
#include <ROL_Objective.hpp>
#include <ROL_Step.hpp>
#include <ROL_Types.hpp>
#include <ROL_UpdateType.hpp>
#include <ROL_Vector.hpp>
#include <Teuchos_ENull.hpp>
#include <Teuchos_RCPDecl.hpp>
#include <Teuchos_RCPNode.hpp>
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

// ROL::Step file:ROL_Step.hpp line:68
struct PyCallBack_ROL_Step_double_t : public ROL::Step<double> {
	using ROL::Step<double>::Step;

	void initialize(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::BoundConstraint<double> & a3, struct ROL::AlgorithmState<double> & a4) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::Step<double> *>(this), "initialize");
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
	void initialize(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2, class ROL::Objective<double> & a3, class ROL::BoundConstraint<double> & a4, struct ROL::AlgorithmState<double> & a5) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::Step<double> *>(this), "initialize");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4, a5);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Step::initialize(a0, a1, a2, a3, a4, a5);
	}
	void initialize(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Vector<double> & a2, const class ROL::Vector<double> & a3, class ROL::Objective<double> & a4, class ROL::Constraint<double> & a5, struct ROL::AlgorithmState<double> & a6) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::Step<double> *>(this), "initialize");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::Step<double> *>(this), "initialize");
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
	void compute(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::BoundConstraint<double> & a3, struct ROL::AlgorithmState<double> & a4) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::Step<double> *>(this), "compute");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Step::compute(a0, a1, a2, a3, a4);
	}
	void update(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::BoundConstraint<double> & a3, struct ROL::AlgorithmState<double> & a4) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::Step<double> *>(this), "update");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Step::update(a0, a1, a2, a3, a4);
	}
	void compute(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2, class ROL::Objective<double> & a3, class ROL::Constraint<double> & a4, struct ROL::AlgorithmState<double> & a5) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::Step<double> *>(this), "compute");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::Step<double> *>(this), "update");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::Step<double> *>(this), "compute");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::Step<double> *>(this), "update");
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
	std::string printHeader() const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::Step<double> *>(this), "printHeader");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<std::string>::value) {
				static pybind11::detail::override_caster_t<std::string> caster;
				return pybind11::detail::cast_ref<std::string>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<std::string>(std::move(o));
		}
		return Step::printHeader();
	}
	std::string printName() const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::Step<double> *>(this), "printName");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<std::string>::value) {
				static pybind11::detail::override_caster_t<std::string> caster;
				return pybind11::detail::cast_ref<std::string>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<std::string>(std::move(o));
		}
		return Step::printName();
	}
	std::string print(struct ROL::AlgorithmState<double> & a0, bool a1) const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::Step<double> *>(this), "print");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1);
			if (pybind11::detail::cast_is_temporary_value_reference<std::string>::value) {
				static pybind11::detail::override_caster_t<std::string> caster;
				return pybind11::detail::cast_ref<std::string>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<std::string>(std::move(o));
		}
		return Step::print(a0, a1);
	}
};

// ROL::DescentDirection_U file:ROL_DescentDirection_U.hpp line:58
struct PyCallBack_ROL_DescentDirection_U_double_t : public ROL::DescentDirection_U<double> {
	using ROL::DescentDirection_U<double>::DescentDirection_U;

	void initialize(const class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::DescentDirection_U<double> *>(this), "initialize");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return DescentDirection_U::initialize(a0, a1);
	}
	void compute(class ROL::Vector<double> & a0, double & a1, double & a2, int & a3, int & a4, const class ROL::Vector<double> & a5, const class ROL::Vector<double> & a6, class ROL::Objective<double> & a7) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::DescentDirection_U<double> *>(this), "compute");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4, a5, a6, a7);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		pybind11::pybind11_fail("Tried to call pure virtual function \"DescentDirection_U::compute\"");
	}
	void update(const class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2, const class ROL::Vector<double> & a3, const double a4, const int a5) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::DescentDirection_U<double> *>(this), "update");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4, a5);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return DescentDirection_U::update(a0, a1, a2, a3, a4, a5);
	}
	std::string printName() const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::DescentDirection_U<double> *>(this), "printName");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<std::string>::value) {
				static pybind11::detail::override_caster_t<std::string> caster;
				return pybind11::detail::cast_ref<std::string>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<std::string>(std::move(o));
		}
		return DescentDirection_U::printName();
	}
};

void bind_ROL_Step(std::function< pybind11::module &(std::string const &namespace_) > &M)
{
	{ // ROL::Step file:ROL_Step.hpp line:68
		pybind11::class_<ROL::Step<double>, Teuchos::RCP<ROL::Step<double>>, PyCallBack_ROL_Step_double_t> cl(M("ROL"), "Step_double_t", "", pybind11::module_local());
		cl.def( pybind11::init( [](){ return new ROL::Step<double>(); }, [](){ return new PyCallBack_ROL_Step_double_t(); } ) );
		cl.def( pybind11::init( [](PyCallBack_ROL_Step_double_t const &o){ return new PyCallBack_ROL_Step_double_t(o); } ) );
		cl.def( pybind11::init( [](ROL::Step<double> const &o){ return new ROL::Step<double>(o); } ) );
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
		cl.def("getStepState", (const class Teuchos::RCP<const struct ROL::StepState<double> > (ROL::Step<double>::*)() const) &ROL::Step<double>::getStepState, "C++: ROL::Step<double>::getStepState() const --> const class Teuchos::RCP<const struct ROL::StepState<double> >");
		cl.def("reset", [](ROL::Step<double> &o) -> void { return o.reset(); }, "");
		cl.def("reset", (void (ROL::Step<double>::*)(const double)) &ROL::Step<double>::reset, "C++: ROL::Step<double>::reset(const double) --> void", pybind11::arg("searchSize"));
		cl.def("assign", (class ROL::Step<double> & (ROL::Step<double>::*)(const class ROL::Step<double> &)) &ROL::Step<double>::operator=, "C++: ROL::Step<double>::operator=(const class ROL::Step<double> &) --> class ROL::Step<double> &", pybind11::return_value_policy::automatic, pybind11::arg(""));
	}
	{ // ROL::DescentDirection_U file:ROL_DescentDirection_U.hpp line:58
		pybind11::class_<ROL::DescentDirection_U<double>, Teuchos::RCP<ROL::DescentDirection_U<double>>, PyCallBack_ROL_DescentDirection_U_double_t> cl(M("ROL"), "DescentDirection_U_double_t", "", pybind11::module_local());
		cl.def(pybind11::init<PyCallBack_ROL_DescentDirection_U_double_t const &>());
		cl.def( pybind11::init( [](){ return new PyCallBack_ROL_DescentDirection_U_double_t(); } ) );
		cl.def("initialize", (void (ROL::DescentDirection_U<double>::*)(const class ROL::Vector<double> &, const class ROL::Vector<double> &)) &ROL::DescentDirection_U<double>::initialize, "C++: ROL::DescentDirection_U<double>::initialize(const class ROL::Vector<double> &, const class ROL::Vector<double> &) --> void", pybind11::arg("x"), pybind11::arg("g"));
		cl.def("compute", (void (ROL::DescentDirection_U<double>::*)(class ROL::Vector<double> &, double &, double &, int &, int &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &)) &ROL::DescentDirection_U<double>::compute, "C++: ROL::DescentDirection_U<double>::compute(class ROL::Vector<double> &, double &, double &, int &, int &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &) --> void", pybind11::arg("s"), pybind11::arg("snorm"), pybind11::arg("sdotg"), pybind11::arg("iter"), pybind11::arg("flag"), pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"));
		cl.def("update", (void (ROL::DescentDirection_U<double>::*)(const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, const double, const int)) &ROL::DescentDirection_U<double>::update, "C++: ROL::DescentDirection_U<double>::update(const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, const double, const int) --> void", pybind11::arg("x"), pybind11::arg("s"), pybind11::arg("gold"), pybind11::arg("gnew"), pybind11::arg("snorm"), pybind11::arg("iter"));
		cl.def("printName", (std::string (ROL::DescentDirection_U<double>::*)() const) &ROL::DescentDirection_U<double>::printName, "C++: ROL::DescentDirection_U<double>::printName() const --> std::string");
		cl.def("assign", (class ROL::DescentDirection_U<double> & (ROL::DescentDirection_U<double>::*)(const class ROL::DescentDirection_U<double> &)) &ROL::DescentDirection_U<double>::operator=, "C++: ROL::DescentDirection_U<double>::operator=(const class ROL::DescentDirection_U<double> &) --> class ROL::DescentDirection_U<double> &", pybind11::return_value_policy::automatic, pybind11::arg(""));
	}
}
