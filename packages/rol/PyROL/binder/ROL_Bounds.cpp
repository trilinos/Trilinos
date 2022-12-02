#include <ROL_BoundConstraint.hpp>
#include <ROL_Bounds.hpp>
#include <ROL_Elementwise_Function.hpp>
#include <ROL_Elementwise_Reduce.hpp>
#include <ROL_Vector.hpp>
#include <Teuchos_ENull.hpp>
#include <Teuchos_RCPDecl.hpp>
#include <Teuchos_RCPNode.hpp>
#include <Teuchos_any.hpp>
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

// ROL::Bounds file:ROL_Bounds.hpp line:59
struct PyCallBack_ROL_Bounds_double_t : public ROL::Bounds<double> {
	using ROL::Bounds<double>::Bounds;

	void project(class ROL::Vector<double> & a0) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::Bounds<double> *>(this), "project");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Bounds::project(a0);
	}
	void projectInterior(class ROL::Vector<double> & a0) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::Bounds<double> *>(this), "projectInterior");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Bounds::projectInterior(a0);
	}
	void pruneUpperActive(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, double a2) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::Bounds<double> *>(this), "pruneUpperActive");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Bounds::pruneUpperActive(a0, a1, a2);
	}
	void pruneUpperActive(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2, double a3, double a4) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::Bounds<double> *>(this), "pruneUpperActive");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Bounds::pruneUpperActive(a0, a1, a2, a3, a4);
	}
	void pruneLowerActive(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, double a2) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::Bounds<double> *>(this), "pruneLowerActive");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Bounds::pruneLowerActive(a0, a1, a2);
	}
	void pruneLowerActive(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2, double a3, double a4) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::Bounds<double> *>(this), "pruneLowerActive");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Bounds::pruneLowerActive(a0, a1, a2, a3, a4);
	}
	bool isFeasible(const class ROL::Vector<double> & a0) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::Bounds<double> *>(this), "isFeasible");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return Bounds::isFeasible(a0);
	}
	void applyInverseScalingFunction(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2, const class ROL::Vector<double> & a3) const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::Bounds<double> *>(this), "applyInverseScalingFunction");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Bounds::applyInverseScalingFunction(a0, a1, a2, a3);
	}
	void applyScalingFunctionJacobian(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2, const class ROL::Vector<double> & a3) const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::Bounds<double> *>(this), "applyScalingFunctionJacobian");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Bounds::applyScalingFunctionJacobian(a0, a1, a2, a3);
	}
	const class Teuchos::RCP<const class ROL::Vector<double> > getLowerBound() const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::Bounds<double> *>(this), "getLowerBound");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<const class Teuchos::RCP<const class ROL::Vector<double> >>::value) {
				static pybind11::detail::override_caster_t<const class Teuchos::RCP<const class ROL::Vector<double> >> caster;
				return pybind11::detail::cast_ref<const class Teuchos::RCP<const class ROL::Vector<double> >>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<const class Teuchos::RCP<const class ROL::Vector<double> >>(std::move(o));
		}
		return BoundConstraint::getLowerBound();
	}
	const class Teuchos::RCP<const class ROL::Vector<double> > getUpperBound() const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::Bounds<double> *>(this), "getUpperBound");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<const class Teuchos::RCP<const class ROL::Vector<double> >>::value) {
				static pybind11::detail::override_caster_t<const class Teuchos::RCP<const class ROL::Vector<double> >> caster;
				return pybind11::detail::cast_ref<const class Teuchos::RCP<const class ROL::Vector<double> >>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<const class Teuchos::RCP<const class ROL::Vector<double> >>(std::move(o));
		}
		return BoundConstraint::getUpperBound();
	}
};

void bind_ROL_Bounds(std::function< pybind11::module &(std::string const &namespace_) > &M)
{
	{ // ROL::Bounds file:ROL_Bounds.hpp line:59
		pybind11::class_<ROL::Bounds<double>, Teuchos::RCP<ROL::Bounds<double>>, PyCallBack_ROL_Bounds_double_t, ROL::BoundConstraint<double>> cl(M("ROL"), "Bounds_double_t", "");
		cl.def( pybind11::init( [](const class ROL::Vector<double> & a0){ return new ROL::Bounds<double>(a0); }, [](const class ROL::Vector<double> & a0){ return new PyCallBack_ROL_Bounds_double_t(a0); } ), "doc");
		cl.def( pybind11::init( [](const class ROL::Vector<double> & a0, bool const & a1){ return new ROL::Bounds<double>(a0, a1); }, [](const class ROL::Vector<double> & a0, bool const & a1){ return new PyCallBack_ROL_Bounds_double_t(a0, a1); } ), "doc");
		cl.def( pybind11::init( [](const class ROL::Vector<double> & a0, bool const & a1, double const & a2){ return new ROL::Bounds<double>(a0, a1, a2); }, [](const class ROL::Vector<double> & a0, bool const & a1, double const & a2){ return new PyCallBack_ROL_Bounds_double_t(a0, a1, a2); } ), "doc");
		cl.def( pybind11::init<const class ROL::Vector<double> &, bool, double, double>(), pybind11::arg("x"), pybind11::arg("isLower"), pybind11::arg("scale"), pybind11::arg("feasTol") );

		cl.def( pybind11::init( [](const class Teuchos::RCP<class ROL::Vector<double> > & a0, const class Teuchos::RCP<class ROL::Vector<double> > & a1){ return new ROL::Bounds<double>(a0, a1); }, [](const class Teuchos::RCP<class ROL::Vector<double> > & a0, const class Teuchos::RCP<class ROL::Vector<double> > & a1){ return new PyCallBack_ROL_Bounds_double_t(a0, a1); } ), "doc");
		cl.def( pybind11::init( [](const class Teuchos::RCP<class ROL::Vector<double> > & a0, const class Teuchos::RCP<class ROL::Vector<double> > & a1, const double & a2){ return new ROL::Bounds<double>(a0, a1, a2); }, [](const class Teuchos::RCP<class ROL::Vector<double> > & a0, const class Teuchos::RCP<class ROL::Vector<double> > & a1, const double & a2){ return new PyCallBack_ROL_Bounds_double_t(a0, a1, a2); } ), "doc");
		cl.def( pybind11::init<const class Teuchos::RCP<class ROL::Vector<double> > &, const class Teuchos::RCP<class ROL::Vector<double> > &, const double, const double>(), pybind11::arg("x_lo"), pybind11::arg("x_up"), pybind11::arg("scale"), pybind11::arg("feasTol") );

		cl.def( pybind11::init( [](PyCallBack_ROL_Bounds_double_t const &o){ return new PyCallBack_ROL_Bounds_double_t(o); } ) );
		cl.def( pybind11::init( [](ROL::Bounds<double> const &o){ return new ROL::Bounds<double>(o); } ) );
		cl.def("project", (void (ROL::Bounds<double>::*)(class ROL::Vector<double> &)) &ROL::Bounds<double>::project, "C++: ROL::Bounds<double>::project(class ROL::Vector<double> &) --> void", pybind11::arg("x"));
		cl.def("projectInterior", (void (ROL::Bounds<double>::*)(class ROL::Vector<double> &)) &ROL::Bounds<double>::projectInterior, "C++: ROL::Bounds<double>::projectInterior(class ROL::Vector<double> &) --> void", pybind11::arg("x"));
		cl.def("pruneUpperActive", [](ROL::Bounds<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1) -> void { return o.pruneUpperActive(a0, a1); }, "", pybind11::arg("v"), pybind11::arg("x"));
		cl.def("pruneUpperActive", (void (ROL::Bounds<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, double)) &ROL::Bounds<double>::pruneUpperActive, "C++: ROL::Bounds<double>::pruneUpperActive(class ROL::Vector<double> &, const class ROL::Vector<double> &, double) --> void", pybind11::arg("v"), pybind11::arg("x"), pybind11::arg("eps"));
		cl.def("pruneUpperActive", [](ROL::Bounds<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2) -> void { return o.pruneUpperActive(a0, a1, a2); }, "", pybind11::arg("v"), pybind11::arg("g"), pybind11::arg("x"));
		cl.def("pruneUpperActive", [](ROL::Bounds<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2, double const & a3) -> void { return o.pruneUpperActive(a0, a1, a2, a3); }, "", pybind11::arg("v"), pybind11::arg("g"), pybind11::arg("x"), pybind11::arg("xeps"));
		cl.def("pruneUpperActive", (void (ROL::Bounds<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, double, double)) &ROL::Bounds<double>::pruneUpperActive, "C++: ROL::Bounds<double>::pruneUpperActive(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, double, double) --> void", pybind11::arg("v"), pybind11::arg("g"), pybind11::arg("x"), pybind11::arg("xeps"), pybind11::arg("geps"));
		cl.def("pruneLowerActive", [](ROL::Bounds<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1) -> void { return o.pruneLowerActive(a0, a1); }, "", pybind11::arg("v"), pybind11::arg("x"));
		cl.def("pruneLowerActive", (void (ROL::Bounds<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, double)) &ROL::Bounds<double>::pruneLowerActive, "C++: ROL::Bounds<double>::pruneLowerActive(class ROL::Vector<double> &, const class ROL::Vector<double> &, double) --> void", pybind11::arg("v"), pybind11::arg("x"), pybind11::arg("eps"));
		cl.def("pruneLowerActive", [](ROL::Bounds<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2) -> void { return o.pruneLowerActive(a0, a1, a2); }, "", pybind11::arg("v"), pybind11::arg("g"), pybind11::arg("x"));
		cl.def("pruneLowerActive", [](ROL::Bounds<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2, double const & a3) -> void { return o.pruneLowerActive(a0, a1, a2, a3); }, "", pybind11::arg("v"), pybind11::arg("g"), pybind11::arg("x"), pybind11::arg("xeps"));
		cl.def("pruneLowerActive", (void (ROL::Bounds<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, double, double)) &ROL::Bounds<double>::pruneLowerActive, "C++: ROL::Bounds<double>::pruneLowerActive(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, double, double) --> void", pybind11::arg("v"), pybind11::arg("g"), pybind11::arg("x"), pybind11::arg("xeps"), pybind11::arg("geps"));
		cl.def("isFeasible", (bool (ROL::Bounds<double>::*)(const class ROL::Vector<double> &)) &ROL::Bounds<double>::isFeasible, "C++: ROL::Bounds<double>::isFeasible(const class ROL::Vector<double> &) --> bool", pybind11::arg("v"));
		cl.def("applyInverseScalingFunction", (void (ROL::Bounds<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &) const) &ROL::Bounds<double>::applyInverseScalingFunction, "C++: ROL::Bounds<double>::applyInverseScalingFunction(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &) const --> void", pybind11::arg("dv"), pybind11::arg("v"), pybind11::arg("x"), pybind11::arg("g"));
		cl.def("applyScalingFunctionJacobian", (void (ROL::Bounds<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &) const) &ROL::Bounds<double>::applyScalingFunctionJacobian, "C++: ROL::Bounds<double>::applyScalingFunctionJacobian(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &) const --> void", pybind11::arg("dv"), pybind11::arg("v"), pybind11::arg("x"), pybind11::arg("g"));
		cl.def("project", (void (ROL::BoundConstraint<double>::*)(class ROL::Vector<double> &)) &ROL::BoundConstraint<double>::project, "C++: ROL::BoundConstraint<double>::project(class ROL::Vector<double> &) --> void", pybind11::arg("x"));
		cl.def("projectInterior", (void (ROL::BoundConstraint<double>::*)(class ROL::Vector<double> &)) &ROL::BoundConstraint<double>::projectInterior, "C++: ROL::BoundConstraint<double>::projectInterior(class ROL::Vector<double> &) --> void", pybind11::arg("x"));
		cl.def("pruneUpperActive", [](ROL::BoundConstraint<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1) -> void { return o.pruneUpperActive(a0, a1); }, "", pybind11::arg("v"), pybind11::arg("x"));
		cl.def("pruneUpperActive", (void (ROL::BoundConstraint<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, double)) &ROL::BoundConstraint<double>::pruneUpperActive, "C++: ROL::BoundConstraint<double>::pruneUpperActive(class ROL::Vector<double> &, const class ROL::Vector<double> &, double) --> void", pybind11::arg("v"), pybind11::arg("x"), pybind11::arg("eps"));
		cl.def("pruneUpperActive", [](ROL::BoundConstraint<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2) -> void { return o.pruneUpperActive(a0, a1, a2); }, "", pybind11::arg("v"), pybind11::arg("g"), pybind11::arg("x"));
		cl.def("pruneUpperActive", [](ROL::BoundConstraint<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2, double const & a3) -> void { return o.pruneUpperActive(a0, a1, a2, a3); }, "", pybind11::arg("v"), pybind11::arg("g"), pybind11::arg("x"), pybind11::arg("xeps"));
		cl.def("pruneUpperActive", (void (ROL::BoundConstraint<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, double, double)) &ROL::BoundConstraint<double>::pruneUpperActive, "C++: ROL::BoundConstraint<double>::pruneUpperActive(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, double, double) --> void", pybind11::arg("v"), pybind11::arg("g"), pybind11::arg("x"), pybind11::arg("xeps"), pybind11::arg("geps"));
		cl.def("pruneLowerActive", [](ROL::BoundConstraint<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1) -> void { return o.pruneLowerActive(a0, a1); }, "", pybind11::arg("v"), pybind11::arg("x"));
		cl.def("pruneLowerActive", (void (ROL::BoundConstraint<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, double)) &ROL::BoundConstraint<double>::pruneLowerActive, "C++: ROL::BoundConstraint<double>::pruneLowerActive(class ROL::Vector<double> &, const class ROL::Vector<double> &, double) --> void", pybind11::arg("v"), pybind11::arg("x"), pybind11::arg("eps"));
		cl.def("pruneLowerActive", [](ROL::BoundConstraint<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2) -> void { return o.pruneLowerActive(a0, a1, a2); }, "", pybind11::arg("v"), pybind11::arg("g"), pybind11::arg("x"));
		cl.def("pruneLowerActive", [](ROL::BoundConstraint<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2, double const & a3) -> void { return o.pruneLowerActive(a0, a1, a2, a3); }, "", pybind11::arg("v"), pybind11::arg("g"), pybind11::arg("x"), pybind11::arg("xeps"));
		cl.def("pruneLowerActive", (void (ROL::BoundConstraint<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, double, double)) &ROL::BoundConstraint<double>::pruneLowerActive, "C++: ROL::BoundConstraint<double>::pruneLowerActive(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, double, double) --> void", pybind11::arg("v"), pybind11::arg("g"), pybind11::arg("x"), pybind11::arg("xeps"), pybind11::arg("geps"));
		cl.def("getLowerBound", (const class Teuchos::RCP<const class ROL::Vector<double> > (ROL::BoundConstraint<double>::*)() const) &ROL::BoundConstraint<double>::getLowerBound, "C++: ROL::BoundConstraint<double>::getLowerBound() const --> const class Teuchos::RCP<const class ROL::Vector<double> >");
		cl.def("getUpperBound", (const class Teuchos::RCP<const class ROL::Vector<double> > (ROL::BoundConstraint<double>::*)() const) &ROL::BoundConstraint<double>::getUpperBound, "C++: ROL::BoundConstraint<double>::getUpperBound() const --> const class Teuchos::RCP<const class ROL::Vector<double> >");
		cl.def("isFeasible", (bool (ROL::BoundConstraint<double>::*)(const class ROL::Vector<double> &)) &ROL::BoundConstraint<double>::isFeasible, "C++: ROL::BoundConstraint<double>::isFeasible(const class ROL::Vector<double> &) --> bool", pybind11::arg("v"));
		cl.def("applyInverseScalingFunction", (void (ROL::BoundConstraint<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &) const) &ROL::BoundConstraint<double>::applyInverseScalingFunction, "C++: ROL::BoundConstraint<double>::applyInverseScalingFunction(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &) const --> void", pybind11::arg("dv"), pybind11::arg("v"), pybind11::arg("x"), pybind11::arg("g"));
		cl.def("applyScalingFunctionJacobian", (void (ROL::BoundConstraint<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &) const) &ROL::BoundConstraint<double>::applyScalingFunctionJacobian, "C++: ROL::BoundConstraint<double>::applyScalingFunctionJacobian(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &) const --> void", pybind11::arg("dv"), pybind11::arg("v"), pybind11::arg("x"), pybind11::arg("g"));
		cl.def("activateLower", (void (ROL::BoundConstraint<double>::*)()) &ROL::BoundConstraint<double>::activateLower, "C++: ROL::BoundConstraint<double>::activateLower() --> void");
		cl.def("activateUpper", (void (ROL::BoundConstraint<double>::*)()) &ROL::BoundConstraint<double>::activateUpper, "C++: ROL::BoundConstraint<double>::activateUpper() --> void");
		cl.def("activate", (void (ROL::BoundConstraint<double>::*)()) &ROL::BoundConstraint<double>::activate, "C++: ROL::BoundConstraint<double>::activate() --> void");
		cl.def("deactivateLower", (void (ROL::BoundConstraint<double>::*)()) &ROL::BoundConstraint<double>::deactivateLower, "C++: ROL::BoundConstraint<double>::deactivateLower() --> void");
		cl.def("deactivateUpper", (void (ROL::BoundConstraint<double>::*)()) &ROL::BoundConstraint<double>::deactivateUpper, "C++: ROL::BoundConstraint<double>::deactivateUpper() --> void");
		cl.def("deactivate", (void (ROL::BoundConstraint<double>::*)()) &ROL::BoundConstraint<double>::deactivate, "C++: ROL::BoundConstraint<double>::deactivate() --> void");
		cl.def("isLowerActivated", (bool (ROL::BoundConstraint<double>::*)() const) &ROL::BoundConstraint<double>::isLowerActivated, "C++: ROL::BoundConstraint<double>::isLowerActivated() const --> bool");
		cl.def("isUpperActivated", (bool (ROL::BoundConstraint<double>::*)() const) &ROL::BoundConstraint<double>::isUpperActivated, "C++: ROL::BoundConstraint<double>::isUpperActivated() const --> bool");
		cl.def("isActivated", (bool (ROL::BoundConstraint<double>::*)() const) &ROL::BoundConstraint<double>::isActivated, "C++: ROL::BoundConstraint<double>::isActivated() const --> bool");
		cl.def("pruneActive", [](ROL::BoundConstraint<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1) -> void { return o.pruneActive(a0, a1); }, "", pybind11::arg("v"), pybind11::arg("x"));
		cl.def("pruneActive", (void (ROL::BoundConstraint<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, double)) &ROL::BoundConstraint<double>::pruneActive, "C++: ROL::BoundConstraint<double>::pruneActive(class ROL::Vector<double> &, const class ROL::Vector<double> &, double) --> void", pybind11::arg("v"), pybind11::arg("x"), pybind11::arg("eps"));
		cl.def("pruneActive", [](ROL::BoundConstraint<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2) -> void { return o.pruneActive(a0, a1, a2); }, "", pybind11::arg("v"), pybind11::arg("g"), pybind11::arg("x"));
		cl.def("pruneActive", [](ROL::BoundConstraint<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2, double const & a3) -> void { return o.pruneActive(a0, a1, a2, a3); }, "", pybind11::arg("v"), pybind11::arg("g"), pybind11::arg("x"), pybind11::arg("xeps"));
		cl.def("pruneActive", (void (ROL::BoundConstraint<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, double, double)) &ROL::BoundConstraint<double>::pruneActive, "C++: ROL::BoundConstraint<double>::pruneActive(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, double, double) --> void", pybind11::arg("v"), pybind11::arg("g"), pybind11::arg("x"), pybind11::arg("xeps"), pybind11::arg("geps"));
		cl.def("pruneLowerInactive", [](ROL::BoundConstraint<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1) -> void { return o.pruneLowerInactive(a0, a1); }, "", pybind11::arg("v"), pybind11::arg("x"));
		cl.def("pruneLowerInactive", (void (ROL::BoundConstraint<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, double)) &ROL::BoundConstraint<double>::pruneLowerInactive, "C++: ROL::BoundConstraint<double>::pruneLowerInactive(class ROL::Vector<double> &, const class ROL::Vector<double> &, double) --> void", pybind11::arg("v"), pybind11::arg("x"), pybind11::arg("eps"));
		cl.def("pruneUpperInactive", [](ROL::BoundConstraint<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1) -> void { return o.pruneUpperInactive(a0, a1); }, "", pybind11::arg("v"), pybind11::arg("x"));
		cl.def("pruneUpperInactive", (void (ROL::BoundConstraint<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, double)) &ROL::BoundConstraint<double>::pruneUpperInactive, "C++: ROL::BoundConstraint<double>::pruneUpperInactive(class ROL::Vector<double> &, const class ROL::Vector<double> &, double) --> void", pybind11::arg("v"), pybind11::arg("x"), pybind11::arg("eps"));
		cl.def("pruneLowerInactive", [](ROL::BoundConstraint<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2) -> void { return o.pruneLowerInactive(a0, a1, a2); }, "", pybind11::arg("v"), pybind11::arg("g"), pybind11::arg("x"));
		cl.def("pruneLowerInactive", [](ROL::BoundConstraint<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2, double const & a3) -> void { return o.pruneLowerInactive(a0, a1, a2, a3); }, "", pybind11::arg("v"), pybind11::arg("g"), pybind11::arg("x"), pybind11::arg("xeps"));
		cl.def("pruneLowerInactive", (void (ROL::BoundConstraint<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, double, double)) &ROL::BoundConstraint<double>::pruneLowerInactive, "C++: ROL::BoundConstraint<double>::pruneLowerInactive(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, double, double) --> void", pybind11::arg("v"), pybind11::arg("g"), pybind11::arg("x"), pybind11::arg("xeps"), pybind11::arg("geps"));
		cl.def("pruneUpperInactive", [](ROL::BoundConstraint<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2) -> void { return o.pruneUpperInactive(a0, a1, a2); }, "", pybind11::arg("v"), pybind11::arg("g"), pybind11::arg("x"));
		cl.def("pruneUpperInactive", [](ROL::BoundConstraint<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2, double const & a3) -> void { return o.pruneUpperInactive(a0, a1, a2, a3); }, "", pybind11::arg("v"), pybind11::arg("g"), pybind11::arg("x"), pybind11::arg("xeps"));
		cl.def("pruneUpperInactive", (void (ROL::BoundConstraint<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, double, double)) &ROL::BoundConstraint<double>::pruneUpperInactive, "C++: ROL::BoundConstraint<double>::pruneUpperInactive(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, double, double) --> void", pybind11::arg("v"), pybind11::arg("g"), pybind11::arg("x"), pybind11::arg("xeps"), pybind11::arg("geps"));
		cl.def("pruneInactive", [](ROL::BoundConstraint<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1) -> void { return o.pruneInactive(a0, a1); }, "", pybind11::arg("v"), pybind11::arg("x"));
		cl.def("pruneInactive", (void (ROL::BoundConstraint<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, double)) &ROL::BoundConstraint<double>::pruneInactive, "C++: ROL::BoundConstraint<double>::pruneInactive(class ROL::Vector<double> &, const class ROL::Vector<double> &, double) --> void", pybind11::arg("v"), pybind11::arg("x"), pybind11::arg("eps"));
		cl.def("pruneInactive", [](ROL::BoundConstraint<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2) -> void { return o.pruneInactive(a0, a1, a2); }, "", pybind11::arg("v"), pybind11::arg("g"), pybind11::arg("x"));
		cl.def("pruneInactive", [](ROL::BoundConstraint<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2, double const & a3) -> void { return o.pruneInactive(a0, a1, a2, a3); }, "", pybind11::arg("v"), pybind11::arg("g"), pybind11::arg("x"), pybind11::arg("xeps"));
		cl.def("pruneInactive", (void (ROL::BoundConstraint<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, double, double)) &ROL::BoundConstraint<double>::pruneInactive, "C++: ROL::BoundConstraint<double>::pruneInactive(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, double, double) --> void", pybind11::arg("v"), pybind11::arg("g"), pybind11::arg("x"), pybind11::arg("xeps"), pybind11::arg("geps"));
		cl.def("computeProjectedGradient", (void (ROL::BoundConstraint<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &)) &ROL::BoundConstraint<double>::computeProjectedGradient, "C++: ROL::BoundConstraint<double>::computeProjectedGradient(class ROL::Vector<double> &, const class ROL::Vector<double> &) --> void", pybind11::arg("g"), pybind11::arg("x"));
		cl.def("computeProjectedStep", (void (ROL::BoundConstraint<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &)) &ROL::BoundConstraint<double>::computeProjectedStep, "C++: ROL::BoundConstraint<double>::computeProjectedStep(class ROL::Vector<double> &, const class ROL::Vector<double> &) --> void", pybind11::arg("v"), pybind11::arg("x"));
		cl.def("assign", (class ROL::BoundConstraint<double> & (ROL::BoundConstraint<double>::*)(const class ROL::BoundConstraint<double> &)) &ROL::BoundConstraint<double>::operator=, "C++: ROL::BoundConstraint<double>::operator=(const class ROL::BoundConstraint<double> &) --> class ROL::BoundConstraint<double> &", pybind11::return_value_policy::automatic, pybind11::arg(""));
	}
}
