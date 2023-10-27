#include <ROL_BoundConstraint.hpp>
#include <ROL_BrentsProjection.hpp>
#include <ROL_Constraint.hpp>
#include <ROL_Elementwise_Function.hpp>
#include <ROL_Elementwise_Reduce.hpp>
#include <ROL_Krylov.hpp>
#include <ROL_KrylovFactory.hpp>
#include <ROL_LinearOperator.hpp>
#include <ROL_PolyhedralProjection.hpp>
#include <ROL_PolyhedralProjectionFactory.hpp>
#include <ROL_RiddersProjection.hpp>
#include <ROL_SemismoothNewtonProjection.hpp>
#include <ROL_UpdateType.hpp>
#include <ROL_Vector.hpp>
#include <Teuchos_ENull.hpp>
#include <Teuchos_FilteredIterator.hpp>
#include <Teuchos_ParameterEntry.hpp>
#include <Teuchos_ParameterEntryValidator.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_ParameterListModifier.hpp>
#include <Teuchos_PtrDecl.hpp>
#include <Teuchos_RCPDecl.hpp>
#include <Teuchos_RCPNode.hpp>
#include <Teuchos_StringIndexedOrderedValueObjectContainer.hpp>
#include <Teuchos_any.hpp>
#include <cwchar>
#include <deque>
#include <ios>
#include <iterator>
#include <locale>
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

// ROL::SemismoothNewtonProjection file:ROL_SemismoothNewtonProjection.hpp line:55
struct PyCallBack_ROL_SemismoothNewtonProjection_double_t : public ROL::SemismoothNewtonProjection<double> {
	using ROL::SemismoothNewtonProjection<double>::SemismoothNewtonProjection;

	void project(class ROL::Vector<double> & a0, std::ostream & a1) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::SemismoothNewtonProjection<double> *>(this), "project");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return SemismoothNewtonProjection::project(a0, a1);
	}
};

// ROL::RiddersProjection file:ROL_RiddersProjection.hpp line:54
struct PyCallBack_ROL_RiddersProjection_double_t : public ROL::RiddersProjection<double> {
	using ROL::RiddersProjection<double>::RiddersProjection;

	void project(class ROL::Vector<double> & a0, std::ostream & a1) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::RiddersProjection<double> *>(this), "project");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return RiddersProjection::project(a0, a1);
	}
};

// ROL::BrentsProjection file:ROL_BrentsProjection.hpp line:54
struct PyCallBack_ROL_BrentsProjection_double_t : public ROL::BrentsProjection<double> {
	using ROL::BrentsProjection<double>::BrentsProjection;

	void project(class ROL::Vector<double> & a0, std::ostream & a1) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::BrentsProjection<double> *>(this), "project");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return BrentsProjection::project(a0, a1);
	}
};

void bind_ROL_KrylovFactory(std::function< pybind11::module &(std::string const &namespace_) > &M)
{
	// ROL::EKrylov file:ROL_KrylovFactory.hpp line:66
	pybind11::enum_<ROL::EKrylov>(M("ROL"), "EKrylov", pybind11::arithmetic(), "Enumeration of Krylov methods.\n\n      \n    CG          Conjugate Gradient Method\n      \n\n    CR          Conjugate Residual Method\n      \n\n    GMRES       Generalized Minimum Residual Method\n      \n\n    MINRES      Minimum Residual Method\n      \n\n    USERDEFINED User defined Krylov method\n      \n\n    LAST        Dummy type", pybind11::module_local())
		.value("KRYLOV_CG", ROL::KRYLOV_CG)
		.value("KRYLOV_CR", ROL::KRYLOV_CR)
		.value("KRYLOV_GMRES", ROL::KRYLOV_GMRES)
		.value("KRYLOV_MINRES", ROL::KRYLOV_MINRES)
		.value("KRYLOV_USERDEFINED", ROL::KRYLOV_USERDEFINED)
		.value("KRYLOV_LAST", ROL::KRYLOV_LAST)
		.export_values();

;

	// ROL::EKrylovToString(enum ROL::EKrylov) file:ROL_KrylovFactory.hpp line:75
	M("ROL").def("EKrylovToString", (std::string (*)(enum ROL::EKrylov)) &ROL::EKrylovToString, "C++: ROL::EKrylovToString(enum ROL::EKrylov) --> std::string", pybind11::arg("type"));

	// ROL::isValidKrylov(enum ROL::EKrylov) file:ROL_KrylovFactory.hpp line:94
	M("ROL").def("isValidKrylov", (int (*)(enum ROL::EKrylov)) &ROL::isValidKrylov, "Verifies validity of a Krylov enum.\n\n      \n  [in]  - enum of the Krylov\n      \n\n 1 if the argument is a valid Secant; 0 otherwise.\n\nC++: ROL::isValidKrylov(enum ROL::EKrylov) --> int", pybind11::arg("type"));

	// ROL::StringToEKrylov(std::string) file:ROL_KrylovFactory.hpp line:121
	M("ROL").def("StringToEKrylov", (enum ROL::EKrylov (*)(std::string)) &ROL::StringToEKrylov, "C++: ROL::StringToEKrylov(std::string) --> enum ROL::EKrylov", pybind11::arg("s"));

	// ROL::KrylovFactory(class Teuchos::ParameterList &) file:ROL_KrylovFactory.hpp line:132
	M("ROL").def("KrylovFactory", (class Teuchos::RCP<class ROL::Krylov<double> > (*)(class Teuchos::ParameterList &)) &ROL::KrylovFactory<double>, "C++: ROL::KrylovFactory(class Teuchos::ParameterList &) --> class Teuchos::RCP<class ROL::Krylov<double> >", pybind11::arg("parlist"));

	{ // ROL::SemismoothNewtonProjection file:ROL_SemismoothNewtonProjection.hpp line:55
		pybind11::class_<ROL::SemismoothNewtonProjection<double>, Teuchos::RCP<ROL::SemismoothNewtonProjection<double>>, PyCallBack_ROL_SemismoothNewtonProjection_double_t, ROL::PolyhedralProjection<double>> cl(M("ROL"), "SemismoothNewtonProjection_double_t", "", pybind11::module_local());
		cl.def( pybind11::init<const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class Teuchos::RCP<class ROL::BoundConstraint<double> > &, const class Teuchos::RCP<class ROL::Constraint<double> > &, const class ROL::Vector<double> &, const class ROL::Vector<double> &>(), pybind11::arg("xprim"), pybind11::arg("xdual"), pybind11::arg("bnd"), pybind11::arg("con"), pybind11::arg("mul"), pybind11::arg("res") );

		cl.def( pybind11::init<const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class Teuchos::RCP<class ROL::BoundConstraint<double> > &, const class Teuchos::RCP<class ROL::Constraint<double> > &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, class Teuchos::ParameterList &>(), pybind11::arg("xprim"), pybind11::arg("xdual"), pybind11::arg("bnd"), pybind11::arg("con"), pybind11::arg("mul"), pybind11::arg("res"), pybind11::arg("list") );

		cl.def( pybind11::init( [](PyCallBack_ROL_SemismoothNewtonProjection_double_t const &o){ return new PyCallBack_ROL_SemismoothNewtonProjection_double_t(o); } ) );
		cl.def( pybind11::init( [](ROL::SemismoothNewtonProjection<double> const &o){ return new ROL::SemismoothNewtonProjection<double>(o); } ) );
		cl.def("project", [](ROL::SemismoothNewtonProjection<double> &o, class ROL::Vector<double> & a0) -> void { return o.project(a0); }, "", pybind11::arg("x"));
		cl.def("project", (void (ROL::SemismoothNewtonProjection<double>::*)(class ROL::Vector<double> &, std::ostream &)) &ROL::SemismoothNewtonProjection<double>::project, "C++: ROL::SemismoothNewtonProjection<double>::project(class ROL::Vector<double> &, std::ostream &) --> void", pybind11::arg("x"), pybind11::arg("stream"));
		cl.def("project", [](ROL::PolyhedralProjection<double> &o, class ROL::Vector<double> & a0) -> void { return o.project(a0); }, "", pybind11::arg("x"));
		cl.def("project", (void (ROL::PolyhedralProjection<double>::*)(class ROL::Vector<double> &, std::ostream &)) &ROL::PolyhedralProjection<double>::project, "C++: ROL::PolyhedralProjection<double>::project(class ROL::Vector<double> &, std::ostream &) --> void", pybind11::arg("x"), pybind11::arg("stream"));
		cl.def("getLinearConstraint", (const class Teuchos::RCP<class ROL::Constraint<double> > (ROL::PolyhedralProjection<double>::*)() const) &ROL::PolyhedralProjection<double>::getLinearConstraint, "C++: ROL::PolyhedralProjection<double>::getLinearConstraint() const --> const class Teuchos::RCP<class ROL::Constraint<double> >");
		cl.def("getBoundConstraint", (const class Teuchos::RCP<class ROL::BoundConstraint<double> > (ROL::PolyhedralProjection<double>::*)() const) &ROL::PolyhedralProjection<double>::getBoundConstraint, "C++: ROL::PolyhedralProjection<double>::getBoundConstraint() const --> const class Teuchos::RCP<class ROL::BoundConstraint<double> >");
		cl.def("getMultiplier", (const class Teuchos::RCP<class ROL::Vector<double> > (ROL::PolyhedralProjection<double>::*)() const) &ROL::PolyhedralProjection<double>::getMultiplier, "C++: ROL::PolyhedralProjection<double>::getMultiplier() const --> const class Teuchos::RCP<class ROL::Vector<double> >");
		cl.def("getResidual", (const class Teuchos::RCP<class ROL::Vector<double> > (ROL::PolyhedralProjection<double>::*)() const) &ROL::PolyhedralProjection<double>::getResidual, "C++: ROL::PolyhedralProjection<double>::getResidual() const --> const class Teuchos::RCP<class ROL::Vector<double> >");
	}
	{ // ROL::RiddersProjection file:ROL_RiddersProjection.hpp line:54
		pybind11::class_<ROL::RiddersProjection<double>, Teuchos::RCP<ROL::RiddersProjection<double>>, PyCallBack_ROL_RiddersProjection_double_t, ROL::PolyhedralProjection<double>> cl(M("ROL"), "RiddersProjection_double_t", "", pybind11::module_local());
		cl.def( pybind11::init<const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class Teuchos::RCP<class ROL::BoundConstraint<double> > &, const class Teuchos::RCP<class ROL::Constraint<double> > &, const class ROL::Vector<double> &, const class ROL::Vector<double> &>(), pybind11::arg("xprim"), pybind11::arg("xdual"), pybind11::arg("bnd"), pybind11::arg("con"), pybind11::arg("mul"), pybind11::arg("res") );

		cl.def( pybind11::init<const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class Teuchos::RCP<class ROL::BoundConstraint<double> > &, const class Teuchos::RCP<class ROL::Constraint<double> > &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, class Teuchos::ParameterList &>(), pybind11::arg("xprim"), pybind11::arg("xdual"), pybind11::arg("bnd"), pybind11::arg("con"), pybind11::arg("mul"), pybind11::arg("res"), pybind11::arg("list") );

		cl.def( pybind11::init( [](PyCallBack_ROL_RiddersProjection_double_t const &o){ return new PyCallBack_ROL_RiddersProjection_double_t(o); } ) );
		cl.def( pybind11::init( [](ROL::RiddersProjection<double> const &o){ return new ROL::RiddersProjection<double>(o); } ) );
		cl.def("project", [](ROL::RiddersProjection<double> &o, class ROL::Vector<double> & a0) -> void { return o.project(a0); }, "", pybind11::arg("x"));
		cl.def("project", (void (ROL::RiddersProjection<double>::*)(class ROL::Vector<double> &, std::ostream &)) &ROL::RiddersProjection<double>::project, "C++: ROL::RiddersProjection<double>::project(class ROL::Vector<double> &, std::ostream &) --> void", pybind11::arg("x"), pybind11::arg("stream"));
		cl.def("project", [](ROL::PolyhedralProjection<double> &o, class ROL::Vector<double> & a0) -> void { return o.project(a0); }, "", pybind11::arg("x"));
		cl.def("project", (void (ROL::PolyhedralProjection<double>::*)(class ROL::Vector<double> &, std::ostream &)) &ROL::PolyhedralProjection<double>::project, "C++: ROL::PolyhedralProjection<double>::project(class ROL::Vector<double> &, std::ostream &) --> void", pybind11::arg("x"), pybind11::arg("stream"));
		cl.def("getLinearConstraint", (const class Teuchos::RCP<class ROL::Constraint<double> > (ROL::PolyhedralProjection<double>::*)() const) &ROL::PolyhedralProjection<double>::getLinearConstraint, "C++: ROL::PolyhedralProjection<double>::getLinearConstraint() const --> const class Teuchos::RCP<class ROL::Constraint<double> >");
		cl.def("getBoundConstraint", (const class Teuchos::RCP<class ROL::BoundConstraint<double> > (ROL::PolyhedralProjection<double>::*)() const) &ROL::PolyhedralProjection<double>::getBoundConstraint, "C++: ROL::PolyhedralProjection<double>::getBoundConstraint() const --> const class Teuchos::RCP<class ROL::BoundConstraint<double> >");
		cl.def("getMultiplier", (const class Teuchos::RCP<class ROL::Vector<double> > (ROL::PolyhedralProjection<double>::*)() const) &ROL::PolyhedralProjection<double>::getMultiplier, "C++: ROL::PolyhedralProjection<double>::getMultiplier() const --> const class Teuchos::RCP<class ROL::Vector<double> >");
		cl.def("getResidual", (const class Teuchos::RCP<class ROL::Vector<double> > (ROL::PolyhedralProjection<double>::*)() const) &ROL::PolyhedralProjection<double>::getResidual, "C++: ROL::PolyhedralProjection<double>::getResidual() const --> const class Teuchos::RCP<class ROL::Vector<double> >");
	}
	{ // ROL::BrentsProjection file:ROL_BrentsProjection.hpp line:54
		pybind11::class_<ROL::BrentsProjection<double>, Teuchos::RCP<ROL::BrentsProjection<double>>, PyCallBack_ROL_BrentsProjection_double_t, ROL::PolyhedralProjection<double>> cl(M("ROL"), "BrentsProjection_double_t", "", pybind11::module_local());
		cl.def( pybind11::init<const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class Teuchos::RCP<class ROL::BoundConstraint<double> > &, const class Teuchos::RCP<class ROL::Constraint<double> > &, const class ROL::Vector<double> &, const class ROL::Vector<double> &>(), pybind11::arg("xprim"), pybind11::arg("xdual"), pybind11::arg("bnd"), pybind11::arg("con"), pybind11::arg("mul"), pybind11::arg("res") );

		cl.def( pybind11::init<const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class Teuchos::RCP<class ROL::BoundConstraint<double> > &, const class Teuchos::RCP<class ROL::Constraint<double> > &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, class Teuchos::ParameterList &>(), pybind11::arg("xprim"), pybind11::arg("xdual"), pybind11::arg("bnd"), pybind11::arg("con"), pybind11::arg("mul"), pybind11::arg("res"), pybind11::arg("list") );

		cl.def( pybind11::init( [](PyCallBack_ROL_BrentsProjection_double_t const &o){ return new PyCallBack_ROL_BrentsProjection_double_t(o); } ) );
		cl.def( pybind11::init( [](ROL::BrentsProjection<double> const &o){ return new ROL::BrentsProjection<double>(o); } ) );
		cl.def("project", [](ROL::BrentsProjection<double> &o, class ROL::Vector<double> & a0) -> void { return o.project(a0); }, "", pybind11::arg("x"));
		cl.def("project", (void (ROL::BrentsProjection<double>::*)(class ROL::Vector<double> &, std::ostream &)) &ROL::BrentsProjection<double>::project, "C++: ROL::BrentsProjection<double>::project(class ROL::Vector<double> &, std::ostream &) --> void", pybind11::arg("x"), pybind11::arg("stream"));
		cl.def("project", [](ROL::PolyhedralProjection<double> &o, class ROL::Vector<double> & a0) -> void { return o.project(a0); }, "", pybind11::arg("x"));
		cl.def("project", (void (ROL::PolyhedralProjection<double>::*)(class ROL::Vector<double> &, std::ostream &)) &ROL::PolyhedralProjection<double>::project, "C++: ROL::PolyhedralProjection<double>::project(class ROL::Vector<double> &, std::ostream &) --> void", pybind11::arg("x"), pybind11::arg("stream"));
		cl.def("getLinearConstraint", (const class Teuchos::RCP<class ROL::Constraint<double> > (ROL::PolyhedralProjection<double>::*)() const) &ROL::PolyhedralProjection<double>::getLinearConstraint, "C++: ROL::PolyhedralProjection<double>::getLinearConstraint() const --> const class Teuchos::RCP<class ROL::Constraint<double> >");
		cl.def("getBoundConstraint", (const class Teuchos::RCP<class ROL::BoundConstraint<double> > (ROL::PolyhedralProjection<double>::*)() const) &ROL::PolyhedralProjection<double>::getBoundConstraint, "C++: ROL::PolyhedralProjection<double>::getBoundConstraint() const --> const class Teuchos::RCP<class ROL::BoundConstraint<double> >");
		cl.def("getMultiplier", (const class Teuchos::RCP<class ROL::Vector<double> > (ROL::PolyhedralProjection<double>::*)() const) &ROL::PolyhedralProjection<double>::getMultiplier, "C++: ROL::PolyhedralProjection<double>::getMultiplier() const --> const class Teuchos::RCP<class ROL::Vector<double> >");
		cl.def("getResidual", (const class Teuchos::RCP<class ROL::Vector<double> > (ROL::PolyhedralProjection<double>::*)() const) &ROL::PolyhedralProjection<double>::getResidual, "C++: ROL::PolyhedralProjection<double>::getResidual() const --> const class Teuchos::RCP<class ROL::Vector<double> >");
	}
	// ROL::EPolyProjAlgo file:ROL_PolyhedralProjectionFactory.hpp line:65
	pybind11::enum_<ROL::EPolyProjAlgo>(M("ROL"), "EPolyProjAlgo", pybind11::arithmetic(), "Enumeration of polyhedral projecdtion algorithm types.\n\n    \n    PPA_DAIFLETCHER describe\n    \n\n    PPA_DYKSTRA     describe\n    \n\n    PPA_NEWTON      describe\n    \n\n    PPA_RIDDERS     describe", pybind11::module_local())
		.value("PPA_DAIFLETCHER", ROL::PPA_DAIFLETCHER)
		.value("PPA_DYKSTRA", ROL::PPA_DYKSTRA)
		.value("PPA_DOUGLASRACHFORD", ROL::PPA_DOUGLASRACHFORD)
		.value("PPA_NEWTON", ROL::PPA_NEWTON)
		.value("PPA_RIDDERS", ROL::PPA_RIDDERS)
		.value("PPA_BRENTS", ROL::PPA_BRENTS)
		.value("PPA_LAST", ROL::PPA_LAST)
		.export_values();

;

	// ROL::EPolyProjAlgoToString(enum ROL::EPolyProjAlgo) file:ROL_PolyhedralProjectionFactory.hpp line:75
	M("ROL").def("EPolyProjAlgoToString", (std::string (*)(enum ROL::EPolyProjAlgo)) &ROL::EPolyProjAlgoToString, "C++: ROL::EPolyProjAlgoToString(enum ROL::EPolyProjAlgo) --> std::string", pybind11::arg("alg"));

	// ROL::isValidPolyProjAlgo(enum ROL::EPolyProjAlgo) file:ROL_PolyhedralProjectionFactory.hpp line:95
	M("ROL").def("isValidPolyProjAlgo", (int (*)(enum ROL::EPolyProjAlgo)) &ROL::isValidPolyProjAlgo, "Verifies validity of a PolyProjAlgo enum.\n\n    \n  [in]  - enum of the PolyProjAlgo\n    \n\n 1 if the argument is a valid PolyProjAlgo; 0 otherwise.\n\nC++: ROL::isValidPolyProjAlgo(enum ROL::EPolyProjAlgo) --> int", pybind11::arg("alg"));

	// ROL::StringToEPolyProjAlgo(std::string) file:ROL_PolyhedralProjectionFactory.hpp line:126
	M("ROL").def("StringToEPolyProjAlgo", (enum ROL::EPolyProjAlgo (*)(std::string)) &ROL::StringToEPolyProjAlgo, "C++: ROL::StringToEPolyProjAlgo(std::string) --> enum ROL::EPolyProjAlgo", pybind11::arg("s"));

	// ROL::PolyhedralProjectionFactory(const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class Teuchos::RCP<class ROL::BoundConstraint<double> > &, const class Teuchos::RCP<class ROL::Constraint<double> > &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, class Teuchos::ParameterList &) file:ROL_PolyhedralProjectionFactory.hpp line:137
	M("ROL").def("PolyhedralProjectionFactory", (class Teuchos::RCP<class ROL::PolyhedralProjection<double> > (*)(const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class Teuchos::RCP<class ROL::BoundConstraint<double> > &, const class Teuchos::RCP<class ROL::Constraint<double> > &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, class Teuchos::ParameterList &)) &ROL::PolyhedralProjectionFactory<double>, "C++: ROL::PolyhedralProjectionFactory(const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class Teuchos::RCP<class ROL::BoundConstraint<double> > &, const class Teuchos::RCP<class ROL::Constraint<double> > &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, class Teuchos::ParameterList &) --> class Teuchos::RCP<class ROL::PolyhedralProjection<double> >", pybind11::arg("xprim"), pybind11::arg("xdual"), pybind11::arg("bnd"), pybind11::arg("con"), pybind11::arg("mul"), pybind11::arg("res"), pybind11::arg("list"));

}
