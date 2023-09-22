#include <ROL_BoundConstraint.hpp>
#include <ROL_Constraint.hpp>
#include <ROL_Elementwise_Function.hpp>
#include <ROL_Elementwise_Reduce.hpp>
#include <ROL_Objective.hpp>
#include <ROL_PolyhedralProjection.hpp>
#include <ROL_Problem.hpp>
#include <ROL_StatusTest.hpp>
#include <ROL_TypeE_Algorithm.hpp>
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
#include <locale>
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
#include <Teuchos_RCP.hpp>


#ifndef BINDER_PYBIND11_TYPE_CASTER
	#define BINDER_PYBIND11_TYPE_CASTER
	PYBIND11_DECLARE_HOLDER_TYPE(T, Teuchos::RCP<T>)
	PYBIND11_DECLARE_HOLDER_TYPE(T, T*)
	PYBIND11_MAKE_OPAQUE(Teuchos::RCP<void>)
#endif

// ROL::TypeE::Algorithm file:ROL_TypeE_Algorithm.hpp line:88
struct PyCallBack_ROL_TypeE_Algorithm_double_t : public ROL::TypeE::Algorithm<double> {
	using ROL::TypeE::Algorithm<double>::Algorithm;

	void run(class ROL::Problem<double> & a0, std::ostream & a1) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TypeE::Algorithm<double> *>(this), "run");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Algorithm::run(a0, a1);
	}
	void run(class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::Constraint<double> & a2, class ROL::Vector<double> & a3, std::ostream & a4) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TypeE::Algorithm<double> *>(this), "run");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Algorithm::run(a0, a1, a2, a3, a4);
	}
	void run(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4, const class ROL::Vector<double> & a5, std::ostream & a6) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TypeE::Algorithm<double> *>(this), "run");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4, a5, a6);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		pybind11::pybind11_fail("Tried to call pure virtual function \"Algorithm::run\"");
	}
	void run(class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::Constraint<double> & a2, class ROL::Vector<double> & a3, class ROL::Constraint<double> & a4, class ROL::Vector<double> & a5, std::ostream & a6) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TypeE::Algorithm<double> *>(this), "run");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4, a5, a6);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Algorithm::run(a0, a1, a2, a3, a4, a5, a6);
	}
	void run(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4, const class ROL::Vector<double> & a5, class ROL::Constraint<double> & a6, class ROL::Vector<double> & a7, const class ROL::Vector<double> & a8, std::ostream & a9) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TypeE::Algorithm<double> *>(this), "run");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Algorithm::run(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9);
	}
	void writeHeader(std::ostream & a0) const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TypeE::Algorithm<double> *>(this), "writeHeader");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Algorithm::writeHeader(a0);
	}
	void writeName(std::ostream & a0) const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TypeE::Algorithm<double> *>(this), "writeName");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Algorithm::writeName(a0);
	}
	void writeOutput(std::ostream & a0, const bool a1) const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TypeE::Algorithm<double> *>(this), "writeOutput");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Algorithm::writeOutput(a0, a1);
	}
	void writeExitStatus(std::ostream & a0) const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TypeE::Algorithm<double> *>(this), "writeExitStatus");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Algorithm::writeExitStatus(a0);
	}
};

void bind_ROL_TypeE_Algorithm(std::function< pybind11::module &(std::string const &namespace_) > &M)
{
	{ // ROL::TypeE::AlgorithmState file:ROL_TypeE_Algorithm.hpp line:60
		pybind11::class_<ROL::TypeE::AlgorithmState<double>, Teuchos::RCP<ROL::TypeE::AlgorithmState<double>>, ROL::AlgorithmState<double>> cl(M("ROL::TypeE"), "AlgorithmState_double_t", "", pybind11::module_local());
		cl.def( pybind11::init( [](){ return new ROL::TypeE::AlgorithmState<double>(); } ) );
		cl.def( pybind11::init( [](ROL::TypeE::AlgorithmState<double> const &o){ return new ROL::TypeE::AlgorithmState<double>(o); } ) );
		cl.def_readwrite("searchSize", &ROL::TypeE::AlgorithmState<double>::searchSize);
		cl.def_readwrite("stepVec", &ROL::TypeE::AlgorithmState<double>::stepVec);
		cl.def_readwrite("gradientVec", &ROL::TypeE::AlgorithmState<double>::gradientVec);
		cl.def_readwrite("constraintVec", &ROL::TypeE::AlgorithmState<double>::constraintVec);
		cl.def("reset", (void (ROL::TypeE::AlgorithmState<double>::*)()) &ROL::TypeE::AlgorithmState<double>::reset, "C++: ROL::TypeE::AlgorithmState<double>::reset() --> void");
		cl.def("assign", (struct ROL::TypeE::AlgorithmState<double> & (ROL::TypeE::AlgorithmState<double>::*)(const struct ROL::TypeE::AlgorithmState<double> &)) &ROL::TypeE::AlgorithmState<double>::operator=, "C++: ROL::TypeE::AlgorithmState<double>::operator=(const struct ROL::TypeE::AlgorithmState<double> &) --> struct ROL::TypeE::AlgorithmState<double> &", pybind11::return_value_policy::automatic, pybind11::arg(""));
		cl.def_readwrite("iter", &ROL::AlgorithmState<double>::iter);
		cl.def_readwrite("minIter", &ROL::AlgorithmState<double>::minIter);
		cl.def_readwrite("nfval", &ROL::AlgorithmState<double>::nfval);
		cl.def_readwrite("ncval", &ROL::AlgorithmState<double>::ncval);
		cl.def_readwrite("ngrad", &ROL::AlgorithmState<double>::ngrad);
		cl.def_readwrite("value", &ROL::AlgorithmState<double>::value);
		cl.def_readwrite("minValue", &ROL::AlgorithmState<double>::minValue);
		cl.def_readwrite("gnorm", &ROL::AlgorithmState<double>::gnorm);
		cl.def_readwrite("cnorm", &ROL::AlgorithmState<double>::cnorm);
		cl.def_readwrite("snorm", &ROL::AlgorithmState<double>::snorm);
		cl.def_readwrite("aggregateGradientNorm", &ROL::AlgorithmState<double>::aggregateGradientNorm);
		cl.def_readwrite("aggregateModelError", &ROL::AlgorithmState<double>::aggregateModelError);
		cl.def_readwrite("flag", &ROL::AlgorithmState<double>::flag);
		cl.def_readwrite("iterateVec", &ROL::AlgorithmState<double>::iterateVec);
		cl.def_readwrite("lagmultVec", &ROL::AlgorithmState<double>::lagmultVec);
		cl.def_readwrite("minIterVec", &ROL::AlgorithmState<double>::minIterVec);
		cl.def_readwrite("statusFlag", &ROL::AlgorithmState<double>::statusFlag);
		cl.def("reset", (void (ROL::AlgorithmState<double>::*)()) &ROL::AlgorithmState<double>::reset, "C++: ROL::AlgorithmState<double>::reset() --> void");
		cl.def("assign", (struct ROL::AlgorithmState<double> & (ROL::AlgorithmState<double>::*)(const struct ROL::AlgorithmState<double> &)) &ROL::AlgorithmState<double>::operator=, "C++: ROL::AlgorithmState<double>::operator=(const struct ROL::AlgorithmState<double> &) --> struct ROL::AlgorithmState<double> &", pybind11::return_value_policy::automatic, pybind11::arg(""));
	}
	{ // ROL::TypeE::Algorithm file:ROL_TypeE_Algorithm.hpp line:88
		pybind11::class_<ROL::TypeE::Algorithm<double>, Teuchos::RCP<ROL::TypeE::Algorithm<double>>, PyCallBack_ROL_TypeE_Algorithm_double_t> cl(M("ROL::TypeE"), "Algorithm_double_t", "", pybind11::module_local());
		cl.def( pybind11::init( [](){ return new PyCallBack_ROL_TypeE_Algorithm_double_t(); } ) );
		cl.def(pybind11::init<PyCallBack_ROL_TypeE_Algorithm_double_t const &>());
		cl.def("setStatusTest", [](ROL::TypeE::Algorithm<double> &o, const class Teuchos::RCP<class ROL::StatusTest<double> > & a0) -> void { return o.setStatusTest(a0); }, "", pybind11::arg("status"));
		cl.def("setStatusTest", (void (ROL::TypeE::Algorithm<double>::*)(const class Teuchos::RCP<class ROL::StatusTest<double> > &, bool)) &ROL::TypeE::Algorithm<double>::setStatusTest, "C++: ROL::TypeE::Algorithm<double>::setStatusTest(const class Teuchos::RCP<class ROL::StatusTest<double> > &, bool) --> void", pybind11::arg("status"), pybind11::arg("combineStatus"));
		cl.def("run", [](ROL::TypeE::Algorithm<double> &o, class ROL::Problem<double> & a0) -> void { return o.run(a0); }, "", pybind11::arg("problem"));
		cl.def("run", (void (ROL::TypeE::Algorithm<double>::*)(class ROL::Problem<double> &, std::ostream &)) &ROL::TypeE::Algorithm<double>::run, "C++: ROL::TypeE::Algorithm<double>::run(class ROL::Problem<double> &, std::ostream &) --> void", pybind11::arg("problem"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeE::Algorithm<double> &o, class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::Constraint<double> & a2, class ROL::Vector<double> & a3) -> void { return o.run(a0, a1, a2, a3); }, "", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("econ"), pybind11::arg("emul"));
		cl.def("run", (void (ROL::TypeE::Algorithm<double>::*)(class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, std::ostream &)) &ROL::TypeE::Algorithm<double>::run, "C++: ROL::TypeE::Algorithm<double>::run(class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, std::ostream &) --> void", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeE::Algorithm<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4, const class ROL::Vector<double> & a5) -> void { return o.run(a0, a1, a2, a3, a4, a5); }, "", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("eres"));
		cl.def("run", (void (ROL::TypeE::Algorithm<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, const class ROL::Vector<double> &, std::ostream &)) &ROL::TypeE::Algorithm<double>::run, "C++: ROL::TypeE::Algorithm<double>::run(class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, const class ROL::Vector<double> &, std::ostream &) --> void", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("eres"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeE::Algorithm<double> &o, class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::Constraint<double> & a2, class ROL::Vector<double> & a3, class ROL::Constraint<double> & a4, class ROL::Vector<double> & a5) -> void { return o.run(a0, a1, a2, a3, a4, a5); }, "", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"));
		cl.def("run", (void (ROL::TypeE::Algorithm<double>::*)(class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, std::ostream &)) &ROL::TypeE::Algorithm<double>::run, "C++: ROL::TypeE::Algorithm<double>::run(class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, std::ostream &) --> void", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeE::Algorithm<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4, const class ROL::Vector<double> & a5, class ROL::Constraint<double> & a6, class ROL::Vector<double> & a7, const class ROL::Vector<double> & a8) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6, a7, a8); }, "", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("eres"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"), pybind11::arg("linear_eres"));
		cl.def("run", (void (ROL::TypeE::Algorithm<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, const class ROL::Vector<double> &, std::ostream &)) &ROL::TypeE::Algorithm<double>::run, "C++: ROL::TypeE::Algorithm<double>::run(class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, const class ROL::Vector<double> &, std::ostream &) --> void", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("eres"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"), pybind11::arg("linear_eres"), pybind11::arg("outStream"));
		cl.def("writeHeader", (void (ROL::TypeE::Algorithm<double>::*)(std::ostream &) const) &ROL::TypeE::Algorithm<double>::writeHeader, "C++: ROL::TypeE::Algorithm<double>::writeHeader(std::ostream &) const --> void", pybind11::arg("os"));
		cl.def("writeName", (void (ROL::TypeE::Algorithm<double>::*)(std::ostream &) const) &ROL::TypeE::Algorithm<double>::writeName, "C++: ROL::TypeE::Algorithm<double>::writeName(std::ostream &) const --> void", pybind11::arg("os"));
		cl.def("writeOutput", [](ROL::TypeE::Algorithm<double> const &o, std::ostream & a0) -> void { return o.writeOutput(a0); }, "", pybind11::arg("os"));
		cl.def("writeOutput", (void (ROL::TypeE::Algorithm<double>::*)(std::ostream &, const bool) const) &ROL::TypeE::Algorithm<double>::writeOutput, "C++: ROL::TypeE::Algorithm<double>::writeOutput(std::ostream &, const bool) const --> void", pybind11::arg("os"), pybind11::arg("write_header"));
		cl.def("writeExitStatus", (void (ROL::TypeE::Algorithm<double>::*)(std::ostream &) const) &ROL::TypeE::Algorithm<double>::writeExitStatus, "C++: ROL::TypeE::Algorithm<double>::writeExitStatus(std::ostream &) const --> void", pybind11::arg("os"));
		cl.def("getState", (class Teuchos::RCP<const struct ROL::TypeE::AlgorithmState<double> > (ROL::TypeE::Algorithm<double>::*)() const) &ROL::TypeE::Algorithm<double>::getState, "C++: ROL::TypeE::Algorithm<double>::getState() const --> class Teuchos::RCP<const struct ROL::TypeE::AlgorithmState<double> >");
		cl.def("reset", (void (ROL::TypeE::Algorithm<double>::*)()) &ROL::TypeE::Algorithm<double>::reset, "C++: ROL::TypeE::Algorithm<double>::reset() --> void");
	}
}
