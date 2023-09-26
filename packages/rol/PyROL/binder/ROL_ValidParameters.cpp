#include <ROL_Bundle_U.hpp>
#include <ROL_Elementwise_Function.hpp>
#include <ROL_Elementwise_Reduce.hpp>
#include <ROL_LineSearch_U.hpp>
#include <ROL_LineSearch_U_Types.hpp>
#include <ROL_Objective.hpp>
#include <ROL_ScalarFunction.hpp>
#include <ROL_UpdateType.hpp>
#include <ROL_ValidParameters.hpp>
#include <ROL_Vector.hpp>
#include <Teuchos_ENull.hpp>
#include <Teuchos_FilteredIterator.hpp>
#include <Teuchos_ParameterEntry.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_ParameterListModifier.hpp>
#include <Teuchos_RCPDecl.hpp>
#include <Teuchos_RCPNode.hpp>
#include <Teuchos_StringIndexedOrderedValueObjectContainer.hpp>
#include <Teuchos_any.hpp>
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

// ROL::Bundle_U file:ROL_Bundle_U.hpp line:59
struct PyCallBack_ROL_Bundle_U_double_t : public ROL::Bundle_U<double> {
	using ROL::Bundle_U<double>::Bundle_U;

	void initialize(const class ROL::Vector<double> & a0) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::Bundle_U<double> *>(this), "initialize");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Bundle_U::initialize(a0);
	}
	unsigned int solveDual(const double a0, const unsigned int a1, const double a2) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::Bundle_U<double> *>(this), "solveDual");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2);
			if (pybind11::detail::cast_is_temporary_value_reference<unsigned int>::value) {
				static pybind11::detail::override_caster_t<unsigned int> caster;
				return pybind11::detail::cast_ref<unsigned int>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<unsigned int>(std::move(o));
		}
		pybind11::pybind11_fail("Tried to call pure virtual function \"Bundle_U::solveDual\"");
	}
};

// ROL::ScalarFunction file:ROL_ScalarFunction.hpp line:56
struct PyCallBack_ROL_ScalarFunction_double_t : public ROL::ScalarFunction<double> {
	using ROL::ScalarFunction<double>::ScalarFunction;

	double value(const double a0) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::ScalarFunction<double> *>(this), "value");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		pybind11::pybind11_fail("Tried to call pure virtual function \"ScalarFunction::value\"");
	}
	double deriv(const double a0) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::ScalarFunction<double> *>(this), "deriv");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return ScalarFunction::deriv(a0);
	}
};

// ROL::LineSearch_U file:ROL_LineSearch_U.hpp line:61
struct PyCallBack_ROL_LineSearch_U_double_t : public ROL::LineSearch_U<double> {
	using ROL::LineSearch_U<double>::LineSearch_U;

	void initialize(const class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::LineSearch_U<double> *>(this), "initialize");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return LineSearch_U::initialize(a0, a1);
	}
	void run(double & a0, double & a1, int & a2, int & a3, const double & a4, const class ROL::Vector<double> & a5, const class ROL::Vector<double> & a6, class ROL::Objective<double> & a7) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::LineSearch_U<double> *>(this), "run");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4, a5, a6, a7);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		pybind11::pybind11_fail("Tried to call pure virtual function \"LineSearch_U::run\"");
	}
	bool status(const enum ROL::ELineSearchU a0, int & a1, int & a2, const double a3, const double a4, const double a5, const double a6, const class ROL::Vector<double> & a7, const class ROL::Vector<double> & a8, class ROL::Objective<double> & a9) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::LineSearch_U<double> *>(this), "status");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return LineSearch_U::status(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9);
	}
	double getInitialAlpha(int & a0, int & a1, const double a2, const double a3, const class ROL::Vector<double> & a4, const class ROL::Vector<double> & a5, class ROL::Objective<double> & a6) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::LineSearch_U<double> *>(this), "getInitialAlpha");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4, a5, a6);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return LineSearch_U::getInitialAlpha(a0, a1, a2, a3, a4, a5, a6);
	}
};

void bind_ROL_ValidParameters(std::function< pybind11::module &(std::string const &namespace_) > &M)
{
	// ROL::getValidROLParameters() file:ROL_ValidParameters.hpp line:55
	M("ROL").def("getValidROLParameters", (class Teuchos::RCP<const class Teuchos::ParameterList> (*)()) &ROL::getValidROLParameters, "C++: ROL::getValidROLParameters() --> class Teuchos::RCP<const class Teuchos::ParameterList>");

	// ROL::getValidSOLParameters() file:ROL_ValidParameters.hpp line:290
	M("ROL").def("getValidSOLParameters", (class Teuchos::RCP<const class Teuchos::ParameterList> (*)()) &ROL::getValidSOLParameters, "C++: ROL::getValidSOLParameters() --> class Teuchos::RCP<const class Teuchos::ParameterList>");

	{ // ROL::Bundle_U file:ROL_Bundle_U.hpp line:59
		pybind11::class_<ROL::Bundle_U<double>, Teuchos::RCP<ROL::Bundle_U<double>>, PyCallBack_ROL_Bundle_U_double_t> cl(M("ROL"), "Bundle_U_double_t", "", pybind11::module_local());
		cl.def( pybind11::init( [](){ return new PyCallBack_ROL_Bundle_U_double_t(); } ), "doc");
		cl.def( pybind11::init( [](const unsigned int & a0){ return new PyCallBack_ROL_Bundle_U_double_t(a0); } ), "doc");
		cl.def( pybind11::init( [](const unsigned int & a0, const double & a1){ return new PyCallBack_ROL_Bundle_U_double_t(a0, a1); } ), "doc");
		cl.def( pybind11::init( [](const unsigned int & a0, const double & a1, const double & a2){ return new PyCallBack_ROL_Bundle_U_double_t(a0, a1, a2); } ), "doc");
		cl.def( pybind11::init<const unsigned int, const double, const double, const unsigned int>(), pybind11::arg("maxSize"), pybind11::arg("coeff"), pybind11::arg("omega"), pybind11::arg("remSize") );

		cl.def(pybind11::init<PyCallBack_ROL_Bundle_U_double_t const &>());
		cl.def("initialize", (void (ROL::Bundle_U<double>::*)(const class ROL::Vector<double> &)) &ROL::Bundle_U<double>::initialize, "C++: ROL::Bundle_U<double>::initialize(const class ROL::Vector<double> &) --> void", pybind11::arg("g"));
		cl.def("solveDual", [](ROL::Bundle_U<double> &o, const double & a0) -> unsigned int { return o.solveDual(a0); }, "", pybind11::arg("t"));
		cl.def("solveDual", [](ROL::Bundle_U<double> &o, const double & a0, const unsigned int & a1) -> unsigned int { return o.solveDual(a0, a1); }, "", pybind11::arg("t"), pybind11::arg("maxit"));
		cl.def("solveDual", (unsigned int (ROL::Bundle_U<double>::*)(const double, const unsigned int, const double)) &ROL::Bundle_U<double>::solveDual, "C++: ROL::Bundle_U<double>::solveDual(const double, const unsigned int, const double) --> unsigned int", pybind11::arg("t"), pybind11::arg("maxit"), pybind11::arg("tol"));
		cl.def("linearizationError", (const double (ROL::Bundle_U<double>::*)(const unsigned int) const) &ROL::Bundle_U<double>::linearizationError, "C++: ROL::Bundle_U<double>::linearizationError(const unsigned int) const --> const double", pybind11::arg("i"));
		cl.def("distanceMeasure", (const double (ROL::Bundle_U<double>::*)(const unsigned int) const) &ROL::Bundle_U<double>::distanceMeasure, "C++: ROL::Bundle_U<double>::distanceMeasure(const unsigned int) const --> const double", pybind11::arg("i"));
		cl.def("subgradient", (const class ROL::Vector<double> & (ROL::Bundle_U<double>::*)(const unsigned int) const) &ROL::Bundle_U<double>::subgradient, "C++: ROL::Bundle_U<double>::subgradient(const unsigned int) const --> const class ROL::Vector<double> &", pybind11::return_value_policy::automatic, pybind11::arg("i"));
		cl.def("getDualVariable", (const double (ROL::Bundle_U<double>::*)(const unsigned int) const) &ROL::Bundle_U<double>::getDualVariable, "C++: ROL::Bundle_U<double>::getDualVariable(const unsigned int) const --> const double", pybind11::arg("i"));
		cl.def("setDualVariable", (void (ROL::Bundle_U<double>::*)(const unsigned int, const double)) &ROL::Bundle_U<double>::setDualVariable, "C++: ROL::Bundle_U<double>::setDualVariable(const unsigned int, const double) --> void", pybind11::arg("i"), pybind11::arg("val"));
		cl.def("resetDualVariables", (void (ROL::Bundle_U<double>::*)()) &ROL::Bundle_U<double>::resetDualVariables, "C++: ROL::Bundle_U<double>::resetDualVariables() --> void");
		cl.def("computeAlpha", (const double (ROL::Bundle_U<double>::*)(const double, const double) const) &ROL::Bundle_U<double>::computeAlpha, "C++: ROL::Bundle_U<double>::computeAlpha(const double, const double) const --> const double", pybind11::arg("dm"), pybind11::arg("le"));
		cl.def("alpha", (const double (ROL::Bundle_U<double>::*)(const unsigned int) const) &ROL::Bundle_U<double>::alpha, "C++: ROL::Bundle_U<double>::alpha(const unsigned int) const --> const double", pybind11::arg("i"));
		cl.def("size", (unsigned int (ROL::Bundle_U<double>::*)() const) &ROL::Bundle_U<double>::size, "C++: ROL::Bundle_U<double>::size() const --> unsigned int");
		cl.def("aggregate", (void (ROL::Bundle_U<double>::*)(class ROL::Vector<double> &, double &, double &) const) &ROL::Bundle_U<double>::aggregate, "C++: ROL::Bundle_U<double>::aggregate(class ROL::Vector<double> &, double &, double &) const --> void", pybind11::arg("aggSubGrad"), pybind11::arg("aggLinErr"), pybind11::arg("aggDistMeas"));
		cl.def("reset", (void (ROL::Bundle_U<double>::*)(const class ROL::Vector<double> &, const double, const double)) &ROL::Bundle_U<double>::reset, "C++: ROL::Bundle_U<double>::reset(const class ROL::Vector<double> &, const double, const double) --> void", pybind11::arg("g"), pybind11::arg("le"), pybind11::arg("dm"));
		cl.def("update", (void (ROL::Bundle_U<double>::*)(const bool, const double, const double, const class ROL::Vector<double> &, const class ROL::Vector<double> &)) &ROL::Bundle_U<double>::update, "C++: ROL::Bundle_U<double>::update(const bool, const double, const double, const class ROL::Vector<double> &, const class ROL::Vector<double> &) --> void", pybind11::arg("flag"), pybind11::arg("linErr"), pybind11::arg("distMeas"), pybind11::arg("g"), pybind11::arg("s"));
		cl.def("assign", (class ROL::Bundle_U<double> & (ROL::Bundle_U<double>::*)(const class ROL::Bundle_U<double> &)) &ROL::Bundle_U<double>::operator=, "C++: ROL::Bundle_U<double>::operator=(const class ROL::Bundle_U<double> &) --> class ROL::Bundle_U<double> &", pybind11::return_value_policy::automatic, pybind11::arg(""));
	}
	{ // ROL::ScalarFunction file:ROL_ScalarFunction.hpp line:56
		pybind11::class_<ROL::ScalarFunction<double>, Teuchos::RCP<ROL::ScalarFunction<double>>, PyCallBack_ROL_ScalarFunction_double_t> cl(M("ROL"), "ScalarFunction_double_t", "", pybind11::module_local());
		cl.def(pybind11::init<PyCallBack_ROL_ScalarFunction_double_t const &>());
		cl.def( pybind11::init( [](){ return new PyCallBack_ROL_ScalarFunction_double_t(); } ) );
		cl.def("value", (double (ROL::ScalarFunction<double>::*)(const double)) &ROL::ScalarFunction<double>::value, "C++: ROL::ScalarFunction<double>::value(const double) --> double", pybind11::arg("alpha"));
		cl.def("deriv", (double (ROL::ScalarFunction<double>::*)(const double)) &ROL::ScalarFunction<double>::deriv, "C++: ROL::ScalarFunction<double>::deriv(const double) --> double", pybind11::arg("alpha"));
		cl.def("assign", (class ROL::ScalarFunction<double> & (ROL::ScalarFunction<double>::*)(const class ROL::ScalarFunction<double> &)) &ROL::ScalarFunction<double>::operator=, "C++: ROL::ScalarFunction<double>::operator=(const class ROL::ScalarFunction<double> &) --> class ROL::ScalarFunction<double> &", pybind11::return_value_policy::automatic, pybind11::arg(""));
	}
	// ROL::EDescentU file:ROL_LineSearch_U_Types.hpp line:64
	pybind11::enum_<ROL::EDescentU>(M("ROL"), "EDescentU", pybind11::arithmetic(), "Enumeration of descent direction types.\n\n      \n    DESCENT_U_STEEPEST        describe\n      \n\n    DESCENT_U_NONLINEARCG     describe\n      \n\n    DESCENT_U_SECANT          describe\n      \n\n    DESCENT_U_NEWTON          describe \n      \n\n    DESCENT_U_NEWTONKRYLOV    describe\n      \n\n    DESCENT_U_SECANTPRECOND   describe", pybind11::module_local())
		.value("DESCENT_U_STEEPEST", ROL::DESCENT_U_STEEPEST)
		.value("DESCENT_U_NONLINEARCG", ROL::DESCENT_U_NONLINEARCG)
		.value("DESCENT_U_SECANT", ROL::DESCENT_U_SECANT)
		.value("DESCENT_U_NEWTON", ROL::DESCENT_U_NEWTON)
		.value("DESCENT_U_NEWTONKRYLOV", ROL::DESCENT_U_NEWTONKRYLOV)
		.value("DESCENT_U_USERDEFINED", ROL::DESCENT_U_USERDEFINED)
		.value("DESCENT_U_LAST", ROL::DESCENT_U_LAST)
		.export_values();

;

	// ROL::EDescentUToString(enum ROL::EDescentU) file:ROL_LineSearch_U_Types.hpp line:74
	M("ROL").def("EDescentUToString", (std::string (*)(enum ROL::EDescentU)) &ROL::EDescentUToString, "C++: ROL::EDescentUToString(enum ROL::EDescentU) --> std::string", pybind11::arg("tr"));

	// ROL::isValidDescentU(enum ROL::EDescentU) file:ROL_LineSearch_U_Types.hpp line:94
	M("ROL").def("isValidDescentU", (int (*)(enum ROL::EDescentU)) &ROL::isValidDescentU, "Verifies validity of a DescentU enum.\n\n      \n  [in]  - enum of the DescentU\n      \n\n 1 if the argument is a valid DescentU; 0 otherwise.\n\nC++: ROL::isValidDescentU(enum ROL::EDescentU) --> int", pybind11::arg("d"));

	// ROL::StringToEDescentU(std::string) file:ROL_LineSearch_U_Types.hpp line:124
	M("ROL").def("StringToEDescentU", (enum ROL::EDescentU (*)(std::string)) &ROL::StringToEDescentU, "C++: ROL::StringToEDescentU(std::string) --> enum ROL::EDescentU", pybind11::arg("s"));

	// ROL::ELineSearchU file:ROL_LineSearch_U_Types.hpp line:144
	pybind11::enum_<ROL::ELineSearchU>(M("ROL"), "ELineSearchU", pybind11::arithmetic(), "Enumeration of line-search types.\n\n      \n    LINESEARCH_U_BACKTRACKING    describe\n      \n\n    LINESEARCH_U_BISECTION       describe\n      \n\n    LINESEARCH_U_GOLDENSECTION   describe\n      \n\n    LINESEARCH_U_CUBICINTERP     describe\n      \n\n    LINESEARCH_U_BRENTS          describe\n      \n\n    LINESEARCH_U_USERDEFINED     describe", pybind11::module_local())
		.value("LINESEARCH_U_ITERATIONSCALING", ROL::LINESEARCH_U_ITERATIONSCALING)
		.value("LINESEARCH_U_PATHBASEDTARGETLEVEL", ROL::LINESEARCH_U_PATHBASEDTARGETLEVEL)
		.value("LINESEARCH_U_BACKTRACKING", ROL::LINESEARCH_U_BACKTRACKING)
		.value("LINESEARCH_U_BISECTION", ROL::LINESEARCH_U_BISECTION)
		.value("LINESEARCH_U_GOLDENSECTION", ROL::LINESEARCH_U_GOLDENSECTION)
		.value("LINESEARCH_U_CUBICINTERP", ROL::LINESEARCH_U_CUBICINTERP)
		.value("LINESEARCH_U_BRENTS", ROL::LINESEARCH_U_BRENTS)
		.value("LINESEARCH_U_USERDEFINED", ROL::LINESEARCH_U_USERDEFINED)
		.value("LINESEARCH_U_LAST", ROL::LINESEARCH_U_LAST)
		.export_values();

;

	// ROL::ELineSearchUToString(enum ROL::ELineSearchU) file:ROL_LineSearch_U_Types.hpp line:156
	M("ROL").def("ELineSearchUToString", (std::string (*)(enum ROL::ELineSearchU)) &ROL::ELineSearchUToString, "C++: ROL::ELineSearchUToString(enum ROL::ELineSearchU) --> std::string", pybind11::arg("ls"));

	// ROL::isValidLineSearchU(enum ROL::ELineSearchU) file:ROL_LineSearch_U_Types.hpp line:178
	M("ROL").def("isValidLineSearchU", (int (*)(enum ROL::ELineSearchU)) &ROL::isValidLineSearchU, "Verifies validity of a LineSearchU enum.\n\n      \n  [in]  - enum of the linesearch\n      \n\n 1 if the argument is a valid linesearch; 0 otherwise.\n\nC++: ROL::isValidLineSearchU(enum ROL::ELineSearchU) --> int", pybind11::arg("ls"));

	// ROL::StringToELineSearchU(std::string) file:ROL_LineSearch_U_Types.hpp line:210
	M("ROL").def("StringToELineSearchU", (enum ROL::ELineSearchU (*)(std::string)) &ROL::StringToELineSearchU, "C++: ROL::StringToELineSearchU(std::string) --> enum ROL::ELineSearchU", pybind11::arg("s"));

	// ROL::ECurvatureConditionU file:ROL_LineSearch_U_Types.hpp line:227
	pybind11::enum_<ROL::ECurvatureConditionU>(M("ROL"), "ECurvatureConditionU", pybind11::arithmetic(), "Enumeration of line-search curvature conditions.\n\n      \n    CURVATURECONDITION_U_WOLFE           describe\n      \n\n    CURVATURECONDITION_U_STRONGWOLFE     describe\n      \n\n    CURVATURECONDITION_U_GOLDSTEIN       describe", pybind11::module_local())
		.value("CURVATURECONDITION_U_WOLFE", ROL::CURVATURECONDITION_U_WOLFE)
		.value("CURVATURECONDITION_U_STRONGWOLFE", ROL::CURVATURECONDITION_U_STRONGWOLFE)
		.value("CURVATURECONDITION_U_GENERALIZEDWOLFE", ROL::CURVATURECONDITION_U_GENERALIZEDWOLFE)
		.value("CURVATURECONDITION_U_APPROXIMATEWOLFE", ROL::CURVATURECONDITION_U_APPROXIMATEWOLFE)
		.value("CURVATURECONDITION_U_GOLDSTEIN", ROL::CURVATURECONDITION_U_GOLDSTEIN)
		.value("CURVATURECONDITION_U_NULL", ROL::CURVATURECONDITION_U_NULL)
		.value("CURVATURECONDITION_U_LAST", ROL::CURVATURECONDITION_U_LAST)
		.export_values();

;

	// ROL::ECurvatureConditionUToString(enum ROL::ECurvatureConditionU) file:ROL_LineSearch_U_Types.hpp line:237
	M("ROL").def("ECurvatureConditionUToString", (std::string (*)(enum ROL::ECurvatureConditionU)) &ROL::ECurvatureConditionUToString, "C++: ROL::ECurvatureConditionUToString(enum ROL::ECurvatureConditionU) --> std::string", pybind11::arg("ls"));

	// ROL::isValidCurvatureConditionU(enum ROL::ECurvatureConditionU) file:ROL_LineSearch_U_Types.hpp line:257
	M("ROL").def("isValidCurvatureConditionU", (int (*)(enum ROL::ECurvatureConditionU)) &ROL::isValidCurvatureConditionU, "Verifies validity of a CurvatureConditionU enum.\n\n      \n  [in]  - enum of the Curvature Conditions\n      \n\n 1 if the argument is a valid curvature condition; 0 otherwise.\n\nC++: ROL::isValidCurvatureConditionU(enum ROL::ECurvatureConditionU) --> int", pybind11::arg("ls"));

	// ROL::StringToECurvatureConditionU(std::string) file:ROL_LineSearch_U_Types.hpp line:287
	M("ROL").def("StringToECurvatureConditionU", (enum ROL::ECurvatureConditionU (*)(std::string)) &ROL::StringToECurvatureConditionU, "C++: ROL::StringToECurvatureConditionU(std::string) --> enum ROL::ECurvatureConditionU", pybind11::arg("s"));

	{ // ROL::LineSearch_U file:ROL_LineSearch_U.hpp line:61
		pybind11::class_<ROL::LineSearch_U<double>, Teuchos::RCP<ROL::LineSearch_U<double>>, PyCallBack_ROL_LineSearch_U_double_t> cl(M("ROL"), "LineSearch_U_double_t", "", pybind11::module_local());
		cl.def( pybind11::init<class Teuchos::ParameterList &>(), pybind11::arg("parlist") );

		cl.def(pybind11::init<PyCallBack_ROL_LineSearch_U_double_t const &>());
		cl.def("initialize", (void (ROL::LineSearch_U<double>::*)(const class ROL::Vector<double> &, const class ROL::Vector<double> &)) &ROL::LineSearch_U<double>::initialize, "C++: ROL::LineSearch_U<double>::initialize(const class ROL::Vector<double> &, const class ROL::Vector<double> &) --> void", pybind11::arg("x"), pybind11::arg("g"));
		cl.def("run", (void (ROL::LineSearch_U<double>::*)(double &, double &, int &, int &, const double &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &)) &ROL::LineSearch_U<double>::run, "C++: ROL::LineSearch_U<double>::run(double &, double &, int &, int &, const double &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &) --> void", pybind11::arg("alpha"), pybind11::arg("fval"), pybind11::arg("ls_neval"), pybind11::arg("ls_ngrad"), pybind11::arg("gs"), pybind11::arg("s"), pybind11::arg("x"), pybind11::arg("obj"));
		cl.def("setMaxitUpdate", (void (ROL::LineSearch_U<double>::*)(double &, double &, const double &)) &ROL::LineSearch_U<double>::setMaxitUpdate, "C++: ROL::LineSearch_U<double>::setMaxitUpdate(double &, double &, const double &) --> void", pybind11::arg("alpha"), pybind11::arg("fnew"), pybind11::arg("fold"));
		cl.def("assign", (class ROL::LineSearch_U<double> & (ROL::LineSearch_U<double>::*)(const class ROL::LineSearch_U<double> &)) &ROL::LineSearch_U<double>::operator=, "C++: ROL::LineSearch_U<double>::operator=(const class ROL::LineSearch_U<double> &) --> class ROL::LineSearch_U<double> &", pybind11::return_value_policy::automatic, pybind11::arg(""));
	}
}
