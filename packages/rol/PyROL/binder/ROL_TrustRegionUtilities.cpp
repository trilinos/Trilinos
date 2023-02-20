#include <ROL_Elementwise_Function.hpp>
#include <ROL_Elementwise_Reduce.hpp>
#include <ROL_Objective.hpp>
#include <ROL_Secant.hpp>
#include <ROL_TrustRegionModel_U.hpp>
#include <ROL_TrustRegionUtilities.hpp>
#include <ROL_TrustRegion_U_Types.hpp>
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

void bind_ROL_TrustRegionUtilities(std::function< pybind11::module &(std::string const &namespace_) > &M)
{
	// ROL::TRUtils::ETRFlag file:ROL_TrustRegionUtilities.hpp line:62
	pybind11::enum_<ROL::TRUtils::ETRFlag>(M("ROL::TRUtils"), "ETRFlag", pybind11::arithmetic(), "Enumation of flags used by trust-region solvers.\n\n    \n SUCCESS        Actual and predicted reductions are positive \n    \n\n POSPREDNEG     Reduction is positive, predicted negative (impossible)\n    \n\n NPOSPREDPOS    Reduction is nonpositive, predicted positive\n    \n\n NPOSPREDNEG    Reduction is nonpositive, predicted negative (impossible)\n    \n\n TRNAN          Actual and/or predicted reduction is NaN", pybind11::module_local())
		.value("SUCCESS", ROL::TRUtils::SUCCESS)
		.value("POSPREDNEG", ROL::TRUtils::POSPREDNEG)
		.value("NPOSPREDPOS", ROL::TRUtils::NPOSPREDPOS)
		.value("NPOSPREDNEG", ROL::TRUtils::NPOSPREDNEG)
		.value("TRNAN", ROL::TRUtils::TRNAN)
		.value("QMINSUFDEC", ROL::TRUtils::QMINSUFDEC)
		.value("UNDEFINED", ROL::TRUtils::UNDEFINED)
		.export_values();

;

	// ROL::TRUtils::ETRFlagToString(enum ROL::TRUtils::ETRFlag) file:ROL_TrustRegionUtilities.hpp line:72
	M("ROL::TRUtils").def("ETRFlagToString", (std::string (*)(enum ROL::TRUtils::ETRFlag)) &ROL::TRUtils::ETRFlagToString, "C++: ROL::TRUtils::ETRFlagToString(enum ROL::TRUtils::ETRFlag) --> std::string", pybind11::arg("trf"));

	// ROL::TRUtils::initialRadius(int &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Vector<double> &, const double, const double, class ROL::Objective<double> &, class ROL::TrustRegionModel_U<double> &, const double, std::ostream &, const bool) file:ROL_TrustRegionUtilities.hpp line:100
	M("ROL::TRUtils").def("initialRadius", [](int & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2, class ROL::Vector<double> & a3, const double & a4, const double & a5, class ROL::Objective<double> & a6, class ROL::TrustRegionModel_U<double> & a7, const double & a8, std::ostream & a9) -> double { return ROL::TRUtils::initialRadius(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9); }, "", pybind11::arg("nfval"), pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("Bg"), pybind11::arg("fx"), pybind11::arg("gnorm"), pybind11::arg("obj"), pybind11::arg("model"), pybind11::arg("delMax"), pybind11::arg("outStream"));
	M("ROL::TRUtils").def("initialRadius", (double (*)(int &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Vector<double> &, const double, const double, class ROL::Objective<double> &, class ROL::TrustRegionModel_U<double> &, const double, std::ostream &, const bool)) &ROL::TRUtils::initialRadius<double>, "C++: ROL::TRUtils::initialRadius(int &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Vector<double> &, const double, const double, class ROL::Objective<double> &, class ROL::TrustRegionModel_U<double> &, const double, std::ostream &, const bool) --> double", pybind11::arg("nfval"), pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("Bg"), pybind11::arg("fx"), pybind11::arg("gnorm"), pybind11::arg("obj"), pybind11::arg("model"), pybind11::arg("delMax"), pybind11::arg("outStream"), pybind11::arg("print"));

	// ROL::TRUtils::analyzeRatio(double &, enum ROL::TRUtils::ETRFlag &, const double, const double, const double, const double, std::ostream &, const bool) file:ROL_TrustRegionUtilities.hpp line:171
	M("ROL::TRUtils").def("analyzeRatio", [](double & a0, enum ROL::TRUtils::ETRFlag & a1, const double & a2, const double & a3, const double & a4, const double & a5) -> void { return ROL::TRUtils::analyzeRatio(a0, a1, a2, a3, a4, a5); }, "", pybind11::arg("rho"), pybind11::arg("flag"), pybind11::arg("fold"), pybind11::arg("ftrial"), pybind11::arg("pRed"), pybind11::arg("epsi"));
	M("ROL::TRUtils").def("analyzeRatio", [](double & a0, enum ROL::TRUtils::ETRFlag & a1, const double & a2, const double & a3, const double & a4, const double & a5, std::ostream & a6) -> void { return ROL::TRUtils::analyzeRatio(a0, a1, a2, a3, a4, a5, a6); }, "", pybind11::arg("rho"), pybind11::arg("flag"), pybind11::arg("fold"), pybind11::arg("ftrial"), pybind11::arg("pRed"), pybind11::arg("epsi"), pybind11::arg("outStream"));
	M("ROL::TRUtils").def("analyzeRatio", (void (*)(double &, enum ROL::TRUtils::ETRFlag &, const double, const double, const double, const double, std::ostream &, const bool)) &ROL::TRUtils::analyzeRatio<double>, "C++: ROL::TRUtils::analyzeRatio(double &, enum ROL::TRUtils::ETRFlag &, const double, const double, const double, const double, std::ostream &, const bool) --> void", pybind11::arg("rho"), pybind11::arg("flag"), pybind11::arg("fold"), pybind11::arg("ftrial"), pybind11::arg("pRed"), pybind11::arg("epsi"), pybind11::arg("outStream"), pybind11::arg("print"));

	// ROL::TRUtils::interpolateRadius(const class ROL::Vector<double> &, const class ROL::Vector<double> &, const double, const double, const double, const double, const double, const double, const double, const double, std::ostream &, const bool) file:ROL_TrustRegionUtilities.hpp line:222
	M("ROL::TRUtils").def("interpolateRadius", [](const class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, const double & a2, const double & a3, const double & a4, const double & a5, const double & a6, const double & a7, const double & a8, const double & a9) -> double { return ROL::TRUtils::interpolateRadius(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9); }, "", pybind11::arg("g"), pybind11::arg("s"), pybind11::arg("snorm"), pybind11::arg("pRed"), pybind11::arg("fold"), pybind11::arg("ftrial"), pybind11::arg("del"), pybind11::arg("gamma0"), pybind11::arg("gamma1"), pybind11::arg("eta2"));
	M("ROL::TRUtils").def("interpolateRadius", [](const class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, const double & a2, const double & a3, const double & a4, const double & a5, const double & a6, const double & a7, const double & a8, const double & a9, std::ostream & a10) -> double { return ROL::TRUtils::interpolateRadius(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10); }, "", pybind11::arg("g"), pybind11::arg("s"), pybind11::arg("snorm"), pybind11::arg("pRed"), pybind11::arg("fold"), pybind11::arg("ftrial"), pybind11::arg("del"), pybind11::arg("gamma0"), pybind11::arg("gamma1"), pybind11::arg("eta2"), pybind11::arg("outStream"));
	M("ROL::TRUtils").def("interpolateRadius", (double (*)(const class ROL::Vector<double> &, const class ROL::Vector<double> &, const double, const double, const double, const double, const double, const double, const double, const double, std::ostream &, const bool)) &ROL::TRUtils::interpolateRadius<double>, "C++: ROL::TRUtils::interpolateRadius(const class ROL::Vector<double> &, const class ROL::Vector<double> &, const double, const double, const double, const double, const double, const double, const double, const double, std::ostream &, const bool) --> double", pybind11::arg("g"), pybind11::arg("s"), pybind11::arg("snorm"), pybind11::arg("pRed"), pybind11::arg("fold"), pybind11::arg("ftrial"), pybind11::arg("del"), pybind11::arg("gamma0"), pybind11::arg("gamma1"), pybind11::arg("eta2"), pybind11::arg("outStream"), pybind11::arg("print"));

}
