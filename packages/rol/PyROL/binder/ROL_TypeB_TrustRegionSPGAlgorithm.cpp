#include <ROL_BoundConstraint.hpp>
#include <ROL_Elementwise_Function.hpp>
#include <ROL_Elementwise_Reduce.hpp>
#include <ROL_Objective.hpp>
#include <ROL_PartitionedVector.hpp>
#include <ROL_Secant.hpp>
#include <ROL_StatusTest.hpp>
#include <ROL_TypeB_Algorithm.hpp>
#include <ROL_TypeB_AlgorithmFactory.hpp>
#include <ROL_TypeB_ColemanLiAlgorithm.hpp>
#include <ROL_TypeB_GradientAlgorithm.hpp>
#include <ROL_TypeB_InteriorPointAlgorithm.hpp>
#include <ROL_TypeB_LinMoreAlgorithm.hpp>
#include <ROL_TypeB_MoreauYosidaAlgorithm.hpp>
#include <ROL_TypeB_PrimalDualActiveSetAlgorithm.hpp>
#include <ROL_TypeB_SpectralGradientAlgorithm.hpp>
#include <ROL_TypeB_TrustRegionSPGAlgorithm.hpp>
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


#ifndef BINDER_PYBIND11_TYPE_CASTER
	#define BINDER_PYBIND11_TYPE_CASTER
	PYBIND11_DECLARE_HOLDER_TYPE(T, std::shared_ptr<T>)
	PYBIND11_DECLARE_HOLDER_TYPE(T, T*)
	PYBIND11_MAKE_OPAQUE(std::shared_ptr<void>)
#endif

// ROL::TypeB::TrustRegionSPGAlgorithm file:ROL_TypeB_TrustRegionSPGAlgorithm.hpp line:59
struct PyCallBack_ROL_TypeB_TrustRegionSPGAlgorithm_double_t : public ROL::TypeB::TrustRegionSPGAlgorithm<double> {
	using ROL::TypeB::TrustRegionSPGAlgorithm<double>::TrustRegionSPGAlgorithm;

};

// ROL::TypeB::ColemanLiAlgorithm file:ROL_TypeB_ColemanLiAlgorithm.hpp line:61
struct PyCallBack_ROL_TypeB_ColemanLiAlgorithm_double_t : public ROL::TypeB::ColemanLiAlgorithm<double> {
	using ROL::TypeB::ColemanLiAlgorithm<double>::ColemanLiAlgorithm;

};

void bind_ROL_TypeB_TrustRegionSPGAlgorithm(std::function< pybind11::module &(std::string const &namespace_) > &M)
{
	{ // ROL::TypeB::TrustRegionSPGAlgorithm file:ROL_TypeB_TrustRegionSPGAlgorithm.hpp line:59
		pybind11::class_<ROL::TypeB::TrustRegionSPGAlgorithm<double>, std::shared_ptr<ROL::TypeB::TrustRegionSPGAlgorithm<double>>, PyCallBack_ROL_TypeB_TrustRegionSPGAlgorithm_double_t, ROL::TypeB::Algorithm<double>> cl(M("ROL::TypeB"), "TrustRegionSPGAlgorithm_double_t", "", pybind11::module_local());
		cl.def( pybind11::init( [](class Teuchos::ParameterList & a0){ return new ROL::TypeB::TrustRegionSPGAlgorithm<double>(a0); }, [](class Teuchos::ParameterList & a0){ return new PyCallBack_ROL_TypeB_TrustRegionSPGAlgorithm_double_t(a0); } ), "doc");
		cl.def( pybind11::init<class Teuchos::ParameterList &, const class std::shared_ptr<class ROL::Secant<double> > &>(), pybind11::arg("list"), pybind11::arg("secant") );

		cl.def( pybind11::init( [](PyCallBack_ROL_TypeB_TrustRegionSPGAlgorithm_double_t const &o){ return new PyCallBack_ROL_TypeB_TrustRegionSPGAlgorithm_double_t(o); } ) );
		cl.def( pybind11::init( [](ROL::TypeB::TrustRegionSPGAlgorithm<double> const &o){ return new ROL::TypeB::TrustRegionSPGAlgorithm<double>(o); } ) );
		cl.def("setStatusTest", [](ROL::TypeB::Algorithm<double> &o, const class std::shared_ptr<class ROL::StatusTest<double> > & a0) -> void { return o.setStatusTest(a0); }, "", pybind11::arg("status"));
		cl.def("setStatusTest", (void (ROL::TypeB::Algorithm<double>::*)(const class std::shared_ptr<class ROL::StatusTest<double> > &, const bool)) &ROL::TypeB::Algorithm<double>::setStatusTest, "C++: ROL::TypeB::Algorithm<double>::setStatusTest(const class std::shared_ptr<class ROL::StatusTest<double> > &, const bool) --> void", pybind11::arg("status"), pybind11::arg("combineStatus"));
		cl.def("getState", (class std::shared_ptr<const struct ROL::TypeB::AlgorithmState<double> > (ROL::TypeB::Algorithm<double>::*)() const) &ROL::TypeB::Algorithm<double>::getState, "C++: ROL::TypeB::Algorithm<double>::getState() const --> class std::shared_ptr<const struct ROL::TypeB::AlgorithmState<double> >");
		cl.def("reset", (void (ROL::TypeB::Algorithm<double>::*)()) &ROL::TypeB::Algorithm<double>::reset, "C++: ROL::TypeB::Algorithm<double>::reset() --> void");
	}
	{ // ROL::TypeB::ColemanLiAlgorithm file:ROL_TypeB_ColemanLiAlgorithm.hpp line:61
		pybind11::class_<ROL::TypeB::ColemanLiAlgorithm<double>, std::shared_ptr<ROL::TypeB::ColemanLiAlgorithm<double>>, PyCallBack_ROL_TypeB_ColemanLiAlgorithm_double_t, ROL::TypeB::Algorithm<double>> cl(M("ROL::TypeB"), "ColemanLiAlgorithm_double_t", "", pybind11::module_local());
		cl.def( pybind11::init( [](class Teuchos::ParameterList & a0){ return new ROL::TypeB::ColemanLiAlgorithm<double>(a0); }, [](class Teuchos::ParameterList & a0){ return new PyCallBack_ROL_TypeB_ColemanLiAlgorithm_double_t(a0); } ), "doc");
		cl.def( pybind11::init<class Teuchos::ParameterList &, const class std::shared_ptr<class ROL::Secant<double> > &>(), pybind11::arg("list"), pybind11::arg("secant") );

		cl.def( pybind11::init( [](PyCallBack_ROL_TypeB_ColemanLiAlgorithm_double_t const &o){ return new PyCallBack_ROL_TypeB_ColemanLiAlgorithm_double_t(o); } ) );
		cl.def( pybind11::init( [](ROL::TypeB::ColemanLiAlgorithm<double> const &o){ return new ROL::TypeB::ColemanLiAlgorithm<double>(o); } ) );
		cl.def("setStatusTest", [](ROL::TypeB::Algorithm<double> &o, const class std::shared_ptr<class ROL::StatusTest<double> > & a0) -> void { return o.setStatusTest(a0); }, "", pybind11::arg("status"));
		cl.def("setStatusTest", (void (ROL::TypeB::Algorithm<double>::*)(const class std::shared_ptr<class ROL::StatusTest<double> > &, const bool)) &ROL::TypeB::Algorithm<double>::setStatusTest, "C++: ROL::TypeB::Algorithm<double>::setStatusTest(const class std::shared_ptr<class ROL::StatusTest<double> > &, const bool) --> void", pybind11::arg("status"), pybind11::arg("combineStatus"));
		cl.def("getState", (class std::shared_ptr<const struct ROL::TypeB::AlgorithmState<double> > (ROL::TypeB::Algorithm<double>::*)() const) &ROL::TypeB::Algorithm<double>::getState, "C++: ROL::TypeB::Algorithm<double>::getState() const --> class std::shared_ptr<const struct ROL::TypeB::AlgorithmState<double> >");
		cl.def("reset", (void (ROL::TypeB::Algorithm<double>::*)()) &ROL::TypeB::Algorithm<double>::reset, "C++: ROL::TypeB::Algorithm<double>::reset() --> void");
	}
	// ROL::TypeB::EAlgorithmB file:ROL_TypeB_AlgorithmFactory.hpp line:74
	pybind11::enum_<ROL::TypeB::EAlgorithmB>(M("ROL::TypeB"), "EAlgorithmB", pybind11::arithmetic(), "Enumeration of bound constrained algorithm types.\n\n    \n    ALGORITHM_B_LINESEARCH          describe\n    \n\n    ALGORITHM_B_TRUSTREGION         describe\n    \n\n    ALGORITHM_B_MOREAUYOSIDA        describe\n    \n\n    ALGORITHM_B_PRIMALDUALACTIVESET describe\n    \n\n    ALGORITHM_B_INTERIORPOINT       describe\n    \n\n    ALGORITHM_B_SPECTRALGRADIENT    describe", pybind11::module_local())
		.value("ALGORITHM_B_LINESEARCH", ROL::TypeB::ALGORITHM_B_LINESEARCH)
		.value("ALGORITHM_B_TRUSTREGION", ROL::TypeB::ALGORITHM_B_TRUSTREGION)
		.value("ALGORITHM_B_MOREAUYOSIDA", ROL::TypeB::ALGORITHM_B_MOREAUYOSIDA)
		.value("ALGORITHM_B_PRIMALDUALACTIVESET", ROL::TypeB::ALGORITHM_B_PRIMALDUALACTIVESET)
		.value("ALGORITHM_B_INTERIORPOINT", ROL::TypeB::ALGORITHM_B_INTERIORPOINT)
		.value("ALGORITHM_B_SPECTRALGRADIENT", ROL::TypeB::ALGORITHM_B_SPECTRALGRADIENT)
		.value("ALGORITHM_B_LAST", ROL::TypeB::ALGORITHM_B_LAST)
		.export_values();

;

	// ROL::TypeB::EAlgorithmBToString(enum ROL::TypeB::EAlgorithmB) file:ROL_TypeB_AlgorithmFactory.hpp line:84
	M("ROL::TypeB").def("EAlgorithmBToString", (std::string (*)(enum ROL::TypeB::EAlgorithmB)) &ROL::TypeB::EAlgorithmBToString, "C++: ROL::TypeB::EAlgorithmBToString(enum ROL::TypeB::EAlgorithmB) --> std::string", pybind11::arg("alg"));

	// ROL::TypeB::isValidAlgorithmB(enum ROL::TypeB::EAlgorithmB) file:ROL_TypeB_AlgorithmFactory.hpp line:104
	M("ROL::TypeB").def("isValidAlgorithmB", (int (*)(enum ROL::TypeB::EAlgorithmB)) &ROL::TypeB::isValidAlgorithmB, "Verifies validity of a AlgorithmB enum.\n\n    \n  [in]  - enum of the AlgorithmB\n    \n\n 1 if the argument is a valid AlgorithmB; 0 otherwise.\n\nC++: ROL::TypeB::isValidAlgorithmB(enum ROL::TypeB::EAlgorithmB) --> int", pybind11::arg("alg"));

	// ROL::TypeB::StringToEAlgorithmB(std::string) file:ROL_TypeB_AlgorithmFactory.hpp line:135
	M("ROL::TypeB").def("StringToEAlgorithmB", (enum ROL::TypeB::EAlgorithmB (*)(std::string)) &ROL::TypeB::StringToEAlgorithmB, "C++: ROL::TypeB::StringToEAlgorithmB(std::string) --> enum ROL::TypeB::EAlgorithmB", pybind11::arg("s"));

	// ROL::TypeB::AlgorithmFactory(class Teuchos::ParameterList &) file:ROL_TypeB_AlgorithmFactory.hpp line:146
	M("ROL::TypeB").def("AlgorithmFactory", (class std::shared_ptr<class ROL::TypeB::Algorithm<double> > (*)(class Teuchos::ParameterList &)) &ROL::TypeB::AlgorithmFactory<double>, "C++: ROL::TypeB::AlgorithmFactory(class Teuchos::ParameterList &) --> class std::shared_ptr<class ROL::TypeB::Algorithm<double> >", pybind11::arg("parlist"));

}
