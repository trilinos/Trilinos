#include <ROL_BoundConstraint.hpp>
#include <ROL_Elementwise_Function.hpp>
#include <ROL_Elementwise_Reduce.hpp>
#include <ROL_PolyhedralProjection.hpp>
#include <ROL_Problem.hpp>
#include <ROL_StatusTest.hpp>
#include <ROL_TypeU_Algorithm.hpp>
#include <ROL_Types.hpp>
#include <ROL_UpdateType.hpp>
#include <Teuchos_ENull.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCPDecl.hpp>
#include <Teuchos_RCPNode.hpp>
#include <cwchar>
#include <ios>
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

void bind_ROL_TypeU_Algorithm(std::function< pybind11::module &(std::string const &namespace_) > &M)
{
	{ // ROL::TypeU::AlgorithmState file:ROL_TypeU_Algorithm.hpp line:60
		pybind11::class_<ROL::TypeU::AlgorithmState<double>, Teuchos::RCP<ROL::TypeU::AlgorithmState<double>>, ROL::AlgorithmState<double>> cl(M("ROL::TypeU"), "AlgorithmState_double_t", "");
		cl.def( pybind11::init( [](){ return new ROL::TypeU::AlgorithmState<double>(); } ) );
		cl.def( pybind11::init( [](ROL::TypeU::AlgorithmState<double> const &o){ return new ROL::TypeU::AlgorithmState<double>(o); } ) );
		cl.def_readwrite("searchSize", &ROL::TypeU::AlgorithmState<double>::searchSize);
		cl.def_readwrite("stepVec", &ROL::TypeU::AlgorithmState<double>::stepVec);
		cl.def_readwrite("gradientVec", &ROL::TypeU::AlgorithmState<double>::gradientVec);
		cl.def("reset", (void (ROL::TypeU::AlgorithmState<double>::*)()) &ROL::TypeU::AlgorithmState<double>::reset, "C++: ROL::TypeU::AlgorithmState<double>::reset() --> void");
		cl.def("assign", (struct ROL::TypeU::AlgorithmState<double> & (ROL::TypeU::AlgorithmState<double>::*)(const struct ROL::TypeU::AlgorithmState<double> &)) &ROL::TypeU::AlgorithmState<double>::operator=, "C++: ROL::TypeU::AlgorithmState<double>::operator=(const struct ROL::TypeU::AlgorithmState<double> &) --> struct ROL::TypeU::AlgorithmState<double> &", pybind11::return_value_policy::automatic, pybind11::arg(""));
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
	{ // ROL::TypeU::Algorithm file:ROL_TypeU_Algorithm.hpp line:83
		pybind11::class_<ROL::TypeU::Algorithm<double>, Teuchos::RCP<ROL::TypeU::Algorithm<double>>> cl(M("ROL::TypeU"), "Algorithm_double_t", "");
		cl.def("setStatusTest", [](ROL::TypeU::Algorithm<double> &o, const class Teuchos::RCP<class ROL::StatusTest<double> > & a0) -> void { return o.setStatusTest(a0); }, "", pybind11::arg("status"));
		cl.def("setStatusTest", (void (ROL::TypeU::Algorithm<double>::*)(const class Teuchos::RCP<class ROL::StatusTest<double> > &, bool)) &ROL::TypeU::Algorithm<double>::setStatusTest, "C++: ROL::TypeU::Algorithm<double>::setStatusTest(const class Teuchos::RCP<class ROL::StatusTest<double> > &, bool) --> void", pybind11::arg("status"), pybind11::arg("combineStatus"));
		cl.def("getState", (class Teuchos::RCP<const struct ROL::TypeU::AlgorithmState<double> > (ROL::TypeU::Algorithm<double>::*)() const) &ROL::TypeU::Algorithm<double>::getState, "C++: ROL::TypeU::Algorithm<double>::getState() const --> class Teuchos::RCP<const struct ROL::TypeU::AlgorithmState<double> >");
		cl.def("reset", (void (ROL::TypeU::Algorithm<double>::*)()) &ROL::TypeU::Algorithm<double>::reset, "C++: ROL::TypeU::Algorithm<double>::reset() --> void");
	}
}
