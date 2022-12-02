#include <ROL_TypeU_Algorithm.hpp>
#include <sstream> // __str__

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

void bind_ROL_TypeU_Algorithm_1(std::function< pybind11::module &(std::string const &namespace_) > &M)
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
}
