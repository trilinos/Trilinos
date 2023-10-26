#include <ROL_TypeU_AlgorithmFactory.hpp>
#include <iterator>
#include <memory>
#include <string>

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

void bind_ROL_TypeU_AlgorithmFactory(std::function< pybind11::module &(std::string const &namespace_) > &M)
{
	// ROL::TypeU::EAlgorithmU file:ROL_TypeU_AlgorithmFactory.hpp line:62
	pybind11::enum_<ROL::TypeU::EAlgorithmU>(M("ROL::TypeU"), "EAlgorithmU", pybind11::arithmetic(), "Enumeration of unconstrained algorithm types.\n\n    \n    ALGORITHM_U_BUNDLE         describe\n    \n\n    ALGORITHM_U_LINESEARCH     describe\n    \n\n    ALGORITHM_U_TRUSTREGION    describe", pybind11::module_local())
		.value("ALGORITHM_U_BUNDLE", ROL::TypeU::ALGORITHM_U_BUNDLE)
		.value("ALGORITHM_U_LINESEARCH", ROL::TypeU::ALGORITHM_U_LINESEARCH)
		.value("ALGORITHM_U_TRUSTREGION", ROL::TypeU::ALGORITHM_U_TRUSTREGION)
		.value("ALGORITHM_U_LAST", ROL::TypeU::ALGORITHM_U_LAST)
		.export_values();

;

	// ROL::TypeU::EAlgorithmUToString(enum ROL::TypeU::EAlgorithmU) file:ROL_TypeU_AlgorithmFactory.hpp line:69
	M("ROL::TypeU").def("EAlgorithmUToString", (std::string (*)(enum ROL::TypeU::EAlgorithmU)) &ROL::TypeU::EAlgorithmUToString, "C++: ROL::TypeU::EAlgorithmUToString(enum ROL::TypeU::EAlgorithmU) --> std::string", pybind11::arg("alg"));

	// ROL::TypeU::isValidAlgorithmU(enum ROL::TypeU::EAlgorithmU) file:ROL_TypeU_AlgorithmFactory.hpp line:86
	M("ROL::TypeU").def("isValidAlgorithmU", (int (*)(enum ROL::TypeU::EAlgorithmU)) &ROL::TypeU::isValidAlgorithmU, "Verifies validity of a AlgorithmU enum.\n\n    \n  [in]  - enum of the AlgorithmU\n    \n\n 1 if the argument is a valid AlgorithmU; 0 otherwise.\n\nC++: ROL::TypeU::isValidAlgorithmU(enum ROL::TypeU::EAlgorithmU) --> int", pybind11::arg("alg"));

	// ROL::TypeU::StringToEAlgorithmU(std::string) file:ROL_TypeU_AlgorithmFactory.hpp line:114
	M("ROL::TypeU").def("StringToEAlgorithmU", (enum ROL::TypeU::EAlgorithmU (*)(std::string)) &ROL::TypeU::StringToEAlgorithmU, "C++: ROL::TypeU::StringToEAlgorithmU(std::string) --> enum ROL::TypeU::EAlgorithmU", pybind11::arg("s"));

}
