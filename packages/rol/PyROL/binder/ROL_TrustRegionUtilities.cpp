#include <ROL_TrustRegionUtilities.hpp>
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

}
