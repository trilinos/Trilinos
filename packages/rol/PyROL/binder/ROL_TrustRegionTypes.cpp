#include <ROL_TrustRegionTypes.hpp>
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

void bind_ROL_TrustRegionTypes(std::function< pybind11::module &(std::string const &namespace_) > &M)
{
	// ROL::ETrustRegion file:ROL_TrustRegionTypes.hpp line:64
	pybind11::enum_<ROL::ETrustRegion>(M("ROL"), "ETrustRegion", pybind11::arithmetic(), "Enumeration of trust-region solver types.\n\n      \n    CAUCHYPOINT     describe\n      \n\n    TRUNCATEDCG     describe\n      \n\n    DOGLEG          describe\n      \n\n    DOUBLEDOGLEG    describe")
		.value("TRUSTREGION_CAUCHYPOINT", ROL::TRUSTREGION_CAUCHYPOINT)
		.value("TRUSTREGION_TRUNCATEDCG", ROL::TRUSTREGION_TRUNCATEDCG)
		.value("TRUSTREGION_DOGLEG", ROL::TRUSTREGION_DOGLEG)
		.value("TRUSTREGION_DOUBLEDOGLEG", ROL::TRUSTREGION_DOUBLEDOGLEG)
		.value("TRUSTREGION_LINMORE", ROL::TRUSTREGION_LINMORE)
		.value("TRUSTREGION_LAST", ROL::TRUSTREGION_LAST)
		.export_values();

;

	// ROL::ETrustRegionToString(enum ROL::ETrustRegion) file:ROL_TrustRegionTypes.hpp line:73
	M("ROL").def("ETrustRegionToString", (std::string (*)(enum ROL::ETrustRegion)) &ROL::ETrustRegionToString, "C++: ROL::ETrustRegionToString(enum ROL::ETrustRegion) --> std::string", pybind11::arg("tr"));

	// ROL::isValidTrustRegion(enum ROL::ETrustRegion) file:ROL_TrustRegionTypes.hpp line:92
	M("ROL").def("isValidTrustRegion", (int (*)(enum ROL::ETrustRegion)) &ROL::isValidTrustRegion, "Verifies validity of a TrustRegion enum.\n\n      \n  [in]  - enum of the TrustRegion\n      \n\n 1 if the argument is a valid TrustRegion; 0 otherwise.\n\nC++: ROL::isValidTrustRegion(enum ROL::ETrustRegion) --> int", pybind11::arg("ls"));

	// ROL::StringToETrustRegion(std::string) file:ROL_TrustRegionTypes.hpp line:121
	M("ROL").def("StringToETrustRegion", (enum ROL::ETrustRegion (*)(std::string)) &ROL::StringToETrustRegion, "C++: ROL::StringToETrustRegion(std::string) --> enum ROL::ETrustRegion", pybind11::arg("s"));

	// ROL::ETrustRegionModel file:ROL_TrustRegionTypes.hpp line:137
	pybind11::enum_<ROL::ETrustRegionModel>(M("ROL"), "ETrustRegionModel", pybind11::arithmetic(), "Enumeration of trust-region model types.\n\n      \n    COLEMANLI   describe\n      \n\n    KELLEYSACHS describe")
		.value("TRUSTREGION_MODEL_COLEMANLI", ROL::TRUSTREGION_MODEL_COLEMANLI)
		.value("TRUSTREGION_MODEL_KELLEYSACHS", ROL::TRUSTREGION_MODEL_KELLEYSACHS)
		.value("TRUSTREGION_MODEL_LINMORE", ROL::TRUSTREGION_MODEL_LINMORE)
		.value("TRUSTREGION_MODEL_LAST", ROL::TRUSTREGION_MODEL_LAST)
		.export_values();

;

	// ROL::ETrustRegionModelToString(enum ROL::ETrustRegionModel) file:ROL_TrustRegionTypes.hpp line:144
	M("ROL").def("ETrustRegionModelToString", (std::string (*)(enum ROL::ETrustRegionModel)) &ROL::ETrustRegionModelToString, "C++: ROL::ETrustRegionModelToString(enum ROL::ETrustRegionModel) --> std::string", pybind11::arg("tr"));

	// ROL::isValidTrustRegionModel(enum ROL::ETrustRegionModel) file:ROL_TrustRegionTypes.hpp line:161
	M("ROL").def("isValidTrustRegionModel", (int (*)(enum ROL::ETrustRegionModel)) &ROL::isValidTrustRegionModel, "Verifies validity of a TrustRegionModel enum.\n\n      \n  [in]  - enum of the TrustRegionModel\n      \n\n 1 if the argument is a valid TrustRegionModel; 0 otherwise.\n\nC++: ROL::isValidTrustRegionModel(enum ROL::ETrustRegionModel) --> int", pybind11::arg("ls"));

	// ROL::StringToETrustRegionModel(std::string) file:ROL_TrustRegionTypes.hpp line:188
	M("ROL").def("StringToETrustRegionModel", (enum ROL::ETrustRegionModel (*)(std::string)) &ROL::StringToETrustRegionModel, "C++: ROL::StringToETrustRegionModel(std::string) --> enum ROL::ETrustRegionModel", pybind11::arg("s"));

	// ROL::isValidTrustRegionSubproblem(enum ROL::ETrustRegion, enum ROL::ETrustRegionModel, bool) file:ROL_TrustRegionTypes.hpp line:198
	M("ROL").def("isValidTrustRegionSubproblem", (bool (*)(enum ROL::ETrustRegion, enum ROL::ETrustRegionModel, bool)) &ROL::isValidTrustRegionSubproblem, "C++: ROL::isValidTrustRegionSubproblem(enum ROL::ETrustRegion, enum ROL::ETrustRegionModel, bool) --> bool", pybind11::arg("etr"), pybind11::arg("etrm"), pybind11::arg("isBnd"));

	// ROL::ETrustRegionFlag file:ROL_TrustRegionTypes.hpp line:232
	pybind11::enum_<ROL::ETrustRegionFlag>(M("ROL"), "ETrustRegionFlag", pybind11::arithmetic(), "Enumation of flags used by trust-region solvers.\n\n      \n TRUSTREGION_FLAG_SUCCESS        Actual and predicted reductions are positive \n      \n\n TRUSTREGION_FLAG_POSPREDNEG     Reduction is positive, predicted negative (impossible)\n      \n\n TRUSTREGION_FLAG_NPOSPREDPOS    Reduction is nonpositive, predicted positive\n      \n\n TRUSTREGION_FLAG_NPOSPREDNEG    Reduction is nonpositive, predicted negative (impossible)\n      \n\n TRUSTREGION_FLAG_QMINSUFDEC     Insufficient decrease of the quadratic model (bound constraint only)\n      \n\n TRUSTREGION_FLAG_NAN            Actual and/or predicted reduction is NaN\n\n  ")
		.value("TRUSTREGION_FLAG_SUCCESS", ROL::TRUSTREGION_FLAG_SUCCESS)
		.value("TRUSTREGION_FLAG_POSPREDNEG", ROL::TRUSTREGION_FLAG_POSPREDNEG)
		.value("TRUSTREGION_FLAG_NPOSPREDPOS", ROL::TRUSTREGION_FLAG_NPOSPREDPOS)
		.value("TRUSTREGION_FLAG_NPOSPREDNEG", ROL::TRUSTREGION_FLAG_NPOSPREDNEG)
		.value("TRUSTREGION_FLAG_QMINSUFDEC", ROL::TRUSTREGION_FLAG_QMINSUFDEC)
		.value("TRUSTREGION_FLAG_NAN", ROL::TRUSTREGION_FLAG_NAN)
		.value("TRUSTREGION_FLAG_UNDEFINED", ROL::TRUSTREGION_FLAG_UNDEFINED)
		.export_values();

;

	// ROL::ETrustRegionFlagToString(enum ROL::ETrustRegionFlag) file:ROL_TrustRegionTypes.hpp line:243
	M("ROL").def("ETrustRegionFlagToString", (std::string (*)(enum ROL::ETrustRegionFlag)) &ROL::ETrustRegionFlagToString, "C++: ROL::ETrustRegionFlagToString(enum ROL::ETrustRegionFlag) --> std::string", pybind11::arg("trf"));

}
