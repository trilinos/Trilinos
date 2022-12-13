#include <map>
#include <algorithm>
#include <functional>
#include <memory>
#include <stdexcept>
#include <string>

#include <pybind11/pybind11.h>

typedef std::function< pybind11::module & (std::string const &) > ModuleGetter;

void bind_std_postypes(std::function< pybind11::module &(std::string const &namespace_) > &M);
void bind_std_typeinfo(std::function< pybind11::module &(std::string const &namespace_) > &M);
void bind_std_locale_classes(std::function< pybind11::module &(std::string const &namespace_) > &M);
void bind_ROL_Elementwise_Function(std::function< pybind11::module &(std::string const &namespace_) > &M);
void bind_std_istream_tcc(std::function< pybind11::module &(std::string const &namespace_) > &M);
void bind_Teuchos_ConstTypeTraits(std::function< pybind11::module &(std::string const &namespace_) > &M);
void bind_ROL_ScalarTraits(std::function< pybind11::module &(std::string const &namespace_) > &M);
void bind_Teuchos_any(std::function< pybind11::module &(std::string const &namespace_) > &M);
void bind_Teuchos_ENull(std::function< pybind11::module &(std::string const &namespace_) > &M);
void bind_Teuchos_Ptr(std::function< pybind11::module &(std::string const &namespace_) > &M);
void bind_Teuchos_RCP(std::function< pybind11::module &(std::string const &namespace_) > &M);
void bind_ROL_Ptr(std::function< pybind11::module &(std::string const &namespace_) > &M);
void bind_ROL_Types(std::function< pybind11::module &(std::string const &namespace_) > &M);
void bind_ROL_Elementwise_Reduce(std::function< pybind11::module &(std::string const &namespace_) > &M);
void bind_ROL_UnaryFunctions(std::function< pybind11::module &(std::string const &namespace_) > &M);
void bind_ROL_UpdateType(std::function< pybind11::module &(std::string const &namespace_) > &M);
void bind_Teuchos_DataAccess(std::function< pybind11::module &(std::string const &namespace_) > &M);
void bind_Teuchos_ParameterListExceptions(std::function< pybind11::module &(std::string const &namespace_) > &M);
void bind_Teuchos_iostream_helpers(std::function< pybind11::module &(std::string const &namespace_) > &M);
void bind_Teuchos_FancyOStream(std::function< pybind11::module &(std::string const &namespace_) > &M);
void bind_Teuchos_Dependency(std::function< pybind11::module &(std::string const &namespace_) > &M);
void bind_ROL_ParameterList(std::function< pybind11::module &(std::string const &namespace_) > &M);
void bind_unknown_unknown(std::function< pybind11::module &(std::string const &namespace_) > &M);
void bind_unknown_unknown_1(std::function< pybind11::module &(std::string const &namespace_) > &M);
void bind_unknown_unknown_2(std::function< pybind11::module &(std::string const &namespace_) > &M);
void bind_unknown_unknown_3(std::function< pybind11::module &(std::string const &namespace_) > &M);
void bind_ROL_Step(std::function< pybind11::module &(std::string const &namespace_) > &M);
void bind_ROL_TypeU_Algorithm(std::function< pybind11::module &(std::string const &namespace_) > &M);
void bind_ROL_ValidParameters(std::function< pybind11::module &(std::string const &namespace_) > &M);
void bind_ROL_TypeU_BundleAlgorithm(std::function< pybind11::module &(std::string const &namespace_) > &M);
void bind_ROL_BundleStatusTest(std::function< pybind11::module &(std::string const &namespace_) > &M);
void bind_ROL_TypeU_LineSearchAlgorithm(std::function< pybind11::module &(std::string const &namespace_) > &M);
void bind_ROL_Gradient_U(std::function< pybind11::module &(std::string const &namespace_) > &M);
void bind_ROL_TrustRegionUtilities(std::function< pybind11::module &(std::string const &namespace_) > &M);
void bind_ROL_TypeU_TrustRegionAlgorithm(std::function< pybind11::module &(std::string const &namespace_) > &M);
void bind_ROL_CauchyPoint_U(std::function< pybind11::module &(std::string const &namespace_) > &M);
void bind_ROL_TypeU_AlgorithmFactory(std::function< pybind11::module &(std::string const &namespace_) > &M);
void bind_ROL_TypeB_Algorithm(std::function< pybind11::module &(std::string const &namespace_) > &M);
void bind_ROL_ReducedLinearConstraint(std::function< pybind11::module &(std::string const &namespace_) > &M);
void bind_ROL_TypeB_LinMoreAlgorithm(std::function< pybind11::module &(std::string const &namespace_) > &M);
void bind_ROL_SingletonVector(std::function< pybind11::module &(std::string const &namespace_) > &M);
void bind_ROL_KrylovFactory(std::function< pybind11::module &(std::string const &namespace_) > &M);
void bind_ROL_TypeU_Algorithm_1(std::function< pybind11::module &(std::string const &namespace_) > &M);
void bind_ROL_ScalarController(std::function< pybind11::module &(std::string const &namespace_) > &M);
void bind_ROL_TypeB_MoreauYosidaAlgorithm(std::function< pybind11::module &(std::string const &namespace_) > &M);
void bind_ROL_InteriorPointObjective(std::function< pybind11::module &(std::string const &namespace_) > &M);
void bind_ROL_TypeB_InteriorPointAlgorithm(std::function< pybind11::module &(std::string const &namespace_) > &M);
void bind_ROL_TypeB_QuasiNewtonAlgorithm(std::function< pybind11::module &(std::string const &namespace_) > &M);
void bind_ROL_PQNObjective(std::function< pybind11::module &(std::string const &namespace_) > &M);
void bind_ROL_TypeB_TrustRegionSPGAlgorithm(std::function< pybind11::module &(std::string const &namespace_) > &M);
void bind_ROL_TypeE_Algorithm(std::function< pybind11::module &(std::string const &namespace_) > &M);
void bind_ROL_ConstraintStatusTest(std::function< pybind11::module &(std::string const &namespace_) > &M);
void bind_ROL_TypeE_AugmentedLagrangianAlgorithm(std::function< pybind11::module &(std::string const &namespace_) > &M);
void bind_ROL_FletcherObjectiveBase(std::function< pybind11::module &(std::string const &namespace_) > &M);
void bind_ROL_TypeE_FletcherAlgorithm(std::function< pybind11::module &(std::string const &namespace_) > &M);
void bind_ROL_ElasticLinearConstraint(std::function< pybind11::module &(std::string const &namespace_) > &M);
void bind_ROL_TypeE_StabilizedLCLAlgorithm(std::function< pybind11::module &(std::string const &namespace_) > &M);
void bind_ROL_BoundConstraint(std::function< pybind11::module &(std::string const &namespace_) > &M);
void bind_ROL_TypeE_AlgorithmFactory(std::function< pybind11::module &(std::string const &namespace_) > &M);
void bind_ROL_TypeG_Algorithm(std::function< pybind11::module &(std::string const &namespace_) > &M);
void bind_ROL_TrustRegionTypes(std::function< pybind11::module &(std::string const &namespace_) > &M);
void bind_ROL_KelleySachsModel(std::function< pybind11::module &(std::string const &namespace_) > &M);
void bind_unknown_unknown_4(std::function< pybind11::module &(std::string const &namespace_) > &M);
void bind_unknown_unknown_5(std::function< pybind11::module &(std::string const &namespace_) > &M);


PYBIND11_MODULE(PyROL, root_module) {
	root_module.doc() = "PyROL module";

	std::map <std::string, pybind11::module> modules;
	ModuleGetter M = [&](std::string const &namespace_) -> pybind11::module & {
		auto it = modules.find(namespace_);
		if( it == modules.end() ) throw std::runtime_error("Attempt to access pybind11::module for namespace " + namespace_ + " before it was created!!!");
		return it->second;
	};

	modules[""] = root_module;

	static std::vector<std::string> const reserved_python_words {"nonlocal", "global", };

	auto mangle_namespace_name(
		[](std::string const &ns) -> std::string {
			if ( std::find(reserved_python_words.begin(), reserved_python_words.end(), ns) == reserved_python_words.end() ) return ns;
			else return ns+'_';
		}
	);

	std::vector< std::pair<std::string, std::string> > sub_modules {
		{"", "ROL"},
		{"ROL", "Elementwise"},
		{"ROL", "Exception"},
		{"ROL", "PyROL"},
		{"ROL", "TRUtils"},
		{"ROL", "TypeB"},
		{"ROL", "TypeE"},
		{"ROL", "TypeG"},
		{"ROL", "TypeU"},
		{"", "Teuchos"},
		{"Teuchos", "Exceptions"},
		{"Teuchos", "PtrPrivateUtilityPack"},
		{"", "std"},
	};
	for(auto &p : sub_modules ) modules[p.first.size() ? p.first+"::"+p.second : p.second] = modules[p.first].def_submodule( mangle_namespace_name(p.second).c_str(), ("Bindings for " + p.first + "::" + p.second + " namespace").c_str() );

	//pybind11::class_<std::shared_ptr<void>>(M(""), "_encapsulated_data_");

	bind_std_postypes(M);
	bind_std_typeinfo(M);
	bind_std_locale_classes(M);
	bind_ROL_Elementwise_Function(M);
	bind_std_istream_tcc(M);
	bind_Teuchos_ConstTypeTraits(M);
	bind_ROL_ScalarTraits(M);
	bind_Teuchos_any(M);
	bind_Teuchos_ENull(M);
	bind_Teuchos_Ptr(M);
	bind_Teuchos_RCP(M);
	bind_ROL_Ptr(M);
	bind_ROL_Types(M);
	bind_ROL_Elementwise_Reduce(M);
	bind_ROL_UnaryFunctions(M);
	bind_ROL_UpdateType(M);
	bind_Teuchos_DataAccess(M);
	bind_Teuchos_ParameterListExceptions(M);
	bind_Teuchos_iostream_helpers(M);
	bind_Teuchos_FancyOStream(M);
	bind_Teuchos_Dependency(M);
	bind_ROL_ParameterList(M);
	bind_unknown_unknown(M);
	bind_unknown_unknown_1(M);
	bind_unknown_unknown_2(M);
	bind_unknown_unknown_3(M);
	bind_ROL_Step(M);
	bind_ROL_TypeU_Algorithm(M);
	bind_ROL_ValidParameters(M);
	bind_ROL_TypeU_BundleAlgorithm(M);
	bind_ROL_BundleStatusTest(M);
	bind_ROL_TypeU_LineSearchAlgorithm(M);
	bind_ROL_Gradient_U(M);
	bind_ROL_TrustRegionUtilities(M);
	bind_ROL_TypeU_TrustRegionAlgorithm(M);
	bind_ROL_CauchyPoint_U(M);
	bind_ROL_TypeU_AlgorithmFactory(M);
	bind_ROL_TypeB_Algorithm(M);
	bind_ROL_ReducedLinearConstraint(M);
	bind_ROL_TypeB_LinMoreAlgorithm(M);
	bind_ROL_SingletonVector(M);
	bind_ROL_KrylovFactory(M);
	bind_ROL_TypeU_Algorithm_1(M);
	bind_ROL_ScalarController(M);
	bind_ROL_TypeB_MoreauYosidaAlgorithm(M);
	bind_ROL_InteriorPointObjective(M);
	bind_ROL_TypeB_InteriorPointAlgorithm(M);
	bind_ROL_TypeB_QuasiNewtonAlgorithm(M);
	bind_ROL_PQNObjective(M);
	bind_ROL_TypeB_TrustRegionSPGAlgorithm(M);
	bind_ROL_TypeE_Algorithm(M);
	bind_ROL_ConstraintStatusTest(M);
	bind_ROL_TypeE_AugmentedLagrangianAlgorithm(M);
	bind_ROL_FletcherObjectiveBase(M);
	bind_ROL_TypeE_FletcherAlgorithm(M);
	bind_ROL_ElasticLinearConstraint(M);
	bind_ROL_TypeE_StabilizedLCLAlgorithm(M);
	bind_ROL_BoundConstraint(M);
	bind_ROL_TypeE_AlgorithmFactory(M);
	bind_ROL_TypeG_Algorithm(M);
	bind_ROL_TrustRegionTypes(M);
	bind_ROL_KelleySachsModel(M);
	bind_unknown_unknown_4(M);
	bind_unknown_unknown_5(M);

}
