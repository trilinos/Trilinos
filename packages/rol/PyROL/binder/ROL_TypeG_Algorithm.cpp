#include <ROL_BoundConstraint.hpp>
#include <ROL_BoundConstraint_Partitioned.hpp>
#include <ROL_Bounds.hpp>
#include <ROL_Constraint.hpp>
#include <ROL_Constraint_Partitioned.hpp>
#include <ROL_Elementwise_Function.hpp>
#include <ROL_Elementwise_Reduce.hpp>
#include <ROL_Objective.hpp>
#include <ROL_PartitionedVector.hpp>
#include <ROL_PolyhedralProjection.hpp>
#include <ROL_Problem.hpp>
#include <ROL_SlacklessObjective.hpp>
#include <ROL_StatusTest.hpp>
#include <ROL_TypeG_Algorithm.hpp>
#include <ROL_TypeG_AlgorithmFactory.hpp>
#include <ROL_TypeG_AugmentedLagrangianAlgorithm.hpp>
#include <ROL_TypeG_InteriorPointAlgorithm.hpp>
#include <ROL_TypeG_MoreauYosidaAlgorithm.hpp>
#include <ROL_TypeG_StabilizedLCLAlgorithm.hpp>
#include <ROL_Types.hpp>
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
#include "PyROL_Smart_Holder.hpp"
#include <string>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>


#ifndef BINDER_PYBIND11_TYPE_CASTER
	#define BINDER_PYBIND11_TYPE_CASTER
	PYBIND11_DECLARE_HOLDER_TYPE(T, std::shared_ptr<T>)
	PYBIND11_DECLARE_HOLDER_TYPE(T, T*)
	PYBIND11_MAKE_OPAQUE(std::shared_ptr<void>)
#endif

// ROL::TypeG::AugmentedLagrangianAlgorithm file:ROL_TypeG_AugmentedLagrangianAlgorithm.hpp line:60
struct PyCallBack_ROL_TypeG_AugmentedLagrangianAlgorithm_double_t : public ROL::TypeG::AugmentedLagrangianAlgorithm<double> {
	using ROL::TypeG::AugmentedLagrangianAlgorithm<double>::AugmentedLagrangianAlgorithm;

};

// ROL::TypeG::MoreauYosidaAlgorithm file:ROL_TypeG_MoreauYosidaAlgorithm.hpp line:58
struct PyCallBack_ROL_TypeG_MoreauYosidaAlgorithm_double_t : public ROL::TypeG::MoreauYosidaAlgorithm<double> {
	using ROL::TypeG::MoreauYosidaAlgorithm<double>::MoreauYosidaAlgorithm;

};

// ROL::TypeG::InteriorPointAlgorithm file:ROL_TypeG_InteriorPointAlgorithm.hpp line:58
struct PyCallBack_ROL_TypeG_InteriorPointAlgorithm_double_t : public ROL::TypeG::InteriorPointAlgorithm<double> {
	using ROL::TypeG::InteriorPointAlgorithm<double>::InteriorPointAlgorithm;

};

// ROL::TypeG::StabilizedLCLAlgorithm file:ROL_TypeG_StabilizedLCLAlgorithm.hpp line:61
struct PyCallBack_ROL_TypeG_StabilizedLCLAlgorithm_double_t : public ROL::TypeG::StabilizedLCLAlgorithm<double> {
	using ROL::TypeG::StabilizedLCLAlgorithm<double>::StabilizedLCLAlgorithm;

};

void bind_ROL_TypeG_Algorithm(std::function< pybind11::module &(std::string const &namespace_) > &M)
{
	{ // ROL::TypeG::AlgorithmState file:ROL_TypeG_Algorithm.hpp line:62
		pybind11::class_<ROL::TypeG::AlgorithmState<double>, std::shared_ptr<ROL::TypeG::AlgorithmState<double>>, ROL::AlgorithmState<double>> cl(M("ROL::TypeG"), "AlgorithmState_double_t", "", pybind11::module_local());
		cl.def( pybind11::init( [](){ return new ROL::TypeG::AlgorithmState<double>(); } ) );
		cl.def( pybind11::init( [](ROL::TypeG::AlgorithmState<double> const &o){ return new ROL::TypeG::AlgorithmState<double>(o); } ) );
		cl.def_readwrite("searchSize", &ROL::TypeG::AlgorithmState<double>::searchSize);
		cl.def_readwrite("stepVec", &ROL::TypeG::AlgorithmState<double>::stepVec);
		cl.def_readwrite("gradientVec", &ROL::TypeG::AlgorithmState<double>::gradientVec);
		cl.def_readwrite("constraintVec", &ROL::TypeG::AlgorithmState<double>::constraintVec);
		cl.def("reset", (void (ROL::TypeG::AlgorithmState<double>::*)()) &ROL::TypeG::AlgorithmState<double>::reset, "C++: ROL::TypeG::AlgorithmState<double>::reset() --> void");
		cl.def("assign", (struct ROL::TypeG::AlgorithmState<double> & (ROL::TypeG::AlgorithmState<double>::*)(const struct ROL::TypeG::AlgorithmState<double> &)) &ROL::TypeG::AlgorithmState<double>::operator=, "C++: ROL::TypeG::AlgorithmState<double>::operator=(const struct ROL::TypeG::AlgorithmState<double> &) --> struct ROL::TypeG::AlgorithmState<double> &", pybind11::return_value_policy::automatic, pybind11::arg(""));
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
	{ // ROL::TypeG::Algorithm file:ROL_TypeG_Algorithm.hpp line:90
		pybind11::class_<ROL::TypeG::Algorithm<double>, std::shared_ptr<ROL::TypeG::Algorithm<double>>> cl(M("ROL::TypeG"), "Algorithm_double_t", "", pybind11::module_local());
		cl.def("setStatusTest", [](ROL::TypeG::Algorithm<double> &o, const class std::shared_ptr<class ROL::StatusTest<double> > & a0) -> void { return o.setStatusTest(a0); }, "", pybind11::arg("status"));
		cl.def("setStatusTest", (void (ROL::TypeG::Algorithm<double>::*)(const class std::shared_ptr<class ROL::StatusTest<double> > &, bool)) &ROL::TypeG::Algorithm<double>::setStatusTest, "C++: ROL::TypeG::Algorithm<double>::setStatusTest(const class std::shared_ptr<class ROL::StatusTest<double> > &, bool) --> void", pybind11::arg("status"), pybind11::arg("combineStatus"));
		cl.def("getState", (class std::shared_ptr<const struct ROL::TypeG::AlgorithmState<double> > (ROL::TypeG::Algorithm<double>::*)() const) &ROL::TypeG::Algorithm<double>::getState, "C++: ROL::TypeG::Algorithm<double>::getState() const --> class std::shared_ptr<const struct ROL::TypeG::AlgorithmState<double> >");
		cl.def("reset", (void (ROL::TypeG::Algorithm<double>::*)()) &ROL::TypeG::Algorithm<double>::reset, "C++: ROL::TypeG::Algorithm<double>::reset() --> void");
	}
	{ // ROL::TypeG::AugmentedLagrangianAlgorithm file:ROL_TypeG_AugmentedLagrangianAlgorithm.hpp line:60
		pybind11::class_<ROL::TypeG::AugmentedLagrangianAlgorithm<double>, std::shared_ptr<ROL::TypeG::AugmentedLagrangianAlgorithm<double>>, PyCallBack_ROL_TypeG_AugmentedLagrangianAlgorithm_double_t, ROL::TypeG::Algorithm<double>> cl(M("ROL::TypeG"), "AugmentedLagrangianAlgorithm_double_t", "", pybind11::module_local());
		cl.def( pybind11::init<class Teuchos::ParameterList &>(), pybind11::arg("list") );

		cl.def( pybind11::init( [](PyCallBack_ROL_TypeG_AugmentedLagrangianAlgorithm_double_t const &o){ return new PyCallBack_ROL_TypeG_AugmentedLagrangianAlgorithm_double_t(o); } ) );
		cl.def( pybind11::init( [](ROL::TypeG::AugmentedLagrangianAlgorithm<double> const &o){ return new ROL::TypeG::AugmentedLagrangianAlgorithm<double>(o); } ) );
		cl.def("setStatusTest", [](ROL::TypeG::Algorithm<double> &o, const class std::shared_ptr<class ROL::StatusTest<double> > & a0) -> void { return o.setStatusTest(a0); }, "", pybind11::arg("status"));
		cl.def("setStatusTest", (void (ROL::TypeG::Algorithm<double>::*)(const class std::shared_ptr<class ROL::StatusTest<double> > &, bool)) &ROL::TypeG::Algorithm<double>::setStatusTest, "C++: ROL::TypeG::Algorithm<double>::setStatusTest(const class std::shared_ptr<class ROL::StatusTest<double> > &, bool) --> void", pybind11::arg("status"), pybind11::arg("combineStatus"));
		cl.def("getState", (class std::shared_ptr<const struct ROL::TypeG::AlgorithmState<double> > (ROL::TypeG::Algorithm<double>::*)() const) &ROL::TypeG::Algorithm<double>::getState, "C++: ROL::TypeG::Algorithm<double>::getState() const --> class std::shared_ptr<const struct ROL::TypeG::AlgorithmState<double> >");
		cl.def("reset", (void (ROL::TypeG::Algorithm<double>::*)()) &ROL::TypeG::Algorithm<double>::reset, "C++: ROL::TypeG::Algorithm<double>::reset() --> void");
	}
	{ // ROL::TypeG::MoreauYosidaAlgorithm file:ROL_TypeG_MoreauYosidaAlgorithm.hpp line:58
		pybind11::class_<ROL::TypeG::MoreauYosidaAlgorithm<double>, std::shared_ptr<ROL::TypeG::MoreauYosidaAlgorithm<double>>, PyCallBack_ROL_TypeG_MoreauYosidaAlgorithm_double_t, ROL::TypeG::Algorithm<double>> cl(M("ROL::TypeG"), "MoreauYosidaAlgorithm_double_t", "", pybind11::module_local());
		cl.def( pybind11::init<class Teuchos::ParameterList &>(), pybind11::arg("list") );

		cl.def( pybind11::init( [](PyCallBack_ROL_TypeG_MoreauYosidaAlgorithm_double_t const &o){ return new PyCallBack_ROL_TypeG_MoreauYosidaAlgorithm_double_t(o); } ) );
		cl.def( pybind11::init( [](ROL::TypeG::MoreauYosidaAlgorithm<double> const &o){ return new ROL::TypeG::MoreauYosidaAlgorithm<double>(o); } ) );
		cl.def("setStatusTest", [](ROL::TypeG::Algorithm<double> &o, const class std::shared_ptr<class ROL::StatusTest<double> > & a0) -> void { return o.setStatusTest(a0); }, "", pybind11::arg("status"));
		cl.def("setStatusTest", (void (ROL::TypeG::Algorithm<double>::*)(const class std::shared_ptr<class ROL::StatusTest<double> > &, bool)) &ROL::TypeG::Algorithm<double>::setStatusTest, "C++: ROL::TypeG::Algorithm<double>::setStatusTest(const class std::shared_ptr<class ROL::StatusTest<double> > &, bool) --> void", pybind11::arg("status"), pybind11::arg("combineStatus"));
		cl.def("getState", (class std::shared_ptr<const struct ROL::TypeG::AlgorithmState<double> > (ROL::TypeG::Algorithm<double>::*)() const) &ROL::TypeG::Algorithm<double>::getState, "C++: ROL::TypeG::Algorithm<double>::getState() const --> class std::shared_ptr<const struct ROL::TypeG::AlgorithmState<double> >");
		cl.def("reset", (void (ROL::TypeG::Algorithm<double>::*)()) &ROL::TypeG::Algorithm<double>::reset, "C++: ROL::TypeG::Algorithm<double>::reset() --> void");
	}
	{ // ROL::TypeG::InteriorPointAlgorithm file:ROL_TypeG_InteriorPointAlgorithm.hpp line:58
		pybind11::class_<ROL::TypeG::InteriorPointAlgorithm<double>, std::shared_ptr<ROL::TypeG::InteriorPointAlgorithm<double>>, PyCallBack_ROL_TypeG_InteriorPointAlgorithm_double_t, ROL::TypeG::Algorithm<double>> cl(M("ROL::TypeG"), "InteriorPointAlgorithm_double_t", "", pybind11::module_local());
		cl.def( pybind11::init<class Teuchos::ParameterList &>(), pybind11::arg("list") );

		cl.def( pybind11::init( [](PyCallBack_ROL_TypeG_InteriorPointAlgorithm_double_t const &o){ return new PyCallBack_ROL_TypeG_InteriorPointAlgorithm_double_t(o); } ) );
		cl.def( pybind11::init( [](ROL::TypeG::InteriorPointAlgorithm<double> const &o){ return new ROL::TypeG::InteriorPointAlgorithm<double>(o); } ) );
		cl.def("setStatusTest", [](ROL::TypeG::Algorithm<double> &o, const class std::shared_ptr<class ROL::StatusTest<double> > & a0) -> void { return o.setStatusTest(a0); }, "", pybind11::arg("status"));
		cl.def("setStatusTest", (void (ROL::TypeG::Algorithm<double>::*)(const class std::shared_ptr<class ROL::StatusTest<double> > &, bool)) &ROL::TypeG::Algorithm<double>::setStatusTest, "C++: ROL::TypeG::Algorithm<double>::setStatusTest(const class std::shared_ptr<class ROL::StatusTest<double> > &, bool) --> void", pybind11::arg("status"), pybind11::arg("combineStatus"));
		cl.def("getState", (class std::shared_ptr<const struct ROL::TypeG::AlgorithmState<double> > (ROL::TypeG::Algorithm<double>::*)() const) &ROL::TypeG::Algorithm<double>::getState, "C++: ROL::TypeG::Algorithm<double>::getState() const --> class std::shared_ptr<const struct ROL::TypeG::AlgorithmState<double> >");
		cl.def("reset", (void (ROL::TypeG::Algorithm<double>::*)()) &ROL::TypeG::Algorithm<double>::reset, "C++: ROL::TypeG::Algorithm<double>::reset() --> void");
	}
	{ // ROL::TypeG::StabilizedLCLAlgorithm file:ROL_TypeG_StabilizedLCLAlgorithm.hpp line:61
		pybind11::class_<ROL::TypeG::StabilizedLCLAlgorithm<double>, std::shared_ptr<ROL::TypeG::StabilizedLCLAlgorithm<double>>, PyCallBack_ROL_TypeG_StabilizedLCLAlgorithm_double_t, ROL::TypeG::Algorithm<double>> cl(M("ROL::TypeG"), "StabilizedLCLAlgorithm_double_t", "", pybind11::module_local());
		cl.def( pybind11::init<class Teuchos::ParameterList &>(), pybind11::arg("list") );

		cl.def( pybind11::init( [](PyCallBack_ROL_TypeG_StabilizedLCLAlgorithm_double_t const &o){ return new PyCallBack_ROL_TypeG_StabilizedLCLAlgorithm_double_t(o); } ) );
		cl.def( pybind11::init( [](ROL::TypeG::StabilizedLCLAlgorithm<double> const &o){ return new ROL::TypeG::StabilizedLCLAlgorithm<double>(o); } ) );
		cl.def("setStatusTest", [](ROL::TypeG::Algorithm<double> &o, const class std::shared_ptr<class ROL::StatusTest<double> > & a0) -> void { return o.setStatusTest(a0); }, "", pybind11::arg("status"));
		cl.def("setStatusTest", (void (ROL::TypeG::Algorithm<double>::*)(const class std::shared_ptr<class ROL::StatusTest<double> > &, bool)) &ROL::TypeG::Algorithm<double>::setStatusTest, "C++: ROL::TypeG::Algorithm<double>::setStatusTest(const class std::shared_ptr<class ROL::StatusTest<double> > &, bool) --> void", pybind11::arg("status"), pybind11::arg("combineStatus"));
		cl.def("getState", (class std::shared_ptr<const struct ROL::TypeG::AlgorithmState<double> > (ROL::TypeG::Algorithm<double>::*)() const) &ROL::TypeG::Algorithm<double>::getState, "C++: ROL::TypeG::Algorithm<double>::getState() const --> class std::shared_ptr<const struct ROL::TypeG::AlgorithmState<double> >");
		cl.def("reset", (void (ROL::TypeG::Algorithm<double>::*)()) &ROL::TypeG::Algorithm<double>::reset, "C++: ROL::TypeG::Algorithm<double>::reset() --> void");
	}
	// ROL::TypeG::EAlgorithmG file:ROL_TypeG_AlgorithmFactory.hpp line:64
	pybind11::enum_<ROL::TypeG::EAlgorithmG>(M("ROL::TypeG"), "EAlgorithmG", pybind11::arithmetic(), "Enumeration of generally constrained algorithm types.\n\n    \n    ALGORITHM_G_AUGMENTEDLAGRANGIAN describe\n    \n\n    ALGORITHM_G_MOREAUYOSIDA        describe\n    \n\n    ALGORITHM_G_INTERIORPOINT       describe\n    \n\n    ALGORITHM_G_STABILIZEDLCL       describe", pybind11::module_local())
		.value("ALGORITHM_G_AUGMENTEDLAGRANGIAN", ROL::TypeG::ALGORITHM_G_AUGMENTEDLAGRANGIAN)
		.value("ALGORITHM_G_MOREAUYOSIDA", ROL::TypeG::ALGORITHM_G_MOREAUYOSIDA)
		.value("ALGORITHM_G_INTERIORPOINT", ROL::TypeG::ALGORITHM_G_INTERIORPOINT)
		.value("ALGORITHM_G_STABILIZEDLCL", ROL::TypeG::ALGORITHM_G_STABILIZEDLCL)
		.value("ALGORITHM_G_LAST", ROL::TypeG::ALGORITHM_G_LAST)
		.export_values();

;

	// ROL::TypeG::EAlgorithmGToString(enum ROL::TypeG::EAlgorithmG) file:ROL_TypeG_AlgorithmFactory.hpp line:72
	M("ROL::TypeG").def("EAlgorithmGToString", (std::string (*)(enum ROL::TypeG::EAlgorithmG)) &ROL::TypeG::EAlgorithmGToString, "C++: ROL::TypeG::EAlgorithmGToString(enum ROL::TypeG::EAlgorithmG) --> std::string", pybind11::arg("alg"));

	// ROL::TypeG::isValidAlgorithmG(enum ROL::TypeG::EAlgorithmG) file:ROL_TypeG_AlgorithmFactory.hpp line:90
	M("ROL::TypeG").def("isValidAlgorithmG", (int (*)(enum ROL::TypeG::EAlgorithmG)) &ROL::TypeG::isValidAlgorithmG, "Verifies validity of a AlgorithmG enum.\n\n    \n  [in]  - enum of the AlgorithmG\n    \n\n 1 if the argument is a valid AlgorithmG; 0 otherwise.\n\nC++: ROL::TypeG::isValidAlgorithmG(enum ROL::TypeG::EAlgorithmG) --> int", pybind11::arg("alg"));

	// ROL::TypeG::StringToEAlgorithmG(std::string) file:ROL_TypeG_AlgorithmFactory.hpp line:119
	M("ROL::TypeG").def("StringToEAlgorithmG", (enum ROL::TypeG::EAlgorithmG (*)(std::string)) &ROL::TypeG::StringToEAlgorithmG, "C++: ROL::TypeG::StringToEAlgorithmG(std::string) --> enum ROL::TypeG::EAlgorithmG", pybind11::arg("s"));

	// ROL::TypeG::AlgorithmFactory(class Teuchos::ParameterList &) file:ROL_TypeG_AlgorithmFactory.hpp line:130
	M("ROL::TypeG").def("AlgorithmFactory", (class std::shared_ptr<class ROL::TypeG::Algorithm<double> > (*)(class Teuchos::ParameterList &)) &ROL::TypeG::AlgorithmFactory<double>, "C++: ROL::TypeG::AlgorithmFactory(class Teuchos::ParameterList &) --> class std::shared_ptr<class ROL::TypeG::Algorithm<double> >", pybind11::arg("parlist"));

}
