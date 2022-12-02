#include <ROL_BoundConstraint.hpp>
#include <ROL_Constraint.hpp>
#include <ROL_Elementwise_Function.hpp>
#include <ROL_Elementwise_Reduce.hpp>
#include <ROL_Krylov.hpp>
#include <ROL_Objective.hpp>
#include <ROL_PolyhedralProjection.hpp>
#include <ROL_Problem.hpp>
#include <ROL_Secant.hpp>
#include <ROL_StatusTest.hpp>
#include <ROL_TypeB_Algorithm.hpp>
#include <ROL_TypeB_GradientAlgorithm.hpp>
#include <ROL_TypeB_NewtonKrylovAlgorithm.hpp>
#include <ROL_Types.hpp>
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

// ROL::TypeB::GradientAlgorithm file:ROL_TypeB_GradientAlgorithm.hpp line:57
struct PyCallBack_ROL_TypeB_GradientAlgorithm_double_t : public ROL::TypeB::GradientAlgorithm<double> {
	using ROL::TypeB::GradientAlgorithm<double>::GradientAlgorithm;

};

// ROL::TypeB::NewtonKrylovAlgorithm file:ROL_TypeB_NewtonKrylovAlgorithm.hpp line:59
struct PyCallBack_ROL_TypeB_NewtonKrylovAlgorithm_double_t : public ROL::TypeB::NewtonKrylovAlgorithm<double> {
	using ROL::TypeB::NewtonKrylovAlgorithm<double>::NewtonKrylovAlgorithm;

};

void bind_ROL_TypeB_Algorithm(std::function< pybind11::module &(std::string const &namespace_) > &M)
{
	{ // ROL::TypeB::AlgorithmState file:ROL_TypeB_Algorithm.hpp line:62
		pybind11::class_<ROL::TypeB::AlgorithmState<double>, Teuchos::RCP<ROL::TypeB::AlgorithmState<double>>, ROL::AlgorithmState<double>> cl(M("ROL::TypeB"), "AlgorithmState_double_t", "");
		cl.def( pybind11::init( [](){ return new ROL::TypeB::AlgorithmState<double>(); } ) );
		cl.def( pybind11::init( [](ROL::TypeB::AlgorithmState<double> const &o){ return new ROL::TypeB::AlgorithmState<double>(o); } ) );
		cl.def_readwrite("searchSize", &ROL::TypeB::AlgorithmState<double>::searchSize);
		cl.def_readwrite("stepVec", &ROL::TypeB::AlgorithmState<double>::stepVec);
		cl.def_readwrite("gradientVec", &ROL::TypeB::AlgorithmState<double>::gradientVec);
		cl.def_readwrite("nproj", &ROL::TypeB::AlgorithmState<double>::nproj);
		cl.def("reset", (void (ROL::TypeB::AlgorithmState<double>::*)()) &ROL::TypeB::AlgorithmState<double>::reset, "C++: ROL::TypeB::AlgorithmState<double>::reset() --> void");
		cl.def("assign", (struct ROL::TypeB::AlgorithmState<double> & (ROL::TypeB::AlgorithmState<double>::*)(const struct ROL::TypeB::AlgorithmState<double> &)) &ROL::TypeB::AlgorithmState<double>::operator=, "C++: ROL::TypeB::AlgorithmState<double>::operator=(const struct ROL::TypeB::AlgorithmState<double> &) --> struct ROL::TypeB::AlgorithmState<double> &", pybind11::return_value_policy::automatic, pybind11::arg(""));
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
	{ // ROL::TypeB::Algorithm file:ROL_TypeB_Algorithm.hpp line:84
		pybind11::class_<ROL::TypeB::Algorithm<double>, Teuchos::RCP<ROL::TypeB::Algorithm<double>>> cl(M("ROL::TypeB"), "Algorithm_double_t", "");
		cl.def("setStatusTest", [](ROL::TypeB::Algorithm<double> &o, const class Teuchos::RCP<class ROL::StatusTest<double> > & a0) -> void { return o.setStatusTest(a0); }, "", pybind11::arg("status"));
		cl.def("setStatusTest", (void (ROL::TypeB::Algorithm<double>::*)(const class Teuchos::RCP<class ROL::StatusTest<double> > &, const bool)) &ROL::TypeB::Algorithm<double>::setStatusTest, "C++: ROL::TypeB::Algorithm<double>::setStatusTest(const class Teuchos::RCP<class ROL::StatusTest<double> > &, const bool) --> void", pybind11::arg("status"), pybind11::arg("combineStatus"));
		cl.def("getState", (class Teuchos::RCP<const struct ROL::TypeB::AlgorithmState<double> > (ROL::TypeB::Algorithm<double>::*)() const) &ROL::TypeB::Algorithm<double>::getState, "C++: ROL::TypeB::Algorithm<double>::getState() const --> class Teuchos::RCP<const struct ROL::TypeB::AlgorithmState<double> >");
		cl.def("reset", (void (ROL::TypeB::Algorithm<double>::*)()) &ROL::TypeB::Algorithm<double>::reset, "C++: ROL::TypeB::Algorithm<double>::reset() --> void");
	}
	{ // ROL::TypeB::GradientAlgorithm file:ROL_TypeB_GradientAlgorithm.hpp line:57
		pybind11::class_<ROL::TypeB::GradientAlgorithm<double>, Teuchos::RCP<ROL::TypeB::GradientAlgorithm<double>>, PyCallBack_ROL_TypeB_GradientAlgorithm_double_t, ROL::TypeB::Algorithm<double>> cl(M("ROL::TypeB"), "GradientAlgorithm_double_t", "");
		cl.def( pybind11::init<class Teuchos::ParameterList &>(), pybind11::arg("list") );

		cl.def( pybind11::init( [](PyCallBack_ROL_TypeB_GradientAlgorithm_double_t const &o){ return new PyCallBack_ROL_TypeB_GradientAlgorithm_double_t(o); } ) );
		cl.def( pybind11::init( [](ROL::TypeB::GradientAlgorithm<double> const &o){ return new ROL::TypeB::GradientAlgorithm<double>(o); } ) );
		cl.def("setStatusTest", [](ROL::TypeB::Algorithm<double> &o, const class Teuchos::RCP<class ROL::StatusTest<double> > & a0) -> void { return o.setStatusTest(a0); }, "", pybind11::arg("status"));
		cl.def("setStatusTest", (void (ROL::TypeB::Algorithm<double>::*)(const class Teuchos::RCP<class ROL::StatusTest<double> > &, const bool)) &ROL::TypeB::Algorithm<double>::setStatusTest, "C++: ROL::TypeB::Algorithm<double>::setStatusTest(const class Teuchos::RCP<class ROL::StatusTest<double> > &, const bool) --> void", pybind11::arg("status"), pybind11::arg("combineStatus"));
		cl.def("getState", (class Teuchos::RCP<const struct ROL::TypeB::AlgorithmState<double> > (ROL::TypeB::Algorithm<double>::*)() const) &ROL::TypeB::Algorithm<double>::getState, "C++: ROL::TypeB::Algorithm<double>::getState() const --> class Teuchos::RCP<const struct ROL::TypeB::AlgorithmState<double> >");
		cl.def("reset", (void (ROL::TypeB::Algorithm<double>::*)()) &ROL::TypeB::Algorithm<double>::reset, "C++: ROL::TypeB::Algorithm<double>::reset() --> void");
	}
	{ // ROL::TypeB::NewtonKrylovAlgorithm file:ROL_TypeB_NewtonKrylovAlgorithm.hpp line:59
		pybind11::class_<ROL::TypeB::NewtonKrylovAlgorithm<double>, Teuchos::RCP<ROL::TypeB::NewtonKrylovAlgorithm<double>>, PyCallBack_ROL_TypeB_NewtonKrylovAlgorithm_double_t, ROL::TypeB::Algorithm<double>> cl(M("ROL::TypeB"), "NewtonKrylovAlgorithm_double_t", "");
		cl.def( pybind11::init( [](class Teuchos::ParameterList & a0){ return new ROL::TypeB::NewtonKrylovAlgorithm<double>(a0); }, [](class Teuchos::ParameterList & a0){ return new PyCallBack_ROL_TypeB_NewtonKrylovAlgorithm_double_t(a0); } ), "doc");
		cl.def( pybind11::init<class Teuchos::ParameterList &, const class Teuchos::RCP<class ROL::Secant<double> > &>(), pybind11::arg("list"), pybind11::arg("secant") );

		cl.def( pybind11::init( [](class Teuchos::ParameterList & a0, const class Teuchos::RCP<class ROL::Krylov<double> > & a1){ return new ROL::TypeB::NewtonKrylovAlgorithm<double>(a0, a1); }, [](class Teuchos::ParameterList & a0, const class Teuchos::RCP<class ROL::Krylov<double> > & a1){ return new PyCallBack_ROL_TypeB_NewtonKrylovAlgorithm_double_t(a0, a1); } ), "doc");
		cl.def( pybind11::init<class Teuchos::ParameterList &, const class Teuchos::RCP<class ROL::Krylov<double> > &, const class Teuchos::RCP<class ROL::Secant<double> > &>(), pybind11::arg("list"), pybind11::arg("krylov"), pybind11::arg("secant") );

		cl.def( pybind11::init( [](PyCallBack_ROL_TypeB_NewtonKrylovAlgorithm_double_t const &o){ return new PyCallBack_ROL_TypeB_NewtonKrylovAlgorithm_double_t(o); } ) );
		cl.def( pybind11::init( [](ROL::TypeB::NewtonKrylovAlgorithm<double> const &o){ return new ROL::TypeB::NewtonKrylovAlgorithm<double>(o); } ) );
		cl.def("setStatusTest", [](ROL::TypeB::Algorithm<double> &o, const class Teuchos::RCP<class ROL::StatusTest<double> > & a0) -> void { return o.setStatusTest(a0); }, "", pybind11::arg("status"));
		cl.def("setStatusTest", (void (ROL::TypeB::Algorithm<double>::*)(const class Teuchos::RCP<class ROL::StatusTest<double> > &, const bool)) &ROL::TypeB::Algorithm<double>::setStatusTest, "C++: ROL::TypeB::Algorithm<double>::setStatusTest(const class Teuchos::RCP<class ROL::StatusTest<double> > &, const bool) --> void", pybind11::arg("status"), pybind11::arg("combineStatus"));
		cl.def("getState", (class Teuchos::RCP<const struct ROL::TypeB::AlgorithmState<double> > (ROL::TypeB::Algorithm<double>::*)() const) &ROL::TypeB::Algorithm<double>::getState, "C++: ROL::TypeB::Algorithm<double>::getState() const --> class Teuchos::RCP<const struct ROL::TypeB::AlgorithmState<double> >");
		cl.def("reset", (void (ROL::TypeB::Algorithm<double>::*)()) &ROL::TypeB::Algorithm<double>::reset, "C++: ROL::TypeB::Algorithm<double>::reset() --> void");
	}
}
