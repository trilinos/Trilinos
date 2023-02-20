#include <ROL_BoundConstraint.hpp>
#include <ROL_BoundConstraint_Partitioned.hpp>
#include <ROL_Bounds.hpp>
#include <ROL_Constraint_Partitioned.hpp>
#include <ROL_PartitionedVector.hpp>
#include <ROL_PolyhedralProjection.hpp>
#include <ROL_Problem.hpp>
#include <ROL_SlacklessObjective.hpp>
#include <ROL_Types.hpp>
#include <Teuchos_FilteredIterator.hpp>
#include <Teuchos_ParameterEntry.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_ParameterListModifier.hpp>
#include <Teuchos_RCPDecl.hpp>
#include <Teuchos_StringIndexedOrderedValueObjectContainer.hpp>
#include <cwchar>
#include <deque>
#include <ios>
#include <iterator>
#include <memory>
#include <ostream>
#include <sstream> // __str__
#include <streambuf>
#include <string>

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

// ROL::Problem file: line:49
struct PyCallBack_ROL_Problem_double_t : public ROL::Problem<double> {
	using ROL::Problem<double>::Problem;

	void edit() override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::Problem<double> *>(this), "edit");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Problem::edit();
	}
};

void bind_unknown_unknown_6(std::function< pybind11::module &(std::string const &namespace_) > &M)
{
	{ // ROL::Problem file: line:49
		pybind11::class_<ROL::Problem<double>, std::shared_ptr<ROL::Problem<double>>, PyCallBack_ROL_Problem_double_t> cl(M("ROL"), "Problem_double_t", "", pybind11::module_local());
		cl.def( pybind11::init( [](const class std::shared_ptr<class ROL::Objective<double> > & a0, const class std::shared_ptr<class ROL::Vector<double> > & a1){ return new ROL::Problem<double>(a0, a1); }, [](const class std::shared_ptr<class ROL::Objective<double> > & a0, const class std::shared_ptr<class ROL::Vector<double> > & a1){ return new PyCallBack_ROL_Problem_double_t(a0, a1); } ), "doc");
		cl.def( pybind11::init<const class std::shared_ptr<class ROL::Objective<double> > &, const class std::shared_ptr<class ROL::Vector<double> > &, const class std::shared_ptr<class ROL::Vector<double> > &>(), pybind11::arg("obj"), pybind11::arg("x"), pybind11::arg("g") );

		cl.def( pybind11::init( [](PyCallBack_ROL_Problem_double_t const &o){ return new PyCallBack_ROL_Problem_double_t(o); } ) );
		cl.def( pybind11::init( [](ROL::Problem<double> const &o){ return new ROL::Problem<double>(o); } ) );
		cl.def("addBoundConstraint", (void (ROL::Problem<double>::*)(const class std::shared_ptr<class ROL::BoundConstraint<double> > &)) &ROL::Problem<double>::addBoundConstraint, "Add a bound constraint.\n\n      \n  bound constraint object\n\nC++: ROL::Problem<double>::addBoundConstraint(const class std::shared_ptr<class ROL::BoundConstraint<double> > &) --> void", pybind11::arg("bnd"));
		cl.def("removeBoundConstraint", (void (ROL::Problem<double>::*)()) &ROL::Problem<double>::removeBoundConstraint, "Remove an existing bound constraint.\n\nC++: ROL::Problem<double>::removeBoundConstraint() --> void");
		cl.def("addConstraint", [](ROL::Problem<double> &o, std::string const & a0, const class std::shared_ptr<class ROL::Constraint<double> > & a1, const class std::shared_ptr<class ROL::Vector<double> > & a2) -> void { return o.addConstraint(a0, a1, a2); }, "", pybind11::arg("name"), pybind11::arg("econ"), pybind11::arg("emul"));
		cl.def("addConstraint", [](ROL::Problem<double> &o, std::string const & a0, const class std::shared_ptr<class ROL::Constraint<double> > & a1, const class std::shared_ptr<class ROL::Vector<double> > & a2, const class std::shared_ptr<class ROL::Vector<double> > & a3) -> void { return o.addConstraint(a0, a1, a2, a3); }, "", pybind11::arg("name"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("eres"));
		cl.def("addConstraint", (void (ROL::Problem<double>::*)(std::string, const class std::shared_ptr<class ROL::Constraint<double> > &, const class std::shared_ptr<class ROL::Vector<double> > &, const class std::shared_ptr<class ROL::Vector<double> > &, bool)) &ROL::Problem<double>::addConstraint, "Add an equality constraint.\n\n      \n   the unique constraint identifier\n      \n\n   constraint object\n      \n\n   dual constraint space vector\n      \n\n   primal constraint space vector\n      \n\n  whether or not to clear constraint container\n\nC++: ROL::Problem<double>::addConstraint(std::string, const class std::shared_ptr<class ROL::Constraint<double> > &, const class std::shared_ptr<class ROL::Vector<double> > &, const class std::shared_ptr<class ROL::Vector<double> > &, bool) --> void", pybind11::arg("name"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("eres"), pybind11::arg("reset"));
		cl.def("addConstraint", [](ROL::Problem<double> &o, std::string const & a0, const class std::shared_ptr<class ROL::Constraint<double> > & a1, const class std::shared_ptr<class ROL::Vector<double> > & a2, const class std::shared_ptr<class ROL::BoundConstraint<double> > & a3) -> void { return o.addConstraint(a0, a1, a2, a3); }, "", pybind11::arg("name"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"));
		cl.def("addConstraint", [](ROL::Problem<double> &o, std::string const & a0, const class std::shared_ptr<class ROL::Constraint<double> > & a1, const class std::shared_ptr<class ROL::Vector<double> > & a2, const class std::shared_ptr<class ROL::BoundConstraint<double> > & a3, const class std::shared_ptr<class ROL::Vector<double> > & a4) -> void { return o.addConstraint(a0, a1, a2, a3, a4); }, "", pybind11::arg("name"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("ires"));
		cl.def("addConstraint", (void (ROL::Problem<double>::*)(std::string, const class std::shared_ptr<class ROL::Constraint<double> > &, const class std::shared_ptr<class ROL::Vector<double> > &, const class std::shared_ptr<class ROL::BoundConstraint<double> > &, const class std::shared_ptr<class ROL::Vector<double> > &, bool)) &ROL::Problem<double>::addConstraint, "Add an inequality constraint.\n\n      \n   the unique constraint identifier\n      \n\n   constraint object\n      \n\n   dual constraint space vector\n      \n\n   bound constraint\n      \n\n   primal constraint space vector\n      \n\n  whether or not to clear constraint container\n\nC++: ROL::Problem<double>::addConstraint(std::string, const class std::shared_ptr<class ROL::Constraint<double> > &, const class std::shared_ptr<class ROL::Vector<double> > &, const class std::shared_ptr<class ROL::BoundConstraint<double> > &, const class std::shared_ptr<class ROL::Vector<double> > &, bool) --> void", pybind11::arg("name"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("ires"), pybind11::arg("reset"));
		cl.def("removeConstraint", (void (ROL::Problem<double>::*)(std::string)) &ROL::Problem<double>::removeConstraint, "Remove an existing constraint.\n\n      \n  the unique constraint identifier\n\nC++: ROL::Problem<double>::removeConstraint(std::string) --> void", pybind11::arg("name"));
		cl.def("addLinearConstraint", [](ROL::Problem<double> &o, std::string const & a0, const class std::shared_ptr<class ROL::Constraint<double> > & a1, const class std::shared_ptr<class ROL::Vector<double> > & a2) -> void { return o.addLinearConstraint(a0, a1, a2); }, "", pybind11::arg("name"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"));
		cl.def("addLinearConstraint", [](ROL::Problem<double> &o, std::string const & a0, const class std::shared_ptr<class ROL::Constraint<double> > & a1, const class std::shared_ptr<class ROL::Vector<double> > & a2, const class std::shared_ptr<class ROL::Vector<double> > & a3) -> void { return o.addLinearConstraint(a0, a1, a2, a3); }, "", pybind11::arg("name"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"), pybind11::arg("linear_eres"));
		cl.def("addLinearConstraint", (void (ROL::Problem<double>::*)(std::string, const class std::shared_ptr<class ROL::Constraint<double> > &, const class std::shared_ptr<class ROL::Vector<double> > &, const class std::shared_ptr<class ROL::Vector<double> > &, bool)) &ROL::Problem<double>::addLinearConstraint, "Add a linear equality constraint.\n\n      \n         the unique constraint identifier\n      \n\n  constraint object\n      \n\n  dual constraint space vector\n      \n\n  primal constraint space vector\n      \n\n        whether or not to clear linear constraint container\n\nC++: ROL::Problem<double>::addLinearConstraint(std::string, const class std::shared_ptr<class ROL::Constraint<double> > &, const class std::shared_ptr<class ROL::Vector<double> > &, const class std::shared_ptr<class ROL::Vector<double> > &, bool) --> void", pybind11::arg("name"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"), pybind11::arg("linear_eres"), pybind11::arg("reset"));
		cl.def("addLinearConstraint", [](ROL::Problem<double> &o, std::string const & a0, const class std::shared_ptr<class ROL::Constraint<double> > & a1, const class std::shared_ptr<class ROL::Vector<double> > & a2, const class std::shared_ptr<class ROL::BoundConstraint<double> > & a3) -> void { return o.addLinearConstraint(a0, a1, a2, a3); }, "", pybind11::arg("name"), pybind11::arg("linear_icon"), pybind11::arg("linear_imul"), pybind11::arg("linear_ibnd"));
		cl.def("addLinearConstraint", [](ROL::Problem<double> &o, std::string const & a0, const class std::shared_ptr<class ROL::Constraint<double> > & a1, const class std::shared_ptr<class ROL::Vector<double> > & a2, const class std::shared_ptr<class ROL::BoundConstraint<double> > & a3, const class std::shared_ptr<class ROL::Vector<double> > & a4) -> void { return o.addLinearConstraint(a0, a1, a2, a3, a4); }, "", pybind11::arg("name"), pybind11::arg("linear_icon"), pybind11::arg("linear_imul"), pybind11::arg("linear_ibnd"), pybind11::arg("linear_ires"));
		cl.def("addLinearConstraint", (void (ROL::Problem<double>::*)(std::string, const class std::shared_ptr<class ROL::Constraint<double> > &, const class std::shared_ptr<class ROL::Vector<double> > &, const class std::shared_ptr<class ROL::BoundConstraint<double> > &, const class std::shared_ptr<class ROL::Vector<double> > &, bool)) &ROL::Problem<double>::addLinearConstraint, "Add a linear inequality constraint.\n\n      \n         the unique constraint identifier\n      \n\n  constraint object\n      \n\n  dual constraint space vector\n      \n\n  bound constraint\n      \n\n  primal constraint space vector\n      \n\n        whether or not to clear linear constraint container\n\nC++: ROL::Problem<double>::addLinearConstraint(std::string, const class std::shared_ptr<class ROL::Constraint<double> > &, const class std::shared_ptr<class ROL::Vector<double> > &, const class std::shared_ptr<class ROL::BoundConstraint<double> > &, const class std::shared_ptr<class ROL::Vector<double> > &, bool) --> void", pybind11::arg("name"), pybind11::arg("linear_icon"), pybind11::arg("linear_imul"), pybind11::arg("linear_ibnd"), pybind11::arg("linear_ires"), pybind11::arg("reset"));
		cl.def("removeLinearConstraint", (void (ROL::Problem<double>::*)(std::string)) &ROL::Problem<double>::removeLinearConstraint, "Remove an existing linear constraint.\n\n      \n  the unique constraint identifier\n\nC++: ROL::Problem<double>::removeLinearConstraint(std::string) --> void", pybind11::arg("name"));
		cl.def("setProjectionAlgorithm", (void (ROL::Problem<double>::*)(class Teuchos::ParameterList &)) &ROL::Problem<double>::setProjectionAlgorithm, "Set polyhedral projection algorithm.\n\n      \n  polyhedral projection algorithm\n\nC++: ROL::Problem<double>::setProjectionAlgorithm(class Teuchos::ParameterList &) --> void", pybind11::arg("list"));
		cl.def("getObjective", (const class std::shared_ptr<class ROL::Objective<double> > & (ROL::Problem<double>::*)()) &ROL::Problem<double>::getObjective, "Get the objective function.\n\nC++: ROL::Problem<double>::getObjective() --> const class std::shared_ptr<class ROL::Objective<double> > &", pybind11::return_value_policy::automatic);
		cl.def("getPrimalOptimizationVector", (const class std::shared_ptr<class ROL::Vector<double> > & (ROL::Problem<double>::*)()) &ROL::Problem<double>::getPrimalOptimizationVector, "Get the primal optimization space vector.\n\nC++: ROL::Problem<double>::getPrimalOptimizationVector() --> const class std::shared_ptr<class ROL::Vector<double> > &", pybind11::return_value_policy::automatic);
		cl.def("getDualOptimizationVector", (const class std::shared_ptr<class ROL::Vector<double> > & (ROL::Problem<double>::*)()) &ROL::Problem<double>::getDualOptimizationVector, "Get the dual optimization space vector.\n\nC++: ROL::Problem<double>::getDualOptimizationVector() --> const class std::shared_ptr<class ROL::Vector<double> > &", pybind11::return_value_policy::automatic);
		cl.def("getBoundConstraint", (const class std::shared_ptr<class ROL::BoundConstraint<double> > & (ROL::Problem<double>::*)()) &ROL::Problem<double>::getBoundConstraint, "Get the bound constraint.\n\nC++: ROL::Problem<double>::getBoundConstraint() --> const class std::shared_ptr<class ROL::BoundConstraint<double> > &", pybind11::return_value_policy::automatic);
		cl.def("getConstraint", (const class std::shared_ptr<class ROL::Constraint<double> > & (ROL::Problem<double>::*)()) &ROL::Problem<double>::getConstraint, "Get the equality constraint.\n\nC++: ROL::Problem<double>::getConstraint() --> const class std::shared_ptr<class ROL::Constraint<double> > &", pybind11::return_value_policy::automatic);
		cl.def("getMultiplierVector", (const class std::shared_ptr<class ROL::Vector<double> > & (ROL::Problem<double>::*)()) &ROL::Problem<double>::getMultiplierVector, "Get the dual constraint space vector.\n\nC++: ROL::Problem<double>::getMultiplierVector() --> const class std::shared_ptr<class ROL::Vector<double> > &", pybind11::return_value_policy::automatic);
		cl.def("getResidualVector", (const class std::shared_ptr<class ROL::Vector<double> > & (ROL::Problem<double>::*)()) &ROL::Problem<double>::getResidualVector, "Get the primal constraint space vector.\n\nC++: ROL::Problem<double>::getResidualVector() --> const class std::shared_ptr<class ROL::Vector<double> > &", pybind11::return_value_policy::automatic);
		cl.def("getPolyhedralProjection", (const class std::shared_ptr<class ROL::PolyhedralProjection<double> > & (ROL::Problem<double>::*)()) &ROL::Problem<double>::getPolyhedralProjection, "Get the polyhedral projection object.  This is a null pointer if\n             no linear constraints and/or bounds are present.\n\nC++: ROL::Problem<double>::getPolyhedralProjection() --> const class std::shared_ptr<class ROL::PolyhedralProjection<double> > &", pybind11::return_value_policy::automatic);
		cl.def("getProblemType", (enum ROL::EProblem (ROL::Problem<double>::*)()) &ROL::Problem<double>::getProblemType, "Get the optimization problem type (U, B, E, or G).\n\nC++: ROL::Problem<double>::getProblemType() --> enum ROL::EProblem");
		cl.def("isFinalized", (bool (ROL::Problem<double>::*)() const) &ROL::Problem<double>::isFinalized, "Indicate whether or no finalize has been called.\n\nC++: ROL::Problem<double>::isFinalized() const --> bool");
		cl.def("edit", (void (ROL::Problem<double>::*)()) &ROL::Problem<double>::edit, "Resume editting optimization problem after finalize has been called.\n\nC++: ROL::Problem<double>::edit() --> void");
		cl.def("finalizeIteration", (void (ROL::Problem<double>::*)()) &ROL::Problem<double>::finalizeIteration, "Transform the optimization variables to the native\n             parameterization after an optimization algorithm has finished.\n\nC++: ROL::Problem<double>::finalizeIteration() --> void");
		cl.def("assign", (class ROL::Problem<double> & (ROL::Problem<double>::*)(const class ROL::Problem<double> &)) &ROL::Problem<double>::operator=, "C++: ROL::Problem<double>::operator=(const class ROL::Problem<double> &) --> class ROL::Problem<double> &", pybind11::return_value_policy::automatic, pybind11::arg(""));
	}
}
