#include <PyROL_Teuchos_Custom.hpp>
#include <PyROL_Teuchos_ETI.hpp>
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
#include <streambuf>
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

void bind_PyROL_Teuchos_ETI(std::function< pybind11::module &(std::string const &namespace_) > &M)
{
	// Teuchos::initiate(class Teuchos::ParameterList) file:PyROL_Teuchos_ETI.hpp line:7
	M("Teuchos").def("initiate", (void (*)(class Teuchos::ParameterList)) &Teuchos::initiate, "C++: Teuchos::initiate(class Teuchos::ParameterList) --> void", pybind11::arg("p"));

}
