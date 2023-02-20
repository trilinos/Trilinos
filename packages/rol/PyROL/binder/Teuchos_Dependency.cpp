#include <PyROL_Teuchos_Custom.hpp>
#include <Teuchos_Dependency.hpp>
#include <Teuchos_DependencySheet.hpp>
#include <Teuchos_Describable.hpp>
#include <Teuchos_ENull.hpp>
#include <Teuchos_FancyOStream.hpp>
#include <Teuchos_FilteredIterator.hpp>
#include <Teuchos_LabeledObject.hpp>
#include <Teuchos_ParameterEntry.hpp>
#include <Teuchos_ParameterEntryValidator.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_ParameterListModifier.hpp>
#include <Teuchos_PtrDecl.hpp>
#include <Teuchos_RCPDecl.hpp>
#include <Teuchos_RCPNode.hpp>
#include <Teuchos_StringIndexedOrderedValueObjectContainer.hpp>
#include <Teuchos_VerbosityLevel.hpp>
#include <Teuchos_XMLParameterListCoreHelpers.hpp>
#include <Teuchos_any.hpp>
#include <cwchar>
#include <deque>
#include <ios>
#include <iterator>
#include <locale>
#include <map>
#include <memory>
#include <ostream>
#include <set>
#include <sstream> // __str__
#include <streambuf>
#include <string>
#include <utility>

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

// Teuchos::Dependency file:Teuchos_Dependency.hpp line:64
struct PyCallBack_Teuchos_Dependency : public Teuchos::Dependency {
	using Teuchos::Dependency::Dependency;

	std::string getTypeAttributeValue() const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Teuchos::Dependency *>(this), "getTypeAttributeValue");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<std::string>::value) {
				static pybind11::detail::override_caster_t<std::string> caster;
				return pybind11::detail::cast_ref<std::string>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<std::string>(std::move(o));
		}
		pybind11::pybind11_fail("Tried to call pure virtual function \"Dependency::getTypeAttributeValue\"");
	}
	void evaluate() override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Teuchos::Dependency *>(this), "evaluate");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		pybind11::pybind11_fail("Tried to call pure virtual function \"Dependency::evaluate\"");
	}
	void print(std::ostream & a0) const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Teuchos::Dependency *>(this), "print");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Dependency::print(a0);
	}
	void validateDep() const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Teuchos::Dependency *>(this), "validateDep");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		pybind11::pybind11_fail("Tried to call pure virtual function \"Dependency::validateDep\"");
	}
	std::string description() const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Teuchos::Dependency *>(this), "description");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<std::string>::value) {
				static pybind11::detail::override_caster_t<std::string> caster;
				return pybind11::detail::cast_ref<std::string>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<std::string>(std::move(o));
		}
		return Describable::description();
	}
	void describe(class Teuchos::basic_FancyOStream<char, struct std::char_traits<char> > & a0, const enum Teuchos::EVerbosityLevel a1) const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Teuchos::Dependency *>(this), "describe");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Describable::describe(a0, a1);
	}
	void setObjectLabel(const std::string & a0) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Teuchos::Dependency *>(this), "setObjectLabel");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return LabeledObject::setObjectLabel(a0);
	}
	std::string getObjectLabel() const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Teuchos::Dependency *>(this), "getObjectLabel");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<std::string>::value) {
				static pybind11::detail::override_caster_t<std::string> caster;
				return pybind11::detail::cast_ref<std::string>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<std::string>(std::move(o));
		}
		return LabeledObject::getObjectLabel();
	}
};

void bind_Teuchos_Dependency(std::function< pybind11::module &(std::string const &namespace_) > &M)
{
	{ // Teuchos::Dependency file:Teuchos_Dependency.hpp line:64
		pybind11::class_<Teuchos::Dependency, Teuchos::RCP<Teuchos::Dependency>, PyCallBack_Teuchos_Dependency> cl(M("Teuchos"), "Dependency", "This class represents a depndency between elements in a Parameter List.\n\n  DependencySheet\n  ParameterList", pybind11::module_local());
		cl.def( pybind11::init<class Teuchos::RCP<const class Teuchos::ParameterEntry>, class Teuchos::RCP<class Teuchos::ParameterEntry>>(), pybind11::arg("dependee"), pybind11::arg("dependent") );

		cl.def(pybind11::init<PyCallBack_Teuchos_Dependency const &>());
		cl.def("getFirstDependee", (class Teuchos::RCP<const class Teuchos::ParameterEntry> (Teuchos::Dependency::*)() const) &Teuchos::Dependency::getFirstDependee, "Gets the first dependee in the dependees list.\n This is a convience function.\n\nC++: Teuchos::Dependency::getFirstDependee() const --> class Teuchos::RCP<const class Teuchos::ParameterEntry>");
		cl.def("getTypeAttributeValue", (std::string (Teuchos::Dependency::*)() const) &Teuchos::Dependency::getTypeAttributeValue, "Returns the string to be used for the value of the\n type attribute when converting the dependency to XML.\n\nC++: Teuchos::Dependency::getTypeAttributeValue() const --> std::string");
		cl.def_static("getXMLTagName", (const std::string & (*)()) &Teuchos::Dependency::getXMLTagName, "Returns the XML tag to use when serializing Dependencies.\n\nC++: Teuchos::Dependency::getXMLTagName() --> const std::string &", pybind11::return_value_policy::automatic);
		cl.def("evaluate", (void (Teuchos::Dependency::*)()) &Teuchos::Dependency::evaluate, "Evaluates the dependency and makes any appropriate changes to the\n dependee based on the dependent.\n\nC++: Teuchos::Dependency::evaluate() --> void");
		cl.def("print", (void (Teuchos::Dependency::*)(std::ostream &) const) &Teuchos::Dependency::print, "prints out information about the dependency. \n\nC++: Teuchos::Dependency::print(std::ostream &) const --> void", pybind11::arg("out"));
		cl.def("assign", (class Teuchos::Dependency & (Teuchos::Dependency::*)(const class Teuchos::Dependency &)) &Teuchos::Dependency::operator=, "C++: Teuchos::Dependency::operator=(const class Teuchos::Dependency &) --> class Teuchos::Dependency &", pybind11::return_value_policy::automatic, pybind11::arg(""));
	}
	{ // Teuchos::DependencySheet file:Teuchos_DependencySheet.hpp line:61
		pybind11::class_<Teuchos::DependencySheet, Teuchos::RCP<Teuchos::DependencySheet>> cl(M("Teuchos"), "DependencySheet", "A Dependency sheet keeps track of dependencies between various\n ParameterEntries", pybind11::module_local());
		cl.def( pybind11::init( [](){ return new Teuchos::DependencySheet(); } ) );
		cl.def( pybind11::init<const std::string &>(), pybind11::arg("name") );

		cl.def( pybind11::init( [](Teuchos::DependencySheet const &o){ return new Teuchos::DependencySheet(o); } ) );
		cl.def("addDependency", (void (Teuchos::DependencySheet::*)(class Teuchos::RCP<class Teuchos::Dependency>)) &Teuchos::DependencySheet::addDependency, "Adds a dependency to the sheet.\n\n \n The dependency to be added.\n\nC++: Teuchos::DependencySheet::addDependency(class Teuchos::RCP<class Teuchos::Dependency>) --> void", pybind11::arg("dependency"));
		cl.def("addDependencies", (void (Teuchos::DependencySheet::*)(class Teuchos::RCP<class Teuchos::DependencySheet>)) &Teuchos::DependencySheet::addDependencies, "Adds a dependencies from another she\n to this sheet.\n\n \n The other sheet from which\n to add dependencies.\n\nC++: Teuchos::DependencySheet::addDependencies(class Teuchos::RCP<class Teuchos::DependencySheet>) --> void", pybind11::arg("otherSheet"));
		cl.def("removeDependency", (void (Teuchos::DependencySheet::*)(class Teuchos::RCP<class Teuchos::Dependency>)) &Teuchos::DependencySheet::removeDependency, "Removes a particular dependency between two parameters.\n\n \n The dependency to be removed.\n \n\n True if the removal was sucessfull, false otherwise.\n\nC++: Teuchos::DependencySheet::removeDependency(class Teuchos::RCP<class Teuchos::Dependency>) --> void", pybind11::arg("dependency"));
		cl.def("setName", (void (Teuchos::DependencySheet::*)(const std::string)) &Teuchos::DependencySheet::setName, "sets the name of the dependency sheet\n\nC++: Teuchos::DependencySheet::setName(const std::string) --> void", pybind11::arg("newName"));
		cl.def("hasDependents", (bool (Teuchos::DependencySheet::*)(class Teuchos::RCP<const class Teuchos::ParameterEntry>) const) &Teuchos::DependencySheet::hasDependents, "Determines whether or not a parameter is depended upon by any another\n parameters or parameter lists.\n\n \n The paramteter to be checked for dependents.\n \n\n True if the parameter you're checking has other dependents, false otherwise.\n\nC++: Teuchos::DependencySheet::hasDependents(class Teuchos::RCP<const class Teuchos::ParameterEntry>) const --> bool", pybind11::arg("dependee"));
		cl.def("getName", (const std::string & (Teuchos::DependencySheet::*)() const) &Teuchos::DependencySheet::getName, "Gets the name of the dependency sheet.\n\nC++: Teuchos::DependencySheet::getName() const --> const std::string &", pybind11::return_value_policy::automatic);
		cl.def("empty", (bool (Teuchos::DependencySheet::*)() const) &Teuchos::DependencySheet::empty, "Determines whether or not this dependency sheet has any dependencies.\n\nC++: Teuchos::DependencySheet::empty() const --> bool");
		cl.def("size", (unsigned long (Teuchos::DependencySheet::*)()) &Teuchos::DependencySheet::size, "Returns the number of Dependencies in this\n DependencySheet.\n\n \n The number of Depenedencies in this\n DependencySheet.\n\nC++: Teuchos::DependencySheet::size() --> unsigned long");
		cl.def("printDeps", (void (Teuchos::DependencySheet::*)(std::ostream &) const) &Teuchos::DependencySheet::printDeps, "Prints out a list of the dependencies in the DependencySheet\n\nC++: Teuchos::DependencySheet::printDeps(std::ostream &) const --> void", pybind11::arg("out"));
		cl.def_static("getNameAttributeName", (const std::string & (*)()) &Teuchos::DependencySheet::getNameAttributeName, "When serializing to XML, this string should be used as the name\n of the name attribute \n\nC++: Teuchos::DependencySheet::getNameAttributeName() --> const std::string &", pybind11::return_value_policy::automatic);
	}
	// Teuchos::updateParametersFromXmlFile(const std::string &, const class Teuchos::Ptr<class Teuchos::ParameterList> &) file:Teuchos_XMLParameterListCoreHelpers.hpp line:71
	M("Teuchos").def("updateParametersFromXmlFile", (void (*)(const std::string &, const class Teuchos::Ptr<class Teuchos::ParameterList> &)) &Teuchos::updateParametersFromXmlFile, "Reads XML parameters from a file and updates those already in the\n given parameter list.\n\n \n [in] The file name containing XML parameter list\n specification.\n\n \n [in/out] On input, *paramList may be empty or\n contain some parameters and sublists. On output, parameters and sublist\n from the file xmlFileName will be set or overide those in\n *paramList.\n\n \n\n \n\nC++: Teuchos::updateParametersFromXmlFile(const std::string &, const class Teuchos::Ptr<class Teuchos::ParameterList> &) --> void", pybind11::arg("xmlFileName"), pybind11::arg("paramList"));

	// Teuchos::getParametersFromXmlFile(const std::string &) file:Teuchos_XMLParameterListCoreHelpers.hpp line:86
	M("Teuchos").def("getParametersFromXmlFile", (class Teuchos::RCP<class Teuchos::ParameterList> (*)(const std::string &)) &Teuchos::getParametersFromXmlFile, "Reads XML parameters from a file and return them in a new parameter\n list.\n\n \n [in] The file name containing XML parameter list\n specification.\n\n \n\n \n\nC++: Teuchos::getParametersFromXmlFile(const std::string &) --> class Teuchos::RCP<class Teuchos::ParameterList>", pybind11::arg("xmlFileName"));

	// Teuchos::getParametersFromXmlFile(const std::string &, class Teuchos::RCP<class Teuchos::DependencySheet>) file:Teuchos_XMLParameterListCoreHelpers.hpp line:101
	M("Teuchos").def("getParametersFromXmlFile", (class Teuchos::RCP<class Teuchos::ParameterList> (*)(const std::string &, class Teuchos::RCP<class Teuchos::DependencySheet>)) &Teuchos::getParametersFromXmlFile, "Reads XML parameters from a file and return them in a new parameter\n list.\n\n \n [in] The file name containing XML parameter list\n specification.\n\n \n [out] The Dependency Sheet into which Dependencies should be\n placed.\n\n \n\n \n\nC++: Teuchos::getParametersFromXmlFile(const std::string &, class Teuchos::RCP<class Teuchos::DependencySheet>) --> class Teuchos::RCP<class Teuchos::ParameterList>", pybind11::arg("xmlFileName"), pybind11::arg("depSheet"));

	// Teuchos::updateParametersFromXmlString(const std::string &, const class Teuchos::Ptr<class Teuchos::ParameterList> &, bool) file:Teuchos_XMLParameterListCoreHelpers.hpp line:122
	M("Teuchos").def("updateParametersFromXmlString", [](const std::string & a0, const class Teuchos::Ptr<class Teuchos::ParameterList> & a1) -> void { return Teuchos::updateParametersFromXmlString(a0, a1); }, "", pybind11::arg("xmlStr"), pybind11::arg("paramList"));
	M("Teuchos").def("updateParametersFromXmlString", (void (*)(const std::string &, const class Teuchos::Ptr<class Teuchos::ParameterList> &, bool)) &Teuchos::updateParametersFromXmlString, "Reads XML parameters from a std::string and updates those already in the\n given parameter list.\n\n \n [in] String containing XML parameter list specification.\n\n \n [in/out] On input, *paramList may be empty or\n contain some parameters and sublists. On output, parameters and sublist\n from the file xmlStr will be set or override (or not) those in\n *paramList depending on the overwrite parameter.\n\n \n [in] If true, parameters and sublists in the xmlStr \n will override those in paramList.  If false, any value set in \n paramList will be kept, only values not set will be updated.\n\n \n\n \n\nC++: Teuchos::updateParametersFromXmlString(const std::string &, const class Teuchos::Ptr<class Teuchos::ParameterList> &, bool) --> void", pybind11::arg("xmlStr"), pybind11::arg("paramList"), pybind11::arg("overwrite"));

	// Teuchos::getParametersFromXmlString(const std::string &) file:Teuchos_XMLParameterListCoreHelpers.hpp line:137
	M("Teuchos").def("getParametersFromXmlString", (class Teuchos::RCP<class Teuchos::ParameterList> (*)(const std::string &)) &Teuchos::getParametersFromXmlString, "Reads XML parameters from a std::string and return them in a new\n parameter list.\n\n \n [in] String containing XML parameter list specification.\n\n \n\n \n\nC++: Teuchos::getParametersFromXmlString(const std::string &) --> class Teuchos::RCP<class Teuchos::ParameterList>", pybind11::arg("xmlStr"));

	// Teuchos::getParametersFromXmlString(const std::string &, class Teuchos::RCP<class Teuchos::DependencySheet>) file:Teuchos_XMLParameterListCoreHelpers.hpp line:150
	M("Teuchos").def("getParametersFromXmlString", (class Teuchos::RCP<class Teuchos::ParameterList> (*)(const std::string &, class Teuchos::RCP<class Teuchos::DependencySheet>)) &Teuchos::getParametersFromXmlString, "Reads XML parameters from a std::string and return them in a new\n parameter list.\n\n \n [in] String containing XML parameter list specification.\n \n\n [in] The Dependency Sheet into which Dependencies should be\n placed.\n\n \n\n \n\nC++: Teuchos::getParametersFromXmlString(const std::string &, class Teuchos::RCP<class Teuchos::DependencySheet>) --> class Teuchos::RCP<class Teuchos::ParameterList>", pybind11::arg("xmlStr"), pybind11::arg("depSheet"));

	// Teuchos::writeParameterListToXmlOStream(const class Teuchos::ParameterList &, std::ostream &, class Teuchos::RCP<const class Teuchos::DependencySheet>) file:Teuchos_XMLParameterListCoreHelpers.hpp line:166
	M("Teuchos").def("writeParameterListToXmlOStream", [](const class Teuchos::ParameterList & a0, std::ostream & a1) -> void { return Teuchos::writeParameterListToXmlOStream(a0, a1); }, "", pybind11::arg("paramList"), pybind11::arg("xmlOut"));
	M("Teuchos").def("writeParameterListToXmlOStream", (void (*)(const class Teuchos::ParameterList &, std::ostream &, class Teuchos::RCP<const class Teuchos::DependencySheet>)) &Teuchos::writeParameterListToXmlOStream, "Write parameters and sublists in XML format to an std::ostream.\n\n \n [in] Contains the parameters and sublists that will be\n written to file.\n\n \n [in] The stream that will get the XML output.\n\n \n [in] The Dependency Sheet which should be written out.\n\n \n\n \n\nC++: Teuchos::writeParameterListToXmlOStream(const class Teuchos::ParameterList &, std::ostream &, class Teuchos::RCP<const class Teuchos::DependencySheet>) --> void", pybind11::arg("paramList"), pybind11::arg("xmlOut"), pybind11::arg("depSheet"));

	// Teuchos::writeParameterListToXmlFile(const class Teuchos::ParameterList &, const std::string &, class Teuchos::RCP<const class Teuchos::DependencySheet>) file:Teuchos_XMLParameterListCoreHelpers.hpp line:186
	M("Teuchos").def("writeParameterListToXmlFile", [](const class Teuchos::ParameterList & a0, const std::string & a1) -> void { return Teuchos::writeParameterListToXmlFile(a0, a1); }, "", pybind11::arg("paramList"), pybind11::arg("xmlFileName"));
	M("Teuchos").def("writeParameterListToXmlFile", (void (*)(const class Teuchos::ParameterList &, const std::string &, class Teuchos::RCP<const class Teuchos::DependencySheet>)) &Teuchos::writeParameterListToXmlFile, "Write parameters and sublist to an XML file.\n\n \n [in] Contains the parameters and sublists that will be\n written to file.\n\n \n [in] The file name that will be create to contain the\n XML version of the parameter list specification.\n\n \n [in] The Dependency Sheet which should be written out.\n\n \n\n \n\nC++: Teuchos::writeParameterListToXmlFile(const class Teuchos::ParameterList &, const std::string &, class Teuchos::RCP<const class Teuchos::DependencySheet>) --> void", pybind11::arg("paramList"), pybind11::arg("xmlFileName"), pybind11::arg("depSheet"));

}
