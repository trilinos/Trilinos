#include <Teuchos_CompObject.hpp>
#include <Teuchos_DataAccess.hpp>
#include <Teuchos_Flops.hpp>
#include <Teuchos_Object.hpp>
#include <cwchar>
#include <ios>
#include <iterator>
#include <locale>
#include <memory>
#include <ostream>
#include <sstream> // __str__
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

// Teuchos::Object file:Teuchos_Object.hpp line:68
struct PyCallBack_Teuchos_Object : public Teuchos::Object {
	using Teuchos::Object::Object;

	void setLabel(const char * a0) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Teuchos::Object *>(this), "setLabel");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Object::setLabel(a0);
	}
	const char * label() const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Teuchos::Object *>(this), "label");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<const char *>::value) {
				static pybind11::detail::override_caster_t<const char *> caster;
				return pybind11::detail::cast_ref<const char *>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<const char *>(std::move(o));
		}
		return Object::label();
	}
	int reportError(const std::string a0, int a1) const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Teuchos::Object *>(this), "reportError");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1);
			if (pybind11::detail::cast_is_temporary_value_reference<int>::value) {
				static pybind11::detail::override_caster_t<int> caster;
				return pybind11::detail::cast_ref<int>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<int>(std::move(o));
		}
		return Object::reportError(a0, a1);
	}
};

void bind_Teuchos_DataAccess(std::function< pybind11::module &(std::string const &namespace_) > &M)
{
	// Teuchos::DataAccess file:Teuchos_DataAccess.hpp line:60
	pybind11::enum_<Teuchos::DataAccess>(M("Teuchos"), "DataAccess", pybind11::arithmetic(), "If set to Copy, user data will be copied at construction.\n      If set to View, user data will be encapsulated and used throughout\n      the life of the object.")
		.value("Copy", Teuchos::Copy)
		.value("View", Teuchos::View)
		.export_values();

;

	{ // Teuchos::Object file:Teuchos_Object.hpp line:68
		pybind11::class_<Teuchos::Object, Teuchos::RCP<Teuchos::Object>, PyCallBack_Teuchos_Object> cl(M("Teuchos"), "Object", "");
		cl.def( pybind11::init( [](){ return new Teuchos::Object(); }, [](){ return new PyCallBack_Teuchos_Object(); } ), "doc");
		cl.def( pybind11::init<int>(), pybind11::arg("tracebackModeIn") );

		cl.def( pybind11::init( [](const char * a0){ return new Teuchos::Object(a0); }, [](const char * a0){ return new PyCallBack_Teuchos_Object(a0); } ), "doc");
		cl.def( pybind11::init<const char *, int>(), pybind11::arg("label"), pybind11::arg("tracebackModeIn") );

		cl.def( pybind11::init( [](const std::string & a0){ return new Teuchos::Object(a0); }, [](const std::string & a0){ return new PyCallBack_Teuchos_Object(a0); } ), "doc");
		cl.def( pybind11::init<const std::string &, int>(), pybind11::arg("label"), pybind11::arg("tracebackModeIn") );

		cl.def( pybind11::init( [](PyCallBack_Teuchos_Object const &o){ return new PyCallBack_Teuchos_Object(o); } ) );
		cl.def( pybind11::init( [](Teuchos::Object const &o){ return new Teuchos::Object(o); } ) );
		cl.def("setLabel", (void (Teuchos::Object::*)(const char *)) &Teuchos::Object::setLabel, "C++: Teuchos::Object::setLabel(const char *) --> void", pybind11::arg("theLabel"));
		cl.def_static("setTracebackMode", (void (*)(int)) &Teuchos::Object::setTracebackMode, "Set the value of the Object error traceback report mode.\n\n TracebackMode controls whether or not traceback information is\n printed when run time integer errors are detected:\n\n <= 0 - No information report\n\n = 1 - Fatal (negative) values are reported\n\n >= 2 - All values (except zero) reported.\n\n \n Default is set to -1 when object is constructed.\n\nC++: Teuchos::Object::setTracebackMode(int) --> void", pybind11::arg("tracebackModeValue"));
		cl.def("label", (const char * (Teuchos::Object::*)() const) &Teuchos::Object::label, "Access the object's label (LEGACY; return std::string instead).\n\nC++: Teuchos::Object::label() const --> const char *", pybind11::return_value_policy::automatic);
		cl.def_static("getTracebackMode", (int (*)()) &Teuchos::Object::getTracebackMode, "Get the value of the Object error traceback report mode.\n\nC++: Teuchos::Object::getTracebackMode() --> int");
		cl.def("reportError", (int (Teuchos::Object::*)(const std::string, int) const) &Teuchos::Object::reportError, "Report an error with this Object.\n\nC++: Teuchos::Object::reportError(const std::string, int) const --> int", pybind11::arg("message"), pybind11::arg("errorCode"));
		cl.def("assign", (class Teuchos::Object & (Teuchos::Object::*)(const class Teuchos::Object &)) &Teuchos::Object::operator=, "C++: Teuchos::Object::operator=(const class Teuchos::Object &) --> class Teuchos::Object &", pybind11::return_value_policy::automatic, pybind11::arg(""));

		cl.def("__str__", [](Teuchos::Object const &o) -> std::string { std::ostringstream s; s << o; return s.str(); } );
	}
	{ // Teuchos::Flops file:Teuchos_Flops.hpp line:66
		pybind11::class_<Teuchos::Flops, Teuchos::RCP<Teuchos::Flops>> cl(M("Teuchos"), "Flops", "");
		cl.def( pybind11::init( [](){ return new Teuchos::Flops(); } ) );
		cl.def( pybind11::init( [](Teuchos::Flops const &o){ return new Teuchos::Flops(o); } ) );
		cl.def("flops", (double (Teuchos::Flops::*)() const) &Teuchos::Flops::flops, "Returns the number of floating point operations with  object and resets the count.\n\nC++: Teuchos::Flops::flops() const --> double");
		cl.def("resetFlops", (void (Teuchos::Flops::*)()) &Teuchos::Flops::resetFlops, "Resets the number of floating point operations to zero for  multi-std::vector.\n\nC++: Teuchos::Flops::resetFlops() --> void");
		cl.def("assign", (class Teuchos::Flops & (Teuchos::Flops::*)(const class Teuchos::Flops &)) &Teuchos::Flops::operator=, "C++: Teuchos::Flops::operator=(const class Teuchos::Flops &) --> class Teuchos::Flops &", pybind11::return_value_policy::automatic, pybind11::arg(""));
	}
	{ // Teuchos::CompObject file:Teuchos_CompObject.hpp line:65
		pybind11::class_<Teuchos::CompObject, Teuchos::RCP<Teuchos::CompObject>> cl(M("Teuchos"), "CompObject", "");
		cl.def( pybind11::init( [](){ return new Teuchos::CompObject(); } ) );
		cl.def( pybind11::init( [](Teuchos::CompObject const &o){ return new Teuchos::CompObject(o); } ) );
		cl.def("setFlopCounter", (void (Teuchos::CompObject::*)(const class Teuchos::Flops &)) &Teuchos::CompObject::setFlopCounter, "Set the internal Teuchos::Flops() pointer.\n\nC++: Teuchos::CompObject::setFlopCounter(const class Teuchos::Flops &) --> void", pybind11::arg("FlopCounter"));
		cl.def("setFlopCounter", (void (Teuchos::CompObject::*)(const class Teuchos::CompObject &)) &Teuchos::CompObject::setFlopCounter, "Set the internal Teuchos::Flops() pointer to the flop counter of another Teuchos::CompObject.\n\nC++: Teuchos::CompObject::setFlopCounter(const class Teuchos::CompObject &) --> void", pybind11::arg("compObject"));
		cl.def("unsetFlopCounter", (void (Teuchos::CompObject::*)()) &Teuchos::CompObject::unsetFlopCounter, "Set the internal Teuchos::Flops() pointer to 0 (no flops counted).\n\nC++: Teuchos::CompObject::unsetFlopCounter() --> void");
		cl.def("getFlopCounter", (class Teuchos::Flops * (Teuchos::CompObject::*)() const) &Teuchos::CompObject::getFlopCounter, "Get the pointer to the Teuchos::Flops() object associated with this object, returns 0 if none.\n\nC++: Teuchos::CompObject::getFlopCounter() const --> class Teuchos::Flops *", pybind11::return_value_policy::automatic);
		cl.def("resetFlops", (void (Teuchos::CompObject::*)() const) &Teuchos::CompObject::resetFlops, "Resets the number of floating point operations to zero for  multi-std::vector.\n\nC++: Teuchos::CompObject::resetFlops() const --> void");
		cl.def("getFlops", (double (Teuchos::CompObject::*)() const) &Teuchos::CompObject::getFlops, "Returns the number of floating point operations with  multi-std::vector.\n\nC++: Teuchos::CompObject::getFlops() const --> double");
		cl.def("updateFlops", (void (Teuchos::CompObject::*)(int) const) &Teuchos::CompObject::updateFlops, "Increment Flop count for  object\n\nC++: Teuchos::CompObject::updateFlops(int) const --> void", pybind11::arg("addflops"));
		cl.def("updateFlops", (void (Teuchos::CompObject::*)(long) const) &Teuchos::CompObject::updateFlops, "Increment Flop count for  object\n\nC++: Teuchos::CompObject::updateFlops(long) const --> void", pybind11::arg("addflops"));
		cl.def("updateFlops", (void (Teuchos::CompObject::*)(double) const) &Teuchos::CompObject::updateFlops, "Increment Flop count for  object\n\nC++: Teuchos::CompObject::updateFlops(double) const --> void", pybind11::arg("addflops"));
		cl.def("updateFlops", (void (Teuchos::CompObject::*)(float) const) &Teuchos::CompObject::updateFlops, "Increment Flop count for  object\n\nC++: Teuchos::CompObject::updateFlops(float) const --> void", pybind11::arg("addflops"));
		cl.def("assign", (class Teuchos::CompObject & (Teuchos::CompObject::*)(const class Teuchos::CompObject &)) &Teuchos::CompObject::operator=, "C++: Teuchos::CompObject::operator=(const class Teuchos::CompObject &) --> class Teuchos::CompObject &", pybind11::return_value_policy::automatic, pybind11::arg(""));
	}
}
