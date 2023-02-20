#include <PyROL_Teuchos_Custom.hpp>
#include <Teuchos_BLAS_types.hpp>
#include <Teuchos_CompObject.hpp>
#include <Teuchos_DataAccess.hpp>
#include <Teuchos_Flops.hpp>
#include <Teuchos_Object.hpp>
#include <Teuchos_OrdinalTraits.hpp>
#include <Teuchos_ScalarTraits.hpp>
#include <Teuchos_SerialDenseMatrix.hpp>
#include <Teuchos_SerialDenseVector.hpp>
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
	void print(std::ostream & a0) const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Teuchos::Object *>(this), "print");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Object::print(a0);
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

// Teuchos::SerialDenseMatrix file:Teuchos_SerialDenseMatrix.hpp line:67
struct PyCallBack_Teuchos_SerialDenseMatrix_int_double_t : public Teuchos::SerialDenseMatrix<int,double> {
	using Teuchos::SerialDenseMatrix<int,double>::SerialDenseMatrix;

	std::ostream & print(std::ostream & a0) const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Teuchos::SerialDenseMatrix<int,double> *>(this), "print");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<std::ostream &>::value) {
				static pybind11::detail::override_caster_t<std::ostream &> caster;
				return pybind11::detail::cast_ref<std::ostream &>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<std::ostream &>(std::move(o));
		}
		return SerialDenseMatrix::print(a0);
	}
};

// Teuchos::SerialDenseVector file:Teuchos_SerialDenseVector.hpp line:60
struct PyCallBack_Teuchos_SerialDenseVector_int_double_t : public Teuchos::SerialDenseVector<int,double> {
	using Teuchos::SerialDenseVector<int,double>::SerialDenseVector;

	std::ostream & print(std::ostream & a0) const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Teuchos::SerialDenseVector<int,double> *>(this), "print");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<std::ostream &>::value) {
				static pybind11::detail::override_caster_t<std::ostream &> caster;
				return pybind11::detail::cast_ref<std::ostream &>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<std::ostream &>(std::move(o));
		}
		return SerialDenseVector::print(a0);
	}
};

void bind_Teuchos_DataAccess(std::function< pybind11::module &(std::string const &namespace_) > &M)
{
	// Teuchos::DataAccess file:Teuchos_DataAccess.hpp line:60
	pybind11::enum_<Teuchos::DataAccess>(M("Teuchos"), "DataAccess", pybind11::arithmetic(), "If set to Copy, user data will be copied at construction.\n      If set to View, user data will be encapsulated and used throughout\n      the life of the object.", pybind11::module_local())
		.value("Copy", Teuchos::Copy)
		.value("View", Teuchos::View)
		.export_values();

;

	{ // Teuchos::Object file:Teuchos_Object.hpp line:68
		pybind11::class_<Teuchos::Object, Teuchos::RCP<Teuchos::Object>, PyCallBack_Teuchos_Object> cl(M("Teuchos"), "Object", "", pybind11::module_local());
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
		cl.def("print", (void (Teuchos::Object::*)(std::ostream &) const) &Teuchos::Object::print, "Print the object to the given output stream.\n\nC++: Teuchos::Object::print(std::ostream &) const --> void", pybind11::arg("os"));
		cl.def("reportError", (int (Teuchos::Object::*)(const std::string, int) const) &Teuchos::Object::reportError, "Report an error with this Object.\n\nC++: Teuchos::Object::reportError(const std::string, int) const --> int", pybind11::arg("message"), pybind11::arg("errorCode"));
		cl.def("assign", (class Teuchos::Object & (Teuchos::Object::*)(const class Teuchos::Object &)) &Teuchos::Object::operator=, "C++: Teuchos::Object::operator=(const class Teuchos::Object &) --> class Teuchos::Object &", pybind11::return_value_policy::automatic, pybind11::arg(""));

		cl.def("__str__", [](Teuchos::Object const &o) -> std::string { std::ostringstream s; Teuchos::operator<<(s, o); return s.str(); } );
	}
	{ // Teuchos::Flops file:Teuchos_Flops.hpp line:66
		pybind11::class_<Teuchos::Flops, Teuchos::RCP<Teuchos::Flops>> cl(M("Teuchos"), "Flops", "", pybind11::module_local());
		cl.def( pybind11::init( [](){ return new Teuchos::Flops(); } ) );
		cl.def( pybind11::init( [](Teuchos::Flops const &o){ return new Teuchos::Flops(o); } ) );
		cl.def("flops", (double (Teuchos::Flops::*)() const) &Teuchos::Flops::flops, "Returns the number of floating point operations with  object and resets the count.\n\nC++: Teuchos::Flops::flops() const --> double");
		cl.def("resetFlops", (void (Teuchos::Flops::*)()) &Teuchos::Flops::resetFlops, "Resets the number of floating point operations to zero for  multi-std::vector.\n\nC++: Teuchos::Flops::resetFlops() --> void");
		cl.def("assign", (class Teuchos::Flops & (Teuchos::Flops::*)(const class Teuchos::Flops &)) &Teuchos::Flops::operator=, "C++: Teuchos::Flops::operator=(const class Teuchos::Flops &) --> class Teuchos::Flops &", pybind11::return_value_policy::automatic, pybind11::arg(""));
	}
	{ // Teuchos::CompObject file:Teuchos_CompObject.hpp line:65
		pybind11::class_<Teuchos::CompObject, Teuchos::RCP<Teuchos::CompObject>> cl(M("Teuchos"), "CompObject", "", pybind11::module_local());
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
	// Teuchos::throwScalarTraitsNanInfError(const std::string &) file:Teuchos_ScalarTraits.hpp line:123
	M("Teuchos").def("throwScalarTraitsNanInfError", (void (*)(const std::string &)) &Teuchos::throwScalarTraitsNanInfError, "C++: Teuchos::throwScalarTraitsNanInfError(const std::string &) --> void", pybind11::arg("errMsg"));

	// Teuchos::generic_real_isnaninf(const float &) file:Teuchos_ScalarTraits.hpp line:127
	M("Teuchos").def("generic_real_isnaninf", (bool (*)(const float &)) &Teuchos::generic_real_isnaninf<float>, "C++: Teuchos::generic_real_isnaninf(const float &) --> bool", pybind11::arg("x"));

	// Teuchos::generic_real_isnaninf(const double &) file:Teuchos_ScalarTraits.hpp line:127
	M("Teuchos").def("generic_real_isnaninf", (bool (*)(const double &)) &Teuchos::generic_real_isnaninf<double>, "C++: Teuchos::generic_real_isnaninf(const double &) --> bool", pybind11::arg("x"));

	{ // Teuchos::ScalarTraits file:Teuchos_ScalarTraits.hpp line:158
		pybind11::class_<Teuchos::ScalarTraits<char>, Teuchos::RCP<Teuchos::ScalarTraits<char>>> cl(M("Teuchos"), "ScalarTraits_char_t", "", pybind11::module_local());
		cl.def( pybind11::init( [](){ return new Teuchos::ScalarTraits<char>(); } ) );
		cl.def_static("magnitude", (char (*)(char)) &Teuchos::ScalarTraits<char>::magnitude, "C++: Teuchos::ScalarTraits<char>::magnitude(char) --> char", pybind11::arg("a"));
		cl.def_static("zero", (char (*)()) &Teuchos::ScalarTraits<char>::zero, "C++: Teuchos::ScalarTraits<char>::zero() --> char");
		cl.def_static("one", (char (*)()) &Teuchos::ScalarTraits<char>::one, "C++: Teuchos::ScalarTraits<char>::one() --> char");
		cl.def_static("conjugate", (char (*)(char)) &Teuchos::ScalarTraits<char>::conjugate, "C++: Teuchos::ScalarTraits<char>::conjugate(char) --> char", pybind11::arg("x"));
		cl.def_static("real", (char (*)(char)) &Teuchos::ScalarTraits<char>::real, "C++: Teuchos::ScalarTraits<char>::real(char) --> char", pybind11::arg("x"));
		cl.def_static("imag", (char (*)(char)) &Teuchos::ScalarTraits<char>::imag, "C++: Teuchos::ScalarTraits<char>::imag(char) --> char", pybind11::arg(""));
		cl.def_static("isnaninf", (bool (*)(char)) &Teuchos::ScalarTraits<char>::isnaninf, "C++: Teuchos::ScalarTraits<char>::isnaninf(char) --> bool", pybind11::arg(""));
		cl.def_static("seedrandom", (void (*)(unsigned int)) &Teuchos::ScalarTraits<char>::seedrandom, "C++: Teuchos::ScalarTraits<char>::seedrandom(unsigned int) --> void", pybind11::arg("s"));
		cl.def_static("random", (char (*)()) &Teuchos::ScalarTraits<char>::random, "C++: Teuchos::ScalarTraits<char>::random() --> char");
		cl.def_static("name", (std::string (*)()) &Teuchos::ScalarTraits<char>::name, "C++: Teuchos::ScalarTraits<char>::name() --> std::string");
		cl.def_static("squareroot", (char (*)(char)) &Teuchos::ScalarTraits<char>::squareroot, "C++: Teuchos::ScalarTraits<char>::squareroot(char) --> char", pybind11::arg("x"));
		cl.def_static("pow", (char (*)(char, char)) &Teuchos::ScalarTraits<char>::pow, "C++: Teuchos::ScalarTraits<char>::pow(char, char) --> char", pybind11::arg("x"), pybind11::arg("y"));
		cl.def_static("log", (char (*)(char)) &Teuchos::ScalarTraits<char>::log, "C++: Teuchos::ScalarTraits<char>::log(char) --> char", pybind11::arg("x"));
		cl.def_static("log10", (char (*)(char)) &Teuchos::ScalarTraits<char>::log10, "C++: Teuchos::ScalarTraits<char>::log10(char) --> char", pybind11::arg("x"));
	}
	{ // Teuchos::ScalarTraits file:Teuchos_ScalarTraits.hpp line:195
		pybind11::class_<Teuchos::ScalarTraits<short>, Teuchos::RCP<Teuchos::ScalarTraits<short>>> cl(M("Teuchos"), "ScalarTraits_short_t", "", pybind11::module_local());
		cl.def( pybind11::init( [](){ return new Teuchos::ScalarTraits<short>(); } ) );
		cl.def_static("magnitude", (short (*)(short)) &Teuchos::ScalarTraits<short>::magnitude, "C++: Teuchos::ScalarTraits<short>::magnitude(short) --> short", pybind11::arg("a"));
		cl.def_static("zero", (short (*)()) &Teuchos::ScalarTraits<short>::zero, "C++: Teuchos::ScalarTraits<short>::zero() --> short");
		cl.def_static("one", (short (*)()) &Teuchos::ScalarTraits<short>::one, "C++: Teuchos::ScalarTraits<short>::one() --> short");
		cl.def_static("conjugate", (short (*)(short)) &Teuchos::ScalarTraits<short>::conjugate, "C++: Teuchos::ScalarTraits<short>::conjugate(short) --> short", pybind11::arg("x"));
		cl.def_static("real", (short (*)(short)) &Teuchos::ScalarTraits<short>::real, "C++: Teuchos::ScalarTraits<short>::real(short) --> short", pybind11::arg("x"));
		cl.def_static("imag", (short (*)(short)) &Teuchos::ScalarTraits<short>::imag, "C++: Teuchos::ScalarTraits<short>::imag(short) --> short", pybind11::arg(""));
		cl.def_static("isnaninf", (bool (*)(short)) &Teuchos::ScalarTraits<short>::isnaninf, "C++: Teuchos::ScalarTraits<short>::isnaninf(short) --> bool", pybind11::arg(""));
		cl.def_static("seedrandom", (void (*)(unsigned int)) &Teuchos::ScalarTraits<short>::seedrandom, "C++: Teuchos::ScalarTraits<short>::seedrandom(unsigned int) --> void", pybind11::arg("s"));
		cl.def_static("random", (short (*)()) &Teuchos::ScalarTraits<short>::random, "C++: Teuchos::ScalarTraits<short>::random() --> short");
		cl.def_static("name", (std::string (*)()) &Teuchos::ScalarTraits<short>::name, "C++: Teuchos::ScalarTraits<short>::name() --> std::string");
		cl.def_static("squareroot", (short (*)(short)) &Teuchos::ScalarTraits<short>::squareroot, "C++: Teuchos::ScalarTraits<short>::squareroot(short) --> short", pybind11::arg("x"));
		cl.def_static("pow", (short (*)(short, short)) &Teuchos::ScalarTraits<short>::pow, "C++: Teuchos::ScalarTraits<short>::pow(short, short) --> short", pybind11::arg("x"), pybind11::arg("y"));
		cl.def_static("log", (short (*)(short)) &Teuchos::ScalarTraits<short>::log, "C++: Teuchos::ScalarTraits<short>::log(short) --> short", pybind11::arg("x"));
		cl.def_static("log10", (short (*)(short)) &Teuchos::ScalarTraits<short>::log10, "C++: Teuchos::ScalarTraits<short>::log10(short) --> short", pybind11::arg("x"));
	}
	{ // Teuchos::ScalarTraits file:Teuchos_ScalarTraits.hpp line:231
		pybind11::class_<Teuchos::ScalarTraits<unsigned short>, Teuchos::RCP<Teuchos::ScalarTraits<unsigned short>>> cl(M("Teuchos"), "ScalarTraits_unsigned_short_t", "", pybind11::module_local());
		cl.def( pybind11::init( [](){ return new Teuchos::ScalarTraits<unsigned short>(); } ) );
		cl.def_static("magnitude", (unsigned short (*)(unsigned short)) &Teuchos::ScalarTraits<unsigned short>::magnitude, "C++: Teuchos::ScalarTraits<unsigned short>::magnitude(unsigned short) --> unsigned short", pybind11::arg("a"));
		cl.def_static("zero", (unsigned short (*)()) &Teuchos::ScalarTraits<unsigned short>::zero, "C++: Teuchos::ScalarTraits<unsigned short>::zero() --> unsigned short");
		cl.def_static("one", (unsigned short (*)()) &Teuchos::ScalarTraits<unsigned short>::one, "C++: Teuchos::ScalarTraits<unsigned short>::one() --> unsigned short");
		cl.def_static("conjugate", (unsigned short (*)(unsigned short)) &Teuchos::ScalarTraits<unsigned short>::conjugate, "C++: Teuchos::ScalarTraits<unsigned short>::conjugate(unsigned short) --> unsigned short", pybind11::arg("x"));
		cl.def_static("real", (unsigned short (*)(unsigned short)) &Teuchos::ScalarTraits<unsigned short>::real, "C++: Teuchos::ScalarTraits<unsigned short>::real(unsigned short) --> unsigned short", pybind11::arg("x"));
		cl.def_static("imag", (unsigned short (*)(unsigned short)) &Teuchos::ScalarTraits<unsigned short>::imag, "C++: Teuchos::ScalarTraits<unsigned short>::imag(unsigned short) --> unsigned short", pybind11::arg(""));
		cl.def_static("isnaninf", (bool (*)(unsigned short)) &Teuchos::ScalarTraits<unsigned short>::isnaninf, "C++: Teuchos::ScalarTraits<unsigned short>::isnaninf(unsigned short) --> bool", pybind11::arg(""));
		cl.def_static("seedrandom", (void (*)(unsigned int)) &Teuchos::ScalarTraits<unsigned short>::seedrandom, "C++: Teuchos::ScalarTraits<unsigned short>::seedrandom(unsigned int) --> void", pybind11::arg("s"));
		cl.def_static("random", (unsigned short (*)()) &Teuchos::ScalarTraits<unsigned short>::random, "C++: Teuchos::ScalarTraits<unsigned short>::random() --> unsigned short");
		cl.def_static("name", (std::string (*)()) &Teuchos::ScalarTraits<unsigned short>::name, "C++: Teuchos::ScalarTraits<unsigned short>::name() --> std::string");
		cl.def_static("squareroot", (unsigned short (*)(unsigned short)) &Teuchos::ScalarTraits<unsigned short>::squareroot, "C++: Teuchos::ScalarTraits<unsigned short>::squareroot(unsigned short) --> unsigned short", pybind11::arg("x"));
		cl.def_static("pow", (unsigned short (*)(unsigned short, unsigned short)) &Teuchos::ScalarTraits<unsigned short>::pow, "C++: Teuchos::ScalarTraits<unsigned short>::pow(unsigned short, unsigned short) --> unsigned short", pybind11::arg("x"), pybind11::arg("y"));
		cl.def_static("log", (unsigned short (*)(unsigned short)) &Teuchos::ScalarTraits<unsigned short>::log, "C++: Teuchos::ScalarTraits<unsigned short>::log(unsigned short) --> unsigned short", pybind11::arg("x"));
		cl.def_static("log10", (unsigned short (*)(unsigned short)) &Teuchos::ScalarTraits<unsigned short>::log10, "C++: Teuchos::ScalarTraits<unsigned short>::log10(unsigned short) --> unsigned short", pybind11::arg("x"));
	}
	{ // Teuchos::ScalarTraits file:Teuchos_ScalarTraits.hpp line:268
		pybind11::class_<Teuchos::ScalarTraits<int>, Teuchos::RCP<Teuchos::ScalarTraits<int>>> cl(M("Teuchos"), "ScalarTraits_int_t", "", pybind11::module_local());
		cl.def( pybind11::init( [](){ return new Teuchos::ScalarTraits<int>(); } ) );
		cl.def_static("magnitude", (int (*)(int)) &Teuchos::ScalarTraits<int>::magnitude, "C++: Teuchos::ScalarTraits<int>::magnitude(int) --> int", pybind11::arg("a"));
		cl.def_static("zero", (int (*)()) &Teuchos::ScalarTraits<int>::zero, "C++: Teuchos::ScalarTraits<int>::zero() --> int");
		cl.def_static("one", (int (*)()) &Teuchos::ScalarTraits<int>::one, "C++: Teuchos::ScalarTraits<int>::one() --> int");
		cl.def_static("conjugate", (int (*)(int)) &Teuchos::ScalarTraits<int>::conjugate, "C++: Teuchos::ScalarTraits<int>::conjugate(int) --> int", pybind11::arg("x"));
		cl.def_static("real", (int (*)(int)) &Teuchos::ScalarTraits<int>::real, "C++: Teuchos::ScalarTraits<int>::real(int) --> int", pybind11::arg("x"));
		cl.def_static("imag", (int (*)(int)) &Teuchos::ScalarTraits<int>::imag, "C++: Teuchos::ScalarTraits<int>::imag(int) --> int", pybind11::arg(""));
		cl.def_static("isnaninf", (bool (*)(int)) &Teuchos::ScalarTraits<int>::isnaninf, "C++: Teuchos::ScalarTraits<int>::isnaninf(int) --> bool", pybind11::arg(""));
		cl.def_static("seedrandom", (void (*)(unsigned int)) &Teuchos::ScalarTraits<int>::seedrandom, "C++: Teuchos::ScalarTraits<int>::seedrandom(unsigned int) --> void", pybind11::arg("s"));
		cl.def_static("random", (int (*)()) &Teuchos::ScalarTraits<int>::random, "C++: Teuchos::ScalarTraits<int>::random() --> int");
		cl.def_static("name", (std::string (*)()) &Teuchos::ScalarTraits<int>::name, "C++: Teuchos::ScalarTraits<int>::name() --> std::string");
		cl.def_static("squareroot", (int (*)(int)) &Teuchos::ScalarTraits<int>::squareroot, "C++: Teuchos::ScalarTraits<int>::squareroot(int) --> int", pybind11::arg("x"));
		cl.def_static("pow", (int (*)(int, int)) &Teuchos::ScalarTraits<int>::pow, "C++: Teuchos::ScalarTraits<int>::pow(int, int) --> int", pybind11::arg("x"), pybind11::arg("y"));
		cl.def_static("log", (int (*)(int)) &Teuchos::ScalarTraits<int>::log, "C++: Teuchos::ScalarTraits<int>::log(int) --> int", pybind11::arg("x"));
		cl.def_static("log10", (int (*)(int)) &Teuchos::ScalarTraits<int>::log10, "C++: Teuchos::ScalarTraits<int>::log10(int) --> int", pybind11::arg("x"));
	}
	{ // Teuchos::ScalarTraits file:Teuchos_ScalarTraits.hpp line:305
		pybind11::class_<Teuchos::ScalarTraits<unsigned int>, Teuchos::RCP<Teuchos::ScalarTraits<unsigned int>>> cl(M("Teuchos"), "ScalarTraits_unsigned_int_t", "", pybind11::module_local());
		cl.def( pybind11::init( [](){ return new Teuchos::ScalarTraits<unsigned int>(); } ) );
		cl.def_static("magnitude", (unsigned int (*)(unsigned int)) &Teuchos::ScalarTraits<unsigned int>::magnitude, "C++: Teuchos::ScalarTraits<unsigned int>::magnitude(unsigned int) --> unsigned int", pybind11::arg("a"));
		cl.def_static("zero", (unsigned int (*)()) &Teuchos::ScalarTraits<unsigned int>::zero, "C++: Teuchos::ScalarTraits<unsigned int>::zero() --> unsigned int");
		cl.def_static("one", (unsigned int (*)()) &Teuchos::ScalarTraits<unsigned int>::one, "C++: Teuchos::ScalarTraits<unsigned int>::one() --> unsigned int");
		cl.def_static("conjugate", (unsigned int (*)(unsigned int)) &Teuchos::ScalarTraits<unsigned int>::conjugate, "C++: Teuchos::ScalarTraits<unsigned int>::conjugate(unsigned int) --> unsigned int", pybind11::arg("x"));
		cl.def_static("real", (unsigned int (*)(unsigned int)) &Teuchos::ScalarTraits<unsigned int>::real, "C++: Teuchos::ScalarTraits<unsigned int>::real(unsigned int) --> unsigned int", pybind11::arg("x"));
		cl.def_static("imag", (unsigned int (*)(unsigned int)) &Teuchos::ScalarTraits<unsigned int>::imag, "C++: Teuchos::ScalarTraits<unsigned int>::imag(unsigned int) --> unsigned int", pybind11::arg(""));
		cl.def_static("isnaninf", (bool (*)(unsigned int)) &Teuchos::ScalarTraits<unsigned int>::isnaninf, "C++: Teuchos::ScalarTraits<unsigned int>::isnaninf(unsigned int) --> bool", pybind11::arg(""));
		cl.def_static("seedrandom", (void (*)(unsigned int)) &Teuchos::ScalarTraits<unsigned int>::seedrandom, "C++: Teuchos::ScalarTraits<unsigned int>::seedrandom(unsigned int) --> void", pybind11::arg("s"));
		cl.def_static("random", (unsigned int (*)()) &Teuchos::ScalarTraits<unsigned int>::random, "C++: Teuchos::ScalarTraits<unsigned int>::random() --> unsigned int");
		cl.def_static("name", (std::string (*)()) &Teuchos::ScalarTraits<unsigned int>::name, "C++: Teuchos::ScalarTraits<unsigned int>::name() --> std::string");
		cl.def_static("squareroot", (unsigned int (*)(unsigned int)) &Teuchos::ScalarTraits<unsigned int>::squareroot, "C++: Teuchos::ScalarTraits<unsigned int>::squareroot(unsigned int) --> unsigned int", pybind11::arg("x"));
		cl.def_static("pow", (unsigned int (*)(unsigned int, unsigned int)) &Teuchos::ScalarTraits<unsigned int>::pow, "C++: Teuchos::ScalarTraits<unsigned int>::pow(unsigned int, unsigned int) --> unsigned int", pybind11::arg("x"), pybind11::arg("y"));
		cl.def_static("log", (unsigned int (*)(unsigned int)) &Teuchos::ScalarTraits<unsigned int>::log, "C++: Teuchos::ScalarTraits<unsigned int>::log(unsigned int) --> unsigned int", pybind11::arg("x"));
		cl.def_static("log10", (unsigned int (*)(unsigned int)) &Teuchos::ScalarTraits<unsigned int>::log10, "C++: Teuchos::ScalarTraits<unsigned int>::log10(unsigned int) --> unsigned int", pybind11::arg("x"));
	}
	{ // Teuchos::ScalarTraits file:Teuchos_ScalarTraits.hpp line:342
		pybind11::class_<Teuchos::ScalarTraits<long>, Teuchos::RCP<Teuchos::ScalarTraits<long>>> cl(M("Teuchos"), "ScalarTraits_long_t", "", pybind11::module_local());
		cl.def( pybind11::init( [](){ return new Teuchos::ScalarTraits<long>(); } ) );
		cl.def_static("magnitude", (long (*)(long)) &Teuchos::ScalarTraits<long>::magnitude, "C++: Teuchos::ScalarTraits<long>::magnitude(long) --> long", pybind11::arg("a"));
		cl.def_static("zero", (long (*)()) &Teuchos::ScalarTraits<long>::zero, "C++: Teuchos::ScalarTraits<long>::zero() --> long");
		cl.def_static("one", (long (*)()) &Teuchos::ScalarTraits<long>::one, "C++: Teuchos::ScalarTraits<long>::one() --> long");
		cl.def_static("conjugate", (long (*)(long)) &Teuchos::ScalarTraits<long>::conjugate, "C++: Teuchos::ScalarTraits<long>::conjugate(long) --> long", pybind11::arg("x"));
		cl.def_static("real", (long (*)(long)) &Teuchos::ScalarTraits<long>::real, "C++: Teuchos::ScalarTraits<long>::real(long) --> long", pybind11::arg("x"));
		cl.def_static("imag", (long (*)(long)) &Teuchos::ScalarTraits<long>::imag, "C++: Teuchos::ScalarTraits<long>::imag(long) --> long", pybind11::arg(""));
		cl.def_static("isnaninf", (bool (*)(long)) &Teuchos::ScalarTraits<long>::isnaninf, "C++: Teuchos::ScalarTraits<long>::isnaninf(long) --> bool", pybind11::arg(""));
		cl.def_static("seedrandom", (void (*)(unsigned int)) &Teuchos::ScalarTraits<long>::seedrandom, "C++: Teuchos::ScalarTraits<long>::seedrandom(unsigned int) --> void", pybind11::arg("s"));
		cl.def_static("random", (long (*)()) &Teuchos::ScalarTraits<long>::random, "C++: Teuchos::ScalarTraits<long>::random() --> long");
		cl.def_static("name", (std::string (*)()) &Teuchos::ScalarTraits<long>::name, "C++: Teuchos::ScalarTraits<long>::name() --> std::string");
		cl.def_static("squareroot", (long (*)(long)) &Teuchos::ScalarTraits<long>::squareroot, "C++: Teuchos::ScalarTraits<long>::squareroot(long) --> long", pybind11::arg("x"));
		cl.def_static("pow", (long (*)(long, long)) &Teuchos::ScalarTraits<long>::pow, "C++: Teuchos::ScalarTraits<long>::pow(long, long) --> long", pybind11::arg("x"), pybind11::arg("y"));
		cl.def_static("log", (long (*)(long)) &Teuchos::ScalarTraits<long>::log, "C++: Teuchos::ScalarTraits<long>::log(long) --> long", pybind11::arg("x"));
		cl.def_static("log10", (long (*)(long)) &Teuchos::ScalarTraits<long>::log10, "C++: Teuchos::ScalarTraits<long>::log10(long) --> long", pybind11::arg("x"));
	}
	{ // Teuchos::ScalarTraits file:Teuchos_ScalarTraits.hpp line:381
		pybind11::class_<Teuchos::ScalarTraits<unsigned long>, Teuchos::RCP<Teuchos::ScalarTraits<unsigned long>>> cl(M("Teuchos"), "ScalarTraits_unsigned_long_t", "", pybind11::module_local());
		cl.def( pybind11::init( [](){ return new Teuchos::ScalarTraits<unsigned long>(); } ) );
		cl.def_static("magnitude", (unsigned long (*)(unsigned long)) &Teuchos::ScalarTraits<unsigned long>::magnitude, "C++: Teuchos::ScalarTraits<unsigned long>::magnitude(unsigned long) --> unsigned long", pybind11::arg("a"));
		cl.def_static("zero", (unsigned long (*)()) &Teuchos::ScalarTraits<unsigned long>::zero, "C++: Teuchos::ScalarTraits<unsigned long>::zero() --> unsigned long");
		cl.def_static("one", (unsigned long (*)()) &Teuchos::ScalarTraits<unsigned long>::one, "C++: Teuchos::ScalarTraits<unsigned long>::one() --> unsigned long");
		cl.def_static("conjugate", (unsigned long (*)(unsigned long)) &Teuchos::ScalarTraits<unsigned long>::conjugate, "C++: Teuchos::ScalarTraits<unsigned long>::conjugate(unsigned long) --> unsigned long", pybind11::arg("x"));
		cl.def_static("real", (unsigned long (*)(unsigned long)) &Teuchos::ScalarTraits<unsigned long>::real, "C++: Teuchos::ScalarTraits<unsigned long>::real(unsigned long) --> unsigned long", pybind11::arg("x"));
		cl.def_static("imag", (unsigned long (*)(unsigned long)) &Teuchos::ScalarTraits<unsigned long>::imag, "C++: Teuchos::ScalarTraits<unsigned long>::imag(unsigned long) --> unsigned long", pybind11::arg(""));
		cl.def_static("isnaninf", (bool (*)(unsigned long)) &Teuchos::ScalarTraits<unsigned long>::isnaninf, "C++: Teuchos::ScalarTraits<unsigned long>::isnaninf(unsigned long) --> bool", pybind11::arg(""));
		cl.def_static("seedrandom", (void (*)(unsigned int)) &Teuchos::ScalarTraits<unsigned long>::seedrandom, "C++: Teuchos::ScalarTraits<unsigned long>::seedrandom(unsigned int) --> void", pybind11::arg("s"));
		cl.def_static("random", (unsigned long (*)()) &Teuchos::ScalarTraits<unsigned long>::random, "C++: Teuchos::ScalarTraits<unsigned long>::random() --> unsigned long");
		cl.def_static("name", (std::string (*)()) &Teuchos::ScalarTraits<unsigned long>::name, "C++: Teuchos::ScalarTraits<unsigned long>::name() --> std::string");
		cl.def_static("squareroot", (unsigned long (*)(unsigned long)) &Teuchos::ScalarTraits<unsigned long>::squareroot, "C++: Teuchos::ScalarTraits<unsigned long>::squareroot(unsigned long) --> unsigned long", pybind11::arg("x"));
		cl.def_static("pow", (unsigned long (*)(unsigned long, unsigned long)) &Teuchos::ScalarTraits<unsigned long>::pow, "C++: Teuchos::ScalarTraits<unsigned long>::pow(unsigned long, unsigned long) --> unsigned long", pybind11::arg("x"), pybind11::arg("y"));
		cl.def_static("log", (unsigned long (*)(unsigned long)) &Teuchos::ScalarTraits<unsigned long>::log, "C++: Teuchos::ScalarTraits<unsigned long>::log(unsigned long) --> unsigned long", pybind11::arg("x"));
		cl.def_static("log10", (unsigned long (*)(unsigned long)) &Teuchos::ScalarTraits<unsigned long>::log10, "C++: Teuchos::ScalarTraits<unsigned long>::log10(unsigned long) --> unsigned long", pybind11::arg("x"));
	}
	{ // Teuchos::ScalarTraits file:Teuchos_ScalarTraits.hpp line:420
		pybind11::class_<Teuchos::ScalarTraits<long long>, Teuchos::RCP<Teuchos::ScalarTraits<long long>>> cl(M("Teuchos"), "ScalarTraits_long_long_t", "", pybind11::module_local());
		cl.def( pybind11::init( [](){ return new Teuchos::ScalarTraits<long long>(); } ) );
		cl.def_static("magnitude", (long long (*)(long long)) &Teuchos::ScalarTraits<long long>::magnitude, "C++: Teuchos::ScalarTraits<long long>::magnitude(long long) --> long long", pybind11::arg("a"));
		cl.def_static("zero", (long long (*)()) &Teuchos::ScalarTraits<long long>::zero, "C++: Teuchos::ScalarTraits<long long>::zero() --> long long");
		cl.def_static("one", (long long (*)()) &Teuchos::ScalarTraits<long long>::one, "C++: Teuchos::ScalarTraits<long long>::one() --> long long");
		cl.def_static("conjugate", (long long (*)(long long)) &Teuchos::ScalarTraits<long long>::conjugate, "C++: Teuchos::ScalarTraits<long long>::conjugate(long long) --> long long", pybind11::arg("x"));
		cl.def_static("real", (long long (*)(long long)) &Teuchos::ScalarTraits<long long>::real, "C++: Teuchos::ScalarTraits<long long>::real(long long) --> long long", pybind11::arg("x"));
		cl.def_static("imag", (long long (*)(long long)) &Teuchos::ScalarTraits<long long>::imag, "C++: Teuchos::ScalarTraits<long long>::imag(long long) --> long long", pybind11::arg(""));
		cl.def_static("isnaninf", (bool (*)(long long)) &Teuchos::ScalarTraits<long long>::isnaninf, "C++: Teuchos::ScalarTraits<long long>::isnaninf(long long) --> bool", pybind11::arg(""));
		cl.def_static("seedrandom", (void (*)(unsigned int)) &Teuchos::ScalarTraits<long long>::seedrandom, "C++: Teuchos::ScalarTraits<long long>::seedrandom(unsigned int) --> void", pybind11::arg("s"));
		cl.def_static("random", (long long (*)()) &Teuchos::ScalarTraits<long long>::random, "C++: Teuchos::ScalarTraits<long long>::random() --> long long");
		cl.def_static("name", (std::string (*)()) &Teuchos::ScalarTraits<long long>::name, "C++: Teuchos::ScalarTraits<long long>::name() --> std::string");
		cl.def_static("squareroot", (long long (*)(long long)) &Teuchos::ScalarTraits<long long>::squareroot, "C++: Teuchos::ScalarTraits<long long>::squareroot(long long) --> long long", pybind11::arg("x"));
		cl.def_static("pow", (long long (*)(long long, long long)) &Teuchos::ScalarTraits<long long>::pow, "C++: Teuchos::ScalarTraits<long long>::pow(long long, long long) --> long long", pybind11::arg("x"), pybind11::arg("y"));
		cl.def_static("log", (long long (*)(long long)) &Teuchos::ScalarTraits<long long>::log, "C++: Teuchos::ScalarTraits<long long>::log(long long) --> long long", pybind11::arg("x"));
		cl.def_static("log10", (long long (*)(long long)) &Teuchos::ScalarTraits<long long>::log10, "C++: Teuchos::ScalarTraits<long long>::log10(long long) --> long long", pybind11::arg("x"));
	}
	{ // Teuchos::ScalarTraits file:Teuchos_ScalarTraits.hpp line:458
		pybind11::class_<Teuchos::ScalarTraits<unsigned long long>, Teuchos::RCP<Teuchos::ScalarTraits<unsigned long long>>> cl(M("Teuchos"), "ScalarTraits_unsigned_long_long_t", "", pybind11::module_local());
		cl.def( pybind11::init( [](){ return new Teuchos::ScalarTraits<unsigned long long>(); } ) );
		cl.def_static("magnitude", (unsigned long long (*)(unsigned long long)) &Teuchos::ScalarTraits<unsigned long long>::magnitude, "C++: Teuchos::ScalarTraits<unsigned long long>::magnitude(unsigned long long) --> unsigned long long", pybind11::arg("a"));
		cl.def_static("zero", (unsigned long long (*)()) &Teuchos::ScalarTraits<unsigned long long>::zero, "C++: Teuchos::ScalarTraits<unsigned long long>::zero() --> unsigned long long");
		cl.def_static("one", (unsigned long long (*)()) &Teuchos::ScalarTraits<unsigned long long>::one, "C++: Teuchos::ScalarTraits<unsigned long long>::one() --> unsigned long long");
		cl.def_static("conjugate", (unsigned long long (*)(unsigned long long)) &Teuchos::ScalarTraits<unsigned long long>::conjugate, "C++: Teuchos::ScalarTraits<unsigned long long>::conjugate(unsigned long long) --> unsigned long long", pybind11::arg("x"));
		cl.def_static("real", (unsigned long long (*)(unsigned long long)) &Teuchos::ScalarTraits<unsigned long long>::real, "C++: Teuchos::ScalarTraits<unsigned long long>::real(unsigned long long) --> unsigned long long", pybind11::arg("x"));
		cl.def_static("imag", (unsigned long long (*)(unsigned long long)) &Teuchos::ScalarTraits<unsigned long long>::imag, "C++: Teuchos::ScalarTraits<unsigned long long>::imag(unsigned long long) --> unsigned long long", pybind11::arg(""));
		cl.def_static("isnaninf", (bool (*)(unsigned long long)) &Teuchos::ScalarTraits<unsigned long long>::isnaninf, "C++: Teuchos::ScalarTraits<unsigned long long>::isnaninf(unsigned long long) --> bool", pybind11::arg(""));
		cl.def_static("seedrandom", (void (*)(unsigned int)) &Teuchos::ScalarTraits<unsigned long long>::seedrandom, "C++: Teuchos::ScalarTraits<unsigned long long>::seedrandom(unsigned int) --> void", pybind11::arg("s"));
		cl.def_static("random", (unsigned long long (*)()) &Teuchos::ScalarTraits<unsigned long long>::random, "C++: Teuchos::ScalarTraits<unsigned long long>::random() --> unsigned long long");
		cl.def_static("name", (std::string (*)()) &Teuchos::ScalarTraits<unsigned long long>::name, "C++: Teuchos::ScalarTraits<unsigned long long>::name() --> std::string");
		cl.def_static("squareroot", (unsigned long long (*)(unsigned long long)) &Teuchos::ScalarTraits<unsigned long long>::squareroot, "C++: Teuchos::ScalarTraits<unsigned long long>::squareroot(unsigned long long) --> unsigned long long", pybind11::arg("x"));
		cl.def_static("pow", (unsigned long long (*)(unsigned long long, unsigned long long)) &Teuchos::ScalarTraits<unsigned long long>::pow, "C++: Teuchos::ScalarTraits<unsigned long long>::pow(unsigned long long, unsigned long long) --> unsigned long long", pybind11::arg("x"), pybind11::arg("y"));
		cl.def_static("log", (unsigned long long (*)(unsigned long long)) &Teuchos::ScalarTraits<unsigned long long>::log, "C++: Teuchos::ScalarTraits<unsigned long long>::log(unsigned long long) --> unsigned long long", pybind11::arg("x"));
		cl.def_static("log10", (unsigned long long (*)(unsigned long long)) &Teuchos::ScalarTraits<unsigned long long>::log10, "C++: Teuchos::ScalarTraits<unsigned long long>::log10(unsigned long long) --> unsigned long long", pybind11::arg("x"));
	}
	{ // Teuchos::ScalarTraits file:Teuchos_ScalarTraits.hpp line:581
		pybind11::class_<Teuchos::ScalarTraits<float>, Teuchos::RCP<Teuchos::ScalarTraits<float>>> cl(M("Teuchos"), "ScalarTraits_float_t", "", pybind11::module_local());
		cl.def( pybind11::init( [](){ return new Teuchos::ScalarTraits<float>(); } ) );
		cl.def_static("eps", (float (*)()) &Teuchos::ScalarTraits<float>::eps, "C++: Teuchos::ScalarTraits<float>::eps() --> float");
		cl.def_static("sfmin", (float (*)()) &Teuchos::ScalarTraits<float>::sfmin, "C++: Teuchos::ScalarTraits<float>::sfmin() --> float");
		cl.def_static("base", (float (*)()) &Teuchos::ScalarTraits<float>::base, "C++: Teuchos::ScalarTraits<float>::base() --> float");
		cl.def_static("prec", (float (*)()) &Teuchos::ScalarTraits<float>::prec, "C++: Teuchos::ScalarTraits<float>::prec() --> float");
		cl.def_static("t", (float (*)()) &Teuchos::ScalarTraits<float>::t, "C++: Teuchos::ScalarTraits<float>::t() --> float");
		cl.def_static("rnd", (float (*)()) &Teuchos::ScalarTraits<float>::rnd, "C++: Teuchos::ScalarTraits<float>::rnd() --> float");
		cl.def_static("emin", (float (*)()) &Teuchos::ScalarTraits<float>::emin, "C++: Teuchos::ScalarTraits<float>::emin() --> float");
		cl.def_static("rmin", (float (*)()) &Teuchos::ScalarTraits<float>::rmin, "C++: Teuchos::ScalarTraits<float>::rmin() --> float");
		cl.def_static("emax", (float (*)()) &Teuchos::ScalarTraits<float>::emax, "C++: Teuchos::ScalarTraits<float>::emax() --> float");
		cl.def_static("rmax", (float (*)()) &Teuchos::ScalarTraits<float>::rmax, "C++: Teuchos::ScalarTraits<float>::rmax() --> float");
		cl.def_static("magnitude", (float (*)(float)) &Teuchos::ScalarTraits<float>::magnitude, "C++: Teuchos::ScalarTraits<float>::magnitude(float) --> float", pybind11::arg("a"));
		cl.def_static("zero", (float (*)()) &Teuchos::ScalarTraits<float>::zero, "C++: Teuchos::ScalarTraits<float>::zero() --> float");
		cl.def_static("one", (float (*)()) &Teuchos::ScalarTraits<float>::one, "C++: Teuchos::ScalarTraits<float>::one() --> float");
		cl.def_static("conjugate", (float (*)(float)) &Teuchos::ScalarTraits<float>::conjugate, "C++: Teuchos::ScalarTraits<float>::conjugate(float) --> float", pybind11::arg("x"));
		cl.def_static("real", (float (*)(float)) &Teuchos::ScalarTraits<float>::real, "C++: Teuchos::ScalarTraits<float>::real(float) --> float", pybind11::arg("x"));
		cl.def_static("imag", (float (*)(float)) &Teuchos::ScalarTraits<float>::imag, "C++: Teuchos::ScalarTraits<float>::imag(float) --> float", pybind11::arg(""));
		cl.def_static("nan", (float (*)()) &Teuchos::ScalarTraits<float>::nan, "C++: Teuchos::ScalarTraits<float>::nan() --> float");
		cl.def_static("isnaninf", (bool (*)(float)) &Teuchos::ScalarTraits<float>::isnaninf, "C++: Teuchos::ScalarTraits<float>::isnaninf(float) --> bool", pybind11::arg("x"));
		cl.def_static("seedrandom", (void (*)(unsigned int)) &Teuchos::ScalarTraits<float>::seedrandom, "C++: Teuchos::ScalarTraits<float>::seedrandom(unsigned int) --> void", pybind11::arg("s"));
		cl.def_static("random", (float (*)()) &Teuchos::ScalarTraits<float>::random, "C++: Teuchos::ScalarTraits<float>::random() --> float");
		cl.def_static("name", (std::string (*)()) &Teuchos::ScalarTraits<float>::name, "C++: Teuchos::ScalarTraits<float>::name() --> std::string");
		cl.def_static("squareroot", (float (*)(float)) &Teuchos::ScalarTraits<float>::squareroot, "C++: Teuchos::ScalarTraits<float>::squareroot(float) --> float", pybind11::arg("x"));
		cl.def_static("pow", (float (*)(float, float)) &Teuchos::ScalarTraits<float>::pow, "C++: Teuchos::ScalarTraits<float>::pow(float, float) --> float", pybind11::arg("x"), pybind11::arg("y"));
		cl.def_static("pi", (float (*)()) &Teuchos::ScalarTraits<float>::pi, "C++: Teuchos::ScalarTraits<float>::pi() --> float");
		cl.def_static("log", (float (*)(float)) &Teuchos::ScalarTraits<float>::log, "C++: Teuchos::ScalarTraits<float>::log(float) --> float", pybind11::arg("x"));
		cl.def_static("log10", (float (*)(float)) &Teuchos::ScalarTraits<float>::log10, "C++: Teuchos::ScalarTraits<float>::log10(float) --> float", pybind11::arg("x"));
	}
	{ // Teuchos::ScalarTraits file:Teuchos_ScalarTraits.hpp line:679
		pybind11::class_<Teuchos::ScalarTraits<double>, Teuchos::RCP<Teuchos::ScalarTraits<double>>> cl(M("Teuchos"), "ScalarTraits_double_t", "", pybind11::module_local());
		cl.def( pybind11::init( [](){ return new Teuchos::ScalarTraits<double>(); } ) );
		cl.def_static("eps", (double (*)()) &Teuchos::ScalarTraits<double>::eps, "C++: Teuchos::ScalarTraits<double>::eps() --> double");
		cl.def_static("sfmin", (double (*)()) &Teuchos::ScalarTraits<double>::sfmin, "C++: Teuchos::ScalarTraits<double>::sfmin() --> double");
		cl.def_static("base", (double (*)()) &Teuchos::ScalarTraits<double>::base, "C++: Teuchos::ScalarTraits<double>::base() --> double");
		cl.def_static("prec", (double (*)()) &Teuchos::ScalarTraits<double>::prec, "C++: Teuchos::ScalarTraits<double>::prec() --> double");
		cl.def_static("t", (double (*)()) &Teuchos::ScalarTraits<double>::t, "C++: Teuchos::ScalarTraits<double>::t() --> double");
		cl.def_static("rnd", (double (*)()) &Teuchos::ScalarTraits<double>::rnd, "C++: Teuchos::ScalarTraits<double>::rnd() --> double");
		cl.def_static("emin", (double (*)()) &Teuchos::ScalarTraits<double>::emin, "C++: Teuchos::ScalarTraits<double>::emin() --> double");
		cl.def_static("rmin", (double (*)()) &Teuchos::ScalarTraits<double>::rmin, "C++: Teuchos::ScalarTraits<double>::rmin() --> double");
		cl.def_static("emax", (double (*)()) &Teuchos::ScalarTraits<double>::emax, "C++: Teuchos::ScalarTraits<double>::emax() --> double");
		cl.def_static("rmax", (double (*)()) &Teuchos::ScalarTraits<double>::rmax, "C++: Teuchos::ScalarTraits<double>::rmax() --> double");
		cl.def_static("magnitude", (double (*)(double)) &Teuchos::ScalarTraits<double>::magnitude, "C++: Teuchos::ScalarTraits<double>::magnitude(double) --> double", pybind11::arg("a"));
		cl.def_static("zero", (double (*)()) &Teuchos::ScalarTraits<double>::zero, "C++: Teuchos::ScalarTraits<double>::zero() --> double");
		cl.def_static("one", (double (*)()) &Teuchos::ScalarTraits<double>::one, "C++: Teuchos::ScalarTraits<double>::one() --> double");
		cl.def_static("conjugate", (double (*)(double)) &Teuchos::ScalarTraits<double>::conjugate, "C++: Teuchos::ScalarTraits<double>::conjugate(double) --> double", pybind11::arg("x"));
		cl.def_static("real", (double (*)(double)) &Teuchos::ScalarTraits<double>::real, "C++: Teuchos::ScalarTraits<double>::real(double) --> double", pybind11::arg("x"));
		cl.def_static("imag", (double (*)(double)) &Teuchos::ScalarTraits<double>::imag, "C++: Teuchos::ScalarTraits<double>::imag(double) --> double", pybind11::arg(""));
		cl.def_static("nan", (double (*)()) &Teuchos::ScalarTraits<double>::nan, "C++: Teuchos::ScalarTraits<double>::nan() --> double");
		cl.def_static("isnaninf", (bool (*)(double)) &Teuchos::ScalarTraits<double>::isnaninf, "C++: Teuchos::ScalarTraits<double>::isnaninf(double) --> bool", pybind11::arg("x"));
		cl.def_static("seedrandom", (void (*)(unsigned int)) &Teuchos::ScalarTraits<double>::seedrandom, "C++: Teuchos::ScalarTraits<double>::seedrandom(unsigned int) --> void", pybind11::arg("s"));
		cl.def_static("random", (double (*)()) &Teuchos::ScalarTraits<double>::random, "C++: Teuchos::ScalarTraits<double>::random() --> double");
		cl.def_static("name", (std::string (*)()) &Teuchos::ScalarTraits<double>::name, "C++: Teuchos::ScalarTraits<double>::name() --> std::string");
		cl.def_static("squareroot", (double (*)(double)) &Teuchos::ScalarTraits<double>::squareroot, "C++: Teuchos::ScalarTraits<double>::squareroot(double) --> double", pybind11::arg("x"));
		cl.def_static("pow", (double (*)(double, double)) &Teuchos::ScalarTraits<double>::pow, "C++: Teuchos::ScalarTraits<double>::pow(double, double) --> double", pybind11::arg("x"), pybind11::arg("y"));
		cl.def_static("pi", (double (*)()) &Teuchos::ScalarTraits<double>::pi, "C++: Teuchos::ScalarTraits<double>::pi() --> double");
		cl.def_static("log", (double (*)(double)) &Teuchos::ScalarTraits<double>::log, "C++: Teuchos::ScalarTraits<double>::log(double) --> double", pybind11::arg("x"));
		cl.def_static("log10", (double (*)(double)) &Teuchos::ScalarTraits<double>::log10, "C++: Teuchos::ScalarTraits<double>::log10(double) --> double", pybind11::arg("x"));
	}
	{ // Teuchos::OrdinalTraits file:Teuchos_OrdinalTraits.hpp line:108
		pybind11::class_<Teuchos::OrdinalTraits<char>, Teuchos::RCP<Teuchos::OrdinalTraits<char>>> cl(M("Teuchos"), "OrdinalTraits_char_t", "", pybind11::module_local());
		cl.def( pybind11::init( [](){ return new Teuchos::OrdinalTraits<char>(); } ) );
		cl.def_static("zero", (char (*)()) &Teuchos::OrdinalTraits<char>::zero, "C++: Teuchos::OrdinalTraits<char>::zero() --> char");
		cl.def_static("one", (char (*)()) &Teuchos::OrdinalTraits<char>::one, "C++: Teuchos::OrdinalTraits<char>::one() --> char");
		cl.def_static("invalid", (char (*)()) &Teuchos::OrdinalTraits<char>::invalid, "C++: Teuchos::OrdinalTraits<char>::invalid() --> char");
		cl.def_static("max", (char (*)()) &Teuchos::OrdinalTraits<char>::max, "C++: Teuchos::OrdinalTraits<char>::max() --> char");
		cl.def_static("name", (std::string (*)()) &Teuchos::OrdinalTraits<char>::name, "C++: Teuchos::OrdinalTraits<char>::name() --> std::string");
	}
	{ // Teuchos::OrdinalTraits file:Teuchos_OrdinalTraits.hpp line:118
		pybind11::class_<Teuchos::OrdinalTraits<short>, Teuchos::RCP<Teuchos::OrdinalTraits<short>>> cl(M("Teuchos"), "OrdinalTraits_short_t", "", pybind11::module_local());
		cl.def( pybind11::init( [](){ return new Teuchos::OrdinalTraits<short>(); } ) );
		cl.def_static("zero", (short (*)()) &Teuchos::OrdinalTraits<short>::zero, "C++: Teuchos::OrdinalTraits<short>::zero() --> short");
		cl.def_static("one", (short (*)()) &Teuchos::OrdinalTraits<short>::one, "C++: Teuchos::OrdinalTraits<short>::one() --> short");
		cl.def_static("invalid", (short (*)()) &Teuchos::OrdinalTraits<short>::invalid, "C++: Teuchos::OrdinalTraits<short>::invalid() --> short");
		cl.def_static("max", (short (*)()) &Teuchos::OrdinalTraits<short>::max, "C++: Teuchos::OrdinalTraits<short>::max() --> short");
		cl.def_static("name", (std::string (*)()) &Teuchos::OrdinalTraits<short>::name, "C++: Teuchos::OrdinalTraits<short>::name() --> std::string");
	}
	{ // Teuchos::OrdinalTraits file:Teuchos_OrdinalTraits.hpp line:128
		pybind11::class_<Teuchos::OrdinalTraits<int>, Teuchos::RCP<Teuchos::OrdinalTraits<int>>> cl(M("Teuchos"), "OrdinalTraits_int_t", "", pybind11::module_local());
		cl.def( pybind11::init( [](){ return new Teuchos::OrdinalTraits<int>(); } ) );
		cl.def_static("zero", (int (*)()) &Teuchos::OrdinalTraits<int>::zero, "C++: Teuchos::OrdinalTraits<int>::zero() --> int");
		cl.def_static("one", (int (*)()) &Teuchos::OrdinalTraits<int>::one, "C++: Teuchos::OrdinalTraits<int>::one() --> int");
		cl.def_static("invalid", (int (*)()) &Teuchos::OrdinalTraits<int>::invalid, "C++: Teuchos::OrdinalTraits<int>::invalid() --> int");
		cl.def_static("max", (int (*)()) &Teuchos::OrdinalTraits<int>::max, "C++: Teuchos::OrdinalTraits<int>::max() --> int");
		cl.def_static("name", (std::string (*)()) &Teuchos::OrdinalTraits<int>::name, "C++: Teuchos::OrdinalTraits<int>::name() --> std::string");
	}
	{ // Teuchos::OrdinalTraits file:Teuchos_OrdinalTraits.hpp line:138
		pybind11::class_<Teuchos::OrdinalTraits<unsigned int>, Teuchos::RCP<Teuchos::OrdinalTraits<unsigned int>>> cl(M("Teuchos"), "OrdinalTraits_unsigned_int_t", "", pybind11::module_local());
		cl.def( pybind11::init( [](){ return new Teuchos::OrdinalTraits<unsigned int>(); } ) );
		cl.def_static("zero", (unsigned int (*)()) &Teuchos::OrdinalTraits<unsigned int>::zero, "C++: Teuchos::OrdinalTraits<unsigned int>::zero() --> unsigned int");
		cl.def_static("one", (unsigned int (*)()) &Teuchos::OrdinalTraits<unsigned int>::one, "C++: Teuchos::OrdinalTraits<unsigned int>::one() --> unsigned int");
		cl.def_static("invalid", (unsigned int (*)()) &Teuchos::OrdinalTraits<unsigned int>::invalid, "C++: Teuchos::OrdinalTraits<unsigned int>::invalid() --> unsigned int");
		cl.def_static("max", (unsigned int (*)()) &Teuchos::OrdinalTraits<unsigned int>::max, "C++: Teuchos::OrdinalTraits<unsigned int>::max() --> unsigned int");
		cl.def_static("name", (std::string (*)()) &Teuchos::OrdinalTraits<unsigned int>::name, "C++: Teuchos::OrdinalTraits<unsigned int>::name() --> std::string");
	}
	{ // Teuchos::OrdinalTraits file:Teuchos_OrdinalTraits.hpp line:148
		pybind11::class_<Teuchos::OrdinalTraits<long>, Teuchos::RCP<Teuchos::OrdinalTraits<long>>> cl(M("Teuchos"), "OrdinalTraits_long_t", "", pybind11::module_local());
		cl.def( pybind11::init( [](){ return new Teuchos::OrdinalTraits<long>(); } ) );
		cl.def_static("zero", (long (*)()) &Teuchos::OrdinalTraits<long>::zero, "C++: Teuchos::OrdinalTraits<long>::zero() --> long");
		cl.def_static("one", (long (*)()) &Teuchos::OrdinalTraits<long>::one, "C++: Teuchos::OrdinalTraits<long>::one() --> long");
		cl.def_static("invalid", (long (*)()) &Teuchos::OrdinalTraits<long>::invalid, "C++: Teuchos::OrdinalTraits<long>::invalid() --> long");
		cl.def_static("max", (long (*)()) &Teuchos::OrdinalTraits<long>::max, "C++: Teuchos::OrdinalTraits<long>::max() --> long");
		cl.def_static("name", (std::string (*)()) &Teuchos::OrdinalTraits<long>::name, "C++: Teuchos::OrdinalTraits<long>::name() --> std::string");
	}
	{ // Teuchos::OrdinalTraits file:Teuchos_OrdinalTraits.hpp line:158
		pybind11::class_<Teuchos::OrdinalTraits<unsigned long>, Teuchos::RCP<Teuchos::OrdinalTraits<unsigned long>>> cl(M("Teuchos"), "OrdinalTraits_unsigned_long_t", "", pybind11::module_local());
		cl.def( pybind11::init( [](){ return new Teuchos::OrdinalTraits<unsigned long>(); } ) );
		cl.def_static("zero", (unsigned long (*)()) &Teuchos::OrdinalTraits<unsigned long>::zero, "C++: Teuchos::OrdinalTraits<unsigned long>::zero() --> unsigned long");
		cl.def_static("one", (unsigned long (*)()) &Teuchos::OrdinalTraits<unsigned long>::one, "C++: Teuchos::OrdinalTraits<unsigned long>::one() --> unsigned long");
		cl.def_static("invalid", (unsigned long (*)()) &Teuchos::OrdinalTraits<unsigned long>::invalid, "C++: Teuchos::OrdinalTraits<unsigned long>::invalid() --> unsigned long");
		cl.def_static("max", (unsigned long (*)()) &Teuchos::OrdinalTraits<unsigned long>::max, "C++: Teuchos::OrdinalTraits<unsigned long>::max() --> unsigned long");
		cl.def_static("name", (std::string (*)()) &Teuchos::OrdinalTraits<unsigned long>::name, "C++: Teuchos::OrdinalTraits<unsigned long>::name() --> std::string");
	}
	{ // Teuchos::OrdinalTraits file:Teuchos_OrdinalTraits.hpp line:168
		pybind11::class_<Teuchos::OrdinalTraits<long long>, Teuchos::RCP<Teuchos::OrdinalTraits<long long>>> cl(M("Teuchos"), "OrdinalTraits_long_long_t", "", pybind11::module_local());
		cl.def( pybind11::init( [](){ return new Teuchos::OrdinalTraits<long long>(); } ) );
		cl.def_static("zero", (long long (*)()) &Teuchos::OrdinalTraits<long long>::zero, "C++: Teuchos::OrdinalTraits<long long>::zero() --> long long");
		cl.def_static("one", (long long (*)()) &Teuchos::OrdinalTraits<long long>::one, "C++: Teuchos::OrdinalTraits<long long>::one() --> long long");
		cl.def_static("invalid", (long long (*)()) &Teuchos::OrdinalTraits<long long>::invalid, "C++: Teuchos::OrdinalTraits<long long>::invalid() --> long long");
		cl.def_static("max", (long long (*)()) &Teuchos::OrdinalTraits<long long>::max, "C++: Teuchos::OrdinalTraits<long long>::max() --> long long");
		cl.def_static("name", (std::string (*)()) &Teuchos::OrdinalTraits<long long>::name, "C++: Teuchos::OrdinalTraits<long long>::name() --> std::string");
	}
	{ // Teuchos::OrdinalTraits file:Teuchos_OrdinalTraits.hpp line:178
		pybind11::class_<Teuchos::OrdinalTraits<unsigned long long>, Teuchos::RCP<Teuchos::OrdinalTraits<unsigned long long>>> cl(M("Teuchos"), "OrdinalTraits_unsigned_long_long_t", "", pybind11::module_local());
		cl.def( pybind11::init( [](){ return new Teuchos::OrdinalTraits<unsigned long long>(); } ) );
		cl.def_static("zero", (unsigned long long (*)()) &Teuchos::OrdinalTraits<unsigned long long>::zero, "C++: Teuchos::OrdinalTraits<unsigned long long>::zero() --> unsigned long long");
		cl.def_static("one", (unsigned long long (*)()) &Teuchos::OrdinalTraits<unsigned long long>::one, "C++: Teuchos::OrdinalTraits<unsigned long long>::one() --> unsigned long long");
		cl.def_static("invalid", (unsigned long long (*)()) &Teuchos::OrdinalTraits<unsigned long long>::invalid, "C++: Teuchos::OrdinalTraits<unsigned long long>::invalid() --> unsigned long long");
		cl.def_static("max", (unsigned long long (*)()) &Teuchos::OrdinalTraits<unsigned long long>::max, "C++: Teuchos::OrdinalTraits<unsigned long long>::max() --> unsigned long long");
		cl.def_static("name", (std::string (*)()) &Teuchos::OrdinalTraits<unsigned long long>::name, "C++: Teuchos::OrdinalTraits<unsigned long long>::name() --> std::string");
	}
	// Teuchos::ESide file:Teuchos_BLAS_types.hpp line:88
	pybind11::enum_<Teuchos::ESide>(M("Teuchos"), "ESide", pybind11::arithmetic(), "", pybind11::module_local())
		.value("LEFT_SIDE", Teuchos::LEFT_SIDE)
		.value("RIGHT_SIDE", Teuchos::RIGHT_SIDE)
		.export_values();

;

	// Teuchos::ETransp file:Teuchos_BLAS_types.hpp line:93
	pybind11::enum_<Teuchos::ETransp>(M("Teuchos"), "ETransp", pybind11::arithmetic(), "", pybind11::module_local())
		.value("NO_TRANS", Teuchos::NO_TRANS)
		.value("TRANS", Teuchos::TRANS)
		.value("CONJ_TRANS", Teuchos::CONJ_TRANS)
		.export_values();

;

	// Teuchos::EUplo file:Teuchos_BLAS_types.hpp line:99
	pybind11::enum_<Teuchos::EUplo>(M("Teuchos"), "EUplo", pybind11::arithmetic(), "", pybind11::module_local())
		.value("UPPER_TRI", Teuchos::UPPER_TRI)
		.value("LOWER_TRI", Teuchos::LOWER_TRI)
		.value("UNDEF_TRI", Teuchos::UNDEF_TRI)
		.export_values();

;

	// Teuchos::EDiag file:Teuchos_BLAS_types.hpp line:105
	pybind11::enum_<Teuchos::EDiag>(M("Teuchos"), "EDiag", pybind11::arithmetic(), "", pybind11::module_local())
		.value("UNIT_DIAG", Teuchos::UNIT_DIAG)
		.value("NON_UNIT_DIAG", Teuchos::NON_UNIT_DIAG)
		.export_values();

;

	// Teuchos::EType file:Teuchos_BLAS_types.hpp line:110
	pybind11::enum_<Teuchos::EType>(M("Teuchos"), "EType", pybind11::arithmetic(), "", pybind11::module_local())
		.value("FULL", Teuchos::FULL)
		.value("LOWER", Teuchos::LOWER)
		.value("UPPER", Teuchos::UPPER)
		.value("HESSENBERG", Teuchos::HESSENBERG)
		.value("SYM_BAND_L", Teuchos::SYM_BAND_L)
		.value("SYM_BAND_U", Teuchos::SYM_BAND_U)
		.value("BAND", Teuchos::BAND)
		.export_values();

;

	{ // Teuchos::SerialDenseMatrix file:Teuchos_SerialDenseMatrix.hpp line:67
		pybind11::class_<Teuchos::SerialDenseMatrix<int,double>, Teuchos::RCP<Teuchos::SerialDenseMatrix<int,double>>, PyCallBack_Teuchos_SerialDenseMatrix_int_double_t, Teuchos::CompObject> cl(M("Teuchos"), "SerialDenseMatrix_int_double_t", "", pybind11::module_local());
		cl.def( pybind11::init( [](){ return new Teuchos::SerialDenseMatrix<int,double>(); }, [](){ return new PyCallBack_Teuchos_SerialDenseMatrix_int_double_t(); } ) );
		cl.def( pybind11::init( [](int const & a0, int const & a1){ return new Teuchos::SerialDenseMatrix<int,double>(a0, a1); }, [](int const & a0, int const & a1){ return new PyCallBack_Teuchos_SerialDenseMatrix_int_double_t(a0, a1); } ), "doc");
		cl.def( pybind11::init<int, int, bool>(), pybind11::arg("numRows_in"), pybind11::arg("numCols_in"), pybind11::arg("zeroOut") );

		cl.def( pybind11::init<enum Teuchos::DataAccess, double *, int, int, int>(), pybind11::arg("CV"), pybind11::arg("values"), pybind11::arg("stride"), pybind11::arg("numRows"), pybind11::arg("numCols") );

		cl.def( pybind11::init( [](PyCallBack_Teuchos_SerialDenseMatrix_int_double_t const &o){ return new PyCallBack_Teuchos_SerialDenseMatrix_int_double_t(o); } ) );
		cl.def( pybind11::init( [](Teuchos::SerialDenseMatrix<int,double> const &o){ return new Teuchos::SerialDenseMatrix<int,double>(o); } ) );
		cl.def( pybind11::init<enum Teuchos::DataAccess, const class Teuchos::SerialDenseMatrix<int, double> &>(), pybind11::arg("CV"), pybind11::arg("Source") );

		cl.def( pybind11::init( [](enum Teuchos::DataAccess const & a0, const class Teuchos::SerialDenseMatrix<int, double> & a1, int const & a2, int const & a3){ return new Teuchos::SerialDenseMatrix<int,double>(a0, a1, a2, a3); }, [](enum Teuchos::DataAccess const & a0, const class Teuchos::SerialDenseMatrix<int, double> & a1, int const & a2, int const & a3){ return new PyCallBack_Teuchos_SerialDenseMatrix_int_double_t(a0, a1, a2, a3); } ), "doc");
		cl.def( pybind11::init( [](enum Teuchos::DataAccess const & a0, const class Teuchos::SerialDenseMatrix<int, double> & a1, int const & a2, int const & a3, int const & a4){ return new Teuchos::SerialDenseMatrix<int,double>(a0, a1, a2, a3, a4); }, [](enum Teuchos::DataAccess const & a0, const class Teuchos::SerialDenseMatrix<int, double> & a1, int const & a2, int const & a3, int const & a4){ return new PyCallBack_Teuchos_SerialDenseMatrix_int_double_t(a0, a1, a2, a3, a4); } ), "doc");
		cl.def( pybind11::init<enum Teuchos::DataAccess, const class Teuchos::SerialDenseMatrix<int, double> &, int, int, int, int>(), pybind11::arg("CV"), pybind11::arg("Source"), pybind11::arg("numRows_in"), pybind11::arg("numCols_in"), pybind11::arg("startRow"), pybind11::arg("startCol") );

		cl.def("shape", (int (Teuchos::SerialDenseMatrix<int,double>::*)(int, int)) &Teuchos::SerialDenseMatrix<int, double>::shape, "C++: Teuchos::SerialDenseMatrix<int, double>::shape(int, int) --> int", pybind11::arg("numRows_in"), pybind11::arg("numCols_in"));
		cl.def("shapeUninitialized", (int (Teuchos::SerialDenseMatrix<int,double>::*)(int, int)) &Teuchos::SerialDenseMatrix<int, double>::shapeUninitialized, "C++: Teuchos::SerialDenseMatrix<int, double>::shapeUninitialized(int, int) --> int", pybind11::arg("numRows"), pybind11::arg("numCols"));
		cl.def("reshape", (int (Teuchos::SerialDenseMatrix<int,double>::*)(int, int)) &Teuchos::SerialDenseMatrix<int, double>::reshape, "C++: Teuchos::SerialDenseMatrix<int, double>::reshape(int, int) --> int", pybind11::arg("numRows_in"), pybind11::arg("numCols_in"));
		cl.def("assign", (class Teuchos::SerialDenseMatrix<int, double> & (Teuchos::SerialDenseMatrix<int,double>::*)(const class Teuchos::SerialDenseMatrix<int, double> &)) &Teuchos::SerialDenseMatrix<int, double>::operator=, "C++: Teuchos::SerialDenseMatrix<int, double>::operator=(const class Teuchos::SerialDenseMatrix<int, double> &) --> class Teuchos::SerialDenseMatrix<int, double> &", pybind11::return_value_policy::automatic, pybind11::arg("Source"));
		cl.def("assign", (class Teuchos::SerialDenseMatrix<int, double> & (Teuchos::SerialDenseMatrix<int,double>::*)(const class Teuchos::SerialDenseMatrix<int, double> &)) &Teuchos::SerialDenseMatrix<int, double>::assign, "C++: Teuchos::SerialDenseMatrix<int, double>::assign(const class Teuchos::SerialDenseMatrix<int, double> &) --> class Teuchos::SerialDenseMatrix<int, double> &", pybind11::return_value_policy::automatic, pybind11::arg("Source"));
		cl.def("assign", (class Teuchos::SerialDenseMatrix<int, double> & (Teuchos::SerialDenseMatrix<int,double>::*)(const double)) &Teuchos::SerialDenseMatrix<int, double>::operator=, "C++: Teuchos::SerialDenseMatrix<int, double>::operator=(const double) --> class Teuchos::SerialDenseMatrix<int, double> &", pybind11::return_value_policy::automatic, pybind11::arg("value"));
		cl.def("putScalar", [](Teuchos::SerialDenseMatrix<int,double> &o) -> int { return o.putScalar(); }, "");
		cl.def("putScalar", (int (Teuchos::SerialDenseMatrix<int,double>::*)(const double)) &Teuchos::SerialDenseMatrix<int, double>::putScalar, "C++: Teuchos::SerialDenseMatrix<int, double>::putScalar(const double) --> int", pybind11::arg("value_in"));
		cl.def("swap", (void (Teuchos::SerialDenseMatrix<int,double>::*)(class Teuchos::SerialDenseMatrix<int, double> &)) &Teuchos::SerialDenseMatrix<int, double>::swap, "C++: Teuchos::SerialDenseMatrix<int, double>::swap(class Teuchos::SerialDenseMatrix<int, double> &) --> void", pybind11::arg("B"));
		cl.def("random", (int (Teuchos::SerialDenseMatrix<int,double>::*)()) &Teuchos::SerialDenseMatrix<int, double>::random, "C++: Teuchos::SerialDenseMatrix<int, double>::random() --> int");
		cl.def("__call__", (double & (Teuchos::SerialDenseMatrix<int,double>::*)(int, int)) &Teuchos::SerialDenseMatrix<int, double>::operator(), "C++: Teuchos::SerialDenseMatrix<int, double>::operator()(int, int) --> double &", pybind11::return_value_policy::automatic, pybind11::arg("rowIndex"), pybind11::arg("colIndex"));
		cl.def("__getitem__", (double * (Teuchos::SerialDenseMatrix<int,double>::*)(int)) &Teuchos::SerialDenseMatrix<int, double>::operator[], "C++: Teuchos::SerialDenseMatrix<int, double>::operator[](int) --> double *", pybind11::return_value_policy::automatic, pybind11::arg("colIndex"));
		cl.def("values", (double * (Teuchos::SerialDenseMatrix<int,double>::*)() const) &Teuchos::SerialDenseMatrix<int, double>::values, "C++: Teuchos::SerialDenseMatrix<int, double>::values() const --> double *", pybind11::return_value_policy::automatic);
		cl.def("__iadd__", (class Teuchos::SerialDenseMatrix<int, double> & (Teuchos::SerialDenseMatrix<int,double>::*)(const class Teuchos::SerialDenseMatrix<int, double> &)) &Teuchos::SerialDenseMatrix<int, double>::operator+=, "C++: Teuchos::SerialDenseMatrix<int, double>::operator+=(const class Teuchos::SerialDenseMatrix<int, double> &) --> class Teuchos::SerialDenseMatrix<int, double> &", pybind11::return_value_policy::automatic, pybind11::arg("Source"));
		cl.def("__isub__", (class Teuchos::SerialDenseMatrix<int, double> & (Teuchos::SerialDenseMatrix<int,double>::*)(const class Teuchos::SerialDenseMatrix<int, double> &)) &Teuchos::SerialDenseMatrix<int, double>::operator-=, "C++: Teuchos::SerialDenseMatrix<int, double>::operator-=(const class Teuchos::SerialDenseMatrix<int, double> &) --> class Teuchos::SerialDenseMatrix<int, double> &", pybind11::return_value_policy::automatic, pybind11::arg("Source"));
		cl.def("__imul__", (class Teuchos::SerialDenseMatrix<int, double> & (Teuchos::SerialDenseMatrix<int,double>::*)(const double)) &Teuchos::SerialDenseMatrix<int, double>::operator*=, "C++: Teuchos::SerialDenseMatrix<int, double>::operator*=(const double) --> class Teuchos::SerialDenseMatrix<int, double> &", pybind11::return_value_policy::automatic, pybind11::arg("alpha"));
		cl.def("scale", (int (Teuchos::SerialDenseMatrix<int,double>::*)(const double)) &Teuchos::SerialDenseMatrix<int, double>::scale, "C++: Teuchos::SerialDenseMatrix<int, double>::scale(const double) --> int", pybind11::arg("alpha"));
		cl.def("scale", (int (Teuchos::SerialDenseMatrix<int,double>::*)(const class Teuchos::SerialDenseMatrix<int, double> &)) &Teuchos::SerialDenseMatrix<int, double>::scale, "C++: Teuchos::SerialDenseMatrix<int, double>::scale(const class Teuchos::SerialDenseMatrix<int, double> &) --> int", pybind11::arg("A"));
		cl.def("multiply", (int (Teuchos::SerialDenseMatrix<int,double>::*)(enum Teuchos::ETransp, enum Teuchos::ETransp, double, const class Teuchos::SerialDenseMatrix<int, double> &, const class Teuchos::SerialDenseMatrix<int, double> &, double)) &Teuchos::SerialDenseMatrix<int, double>::multiply, "C++: Teuchos::SerialDenseMatrix<int, double>::multiply(enum Teuchos::ETransp, enum Teuchos::ETransp, double, const class Teuchos::SerialDenseMatrix<int, double> &, const class Teuchos::SerialDenseMatrix<int, double> &, double) --> int", pybind11::arg("transa"), pybind11::arg("transb"), pybind11::arg("alpha"), pybind11::arg("A"), pybind11::arg("B"), pybind11::arg("beta"));
		cl.def("__eq__", (bool (Teuchos::SerialDenseMatrix<int,double>::*)(const class Teuchos::SerialDenseMatrix<int, double> &) const) &Teuchos::SerialDenseMatrix<int, double>::operator==, "C++: Teuchos::SerialDenseMatrix<int, double>::operator==(const class Teuchos::SerialDenseMatrix<int, double> &) const --> bool", pybind11::arg("Operand"));
		cl.def("__ne__", (bool (Teuchos::SerialDenseMatrix<int,double>::*)(const class Teuchos::SerialDenseMatrix<int, double> &) const) &Teuchos::SerialDenseMatrix<int, double>::operator!=, "C++: Teuchos::SerialDenseMatrix<int, double>::operator!=(const class Teuchos::SerialDenseMatrix<int, double> &) const --> bool", pybind11::arg("Operand"));
		cl.def("numRows", (int (Teuchos::SerialDenseMatrix<int,double>::*)() const) &Teuchos::SerialDenseMatrix<int, double>::numRows, "C++: Teuchos::SerialDenseMatrix<int, double>::numRows() const --> int");
		cl.def("numCols", (int (Teuchos::SerialDenseMatrix<int,double>::*)() const) &Teuchos::SerialDenseMatrix<int, double>::numCols, "C++: Teuchos::SerialDenseMatrix<int, double>::numCols() const --> int");
		cl.def("stride", (int (Teuchos::SerialDenseMatrix<int,double>::*)() const) &Teuchos::SerialDenseMatrix<int, double>::stride, "C++: Teuchos::SerialDenseMatrix<int, double>::stride() const --> int");
		cl.def("empty", (bool (Teuchos::SerialDenseMatrix<int,double>::*)() const) &Teuchos::SerialDenseMatrix<int, double>::empty, "C++: Teuchos::SerialDenseMatrix<int, double>::empty() const --> bool");
		cl.def("normOne", (double (Teuchos::SerialDenseMatrix<int,double>::*)() const) &Teuchos::SerialDenseMatrix<int, double>::normOne, "C++: Teuchos::SerialDenseMatrix<int, double>::normOne() const --> double");
		cl.def("normInf", (double (Teuchos::SerialDenseMatrix<int,double>::*)() const) &Teuchos::SerialDenseMatrix<int, double>::normInf, "C++: Teuchos::SerialDenseMatrix<int, double>::normInf() const --> double");
		cl.def("normFrobenius", (double (Teuchos::SerialDenseMatrix<int,double>::*)() const) &Teuchos::SerialDenseMatrix<int, double>::normFrobenius, "C++: Teuchos::SerialDenseMatrix<int, double>::normFrobenius() const --> double");
		cl.def("print", (std::ostream & (Teuchos::SerialDenseMatrix<int,double>::*)(std::ostream &) const) &Teuchos::SerialDenseMatrix<int, double>::print, "C++: Teuchos::SerialDenseMatrix<int, double>::print(std::ostream &) const --> std::ostream &", pybind11::return_value_policy::automatic, pybind11::arg("os"));
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
	{ // Teuchos::SerialDenseVector file:Teuchos_SerialDenseVector.hpp line:60
		pybind11::class_<Teuchos::SerialDenseVector<int,double>, Teuchos::RCP<Teuchos::SerialDenseVector<int,double>>, PyCallBack_Teuchos_SerialDenseVector_int_double_t, Teuchos::SerialDenseMatrix<int,double>> cl(M("Teuchos"), "SerialDenseVector_int_double_t", "", pybind11::module_local());
		cl.def( pybind11::init( [](){ return new Teuchos::SerialDenseVector<int,double>(); }, [](){ return new PyCallBack_Teuchos_SerialDenseVector_int_double_t(); } ) );
		cl.def( pybind11::init( [](int const & a0){ return new Teuchos::SerialDenseVector<int,double>(a0); }, [](int const & a0){ return new PyCallBack_Teuchos_SerialDenseVector_int_double_t(a0); } ), "doc");
		cl.def( pybind11::init<int, bool>(), pybind11::arg("length_in"), pybind11::arg("zeroOut") );

		cl.def( pybind11::init<enum Teuchos::DataAccess, double *, int>(), pybind11::arg("CV"), pybind11::arg("values"), pybind11::arg("length") );

		cl.def( pybind11::init( [](PyCallBack_Teuchos_SerialDenseVector_int_double_t const &o){ return new PyCallBack_Teuchos_SerialDenseVector_int_double_t(o); } ) );
		cl.def( pybind11::init( [](Teuchos::SerialDenseVector<int,double> const &o){ return new Teuchos::SerialDenseVector<int,double>(o); } ) );
		cl.def( pybind11::init<enum Teuchos::DataAccess, const class Teuchos::SerialDenseVector<int, double> &>(), pybind11::arg("CV"), pybind11::arg("Source") );

		cl.def("size", (int (Teuchos::SerialDenseVector<int,double>::*)(int)) &Teuchos::SerialDenseVector<int, double>::size, "C++: Teuchos::SerialDenseVector<int, double>::size(int) --> int", pybind11::arg("length_in"));
		cl.def("sizeUninitialized", (int (Teuchos::SerialDenseVector<int,double>::*)(int)) &Teuchos::SerialDenseVector<int, double>::sizeUninitialized, "C++: Teuchos::SerialDenseVector<int, double>::sizeUninitialized(int) --> int", pybind11::arg("length_in"));
		cl.def("resize", (int (Teuchos::SerialDenseVector<int,double>::*)(int)) &Teuchos::SerialDenseVector<int, double>::resize, "C++: Teuchos::SerialDenseVector<int, double>::resize(int) --> int", pybind11::arg("length_in"));
		cl.def("assign", (class Teuchos::SerialDenseVector<int, double> & (Teuchos::SerialDenseVector<int,double>::*)(const double)) &Teuchos::SerialDenseVector<int, double>::operator=, "C++: Teuchos::SerialDenseVector<int, double>::operator=(const double) --> class Teuchos::SerialDenseVector<int, double> &", pybind11::return_value_policy::automatic, pybind11::arg("value"));
		cl.def("__eq__", (bool (Teuchos::SerialDenseVector<int,double>::*)(const class Teuchos::SerialDenseVector<int, double> &) const) &Teuchos::SerialDenseVector<int, double>::operator==, "C++: Teuchos::SerialDenseVector<int, double>::operator==(const class Teuchos::SerialDenseVector<int, double> &) const --> bool", pybind11::arg("Operand"));
		cl.def("__ne__", (bool (Teuchos::SerialDenseVector<int,double>::*)(const class Teuchos::SerialDenseVector<int, double> &) const) &Teuchos::SerialDenseVector<int, double>::operator!=, "C++: Teuchos::SerialDenseVector<int, double>::operator!=(const class Teuchos::SerialDenseVector<int, double> &) const --> bool", pybind11::arg("Operand"));
		cl.def("assign", (class Teuchos::SerialDenseVector<int, double> & (Teuchos::SerialDenseVector<int,double>::*)(const class Teuchos::SerialDenseVector<int, double> &)) &Teuchos::SerialDenseVector<int, double>::operator=, "C++: Teuchos::SerialDenseVector<int, double>::operator=(const class Teuchos::SerialDenseVector<int, double> &) --> class Teuchos::SerialDenseVector<int, double> &", pybind11::return_value_policy::automatic, pybind11::arg("Source"));
		cl.def("__call__", (double & (Teuchos::SerialDenseVector<int,double>::*)(int)) &Teuchos::SerialDenseVector<int, double>::operator(), "C++: Teuchos::SerialDenseVector<int, double>::operator()(int) --> double &", pybind11::return_value_policy::automatic, pybind11::arg("index"));
		cl.def("__getitem__", (double & (Teuchos::SerialDenseVector<int,double>::*)(int)) &Teuchos::SerialDenseVector<int, double>::operator[], "C++: Teuchos::SerialDenseVector<int, double>::operator[](int) --> double &", pybind11::return_value_policy::automatic, pybind11::arg("index"));
		cl.def("dot", (double (Teuchos::SerialDenseVector<int,double>::*)(const class Teuchos::SerialDenseVector<int, double> &) const) &Teuchos::SerialDenseVector<int, double>::dot, "C++: Teuchos::SerialDenseVector<int, double>::dot(const class Teuchos::SerialDenseVector<int, double> &) const --> double", pybind11::arg("x"));
		cl.def("length", (int (Teuchos::SerialDenseVector<int,double>::*)() const) &Teuchos::SerialDenseVector<int, double>::length, "C++: Teuchos::SerialDenseVector<int, double>::length() const --> int");
		cl.def("print", (std::ostream & (Teuchos::SerialDenseVector<int,double>::*)(std::ostream &) const) &Teuchos::SerialDenseVector<int, double>::print, "C++: Teuchos::SerialDenseVector<int, double>::print(std::ostream &) const --> std::ostream &", pybind11::return_value_policy::automatic, pybind11::arg("os"));
		cl.def("shape", (int (Teuchos::SerialDenseMatrix<int,double>::*)(int, int)) &Teuchos::SerialDenseMatrix<int, double>::shape, "C++: Teuchos::SerialDenseMatrix<int, double>::shape(int, int) --> int", pybind11::arg("numRows_in"), pybind11::arg("numCols_in"));
		cl.def("shapeUninitialized", (int (Teuchos::SerialDenseMatrix<int,double>::*)(int, int)) &Teuchos::SerialDenseMatrix<int, double>::shapeUninitialized, "C++: Teuchos::SerialDenseMatrix<int, double>::shapeUninitialized(int, int) --> int", pybind11::arg("numRows"), pybind11::arg("numCols"));
		cl.def("reshape", (int (Teuchos::SerialDenseMatrix<int,double>::*)(int, int)) &Teuchos::SerialDenseMatrix<int, double>::reshape, "C++: Teuchos::SerialDenseMatrix<int, double>::reshape(int, int) --> int", pybind11::arg("numRows_in"), pybind11::arg("numCols_in"));
		cl.def("assign", (class Teuchos::SerialDenseMatrix<int, double> & (Teuchos::SerialDenseMatrix<int,double>::*)(const class Teuchos::SerialDenseMatrix<int, double> &)) &Teuchos::SerialDenseMatrix<int, double>::operator=, "C++: Teuchos::SerialDenseMatrix<int, double>::operator=(const class Teuchos::SerialDenseMatrix<int, double> &) --> class Teuchos::SerialDenseMatrix<int, double> &", pybind11::return_value_policy::automatic, pybind11::arg("Source"));
		cl.def("assign", (class Teuchos::SerialDenseMatrix<int, double> & (Teuchos::SerialDenseMatrix<int,double>::*)(const class Teuchos::SerialDenseMatrix<int, double> &)) &Teuchos::SerialDenseMatrix<int, double>::assign, "C++: Teuchos::SerialDenseMatrix<int, double>::assign(const class Teuchos::SerialDenseMatrix<int, double> &) --> class Teuchos::SerialDenseMatrix<int, double> &", pybind11::return_value_policy::automatic, pybind11::arg("Source"));
		cl.def("assign", (class Teuchos::SerialDenseMatrix<int, double> & (Teuchos::SerialDenseMatrix<int,double>::*)(const double)) &Teuchos::SerialDenseMatrix<int, double>::operator=, "C++: Teuchos::SerialDenseMatrix<int, double>::operator=(const double) --> class Teuchos::SerialDenseMatrix<int, double> &", pybind11::return_value_policy::automatic, pybind11::arg("value"));
		cl.def("putScalar", [](Teuchos::SerialDenseMatrix<int,double> &o) -> int { return o.putScalar(); }, "");
		cl.def("putScalar", (int (Teuchos::SerialDenseMatrix<int,double>::*)(const double)) &Teuchos::SerialDenseMatrix<int, double>::putScalar, "C++: Teuchos::SerialDenseMatrix<int, double>::putScalar(const double) --> int", pybind11::arg("value_in"));
		cl.def("swap", (void (Teuchos::SerialDenseMatrix<int,double>::*)(class Teuchos::SerialDenseMatrix<int, double> &)) &Teuchos::SerialDenseMatrix<int, double>::swap, "C++: Teuchos::SerialDenseMatrix<int, double>::swap(class Teuchos::SerialDenseMatrix<int, double> &) --> void", pybind11::arg("B"));
		cl.def("random", (int (Teuchos::SerialDenseMatrix<int,double>::*)()) &Teuchos::SerialDenseMatrix<int, double>::random, "C++: Teuchos::SerialDenseMatrix<int, double>::random() --> int");
		cl.def("__call__", (double & (Teuchos::SerialDenseMatrix<int,double>::*)(int, int)) &Teuchos::SerialDenseMatrix<int, double>::operator(), "C++: Teuchos::SerialDenseMatrix<int, double>::operator()(int, int) --> double &", pybind11::return_value_policy::automatic, pybind11::arg("rowIndex"), pybind11::arg("colIndex"));
		cl.def("__getitem__", (double * (Teuchos::SerialDenseMatrix<int,double>::*)(int)) &Teuchos::SerialDenseMatrix<int, double>::operator[], "C++: Teuchos::SerialDenseMatrix<int, double>::operator[](int) --> double *", pybind11::return_value_policy::automatic, pybind11::arg("colIndex"));
		cl.def("values", (double * (Teuchos::SerialDenseMatrix<int,double>::*)() const) &Teuchos::SerialDenseMatrix<int, double>::values, "C++: Teuchos::SerialDenseMatrix<int, double>::values() const --> double *", pybind11::return_value_policy::automatic);
		cl.def("__iadd__", (class Teuchos::SerialDenseMatrix<int, double> & (Teuchos::SerialDenseMatrix<int,double>::*)(const class Teuchos::SerialDenseMatrix<int, double> &)) &Teuchos::SerialDenseMatrix<int, double>::operator+=, "C++: Teuchos::SerialDenseMatrix<int, double>::operator+=(const class Teuchos::SerialDenseMatrix<int, double> &) --> class Teuchos::SerialDenseMatrix<int, double> &", pybind11::return_value_policy::automatic, pybind11::arg("Source"));
		cl.def("__isub__", (class Teuchos::SerialDenseMatrix<int, double> & (Teuchos::SerialDenseMatrix<int,double>::*)(const class Teuchos::SerialDenseMatrix<int, double> &)) &Teuchos::SerialDenseMatrix<int, double>::operator-=, "C++: Teuchos::SerialDenseMatrix<int, double>::operator-=(const class Teuchos::SerialDenseMatrix<int, double> &) --> class Teuchos::SerialDenseMatrix<int, double> &", pybind11::return_value_policy::automatic, pybind11::arg("Source"));
		cl.def("__imul__", (class Teuchos::SerialDenseMatrix<int, double> & (Teuchos::SerialDenseMatrix<int,double>::*)(const double)) &Teuchos::SerialDenseMatrix<int, double>::operator*=, "C++: Teuchos::SerialDenseMatrix<int, double>::operator*=(const double) --> class Teuchos::SerialDenseMatrix<int, double> &", pybind11::return_value_policy::automatic, pybind11::arg("alpha"));
		cl.def("scale", (int (Teuchos::SerialDenseMatrix<int,double>::*)(const double)) &Teuchos::SerialDenseMatrix<int, double>::scale, "C++: Teuchos::SerialDenseMatrix<int, double>::scale(const double) --> int", pybind11::arg("alpha"));
		cl.def("scale", (int (Teuchos::SerialDenseMatrix<int,double>::*)(const class Teuchos::SerialDenseMatrix<int, double> &)) &Teuchos::SerialDenseMatrix<int, double>::scale, "C++: Teuchos::SerialDenseMatrix<int, double>::scale(const class Teuchos::SerialDenseMatrix<int, double> &) --> int", pybind11::arg("A"));
		cl.def("multiply", (int (Teuchos::SerialDenseMatrix<int,double>::*)(enum Teuchos::ETransp, enum Teuchos::ETransp, double, const class Teuchos::SerialDenseMatrix<int, double> &, const class Teuchos::SerialDenseMatrix<int, double> &, double)) &Teuchos::SerialDenseMatrix<int, double>::multiply, "C++: Teuchos::SerialDenseMatrix<int, double>::multiply(enum Teuchos::ETransp, enum Teuchos::ETransp, double, const class Teuchos::SerialDenseMatrix<int, double> &, const class Teuchos::SerialDenseMatrix<int, double> &, double) --> int", pybind11::arg("transa"), pybind11::arg("transb"), pybind11::arg("alpha"), pybind11::arg("A"), pybind11::arg("B"), pybind11::arg("beta"));
		cl.def("__eq__", (bool (Teuchos::SerialDenseMatrix<int,double>::*)(const class Teuchos::SerialDenseMatrix<int, double> &) const) &Teuchos::SerialDenseMatrix<int, double>::operator==, "C++: Teuchos::SerialDenseMatrix<int, double>::operator==(const class Teuchos::SerialDenseMatrix<int, double> &) const --> bool", pybind11::arg("Operand"));
		cl.def("__ne__", (bool (Teuchos::SerialDenseMatrix<int,double>::*)(const class Teuchos::SerialDenseMatrix<int, double> &) const) &Teuchos::SerialDenseMatrix<int, double>::operator!=, "C++: Teuchos::SerialDenseMatrix<int, double>::operator!=(const class Teuchos::SerialDenseMatrix<int, double> &) const --> bool", pybind11::arg("Operand"));
		cl.def("numRows", (int (Teuchos::SerialDenseMatrix<int,double>::*)() const) &Teuchos::SerialDenseMatrix<int, double>::numRows, "C++: Teuchos::SerialDenseMatrix<int, double>::numRows() const --> int");
		cl.def("numCols", (int (Teuchos::SerialDenseMatrix<int,double>::*)() const) &Teuchos::SerialDenseMatrix<int, double>::numCols, "C++: Teuchos::SerialDenseMatrix<int, double>::numCols() const --> int");
		cl.def("stride", (int (Teuchos::SerialDenseMatrix<int,double>::*)() const) &Teuchos::SerialDenseMatrix<int, double>::stride, "C++: Teuchos::SerialDenseMatrix<int, double>::stride() const --> int");
		cl.def("empty", (bool (Teuchos::SerialDenseMatrix<int,double>::*)() const) &Teuchos::SerialDenseMatrix<int, double>::empty, "C++: Teuchos::SerialDenseMatrix<int, double>::empty() const --> bool");
		cl.def("normOne", (double (Teuchos::SerialDenseMatrix<int,double>::*)() const) &Teuchos::SerialDenseMatrix<int, double>::normOne, "C++: Teuchos::SerialDenseMatrix<int, double>::normOne() const --> double");
		cl.def("normInf", (double (Teuchos::SerialDenseMatrix<int,double>::*)() const) &Teuchos::SerialDenseMatrix<int, double>::normInf, "C++: Teuchos::SerialDenseMatrix<int, double>::normInf() const --> double");
		cl.def("normFrobenius", (double (Teuchos::SerialDenseMatrix<int,double>::*)() const) &Teuchos::SerialDenseMatrix<int, double>::normFrobenius, "C++: Teuchos::SerialDenseMatrix<int, double>::normFrobenius() const --> double");
		cl.def("print", (std::ostream & (Teuchos::SerialDenseMatrix<int,double>::*)(std::ostream &) const) &Teuchos::SerialDenseMatrix<int, double>::print, "C++: Teuchos::SerialDenseMatrix<int, double>::print(std::ostream &) const --> std::ostream &", pybind11::return_value_policy::automatic, pybind11::arg("os"));
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
