#include <Teuchos_BLAS.hpp>
#include <Teuchos_BLAS_types.hpp>
#include <Teuchos_CompObject.hpp>
#include <Teuchos_DataAccess.hpp>
#include <Teuchos_Flops.hpp>
#include <Teuchos_LAPACK.hpp>
#include <Teuchos_Object.hpp>
#include <Teuchos_OrdinalTraits.hpp>
#include <Teuchos_ScalarTraits.hpp>
#include <Teuchos_SerialDenseMatrix.hpp>
#include <Teuchos_SerialDenseVector.hpp>
#include <complex>
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

// Teuchos::SerialDenseMatrix file:Teuchos_SerialDenseMatrix.hpp line:67
struct PyCallBack_Teuchos_SerialDenseMatrix_int_double_t : public Teuchos::SerialDenseMatrix<int,double> {
	using Teuchos::SerialDenseMatrix<int,double>::SerialDenseMatrix;

};

// Teuchos::SerialDenseVector file:Teuchos_SerialDenseVector.hpp line:60
struct PyCallBack_Teuchos_SerialDenseVector_int_double_t : public Teuchos::SerialDenseVector<int,double> {
	using Teuchos::SerialDenseVector<int,double>::SerialDenseVector;

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
	// Teuchos::throwScalarTraitsNanInfError(const std::string &) file:Teuchos_ScalarTraits.hpp line:123
	M("Teuchos").def("throwScalarTraitsNanInfError", (void (*)(const std::string &)) &Teuchos::throwScalarTraitsNanInfError, "C++: Teuchos::throwScalarTraitsNanInfError(const std::string &) --> void", pybind11::arg("errMsg"));

	// Teuchos::generic_real_isnaninf(const float &) file:Teuchos_ScalarTraits.hpp line:127
	M("Teuchos").def("generic_real_isnaninf", (bool (*)(const float &)) &Teuchos::generic_real_isnaninf<float>, "C++: Teuchos::generic_real_isnaninf(const float &) --> bool", pybind11::arg("x"));

	// Teuchos::generic_real_isnaninf(const double &) file:Teuchos_ScalarTraits.hpp line:127
	M("Teuchos").def("generic_real_isnaninf", (bool (*)(const double &)) &Teuchos::generic_real_isnaninf<double>, "C++: Teuchos::generic_real_isnaninf(const double &) --> bool", pybind11::arg("x"));

	// Teuchos::ESide file:Teuchos_BLAS_types.hpp line:88
	pybind11::enum_<Teuchos::ESide>(M("Teuchos"), "ESide", pybind11::arithmetic(), "")
		.value("LEFT_SIDE", Teuchos::LEFT_SIDE)
		.value("RIGHT_SIDE", Teuchos::RIGHT_SIDE)
		.export_values();

;

	// Teuchos::ETransp file:Teuchos_BLAS_types.hpp line:93
	pybind11::enum_<Teuchos::ETransp>(M("Teuchos"), "ETransp", pybind11::arithmetic(), "")
		.value("NO_TRANS", Teuchos::NO_TRANS)
		.value("TRANS", Teuchos::TRANS)
		.value("CONJ_TRANS", Teuchos::CONJ_TRANS)
		.export_values();

;

	// Teuchos::EUplo file:Teuchos_BLAS_types.hpp line:99
	pybind11::enum_<Teuchos::EUplo>(M("Teuchos"), "EUplo", pybind11::arithmetic(), "")
		.value("UPPER_TRI", Teuchos::UPPER_TRI)
		.value("LOWER_TRI", Teuchos::LOWER_TRI)
		.value("UNDEF_TRI", Teuchos::UNDEF_TRI)
		.export_values();

;

	// Teuchos::EDiag file:Teuchos_BLAS_types.hpp line:105
	pybind11::enum_<Teuchos::EDiag>(M("Teuchos"), "EDiag", pybind11::arithmetic(), "")
		.value("UNIT_DIAG", Teuchos::UNIT_DIAG)
		.value("NON_UNIT_DIAG", Teuchos::NON_UNIT_DIAG)
		.export_values();

;

	// Teuchos::EType file:Teuchos_BLAS_types.hpp line:110
	pybind11::enum_<Teuchos::EType>(M("Teuchos"), "EType", pybind11::arithmetic(), "")
		.value("FULL", Teuchos::FULL)
		.value("LOWER", Teuchos::LOWER)
		.value("UPPER", Teuchos::UPPER)
		.value("HESSENBERG", Teuchos::HESSENBERG)
		.value("SYM_BAND_L", Teuchos::SYM_BAND_L)
		.value("SYM_BAND_U", Teuchos::SYM_BAND_U)
		.value("BAND", Teuchos::BAND)
		.export_values();

;

	{ // Teuchos::BLAS file:Teuchos_BLAS.hpp line:2285
		pybind11::class_<Teuchos::BLAS<int,double>, Teuchos::RCP<Teuchos::BLAS<int,double>>> cl(M("Teuchos"), "BLAS_int_double_t", "");
		cl.def( pybind11::init( [](){ return new Teuchos::BLAS<int,double>(); } ) );
		cl.def( pybind11::init( [](Teuchos::BLAS<int,double> const &o){ return new Teuchos::BLAS<int,double>(o); } ) );
		cl.def("ROTG", (void (Teuchos::BLAS<int,double>::*)(double *, double *, double *, double *) const) &Teuchos::BLAS<int, double>::ROTG, "C++: Teuchos::BLAS<int, double>::ROTG(double *, double *, double *, double *) const --> void", pybind11::arg("da"), pybind11::arg("db"), pybind11::arg("c"), pybind11::arg("s"));
		cl.def("ROT", (void (Teuchos::BLAS<int,double>::*)(const int &, double *, const int &, double *, const int &, double *, double *) const) &Teuchos::BLAS<int, double>::ROT, "C++: Teuchos::BLAS<int, double>::ROT(const int &, double *, const int &, double *, const int &, double *, double *) const --> void", pybind11::arg("n"), pybind11::arg("dx"), pybind11::arg("incx"), pybind11::arg("dy"), pybind11::arg("incy"), pybind11::arg("c"), pybind11::arg("s"));
		cl.def("ASUM", (double (Teuchos::BLAS<int,double>::*)(const int &, const double *, const int &) const) &Teuchos::BLAS<int, double>::ASUM, "C++: Teuchos::BLAS<int, double>::ASUM(const int &, const double *, const int &) const --> double", pybind11::arg("n"), pybind11::arg("x"), pybind11::arg("incx"));
		cl.def("AXPY", (void (Teuchos::BLAS<int,double>::*)(const int &, const double &, const double *, const int &, double *, const int &) const) &Teuchos::BLAS<int, double>::AXPY, "C++: Teuchos::BLAS<int, double>::AXPY(const int &, const double &, const double *, const int &, double *, const int &) const --> void", pybind11::arg("n"), pybind11::arg("alpha"), pybind11::arg("x"), pybind11::arg("incx"), pybind11::arg("y"), pybind11::arg("incy"));
		cl.def("COPY", (void (Teuchos::BLAS<int,double>::*)(const int &, const double *, const int &, double *, const int &) const) &Teuchos::BLAS<int, double>::COPY, "C++: Teuchos::BLAS<int, double>::COPY(const int &, const double *, const int &, double *, const int &) const --> void", pybind11::arg("n"), pybind11::arg("x"), pybind11::arg("incx"), pybind11::arg("y"), pybind11::arg("incy"));
		cl.def("DOT", (double (Teuchos::BLAS<int,double>::*)(const int &, const double *, const int &, const double *, const int &) const) &Teuchos::BLAS<int, double>::DOT, "C++: Teuchos::BLAS<int, double>::DOT(const int &, const double *, const int &, const double *, const int &) const --> double", pybind11::arg("n"), pybind11::arg("x"), pybind11::arg("incx"), pybind11::arg("y"), pybind11::arg("incy"));
		cl.def("NRM2", (double (Teuchos::BLAS<int,double>::*)(const int &, const double *, const int &) const) &Teuchos::BLAS<int, double>::NRM2, "C++: Teuchos::BLAS<int, double>::NRM2(const int &, const double *, const int &) const --> double", pybind11::arg("n"), pybind11::arg("x"), pybind11::arg("incx"));
		cl.def("SCAL", (void (Teuchos::BLAS<int,double>::*)(const int &, const double &, double *, const int &) const) &Teuchos::BLAS<int, double>::SCAL, "C++: Teuchos::BLAS<int, double>::SCAL(const int &, const double &, double *, const int &) const --> void", pybind11::arg("n"), pybind11::arg("alpha"), pybind11::arg("x"), pybind11::arg("incx"));
		cl.def("IAMAX", (int (Teuchos::BLAS<int,double>::*)(const int &, const double *, const int &) const) &Teuchos::BLAS<int, double>::IAMAX, "C++: Teuchos::BLAS<int, double>::IAMAX(const int &, const double *, const int &) const --> int", pybind11::arg("n"), pybind11::arg("x"), pybind11::arg("incx"));
		cl.def("GEMV", (void (Teuchos::BLAS<int,double>::*)(enum Teuchos::ETransp, const int &, const int &, const double &, const double *, const int &, const double *, const int &, const double &, double *, const int &) const) &Teuchos::BLAS<int, double>::GEMV, "C++: Teuchos::BLAS<int, double>::GEMV(enum Teuchos::ETransp, const int &, const int &, const double &, const double *, const int &, const double *, const int &, const double &, double *, const int &) const --> void", pybind11::arg("trans"), pybind11::arg("m"), pybind11::arg("n"), pybind11::arg("alpha"), pybind11::arg("A"), pybind11::arg("lda"), pybind11::arg("x"), pybind11::arg("incx"), pybind11::arg("beta"), pybind11::arg("y"), pybind11::arg("incy"));
		cl.def("TRMV", (void (Teuchos::BLAS<int,double>::*)(enum Teuchos::EUplo, enum Teuchos::ETransp, enum Teuchos::EDiag, const int &, const double *, const int &, double *, const int &) const) &Teuchos::BLAS<int, double>::TRMV, "C++: Teuchos::BLAS<int, double>::TRMV(enum Teuchos::EUplo, enum Teuchos::ETransp, enum Teuchos::EDiag, const int &, const double *, const int &, double *, const int &) const --> void", pybind11::arg("uplo"), pybind11::arg("trans"), pybind11::arg("diag"), pybind11::arg("n"), pybind11::arg("A"), pybind11::arg("lda"), pybind11::arg("x"), pybind11::arg("incx"));
		cl.def("GER", (void (Teuchos::BLAS<int,double>::*)(const int &, const int &, const double &, const double *, const int &, const double *, const int &, double *, const int &) const) &Teuchos::BLAS<int, double>::GER, "C++: Teuchos::BLAS<int, double>::GER(const int &, const int &, const double &, const double *, const int &, const double *, const int &, double *, const int &) const --> void", pybind11::arg("m"), pybind11::arg("n"), pybind11::arg("alpha"), pybind11::arg("x"), pybind11::arg("incx"), pybind11::arg("y"), pybind11::arg("incy"), pybind11::arg("A"), pybind11::arg("lda"));
		cl.def("GEMM", (void (Teuchos::BLAS<int,double>::*)(enum Teuchos::ETransp, enum Teuchos::ETransp, const int &, const int &, const int &, const double &, const double *, const int &, const double *, const int &, const double &, double *, const int &) const) &Teuchos::BLAS<int, double>::GEMM, "C++: Teuchos::BLAS<int, double>::GEMM(enum Teuchos::ETransp, enum Teuchos::ETransp, const int &, const int &, const int &, const double &, const double *, const int &, const double *, const int &, const double &, double *, const int &) const --> void", pybind11::arg("transa"), pybind11::arg("transb"), pybind11::arg("m"), pybind11::arg("n"), pybind11::arg("k"), pybind11::arg("alpha"), pybind11::arg("A"), pybind11::arg("lda"), pybind11::arg("B"), pybind11::arg("ldb"), pybind11::arg("beta"), pybind11::arg("C"), pybind11::arg("ldc"));
		cl.def("SWAP", (void (Teuchos::BLAS<int,double>::*)(const int &, double *const, const int &, double *const, const int &) const) &Teuchos::BLAS<int, double>::SWAP, "C++: Teuchos::BLAS<int, double>::SWAP(const int &, double *const, const int &, double *const, const int &) const --> void", pybind11::arg("n"), pybind11::arg("x"), pybind11::arg("incx"), pybind11::arg("y"), pybind11::arg("incy"));
		cl.def("SYMM", (void (Teuchos::BLAS<int,double>::*)(enum Teuchos::ESide, enum Teuchos::EUplo, const int &, const int &, const double &, const double *, const int &, const double *, const int &, const double &, double *, const int &) const) &Teuchos::BLAS<int, double>::SYMM, "C++: Teuchos::BLAS<int, double>::SYMM(enum Teuchos::ESide, enum Teuchos::EUplo, const int &, const int &, const double &, const double *, const int &, const double *, const int &, const double &, double *, const int &) const --> void", pybind11::arg("side"), pybind11::arg("uplo"), pybind11::arg("m"), pybind11::arg("n"), pybind11::arg("alpha"), pybind11::arg("A"), pybind11::arg("lda"), pybind11::arg("B"), pybind11::arg("ldb"), pybind11::arg("beta"), pybind11::arg("C"), pybind11::arg("ldc"));
		cl.def("SYRK", (void (Teuchos::BLAS<int,double>::*)(enum Teuchos::EUplo, enum Teuchos::ETransp, const int &, const int &, const double &, const double *, const int &, const double &, double *, const int &) const) &Teuchos::BLAS<int, double>::SYRK, "C++: Teuchos::BLAS<int, double>::SYRK(enum Teuchos::EUplo, enum Teuchos::ETransp, const int &, const int &, const double &, const double *, const int &, const double &, double *, const int &) const --> void", pybind11::arg("uplo"), pybind11::arg("trans"), pybind11::arg("n"), pybind11::arg("k"), pybind11::arg("alpha"), pybind11::arg("A"), pybind11::arg("lda"), pybind11::arg("beta"), pybind11::arg("C"), pybind11::arg("ldc"));
		cl.def("HERK", (void (Teuchos::BLAS<int,double>::*)(enum Teuchos::EUplo, enum Teuchos::ETransp, const int &, const int &, const double &, const double *, const int &, const double &, double *, const int &) const) &Teuchos::BLAS<int, double>::HERK, "C++: Teuchos::BLAS<int, double>::HERK(enum Teuchos::EUplo, enum Teuchos::ETransp, const int &, const int &, const double &, const double *, const int &, const double &, double *, const int &) const --> void", pybind11::arg("uplo"), pybind11::arg("trans"), pybind11::arg("n"), pybind11::arg("k"), pybind11::arg("alpha"), pybind11::arg("A"), pybind11::arg("lda"), pybind11::arg("beta"), pybind11::arg("C"), pybind11::arg("ldc"));
		cl.def("TRMM", (void (Teuchos::BLAS<int,double>::*)(enum Teuchos::ESide, enum Teuchos::EUplo, enum Teuchos::ETransp, enum Teuchos::EDiag, const int &, const int &, const double &, const double *, const int &, double *, const int &) const) &Teuchos::BLAS<int, double>::TRMM, "C++: Teuchos::BLAS<int, double>::TRMM(enum Teuchos::ESide, enum Teuchos::EUplo, enum Teuchos::ETransp, enum Teuchos::EDiag, const int &, const int &, const double &, const double *, const int &, double *, const int &) const --> void", pybind11::arg("side"), pybind11::arg("uplo"), pybind11::arg("transa"), pybind11::arg("diag"), pybind11::arg("m"), pybind11::arg("n"), pybind11::arg("alpha"), pybind11::arg("A"), pybind11::arg("lda"), pybind11::arg("B"), pybind11::arg("ldb"));
		cl.def("TRSM", (void (Teuchos::BLAS<int,double>::*)(enum Teuchos::ESide, enum Teuchos::EUplo, enum Teuchos::ETransp, enum Teuchos::EDiag, const int &, const int &, const double &, const double *, const int &, double *, const int &) const) &Teuchos::BLAS<int, double>::TRSM, "C++: Teuchos::BLAS<int, double>::TRSM(enum Teuchos::ESide, enum Teuchos::EUplo, enum Teuchos::ETransp, enum Teuchos::EDiag, const int &, const int &, const double &, const double *, const int &, double *, const int &) const --> void", pybind11::arg("side"), pybind11::arg("uplo"), pybind11::arg("transa"), pybind11::arg("diag"), pybind11::arg("m"), pybind11::arg("n"), pybind11::arg("alpha"), pybind11::arg("A"), pybind11::arg("lda"), pybind11::arg("B"), pybind11::arg("ldb"));
		cl.def("assign", (class Teuchos::BLAS<int, double> & (Teuchos::BLAS<int,double>::*)(const class Teuchos::BLAS<int, double> &)) &Teuchos::BLAS<int, double>::operator=, "C++: Teuchos::BLAS<int, double>::operator=(const class Teuchos::BLAS<int, double> &) --> class Teuchos::BLAS<int, double> &", pybind11::return_value_policy::automatic, pybind11::arg(""));
	}
	{ // Teuchos::SerialDenseMatrix file:Teuchos_SerialDenseMatrix.hpp line:67
		pybind11::class_<Teuchos::SerialDenseMatrix<int,double>, Teuchos::RCP<Teuchos::SerialDenseMatrix<int,double>>, PyCallBack_Teuchos_SerialDenseMatrix_int_double_t, Teuchos::CompObject, Teuchos::BLAS<int,double>> cl(M("Teuchos"), "SerialDenseMatrix_int_double_t", "");
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
		cl.def("ROTG", (void (Teuchos::BLAS<int,double>::*)(double *, double *, double *, double *) const) &Teuchos::BLAS<int, double>::ROTG, "C++: Teuchos::BLAS<int, double>::ROTG(double *, double *, double *, double *) const --> void", pybind11::arg("da"), pybind11::arg("db"), pybind11::arg("c"), pybind11::arg("s"));
		cl.def("ROT", (void (Teuchos::BLAS<int,double>::*)(const int &, double *, const int &, double *, const int &, double *, double *) const) &Teuchos::BLAS<int, double>::ROT, "C++: Teuchos::BLAS<int, double>::ROT(const int &, double *, const int &, double *, const int &, double *, double *) const --> void", pybind11::arg("n"), pybind11::arg("dx"), pybind11::arg("incx"), pybind11::arg("dy"), pybind11::arg("incy"), pybind11::arg("c"), pybind11::arg("s"));
		cl.def("ASUM", (double (Teuchos::BLAS<int,double>::*)(const int &, const double *, const int &) const) &Teuchos::BLAS<int, double>::ASUM, "C++: Teuchos::BLAS<int, double>::ASUM(const int &, const double *, const int &) const --> double", pybind11::arg("n"), pybind11::arg("x"), pybind11::arg("incx"));
		cl.def("AXPY", (void (Teuchos::BLAS<int,double>::*)(const int &, const double &, const double *, const int &, double *, const int &) const) &Teuchos::BLAS<int, double>::AXPY, "C++: Teuchos::BLAS<int, double>::AXPY(const int &, const double &, const double *, const int &, double *, const int &) const --> void", pybind11::arg("n"), pybind11::arg("alpha"), pybind11::arg("x"), pybind11::arg("incx"), pybind11::arg("y"), pybind11::arg("incy"));
		cl.def("COPY", (void (Teuchos::BLAS<int,double>::*)(const int &, const double *, const int &, double *, const int &) const) &Teuchos::BLAS<int, double>::COPY, "C++: Teuchos::BLAS<int, double>::COPY(const int &, const double *, const int &, double *, const int &) const --> void", pybind11::arg("n"), pybind11::arg("x"), pybind11::arg("incx"), pybind11::arg("y"), pybind11::arg("incy"));
		cl.def("DOT", (double (Teuchos::BLAS<int,double>::*)(const int &, const double *, const int &, const double *, const int &) const) &Teuchos::BLAS<int, double>::DOT, "C++: Teuchos::BLAS<int, double>::DOT(const int &, const double *, const int &, const double *, const int &) const --> double", pybind11::arg("n"), pybind11::arg("x"), pybind11::arg("incx"), pybind11::arg("y"), pybind11::arg("incy"));
		cl.def("NRM2", (double (Teuchos::BLAS<int,double>::*)(const int &, const double *, const int &) const) &Teuchos::BLAS<int, double>::NRM2, "C++: Teuchos::BLAS<int, double>::NRM2(const int &, const double *, const int &) const --> double", pybind11::arg("n"), pybind11::arg("x"), pybind11::arg("incx"));
		cl.def("SCAL", (void (Teuchos::BLAS<int,double>::*)(const int &, const double &, double *, const int &) const) &Teuchos::BLAS<int, double>::SCAL, "C++: Teuchos::BLAS<int, double>::SCAL(const int &, const double &, double *, const int &) const --> void", pybind11::arg("n"), pybind11::arg("alpha"), pybind11::arg("x"), pybind11::arg("incx"));
		cl.def("IAMAX", (int (Teuchos::BLAS<int,double>::*)(const int &, const double *, const int &) const) &Teuchos::BLAS<int, double>::IAMAX, "C++: Teuchos::BLAS<int, double>::IAMAX(const int &, const double *, const int &) const --> int", pybind11::arg("n"), pybind11::arg("x"), pybind11::arg("incx"));
		cl.def("GEMV", (void (Teuchos::BLAS<int,double>::*)(enum Teuchos::ETransp, const int &, const int &, const double &, const double *, const int &, const double *, const int &, const double &, double *, const int &) const) &Teuchos::BLAS<int, double>::GEMV, "C++: Teuchos::BLAS<int, double>::GEMV(enum Teuchos::ETransp, const int &, const int &, const double &, const double *, const int &, const double *, const int &, const double &, double *, const int &) const --> void", pybind11::arg("trans"), pybind11::arg("m"), pybind11::arg("n"), pybind11::arg("alpha"), pybind11::arg("A"), pybind11::arg("lda"), pybind11::arg("x"), pybind11::arg("incx"), pybind11::arg("beta"), pybind11::arg("y"), pybind11::arg("incy"));
		cl.def("TRMV", (void (Teuchos::BLAS<int,double>::*)(enum Teuchos::EUplo, enum Teuchos::ETransp, enum Teuchos::EDiag, const int &, const double *, const int &, double *, const int &) const) &Teuchos::BLAS<int, double>::TRMV, "C++: Teuchos::BLAS<int, double>::TRMV(enum Teuchos::EUplo, enum Teuchos::ETransp, enum Teuchos::EDiag, const int &, const double *, const int &, double *, const int &) const --> void", pybind11::arg("uplo"), pybind11::arg("trans"), pybind11::arg("diag"), pybind11::arg("n"), pybind11::arg("A"), pybind11::arg("lda"), pybind11::arg("x"), pybind11::arg("incx"));
		cl.def("GER", (void (Teuchos::BLAS<int,double>::*)(const int &, const int &, const double &, const double *, const int &, const double *, const int &, double *, const int &) const) &Teuchos::BLAS<int, double>::GER, "C++: Teuchos::BLAS<int, double>::GER(const int &, const int &, const double &, const double *, const int &, const double *, const int &, double *, const int &) const --> void", pybind11::arg("m"), pybind11::arg("n"), pybind11::arg("alpha"), pybind11::arg("x"), pybind11::arg("incx"), pybind11::arg("y"), pybind11::arg("incy"), pybind11::arg("A"), pybind11::arg("lda"));
		cl.def("GEMM", (void (Teuchos::BLAS<int,double>::*)(enum Teuchos::ETransp, enum Teuchos::ETransp, const int &, const int &, const int &, const double &, const double *, const int &, const double *, const int &, const double &, double *, const int &) const) &Teuchos::BLAS<int, double>::GEMM, "C++: Teuchos::BLAS<int, double>::GEMM(enum Teuchos::ETransp, enum Teuchos::ETransp, const int &, const int &, const int &, const double &, const double *, const int &, const double *, const int &, const double &, double *, const int &) const --> void", pybind11::arg("transa"), pybind11::arg("transb"), pybind11::arg("m"), pybind11::arg("n"), pybind11::arg("k"), pybind11::arg("alpha"), pybind11::arg("A"), pybind11::arg("lda"), pybind11::arg("B"), pybind11::arg("ldb"), pybind11::arg("beta"), pybind11::arg("C"), pybind11::arg("ldc"));
		cl.def("SWAP", (void (Teuchos::BLAS<int,double>::*)(const int &, double *const, const int &, double *const, const int &) const) &Teuchos::BLAS<int, double>::SWAP, "C++: Teuchos::BLAS<int, double>::SWAP(const int &, double *const, const int &, double *const, const int &) const --> void", pybind11::arg("n"), pybind11::arg("x"), pybind11::arg("incx"), pybind11::arg("y"), pybind11::arg("incy"));
		cl.def("SYMM", (void (Teuchos::BLAS<int,double>::*)(enum Teuchos::ESide, enum Teuchos::EUplo, const int &, const int &, const double &, const double *, const int &, const double *, const int &, const double &, double *, const int &) const) &Teuchos::BLAS<int, double>::SYMM, "C++: Teuchos::BLAS<int, double>::SYMM(enum Teuchos::ESide, enum Teuchos::EUplo, const int &, const int &, const double &, const double *, const int &, const double *, const int &, const double &, double *, const int &) const --> void", pybind11::arg("side"), pybind11::arg("uplo"), pybind11::arg("m"), pybind11::arg("n"), pybind11::arg("alpha"), pybind11::arg("A"), pybind11::arg("lda"), pybind11::arg("B"), pybind11::arg("ldb"), pybind11::arg("beta"), pybind11::arg("C"), pybind11::arg("ldc"));
		cl.def("SYRK", (void (Teuchos::BLAS<int,double>::*)(enum Teuchos::EUplo, enum Teuchos::ETransp, const int &, const int &, const double &, const double *, const int &, const double &, double *, const int &) const) &Teuchos::BLAS<int, double>::SYRK, "C++: Teuchos::BLAS<int, double>::SYRK(enum Teuchos::EUplo, enum Teuchos::ETransp, const int &, const int &, const double &, const double *, const int &, const double &, double *, const int &) const --> void", pybind11::arg("uplo"), pybind11::arg("trans"), pybind11::arg("n"), pybind11::arg("k"), pybind11::arg("alpha"), pybind11::arg("A"), pybind11::arg("lda"), pybind11::arg("beta"), pybind11::arg("C"), pybind11::arg("ldc"));
		cl.def("HERK", (void (Teuchos::BLAS<int,double>::*)(enum Teuchos::EUplo, enum Teuchos::ETransp, const int &, const int &, const double &, const double *, const int &, const double &, double *, const int &) const) &Teuchos::BLAS<int, double>::HERK, "C++: Teuchos::BLAS<int, double>::HERK(enum Teuchos::EUplo, enum Teuchos::ETransp, const int &, const int &, const double &, const double *, const int &, const double &, double *, const int &) const --> void", pybind11::arg("uplo"), pybind11::arg("trans"), pybind11::arg("n"), pybind11::arg("k"), pybind11::arg("alpha"), pybind11::arg("A"), pybind11::arg("lda"), pybind11::arg("beta"), pybind11::arg("C"), pybind11::arg("ldc"));
		cl.def("TRMM", (void (Teuchos::BLAS<int,double>::*)(enum Teuchos::ESide, enum Teuchos::EUplo, enum Teuchos::ETransp, enum Teuchos::EDiag, const int &, const int &, const double &, const double *, const int &, double *, const int &) const) &Teuchos::BLAS<int, double>::TRMM, "C++: Teuchos::BLAS<int, double>::TRMM(enum Teuchos::ESide, enum Teuchos::EUplo, enum Teuchos::ETransp, enum Teuchos::EDiag, const int &, const int &, const double &, const double *, const int &, double *, const int &) const --> void", pybind11::arg("side"), pybind11::arg("uplo"), pybind11::arg("transa"), pybind11::arg("diag"), pybind11::arg("m"), pybind11::arg("n"), pybind11::arg("alpha"), pybind11::arg("A"), pybind11::arg("lda"), pybind11::arg("B"), pybind11::arg("ldb"));
		cl.def("TRSM", (void (Teuchos::BLAS<int,double>::*)(enum Teuchos::ESide, enum Teuchos::EUplo, enum Teuchos::ETransp, enum Teuchos::EDiag, const int &, const int &, const double &, const double *, const int &, double *, const int &) const) &Teuchos::BLAS<int, double>::TRSM, "C++: Teuchos::BLAS<int, double>::TRSM(enum Teuchos::ESide, enum Teuchos::EUplo, enum Teuchos::ETransp, enum Teuchos::EDiag, const int &, const int &, const double &, const double *, const int &, double *, const int &) const --> void", pybind11::arg("side"), pybind11::arg("uplo"), pybind11::arg("transa"), pybind11::arg("diag"), pybind11::arg("m"), pybind11::arg("n"), pybind11::arg("alpha"), pybind11::arg("A"), pybind11::arg("lda"), pybind11::arg("B"), pybind11::arg("ldb"));
		cl.def("assign", (class Teuchos::BLAS<int, double> & (Teuchos::BLAS<int,double>::*)(const class Teuchos::BLAS<int, double> &)) &Teuchos::BLAS<int, double>::operator=, "C++: Teuchos::BLAS<int, double>::operator=(const class Teuchos::BLAS<int, double> &) --> class Teuchos::BLAS<int, double> &", pybind11::return_value_policy::automatic, pybind11::arg(""));
	}
	{ // Teuchos::SerialDenseVector file:Teuchos_SerialDenseVector.hpp line:60
		pybind11::class_<Teuchos::SerialDenseVector<int,double>, Teuchos::RCP<Teuchos::SerialDenseVector<int,double>>, PyCallBack_Teuchos_SerialDenseVector_int_double_t, Teuchos::SerialDenseMatrix<int,double>> cl(M("Teuchos"), "SerialDenseVector_int_double_t", "");
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
		cl.def("ROTG", (void (Teuchos::BLAS<int,double>::*)(double *, double *, double *, double *) const) &Teuchos::BLAS<int, double>::ROTG, "C++: Teuchos::BLAS<int, double>::ROTG(double *, double *, double *, double *) const --> void", pybind11::arg("da"), pybind11::arg("db"), pybind11::arg("c"), pybind11::arg("s"));
		cl.def("ROT", (void (Teuchos::BLAS<int,double>::*)(const int &, double *, const int &, double *, const int &, double *, double *) const) &Teuchos::BLAS<int, double>::ROT, "C++: Teuchos::BLAS<int, double>::ROT(const int &, double *, const int &, double *, const int &, double *, double *) const --> void", pybind11::arg("n"), pybind11::arg("dx"), pybind11::arg("incx"), pybind11::arg("dy"), pybind11::arg("incy"), pybind11::arg("c"), pybind11::arg("s"));
		cl.def("ASUM", (double (Teuchos::BLAS<int,double>::*)(const int &, const double *, const int &) const) &Teuchos::BLAS<int, double>::ASUM, "C++: Teuchos::BLAS<int, double>::ASUM(const int &, const double *, const int &) const --> double", pybind11::arg("n"), pybind11::arg("x"), pybind11::arg("incx"));
		cl.def("AXPY", (void (Teuchos::BLAS<int,double>::*)(const int &, const double &, const double *, const int &, double *, const int &) const) &Teuchos::BLAS<int, double>::AXPY, "C++: Teuchos::BLAS<int, double>::AXPY(const int &, const double &, const double *, const int &, double *, const int &) const --> void", pybind11::arg("n"), pybind11::arg("alpha"), pybind11::arg("x"), pybind11::arg("incx"), pybind11::arg("y"), pybind11::arg("incy"));
		cl.def("COPY", (void (Teuchos::BLAS<int,double>::*)(const int &, const double *, const int &, double *, const int &) const) &Teuchos::BLAS<int, double>::COPY, "C++: Teuchos::BLAS<int, double>::COPY(const int &, const double *, const int &, double *, const int &) const --> void", pybind11::arg("n"), pybind11::arg("x"), pybind11::arg("incx"), pybind11::arg("y"), pybind11::arg("incy"));
		cl.def("DOT", (double (Teuchos::BLAS<int,double>::*)(const int &, const double *, const int &, const double *, const int &) const) &Teuchos::BLAS<int, double>::DOT, "C++: Teuchos::BLAS<int, double>::DOT(const int &, const double *, const int &, const double *, const int &) const --> double", pybind11::arg("n"), pybind11::arg("x"), pybind11::arg("incx"), pybind11::arg("y"), pybind11::arg("incy"));
		cl.def("NRM2", (double (Teuchos::BLAS<int,double>::*)(const int &, const double *, const int &) const) &Teuchos::BLAS<int, double>::NRM2, "C++: Teuchos::BLAS<int, double>::NRM2(const int &, const double *, const int &) const --> double", pybind11::arg("n"), pybind11::arg("x"), pybind11::arg("incx"));
		cl.def("SCAL", (void (Teuchos::BLAS<int,double>::*)(const int &, const double &, double *, const int &) const) &Teuchos::BLAS<int, double>::SCAL, "C++: Teuchos::BLAS<int, double>::SCAL(const int &, const double &, double *, const int &) const --> void", pybind11::arg("n"), pybind11::arg("alpha"), pybind11::arg("x"), pybind11::arg("incx"));
		cl.def("IAMAX", (int (Teuchos::BLAS<int,double>::*)(const int &, const double *, const int &) const) &Teuchos::BLAS<int, double>::IAMAX, "C++: Teuchos::BLAS<int, double>::IAMAX(const int &, const double *, const int &) const --> int", pybind11::arg("n"), pybind11::arg("x"), pybind11::arg("incx"));
		cl.def("GEMV", (void (Teuchos::BLAS<int,double>::*)(enum Teuchos::ETransp, const int &, const int &, const double &, const double *, const int &, const double *, const int &, const double &, double *, const int &) const) &Teuchos::BLAS<int, double>::GEMV, "C++: Teuchos::BLAS<int, double>::GEMV(enum Teuchos::ETransp, const int &, const int &, const double &, const double *, const int &, const double *, const int &, const double &, double *, const int &) const --> void", pybind11::arg("trans"), pybind11::arg("m"), pybind11::arg("n"), pybind11::arg("alpha"), pybind11::arg("A"), pybind11::arg("lda"), pybind11::arg("x"), pybind11::arg("incx"), pybind11::arg("beta"), pybind11::arg("y"), pybind11::arg("incy"));
		cl.def("TRMV", (void (Teuchos::BLAS<int,double>::*)(enum Teuchos::EUplo, enum Teuchos::ETransp, enum Teuchos::EDiag, const int &, const double *, const int &, double *, const int &) const) &Teuchos::BLAS<int, double>::TRMV, "C++: Teuchos::BLAS<int, double>::TRMV(enum Teuchos::EUplo, enum Teuchos::ETransp, enum Teuchos::EDiag, const int &, const double *, const int &, double *, const int &) const --> void", pybind11::arg("uplo"), pybind11::arg("trans"), pybind11::arg("diag"), pybind11::arg("n"), pybind11::arg("A"), pybind11::arg("lda"), pybind11::arg("x"), pybind11::arg("incx"));
		cl.def("GER", (void (Teuchos::BLAS<int,double>::*)(const int &, const int &, const double &, const double *, const int &, const double *, const int &, double *, const int &) const) &Teuchos::BLAS<int, double>::GER, "C++: Teuchos::BLAS<int, double>::GER(const int &, const int &, const double &, const double *, const int &, const double *, const int &, double *, const int &) const --> void", pybind11::arg("m"), pybind11::arg("n"), pybind11::arg("alpha"), pybind11::arg("x"), pybind11::arg("incx"), pybind11::arg("y"), pybind11::arg("incy"), pybind11::arg("A"), pybind11::arg("lda"));
		cl.def("GEMM", (void (Teuchos::BLAS<int,double>::*)(enum Teuchos::ETransp, enum Teuchos::ETransp, const int &, const int &, const int &, const double &, const double *, const int &, const double *, const int &, const double &, double *, const int &) const) &Teuchos::BLAS<int, double>::GEMM, "C++: Teuchos::BLAS<int, double>::GEMM(enum Teuchos::ETransp, enum Teuchos::ETransp, const int &, const int &, const int &, const double &, const double *, const int &, const double *, const int &, const double &, double *, const int &) const --> void", pybind11::arg("transa"), pybind11::arg("transb"), pybind11::arg("m"), pybind11::arg("n"), pybind11::arg("k"), pybind11::arg("alpha"), pybind11::arg("A"), pybind11::arg("lda"), pybind11::arg("B"), pybind11::arg("ldb"), pybind11::arg("beta"), pybind11::arg("C"), pybind11::arg("ldc"));
		cl.def("SWAP", (void (Teuchos::BLAS<int,double>::*)(const int &, double *const, const int &, double *const, const int &) const) &Teuchos::BLAS<int, double>::SWAP, "C++: Teuchos::BLAS<int, double>::SWAP(const int &, double *const, const int &, double *const, const int &) const --> void", pybind11::arg("n"), pybind11::arg("x"), pybind11::arg("incx"), pybind11::arg("y"), pybind11::arg("incy"));
		cl.def("SYMM", (void (Teuchos::BLAS<int,double>::*)(enum Teuchos::ESide, enum Teuchos::EUplo, const int &, const int &, const double &, const double *, const int &, const double *, const int &, const double &, double *, const int &) const) &Teuchos::BLAS<int, double>::SYMM, "C++: Teuchos::BLAS<int, double>::SYMM(enum Teuchos::ESide, enum Teuchos::EUplo, const int &, const int &, const double &, const double *, const int &, const double *, const int &, const double &, double *, const int &) const --> void", pybind11::arg("side"), pybind11::arg("uplo"), pybind11::arg("m"), pybind11::arg("n"), pybind11::arg("alpha"), pybind11::arg("A"), pybind11::arg("lda"), pybind11::arg("B"), pybind11::arg("ldb"), pybind11::arg("beta"), pybind11::arg("C"), pybind11::arg("ldc"));
		cl.def("SYRK", (void (Teuchos::BLAS<int,double>::*)(enum Teuchos::EUplo, enum Teuchos::ETransp, const int &, const int &, const double &, const double *, const int &, const double &, double *, const int &) const) &Teuchos::BLAS<int, double>::SYRK, "C++: Teuchos::BLAS<int, double>::SYRK(enum Teuchos::EUplo, enum Teuchos::ETransp, const int &, const int &, const double &, const double *, const int &, const double &, double *, const int &) const --> void", pybind11::arg("uplo"), pybind11::arg("trans"), pybind11::arg("n"), pybind11::arg("k"), pybind11::arg("alpha"), pybind11::arg("A"), pybind11::arg("lda"), pybind11::arg("beta"), pybind11::arg("C"), pybind11::arg("ldc"));
		cl.def("HERK", (void (Teuchos::BLAS<int,double>::*)(enum Teuchos::EUplo, enum Teuchos::ETransp, const int &, const int &, const double &, const double *, const int &, const double &, double *, const int &) const) &Teuchos::BLAS<int, double>::HERK, "C++: Teuchos::BLAS<int, double>::HERK(enum Teuchos::EUplo, enum Teuchos::ETransp, const int &, const int &, const double &, const double *, const int &, const double &, double *, const int &) const --> void", pybind11::arg("uplo"), pybind11::arg("trans"), pybind11::arg("n"), pybind11::arg("k"), pybind11::arg("alpha"), pybind11::arg("A"), pybind11::arg("lda"), pybind11::arg("beta"), pybind11::arg("C"), pybind11::arg("ldc"));
		cl.def("TRMM", (void (Teuchos::BLAS<int,double>::*)(enum Teuchos::ESide, enum Teuchos::EUplo, enum Teuchos::ETransp, enum Teuchos::EDiag, const int &, const int &, const double &, const double *, const int &, double *, const int &) const) &Teuchos::BLAS<int, double>::TRMM, "C++: Teuchos::BLAS<int, double>::TRMM(enum Teuchos::ESide, enum Teuchos::EUplo, enum Teuchos::ETransp, enum Teuchos::EDiag, const int &, const int &, const double &, const double *, const int &, double *, const int &) const --> void", pybind11::arg("side"), pybind11::arg("uplo"), pybind11::arg("transa"), pybind11::arg("diag"), pybind11::arg("m"), pybind11::arg("n"), pybind11::arg("alpha"), pybind11::arg("A"), pybind11::arg("lda"), pybind11::arg("B"), pybind11::arg("ldb"));
		cl.def("TRSM", (void (Teuchos::BLAS<int,double>::*)(enum Teuchos::ESide, enum Teuchos::EUplo, enum Teuchos::ETransp, enum Teuchos::EDiag, const int &, const int &, const double &, const double *, const int &, double *, const int &) const) &Teuchos::BLAS<int, double>::TRSM, "C++: Teuchos::BLAS<int, double>::TRSM(enum Teuchos::ESide, enum Teuchos::EUplo, enum Teuchos::ETransp, enum Teuchos::EDiag, const int &, const int &, const double &, const double *, const int &, double *, const int &) const --> void", pybind11::arg("side"), pybind11::arg("uplo"), pybind11::arg("transa"), pybind11::arg("diag"), pybind11::arg("m"), pybind11::arg("n"), pybind11::arg("alpha"), pybind11::arg("A"), pybind11::arg("lda"), pybind11::arg("B"), pybind11::arg("ldb"));
		cl.def("assign", (class Teuchos::BLAS<int, double> & (Teuchos::BLAS<int,double>::*)(const class Teuchos::BLAS<int, double> &)) &Teuchos::BLAS<int, double>::operator=, "C++: Teuchos::BLAS<int, double>::operator=(const class Teuchos::BLAS<int, double> &) --> class Teuchos::BLAS<int, double> &", pybind11::return_value_policy::automatic, pybind11::arg(""));
	}
	// Teuchos::convert_Fortran_complex_to_CXX_complex(_Complex double) file:Teuchos_LAPACK.hpp line:95
	M("Teuchos").def("convert_Fortran_complex_to_CXX_complex", (struct std::complex<double> (*)(_Complex double)) &Teuchos::convert_Fortran_complex_to_CXX_complex, "C++: Teuchos::convert_Fortran_complex_to_CXX_complex(_Complex double) --> struct std::complex<double>", pybind11::arg("val"));

	// Teuchos::convert_Fortran_complex_to_CXX_complex(_Complex float) file:Teuchos_LAPACK.hpp line:97
	M("Teuchos").def("convert_Fortran_complex_to_CXX_complex", (struct std::complex<float> (*)(_Complex float)) &Teuchos::convert_Fortran_complex_to_CXX_complex, "C++: Teuchos::convert_Fortran_complex_to_CXX_complex(_Complex float) --> struct std::complex<float>", pybind11::arg("val"));

}
