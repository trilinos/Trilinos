#include <cwchar>
#include <fstream>
#include <ios>
#include <iterator>
#include <locale>
#include <memory>
#include <sstream> // __str__
#include <streambuf>
#include <string>

#include <functional>
#include "PyROL_Smart_Holder.hpp"
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

// std::basic_filebuf file:bits/fstream.tcc line:1086
struct PyCallBack_std_filebuf : public std::filebuf {
	using std::filebuf::basic_filebuf;

	long showmanyc() override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const std::filebuf *>(this), "showmanyc");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<long>::value) {
				static pybind11::detail::override_caster_t<long> caster;
				return pybind11::detail::cast_ref<long>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<long>(std::move(o));
		}
		return basic_filebuf::showmanyc();
	}
	int underflow() override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const std::filebuf *>(this), "underflow");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<int>::value) {
				static pybind11::detail::override_caster_t<int> caster;
				return pybind11::detail::cast_ref<int>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<int>(std::move(o));
		}
		return basic_filebuf::underflow();
	}
	int pbackfail(int a0) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const std::filebuf *>(this), "pbackfail");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<int>::value) {
				static pybind11::detail::override_caster_t<int> caster;
				return pybind11::detail::cast_ref<int>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<int>(std::move(o));
		}
		return basic_filebuf::pbackfail(a0);
	}
	int overflow(int a0) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const std::filebuf *>(this), "overflow");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<int>::value) {
				static pybind11::detail::override_caster_t<int> caster;
				return pybind11::detail::cast_ref<int>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<int>(std::move(o));
		}
		return basic_filebuf::overflow(a0);
	}
	class std::basic_streambuf<char> * setbuf(char * a0, long a1) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const std::filebuf *>(this), "setbuf");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1);
			if (pybind11::detail::cast_is_temporary_value_reference<class std::basic_streambuf<char> *>::value) {
				static pybind11::detail::override_caster_t<class std::basic_streambuf<char> *> caster;
				return pybind11::detail::cast_ref<class std::basic_streambuf<char> *>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<class std::basic_streambuf<char> *>(std::move(o));
		}
		return basic_filebuf::setbuf(a0, a1);
	}
	class std::fpos<__mbstate_t> seekoff(long a0, enum std::_Ios_Seekdir a1, enum std::_Ios_Openmode a2) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const std::filebuf *>(this), "seekoff");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2);
			if (pybind11::detail::cast_is_temporary_value_reference<class std::fpos<__mbstate_t>>::value) {
				static pybind11::detail::override_caster_t<class std::fpos<__mbstate_t>> caster;
				return pybind11::detail::cast_ref<class std::fpos<__mbstate_t>>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<class std::fpos<__mbstate_t>>(std::move(o));
		}
		return basic_filebuf::seekoff(a0, a1, a2);
	}
	class std::fpos<__mbstate_t> seekpos(class std::fpos<__mbstate_t> a0, enum std::_Ios_Openmode a1) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const std::filebuf *>(this), "seekpos");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1);
			if (pybind11::detail::cast_is_temporary_value_reference<class std::fpos<__mbstate_t>>::value) {
				static pybind11::detail::override_caster_t<class std::fpos<__mbstate_t>> caster;
				return pybind11::detail::cast_ref<class std::fpos<__mbstate_t>>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<class std::fpos<__mbstate_t>>(std::move(o));
		}
		return basic_filebuf::seekpos(a0, a1);
	}
	int sync() override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const std::filebuf *>(this), "sync");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<int>::value) {
				static pybind11::detail::override_caster_t<int> caster;
				return pybind11::detail::cast_ref<int>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<int>(std::move(o));
		}
		return basic_filebuf::sync();
	}
	void imbue(const class std::locale & a0) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const std::filebuf *>(this), "imbue");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return basic_filebuf::imbue(a0);
	}
	long xsgetn(char * a0, long a1) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const std::filebuf *>(this), "xsgetn");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1);
			if (pybind11::detail::cast_is_temporary_value_reference<long>::value) {
				static pybind11::detail::override_caster_t<long> caster;
				return pybind11::detail::cast_ref<long>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<long>(std::move(o));
		}
		return basic_filebuf::xsgetn(a0, a1);
	}
	long xsputn(const char * a0, long a1) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const std::filebuf *>(this), "xsputn");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1);
			if (pybind11::detail::cast_is_temporary_value_reference<long>::value) {
				static pybind11::detail::override_caster_t<long> caster;
				return pybind11::detail::cast_ref<long>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<long>(std::move(o));
		}
		return basic_filebuf::xsputn(a0, a1);
	}
	int uflow() override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const std::filebuf *>(this), "uflow");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<int>::value) {
				static pybind11::detail::override_caster_t<int> caster;
				return pybind11::detail::cast_ref<int>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<int>(std::move(o));
		}
		return basic_streambuf::uflow();
	}
};

void bind_std_fstream_tcc(std::function< pybind11::module &(std::string const &namespace_) > &M)
{
	{ // std::basic_filebuf file:bits/fstream.tcc line:1086
		pybind11::class_<std::filebuf, Teuchos::RCP<std::filebuf>, PyCallBack_std_filebuf, std::streambuf> cl(M("std"), "filebuf", "", pybind11::module_local());
		cl.def( pybind11::init( [](){ return new std::filebuf(); }, [](){ return new PyCallBack_std_filebuf(); } ) );
		cl.def("swap", (void (std::filebuf::*)(class std::basic_filebuf<char> &)) &std::basic_filebuf<char, std::char_traits<char> >::swap, "C++: std::basic_filebuf<char, std::char_traits<char> >::swap(class std::basic_filebuf<char> &) --> void", pybind11::arg(""));
		cl.def("is_open", (bool (std::filebuf::*)() const) &std::basic_filebuf<char, std::char_traits<char> >::is_open, "C++: std::basic_filebuf<char, std::char_traits<char> >::is_open() const --> bool");
		cl.def("open", (class std::basic_filebuf<char> * (std::filebuf::*)(const char *, enum std::_Ios_Openmode)) &std::basic_filebuf<char, std::char_traits<char> >::open, "C++: std::basic_filebuf<char, std::char_traits<char> >::open(const char *, enum std::_Ios_Openmode) --> class std::basic_filebuf<char> *", pybind11::return_value_policy::automatic, pybind11::arg("__s"), pybind11::arg("__mode"));
		cl.def("open", (class std::basic_filebuf<char> * (std::filebuf::*)(const std::string &, enum std::_Ios_Openmode)) &std::basic_filebuf<char, std::char_traits<char> >::open, "C++: std::basic_filebuf<char, std::char_traits<char> >::open(const std::string &, enum std::_Ios_Openmode) --> class std::basic_filebuf<char> *", pybind11::return_value_policy::automatic, pybind11::arg("__s"), pybind11::arg("__mode"));
		cl.def("close", (class std::basic_filebuf<char> * (std::filebuf::*)()) &std::basic_filebuf<char, std::char_traits<char> >::close, "C++: std::basic_filebuf<char, std::char_traits<char> >::close() --> class std::basic_filebuf<char> *", pybind11::return_value_policy::automatic);
		cl.def("pubimbue", (class std::locale (std::streambuf::*)(const class std::locale &)) &std::basic_streambuf<char, std::char_traits<char> >::pubimbue, "C++: std::basic_streambuf<char, std::char_traits<char> >::pubimbue(const class std::locale &) --> class std::locale", pybind11::arg("__loc"));
		cl.def("getloc", (class std::locale (std::streambuf::*)() const) &std::basic_streambuf<char, std::char_traits<char> >::getloc, "C++: std::basic_streambuf<char, std::char_traits<char> >::getloc() const --> class std::locale");
		cl.def("pubsetbuf", (class std::basic_streambuf<char> * (std::streambuf::*)(char *, long)) &std::basic_streambuf<char, std::char_traits<char> >::pubsetbuf, "C++: std::basic_streambuf<char, std::char_traits<char> >::pubsetbuf(char *, long) --> class std::basic_streambuf<char> *", pybind11::return_value_policy::automatic, pybind11::arg("__s"), pybind11::arg("__n"));
		cl.def("pubseekoff", [](std::streambuf &o, long const & a0, enum std::_Ios_Seekdir const & a1) -> std::fpos<__mbstate_t> { return o.pubseekoff(a0, a1); }, "", pybind11::arg("__off"), pybind11::arg("__way"));
		cl.def("pubseekoff", (class std::fpos<__mbstate_t> (std::streambuf::*)(long, enum std::_Ios_Seekdir, enum std::_Ios_Openmode)) &std::basic_streambuf<char, std::char_traits<char> >::pubseekoff, "C++: std::basic_streambuf<char, std::char_traits<char> >::pubseekoff(long, enum std::_Ios_Seekdir, enum std::_Ios_Openmode) --> class std::fpos<__mbstate_t>", pybind11::arg("__off"), pybind11::arg("__way"), pybind11::arg("__mode"));
		cl.def("pubseekpos", [](std::streambuf &o, class std::fpos<__mbstate_t> const & a0) -> std::fpos<__mbstate_t> { return o.pubseekpos(a0); }, "", pybind11::arg("__sp"));
		cl.def("pubseekpos", (class std::fpos<__mbstate_t> (std::streambuf::*)(class std::fpos<__mbstate_t>, enum std::_Ios_Openmode)) &std::basic_streambuf<char, std::char_traits<char> >::pubseekpos, "C++: std::basic_streambuf<char, std::char_traits<char> >::pubseekpos(class std::fpos<__mbstate_t>, enum std::_Ios_Openmode) --> class std::fpos<__mbstate_t>", pybind11::arg("__sp"), pybind11::arg("__mode"));
		cl.def("pubsync", (int (std::streambuf::*)()) &std::basic_streambuf<char, std::char_traits<char> >::pubsync, "C++: std::basic_streambuf<char, std::char_traits<char> >::pubsync() --> int");
		cl.def("in_avail", (long (std::streambuf::*)()) &std::basic_streambuf<char, std::char_traits<char> >::in_avail, "C++: std::basic_streambuf<char, std::char_traits<char> >::in_avail() --> long");
		cl.def("snextc", (int (std::streambuf::*)()) &std::basic_streambuf<char, std::char_traits<char> >::snextc, "C++: std::basic_streambuf<char, std::char_traits<char> >::snextc() --> int");
		cl.def("sbumpc", (int (std::streambuf::*)()) &std::basic_streambuf<char, std::char_traits<char> >::sbumpc, "C++: std::basic_streambuf<char, std::char_traits<char> >::sbumpc() --> int");
		cl.def("sgetc", (int (std::streambuf::*)()) &std::basic_streambuf<char, std::char_traits<char> >::sgetc, "C++: std::basic_streambuf<char, std::char_traits<char> >::sgetc() --> int");
		cl.def("sgetn", (long (std::streambuf::*)(char *, long)) &std::basic_streambuf<char, std::char_traits<char> >::sgetn, "C++: std::basic_streambuf<char, std::char_traits<char> >::sgetn(char *, long) --> long", pybind11::arg("__s"), pybind11::arg("__n"));
		cl.def("sputbackc", (int (std::streambuf::*)(char)) &std::basic_streambuf<char, std::char_traits<char> >::sputbackc, "C++: std::basic_streambuf<char, std::char_traits<char> >::sputbackc(char) --> int", pybind11::arg("__c"));
		cl.def("sungetc", (int (std::streambuf::*)()) &std::basic_streambuf<char, std::char_traits<char> >::sungetc, "C++: std::basic_streambuf<char, std::char_traits<char> >::sungetc() --> int");
		cl.def("sputc", (int (std::streambuf::*)(char)) &std::basic_streambuf<char, std::char_traits<char> >::sputc, "C++: std::basic_streambuf<char, std::char_traits<char> >::sputc(char) --> int", pybind11::arg("__c"));
		cl.def("sputn", (long (std::streambuf::*)(const char *, long)) &std::basic_streambuf<char, std::char_traits<char> >::sputn, "C++: std::basic_streambuf<char, std::char_traits<char> >::sputn(const char *, long) --> long", pybind11::arg("__s"), pybind11::arg("__n"));
		cl.def("stossc", (void (std::streambuf::*)()) &std::basic_streambuf<char, std::char_traits<char> >::stossc, "C++: std::basic_streambuf<char, std::char_traits<char> >::stossc() --> void");
		cl.def("__safe_gbump", (void (std::streambuf::*)(long)) &std::basic_streambuf<char, std::char_traits<char> >::__safe_gbump, "C++: std::basic_streambuf<char, std::char_traits<char> >::__safe_gbump(long) --> void", pybind11::arg("__n"));
		cl.def("__safe_pbump", (void (std::streambuf::*)(long)) &std::basic_streambuf<char, std::char_traits<char> >::__safe_pbump, "C++: std::basic_streambuf<char, std::char_traits<char> >::__safe_pbump(long) --> void", pybind11::arg("__n"));
	}
	{ // std::basic_ofstream file:bits/fstream.tcc line:1088
		pybind11::class_<std::ofstream, Teuchos::RCP<std::ofstream>, std::ostream> cl(M("std"), "ofstream", "", pybind11::module_local());
		cl.def( pybind11::init( [](){ return new std::ofstream(); } ) );
		cl.def( pybind11::init( [](const char * a0){ return new std::ofstream(a0); } ), "doc" , pybind11::arg("__s"));
		cl.def( pybind11::init<const char *, enum std::_Ios_Openmode>(), pybind11::arg("__s"), pybind11::arg("__mode") );

		cl.def( pybind11::init( [](const std::string & a0){ return new std::ofstream(a0); } ), "doc" , pybind11::arg("__s"));
		cl.def( pybind11::init<const std::string &, enum std::_Ios_Openmode>(), pybind11::arg("__s"), pybind11::arg("__mode") );

		cl.def("swap", (void (std::ofstream::*)(class std::basic_ofstream<char> &)) &std::basic_ofstream<char, std::char_traits<char> >::swap, "C++: std::basic_ofstream<char, std::char_traits<char> >::swap(class std::basic_ofstream<char> &) --> void", pybind11::arg("__rhs"));
		cl.def("rdbuf", (class std::basic_filebuf<char> * (std::ofstream::*)() const) &std::basic_ofstream<char, std::char_traits<char> >::rdbuf, "C++: std::basic_ofstream<char, std::char_traits<char> >::rdbuf() const --> class std::basic_filebuf<char> *", pybind11::return_value_policy::automatic);
		cl.def("is_open", (bool (std::ofstream::*)()) &std::basic_ofstream<char, std::char_traits<char> >::is_open, "C++: std::basic_ofstream<char, std::char_traits<char> >::is_open() --> bool");
		cl.def("open", [](std::ofstream &o, const char * a0) -> void { return o.open(a0); }, "", pybind11::arg("__s"));
		cl.def("open", (void (std::ofstream::*)(const char *, enum std::_Ios_Openmode)) &std::basic_ofstream<char, std::char_traits<char> >::open, "C++: std::basic_ofstream<char, std::char_traits<char> >::open(const char *, enum std::_Ios_Openmode) --> void", pybind11::arg("__s"), pybind11::arg("__mode"));
		cl.def("open", [](std::ofstream &o, const std::string & a0) -> void { return o.open(a0); }, "", pybind11::arg("__s"));
		cl.def("open", (void (std::ofstream::*)(const std::string &, enum std::_Ios_Openmode)) &std::basic_ofstream<char, std::char_traits<char> >::open, "C++: std::basic_ofstream<char, std::char_traits<char> >::open(const std::string &, enum std::_Ios_Openmode) --> void", pybind11::arg("__s"), pybind11::arg("__mode"));
		cl.def("close", (void (std::ofstream::*)()) &std::basic_ofstream<char, std::char_traits<char> >::close, "C++: std::basic_ofstream<char, std::char_traits<char> >::close() --> void");
		cl.def("__lshift__", (std::ostream & (std::ostream::*)(long)) &std::basic_ostream<char, std::char_traits<char> >::operator<<, "C++: std::basic_ostream<char, std::char_traits<char> >::operator<<(long) --> std::ostream &", pybind11::return_value_policy::automatic, pybind11::arg("__n"));
		cl.def("__lshift__", (std::ostream & (std::ostream::*)(unsigned long)) &std::basic_ostream<char, std::char_traits<char> >::operator<<, "C++: std::basic_ostream<char, std::char_traits<char> >::operator<<(unsigned long) --> std::ostream &", pybind11::return_value_policy::automatic, pybind11::arg("__n"));
		cl.def("__lshift__", (std::ostream & (std::ostream::*)(bool)) &std::basic_ostream<char, std::char_traits<char> >::operator<<, "C++: std::basic_ostream<char, std::char_traits<char> >::operator<<(bool) --> std::ostream &", pybind11::return_value_policy::automatic, pybind11::arg("__n"));
		cl.def("__lshift__", (std::ostream & (std::ostream::*)(short)) &std::basic_ostream<char, std::char_traits<char> >::operator<<, "C++: std::basic_ostream<char, std::char_traits<char> >::operator<<(short) --> std::ostream &", pybind11::return_value_policy::automatic, pybind11::arg("__n"));
		cl.def("__lshift__", (std::ostream & (std::ostream::*)(unsigned short)) &std::basic_ostream<char, std::char_traits<char> >::operator<<, "C++: std::basic_ostream<char, std::char_traits<char> >::operator<<(unsigned short) --> std::ostream &", pybind11::return_value_policy::automatic, pybind11::arg("__n"));
		cl.def("__lshift__", (std::ostream & (std::ostream::*)(int)) &std::basic_ostream<char, std::char_traits<char> >::operator<<, "C++: std::basic_ostream<char, std::char_traits<char> >::operator<<(int) --> std::ostream &", pybind11::return_value_policy::automatic, pybind11::arg("__n"));
		cl.def("__lshift__", (std::ostream & (std::ostream::*)(unsigned int)) &std::basic_ostream<char, std::char_traits<char> >::operator<<, "C++: std::basic_ostream<char, std::char_traits<char> >::operator<<(unsigned int) --> std::ostream &", pybind11::return_value_policy::automatic, pybind11::arg("__n"));
		cl.def("__lshift__", (std::ostream & (std::ostream::*)(long long)) &std::basic_ostream<char, std::char_traits<char> >::operator<<, "C++: std::basic_ostream<char, std::char_traits<char> >::operator<<(long long) --> std::ostream &", pybind11::return_value_policy::automatic, pybind11::arg("__n"));
		cl.def("__lshift__", (std::ostream & (std::ostream::*)(unsigned long long)) &std::basic_ostream<char, std::char_traits<char> >::operator<<, "C++: std::basic_ostream<char, std::char_traits<char> >::operator<<(unsigned long long) --> std::ostream &", pybind11::return_value_policy::automatic, pybind11::arg("__n"));
		cl.def("__lshift__", (std::ostream & (std::ostream::*)(double)) &std::basic_ostream<char, std::char_traits<char> >::operator<<, "C++: std::basic_ostream<char, std::char_traits<char> >::operator<<(double) --> std::ostream &", pybind11::return_value_policy::automatic, pybind11::arg("__f"));
		cl.def("__lshift__", (std::ostream & (std::ostream::*)(float)) &std::basic_ostream<char, std::char_traits<char> >::operator<<, "C++: std::basic_ostream<char, std::char_traits<char> >::operator<<(float) --> std::ostream &", pybind11::return_value_policy::automatic, pybind11::arg("__f"));
		cl.def("__lshift__", (std::ostream & (std::ostream::*)(long double)) &std::basic_ostream<char, std::char_traits<char> >::operator<<, "C++: std::basic_ostream<char, std::char_traits<char> >::operator<<(long double) --> std::ostream &", pybind11::return_value_policy::automatic, pybind11::arg("__f"));
		cl.def("__lshift__", (std::ostream & (std::ostream::*)(const void *)) &std::basic_ostream<char, std::char_traits<char> >::operator<<, "C++: std::basic_ostream<char, std::char_traits<char> >::operator<<(const void *) --> std::ostream &", pybind11::return_value_policy::automatic, pybind11::arg("__p"));
		cl.def("__lshift__", (std::ostream & (std::ostream::*)(class std::basic_streambuf<char> *)) &std::basic_ostream<char, std::char_traits<char> >::operator<<, "C++: std::basic_ostream<char, std::char_traits<char> >::operator<<(class std::basic_streambuf<char> *) --> std::ostream &", pybind11::return_value_policy::automatic, pybind11::arg("__sb"));
		cl.def("put", (std::ostream & (std::ostream::*)(char)) &std::basic_ostream<char, std::char_traits<char> >::put, "C++: std::basic_ostream<char, std::char_traits<char> >::put(char) --> std::ostream &", pybind11::return_value_policy::automatic, pybind11::arg("__c"));
		cl.def("_M_write", (void (std::ostream::*)(const char *, long)) &std::basic_ostream<char, std::char_traits<char> >::_M_write, "C++: std::basic_ostream<char, std::char_traits<char> >::_M_write(const char *, long) --> void", pybind11::arg("__s"), pybind11::arg("__n"));
		cl.def("write", (std::ostream & (std::ostream::*)(const char *, long)) &std::basic_ostream<char, std::char_traits<char> >::write, "C++: std::basic_ostream<char, std::char_traits<char> >::write(const char *, long) --> std::ostream &", pybind11::return_value_policy::automatic, pybind11::arg("__s"), pybind11::arg("__n"));
		cl.def("flush", (std::ostream & (std::ostream::*)()) &std::basic_ostream<char, std::char_traits<char> >::flush, "C++: std::basic_ostream<char, std::char_traits<char> >::flush() --> std::ostream &", pybind11::return_value_policy::automatic);
		cl.def("tellp", (class std::fpos<__mbstate_t> (std::ostream::*)()) &std::basic_ostream<char, std::char_traits<char> >::tellp, "C++: std::basic_ostream<char, std::char_traits<char> >::tellp() --> class std::fpos<__mbstate_t>");
		cl.def("seekp", (std::ostream & (std::ostream::*)(class std::fpos<__mbstate_t>)) &std::basic_ostream<char, std::char_traits<char> >::seekp, "C++: std::basic_ostream<char, std::char_traits<char> >::seekp(class std::fpos<__mbstate_t>) --> std::ostream &", pybind11::return_value_policy::automatic, pybind11::arg(""));
		cl.def("seekp", (std::ostream & (std::ostream::*)(long, enum std::_Ios_Seekdir)) &std::basic_ostream<char, std::char_traits<char> >::seekp, "C++: std::basic_ostream<char, std::char_traits<char> >::seekp(long, enum std::_Ios_Seekdir) --> std::ostream &", pybind11::return_value_policy::automatic, pybind11::arg(""), pybind11::arg(""));
	}
}
