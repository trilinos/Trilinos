#include <complex>
#include <cwchar>
#include <ios>
#include <istream>
#include <iterator>
#include <locale>
#include <memory>
#include <sstream>
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

// std::basic_stringbuf file:bits/sstream.tcc line:291
struct PyCallBack_std_stringbuf : public std::stringbuf {
	using std::stringbuf::basic_stringbuf;

	long showmanyc() override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const std::stringbuf *>(this), "showmanyc");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<long>::value) {
				static pybind11::detail::override_caster_t<long> caster;
				return pybind11::detail::cast_ref<long>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<long>(std::move(o));
		}
		return basic_stringbuf::showmanyc();
	}
	int underflow() override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const std::stringbuf *>(this), "underflow");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<int>::value) {
				static pybind11::detail::override_caster_t<int> caster;
				return pybind11::detail::cast_ref<int>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<int>(std::move(o));
		}
		return basic_stringbuf::underflow();
	}
	int pbackfail(int a0) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const std::stringbuf *>(this), "pbackfail");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<int>::value) {
				static pybind11::detail::override_caster_t<int> caster;
				return pybind11::detail::cast_ref<int>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<int>(std::move(o));
		}
		return basic_stringbuf::pbackfail(a0);
	}
	int overflow(int a0) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const std::stringbuf *>(this), "overflow");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<int>::value) {
				static pybind11::detail::override_caster_t<int> caster;
				return pybind11::detail::cast_ref<int>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<int>(std::move(o));
		}
		return basic_stringbuf::overflow(a0);
	}
	class std::basic_streambuf<char> * setbuf(char * a0, long a1) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const std::stringbuf *>(this), "setbuf");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1);
			if (pybind11::detail::cast_is_temporary_value_reference<class std::basic_streambuf<char> *>::value) {
				static pybind11::detail::override_caster_t<class std::basic_streambuf<char> *> caster;
				return pybind11::detail::cast_ref<class std::basic_streambuf<char> *>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<class std::basic_streambuf<char> *>(std::move(o));
		}
		return basic_stringbuf::setbuf(a0, a1);
	}
	class std::fpos<__mbstate_t> seekoff(long a0, enum std::_Ios_Seekdir a1, enum std::_Ios_Openmode a2) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const std::stringbuf *>(this), "seekoff");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2);
			if (pybind11::detail::cast_is_temporary_value_reference<class std::fpos<__mbstate_t>>::value) {
				static pybind11::detail::override_caster_t<class std::fpos<__mbstate_t>> caster;
				return pybind11::detail::cast_ref<class std::fpos<__mbstate_t>>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<class std::fpos<__mbstate_t>>(std::move(o));
		}
		return basic_stringbuf::seekoff(a0, a1, a2);
	}
	class std::fpos<__mbstate_t> seekpos(class std::fpos<__mbstate_t> a0, enum std::_Ios_Openmode a1) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const std::stringbuf *>(this), "seekpos");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1);
			if (pybind11::detail::cast_is_temporary_value_reference<class std::fpos<__mbstate_t>>::value) {
				static pybind11::detail::override_caster_t<class std::fpos<__mbstate_t>> caster;
				return pybind11::detail::cast_ref<class std::fpos<__mbstate_t>>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<class std::fpos<__mbstate_t>>(std::move(o));
		}
		return basic_stringbuf::seekpos(a0, a1);
	}
	void imbue(const class std::locale & a0) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const std::stringbuf *>(this), "imbue");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return basic_streambuf::imbue(a0);
	}
	int sync() override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const std::stringbuf *>(this), "sync");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<int>::value) {
				static pybind11::detail::override_caster_t<int> caster;
				return pybind11::detail::cast_ref<int>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<int>(std::move(o));
		}
		return basic_streambuf::sync();
	}
	long xsgetn(char * a0, long a1) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const std::stringbuf *>(this), "xsgetn");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1);
			if (pybind11::detail::cast_is_temporary_value_reference<long>::value) {
				static pybind11::detail::override_caster_t<long> caster;
				return pybind11::detail::cast_ref<long>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<long>(std::move(o));
		}
		return basic_streambuf::xsgetn(a0, a1);
	}
	int uflow() override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const std::stringbuf *>(this), "uflow");
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
	long xsputn(const char * a0, long a1) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const std::stringbuf *>(this), "xsputn");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1);
			if (pybind11::detail::cast_is_temporary_value_reference<long>::value) {
				static pybind11::detail::override_caster_t<long> caster;
				return pybind11::detail::cast_ref<long>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<long>(std::move(o));
		}
		return basic_streambuf::xsputn(a0, a1);
	}
};

void bind_std_istream_tcc(std::function< pybind11::module &(std::string const &namespace_) > &M)
{
	{ // std::basic_istream file:bits/istream.tcc line:1048
		pybind11::class_<std::istream, Teuchos::RCP<std::istream>> cl(M("std"), "istream", "");
		cl.def( pybind11::init<class std::basic_streambuf<char> *>(), pybind11::arg("__sb") );

		cl.def("gcount", (long (std::istream::*)() const) &std::basic_istream<char, std::char_traits<char> >::gcount, "C++: std::basic_istream<char, std::char_traits<char> >::gcount() const --> long");
		cl.def("get", (int (std::istream::*)()) &std::basic_istream<char, std::char_traits<char> >::get, "C++: std::basic_istream<char, std::char_traits<char> >::get() --> int");
		cl.def("get", (std::istream & (std::istream::*)(char &)) &std::basic_istream<char, std::char_traits<char> >::get, "C++: std::basic_istream<char, std::char_traits<char> >::get(char &) --> std::istream &", pybind11::return_value_policy::automatic, pybind11::arg("__c"));
		cl.def("get", (std::istream & (std::istream::*)(char *, long, char)) &std::basic_istream<char, std::char_traits<char> >::get, "C++: std::basic_istream<char, std::char_traits<char> >::get(char *, long, char) --> std::istream &", pybind11::return_value_policy::automatic, pybind11::arg("__s"), pybind11::arg("__n"), pybind11::arg("__delim"));
		cl.def("get", (std::istream & (std::istream::*)(char *, long)) &std::basic_istream<char, std::char_traits<char> >::get, "C++: std::basic_istream<char, std::char_traits<char> >::get(char *, long) --> std::istream &", pybind11::return_value_policy::automatic, pybind11::arg("__s"), pybind11::arg("__n"));
		cl.def("get", (std::istream & (std::istream::*)(class std::basic_streambuf<char> &, char)) &std::basic_istream<char, std::char_traits<char> >::get, "C++: std::basic_istream<char, std::char_traits<char> >::get(class std::basic_streambuf<char> &, char) --> std::istream &", pybind11::return_value_policy::automatic, pybind11::arg("__sb"), pybind11::arg("__delim"));
		cl.def("get", (std::istream & (std::istream::*)(class std::basic_streambuf<char> &)) &std::basic_istream<char, std::char_traits<char> >::get, "C++: std::basic_istream<char, std::char_traits<char> >::get(class std::basic_streambuf<char> &) --> std::istream &", pybind11::return_value_policy::automatic, pybind11::arg("__sb"));
		cl.def("getline", (std::istream & (std::istream::*)(char *, long, char)) &std::basic_istream<char, std::char_traits<char> >::getline, "C++: std::basic_istream<char, std::char_traits<char> >::getline(char *, long, char) --> std::istream &", pybind11::return_value_policy::automatic, pybind11::arg("__s"), pybind11::arg("__n"), pybind11::arg("__delim"));
		cl.def("getline", (std::istream & (std::istream::*)(char *, long)) &std::basic_istream<char, std::char_traits<char> >::getline, "C++: std::basic_istream<char, std::char_traits<char> >::getline(char *, long) --> std::istream &", pybind11::return_value_policy::automatic, pybind11::arg("__s"), pybind11::arg("__n"));
		cl.def("ignore", (std::istream & (std::istream::*)(long, int)) &std::basic_istream<char, std::char_traits<char> >::ignore, "C++: std::basic_istream<char, std::char_traits<char> >::ignore(long, int) --> std::istream &", pybind11::return_value_policy::automatic, pybind11::arg("__n"), pybind11::arg("__delim"));
		cl.def("ignore", (std::istream & (std::istream::*)(long)) &std::basic_istream<char, std::char_traits<char> >::ignore, "C++: std::basic_istream<char, std::char_traits<char> >::ignore(long) --> std::istream &", pybind11::return_value_policy::automatic, pybind11::arg("__n"));
		cl.def("ignore", (std::istream & (std::istream::*)()) &std::basic_istream<char, std::char_traits<char> >::ignore, "C++: std::basic_istream<char, std::char_traits<char> >::ignore() --> std::istream &", pybind11::return_value_policy::automatic);
		cl.def("peek", (int (std::istream::*)()) &std::basic_istream<char, std::char_traits<char> >::peek, "C++: std::basic_istream<char, std::char_traits<char> >::peek() --> int");
		cl.def("read", (std::istream & (std::istream::*)(char *, long)) &std::basic_istream<char, std::char_traits<char> >::read, "C++: std::basic_istream<char, std::char_traits<char> >::read(char *, long) --> std::istream &", pybind11::return_value_policy::automatic, pybind11::arg("__s"), pybind11::arg("__n"));
		cl.def("readsome", (long (std::istream::*)(char *, long)) &std::basic_istream<char, std::char_traits<char> >::readsome, "C++: std::basic_istream<char, std::char_traits<char> >::readsome(char *, long) --> long", pybind11::arg("__s"), pybind11::arg("__n"));
		cl.def("putback", (std::istream & (std::istream::*)(char)) &std::basic_istream<char, std::char_traits<char> >::putback, "C++: std::basic_istream<char, std::char_traits<char> >::putback(char) --> std::istream &", pybind11::return_value_policy::automatic, pybind11::arg("__c"));
		cl.def("unget", (std::istream & (std::istream::*)()) &std::basic_istream<char, std::char_traits<char> >::unget, "C++: std::basic_istream<char, std::char_traits<char> >::unget() --> std::istream &", pybind11::return_value_policy::automatic);
		cl.def("sync", (int (std::istream::*)()) &std::basic_istream<char, std::char_traits<char> >::sync, "C++: std::basic_istream<char, std::char_traits<char> >::sync() --> int");
		cl.def("tellg", (class std::fpos<__mbstate_t> (std::istream::*)()) &std::basic_istream<char, std::char_traits<char> >::tellg, "C++: std::basic_istream<char, std::char_traits<char> >::tellg() --> class std::fpos<__mbstate_t>");
		cl.def("seekg", (std::istream & (std::istream::*)(class std::fpos<__mbstate_t>)) &std::basic_istream<char, std::char_traits<char> >::seekg, "C++: std::basic_istream<char, std::char_traits<char> >::seekg(class std::fpos<__mbstate_t>) --> std::istream &", pybind11::return_value_policy::automatic, pybind11::arg(""));
		cl.def("seekg", (std::istream & (std::istream::*)(long, enum std::_Ios_Seekdir)) &std::basic_istream<char, std::char_traits<char> >::seekg, "C++: std::basic_istream<char, std::char_traits<char> >::seekg(long, enum std::_Ios_Seekdir) --> std::istream &", pybind11::return_value_policy::automatic, pybind11::arg(""), pybind11::arg(""));

		{ // std::basic_istream<char, std::char_traits<char> >::sentry file:istream line:107
			auto & enclosing_class = cl;
			pybind11::class_<std::basic_istream<char, std::char_traits<char> >::sentry, Teuchos::RCP<std::basic_istream<char, std::char_traits<char> >::sentry>> cl(enclosing_class, "sentry", "");
			cl.def( pybind11::init( [](std::istream & a0){ return new std::basic_istream<char, std::char_traits<char> >::sentry(a0); } ), "doc" , pybind11::arg("__is"));
			cl.def( pybind11::init<std::istream &, bool>(), pybind11::arg("__is"), pybind11::arg("__noskipws") );

		}

	}
	{ // std::basic_stringbuf file:bits/sstream.tcc line:291
		pybind11::class_<std::stringbuf, Teuchos::RCP<std::stringbuf>, PyCallBack_std_stringbuf, std::streambuf> cl(M("std"), "stringbuf", "");
		cl.def( pybind11::init( [](){ return new std::stringbuf(); }, [](){ return new PyCallBack_std_stringbuf(); } ), "doc");
		cl.def( pybind11::init<enum std::_Ios_Openmode>(), pybind11::arg("__mode") );

		cl.def( pybind11::init( [](const std::string & a0){ return new std::stringbuf(a0); }, [](const std::string & a0){ return new PyCallBack_std_stringbuf(a0); } ), "doc");
		cl.def( pybind11::init<const std::string &, enum std::_Ios_Openmode>(), pybind11::arg("__str"), pybind11::arg("__mode") );

		cl.def("swap", (void (std::stringbuf::*)(class std::basic_stringbuf<char> &)) &std::basic_stringbuf<char, std::char_traits<char>, std::allocator<char> >::swap, "C++: std::basic_stringbuf<char, std::char_traits<char>, std::allocator<char> >::swap(class std::basic_stringbuf<char> &) --> void", pybind11::arg("__rhs"));
		cl.def("str", (std::string (std::stringbuf::*)() const) &std::basic_stringbuf<char, std::char_traits<char>, std::allocator<char> >::str, "C++: std::basic_stringbuf<char, std::char_traits<char>, std::allocator<char> >::str() const --> std::string");
		cl.def("str", (void (std::stringbuf::*)(const std::string &)) &std::basic_stringbuf<char, std::char_traits<char>, std::allocator<char> >::str, "C++: std::basic_stringbuf<char, std::char_traits<char>, std::allocator<char> >::str(const std::string &) --> void", pybind11::arg("__s"));
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
	{ // std::basic_istringstream file:bits/sstream.tcc line:292
		pybind11::class_<std::istringstream, Teuchos::RCP<std::istringstream>, std::istream> cl(M("std"), "istringstream", "");
		cl.def( pybind11::init( [](){ return new std::istringstream(); } ), "doc" );
		cl.def( pybind11::init<enum std::_Ios_Openmode>(), pybind11::arg("__mode") );

		cl.def( pybind11::init( [](const std::string & a0){ return new std::istringstream(a0); } ), "doc" , pybind11::arg("__str"));
		cl.def( pybind11::init<const std::string &, enum std::_Ios_Openmode>(), pybind11::arg("__str"), pybind11::arg("__mode") );

		cl.def("swap", (void (std::istringstream::*)(class std::basic_istringstream<char> &)) &std::basic_istringstream<char, std::char_traits<char>, std::allocator<char> >::swap, "C++: std::basic_istringstream<char, std::char_traits<char>, std::allocator<char> >::swap(class std::basic_istringstream<char> &) --> void", pybind11::arg("__rhs"));
		cl.def("rdbuf", (class std::basic_stringbuf<char> * (std::istringstream::*)() const) &std::basic_istringstream<char, std::char_traits<char>, std::allocator<char> >::rdbuf, "C++: std::basic_istringstream<char, std::char_traits<char>, std::allocator<char> >::rdbuf() const --> class std::basic_stringbuf<char> *", pybind11::return_value_policy::automatic);
		cl.def("str", (std::string (std::istringstream::*)() const) &std::basic_istringstream<char, std::char_traits<char>, std::allocator<char> >::str, "C++: std::basic_istringstream<char, std::char_traits<char>, std::allocator<char> >::str() const --> std::string");
		cl.def("str", (void (std::istringstream::*)(const std::string &)) &std::basic_istringstream<char, std::char_traits<char>, std::allocator<char> >::str, "C++: std::basic_istringstream<char, std::char_traits<char>, std::allocator<char> >::str(const std::string &) --> void", pybind11::arg("__s"));
		cl.def("gcount", (long (std::istream::*)() const) &std::basic_istream<char, std::char_traits<char> >::gcount, "C++: std::basic_istream<char, std::char_traits<char> >::gcount() const --> long");
		cl.def("get", (int (std::istream::*)()) &std::basic_istream<char, std::char_traits<char> >::get, "C++: std::basic_istream<char, std::char_traits<char> >::get() --> int");
		cl.def("get", (std::istream & (std::istream::*)(char &)) &std::basic_istream<char, std::char_traits<char> >::get, "C++: std::basic_istream<char, std::char_traits<char> >::get(char &) --> std::istream &", pybind11::return_value_policy::automatic, pybind11::arg("__c"));
		cl.def("get", (std::istream & (std::istream::*)(char *, long, char)) &std::basic_istream<char, std::char_traits<char> >::get, "C++: std::basic_istream<char, std::char_traits<char> >::get(char *, long, char) --> std::istream &", pybind11::return_value_policy::automatic, pybind11::arg("__s"), pybind11::arg("__n"), pybind11::arg("__delim"));
		cl.def("get", (std::istream & (std::istream::*)(char *, long)) &std::basic_istream<char, std::char_traits<char> >::get, "C++: std::basic_istream<char, std::char_traits<char> >::get(char *, long) --> std::istream &", pybind11::return_value_policy::automatic, pybind11::arg("__s"), pybind11::arg("__n"));
		cl.def("get", (std::istream & (std::istream::*)(class std::basic_streambuf<char> &, char)) &std::basic_istream<char, std::char_traits<char> >::get, "C++: std::basic_istream<char, std::char_traits<char> >::get(class std::basic_streambuf<char> &, char) --> std::istream &", pybind11::return_value_policy::automatic, pybind11::arg("__sb"), pybind11::arg("__delim"));
		cl.def("get", (std::istream & (std::istream::*)(class std::basic_streambuf<char> &)) &std::basic_istream<char, std::char_traits<char> >::get, "C++: std::basic_istream<char, std::char_traits<char> >::get(class std::basic_streambuf<char> &) --> std::istream &", pybind11::return_value_policy::automatic, pybind11::arg("__sb"));
		cl.def("getline", (std::istream & (std::istream::*)(char *, long, char)) &std::basic_istream<char, std::char_traits<char> >::getline, "C++: std::basic_istream<char, std::char_traits<char> >::getline(char *, long, char) --> std::istream &", pybind11::return_value_policy::automatic, pybind11::arg("__s"), pybind11::arg("__n"), pybind11::arg("__delim"));
		cl.def("getline", (std::istream & (std::istream::*)(char *, long)) &std::basic_istream<char, std::char_traits<char> >::getline, "C++: std::basic_istream<char, std::char_traits<char> >::getline(char *, long) --> std::istream &", pybind11::return_value_policy::automatic, pybind11::arg("__s"), pybind11::arg("__n"));
		cl.def("ignore", (std::istream & (std::istream::*)(long, int)) &std::basic_istream<char, std::char_traits<char> >::ignore, "C++: std::basic_istream<char, std::char_traits<char> >::ignore(long, int) --> std::istream &", pybind11::return_value_policy::automatic, pybind11::arg("__n"), pybind11::arg("__delim"));
		cl.def("ignore", (std::istream & (std::istream::*)(long)) &std::basic_istream<char, std::char_traits<char> >::ignore, "C++: std::basic_istream<char, std::char_traits<char> >::ignore(long) --> std::istream &", pybind11::return_value_policy::automatic, pybind11::arg("__n"));
		cl.def("ignore", (std::istream & (std::istream::*)()) &std::basic_istream<char, std::char_traits<char> >::ignore, "C++: std::basic_istream<char, std::char_traits<char> >::ignore() --> std::istream &", pybind11::return_value_policy::automatic);
		cl.def("peek", (int (std::istream::*)()) &std::basic_istream<char, std::char_traits<char> >::peek, "C++: std::basic_istream<char, std::char_traits<char> >::peek() --> int");
		cl.def("read", (std::istream & (std::istream::*)(char *, long)) &std::basic_istream<char, std::char_traits<char> >::read, "C++: std::basic_istream<char, std::char_traits<char> >::read(char *, long) --> std::istream &", pybind11::return_value_policy::automatic, pybind11::arg("__s"), pybind11::arg("__n"));
		cl.def("readsome", (long (std::istream::*)(char *, long)) &std::basic_istream<char, std::char_traits<char> >::readsome, "C++: std::basic_istream<char, std::char_traits<char> >::readsome(char *, long) --> long", pybind11::arg("__s"), pybind11::arg("__n"));
		cl.def("putback", (std::istream & (std::istream::*)(char)) &std::basic_istream<char, std::char_traits<char> >::putback, "C++: std::basic_istream<char, std::char_traits<char> >::putback(char) --> std::istream &", pybind11::return_value_policy::automatic, pybind11::arg("__c"));
		cl.def("unget", (std::istream & (std::istream::*)()) &std::basic_istream<char, std::char_traits<char> >::unget, "C++: std::basic_istream<char, std::char_traits<char> >::unget() --> std::istream &", pybind11::return_value_policy::automatic);
		cl.def("sync", (int (std::istream::*)()) &std::basic_istream<char, std::char_traits<char> >::sync, "C++: std::basic_istream<char, std::char_traits<char> >::sync() --> int");
		cl.def("tellg", (class std::fpos<__mbstate_t> (std::istream::*)()) &std::basic_istream<char, std::char_traits<char> >::tellg, "C++: std::basic_istream<char, std::char_traits<char> >::tellg() --> class std::fpos<__mbstate_t>");
		cl.def("seekg", (std::istream & (std::istream::*)(class std::fpos<__mbstate_t>)) &std::basic_istream<char, std::char_traits<char> >::seekg, "C++: std::basic_istream<char, std::char_traits<char> >::seekg(class std::fpos<__mbstate_t>) --> std::istream &", pybind11::return_value_policy::automatic, pybind11::arg(""));
		cl.def("seekg", (std::istream & (std::istream::*)(long, enum std::_Ios_Seekdir)) &std::basic_istream<char, std::char_traits<char> >::seekg, "C++: std::basic_istream<char, std::char_traits<char> >::seekg(long, enum std::_Ios_Seekdir) --> std::istream &", pybind11::return_value_policy::automatic, pybind11::arg(""), pybind11::arg(""));
	}
	{ // std::basic_ostringstream file:bits/sstream.tcc line:293
		pybind11::class_<std::ostringstream, Teuchos::RCP<std::ostringstream>> cl(M("std"), "ostringstream", "");
		cl.def( pybind11::init( [](){ return new std::ostringstream(); } ), "doc" );
		cl.def( pybind11::init<enum std::_Ios_Openmode>(), pybind11::arg("__mode") );

		cl.def( pybind11::init( [](const std::string & a0){ return new std::ostringstream(a0); } ), "doc" , pybind11::arg("__str"));
		cl.def( pybind11::init<const std::string &, enum std::_Ios_Openmode>(), pybind11::arg("__str"), pybind11::arg("__mode") );

		cl.def("swap", (void (std::ostringstream::*)(class std::basic_ostringstream<char> &)) &std::basic_ostringstream<char, std::char_traits<char>, std::allocator<char> >::swap, "C++: std::basic_ostringstream<char, std::char_traits<char>, std::allocator<char> >::swap(class std::basic_ostringstream<char> &) --> void", pybind11::arg("__rhs"));
		cl.def("rdbuf", (class std::basic_stringbuf<char> * (std::ostringstream::*)() const) &std::basic_ostringstream<char, std::char_traits<char>, std::allocator<char> >::rdbuf, "C++: std::basic_ostringstream<char, std::char_traits<char>, std::allocator<char> >::rdbuf() const --> class std::basic_stringbuf<char> *", pybind11::return_value_policy::automatic);
		cl.def("str", (std::string (std::ostringstream::*)() const) &std::basic_ostringstream<char, std::char_traits<char>, std::allocator<char> >::str, "C++: std::basic_ostringstream<char, std::char_traits<char>, std::allocator<char> >::str() const --> std::string");
		cl.def("str", (void (std::ostringstream::*)(const std::string &)) &std::basic_ostringstream<char, std::char_traits<char>, std::allocator<char> >::str, "C++: std::basic_ostringstream<char, std::char_traits<char>, std::allocator<char> >::str(const std::string &) --> void", pybind11::arg("__s"));
	}
}
