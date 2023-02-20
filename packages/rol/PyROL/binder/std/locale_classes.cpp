#include <cwchar>
#include <ios>
#include <iterator>
#include <locale>
#include <memory>
#include <sstream> // __str__
#include <stdexcept>
#include <streambuf>
#include <string>

#include <functional>
#include <pybind11/smart_holder.h>
#include <string>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>


#ifndef BINDER_PYBIND11_TYPE_CASTER
	#define BINDER_PYBIND11_TYPE_CASTER
	PYBIND11_DECLARE_HOLDER_TYPE(T, std::shared_ptr<T>)
	PYBIND11_DECLARE_HOLDER_TYPE(T, T*)
	PYBIND11_MAKE_OPAQUE(std::shared_ptr<void>)
#endif

// std::logic_error file:stdexcept line:113
struct PyCallBack_std_logic_error : public std::logic_error {
	using std::logic_error::logic_error;

	const char * what() const throw() override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const std::logic_error *>(this), "what");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<const char *>::value) {
				static pybind11::detail::override_caster_t<const char *> caster;
				return pybind11::detail::cast_ref<const char *>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<const char *>(std::move(o));
		}
		return logic_error::what();
	}
};

// std::invalid_argument file:stdexcept line:168
struct PyCallBack_std_invalid_argument : public std::invalid_argument {
	using std::invalid_argument::invalid_argument;

	const char * what() const throw() override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const std::invalid_argument *>(this), "what");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<const char *>::value) {
				static pybind11::detail::override_caster_t<const char *> caster;
				return pybind11::detail::cast_ref<const char *>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<const char *>(std::move(o));
		}
		return logic_error::what();
	}
};

// std::runtime_error file:stdexcept line:219
struct PyCallBack_std_runtime_error : public std::runtime_error {
	using std::runtime_error::runtime_error;

	const char * what() const throw() override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const std::runtime_error *>(this), "what");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<const char *>::value) {
				static pybind11::detail::override_caster_t<const char *> caster;
				return pybind11::detail::cast_ref<const char *>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<const char *>(std::move(o));
		}
		return runtime_error::what();
	}
};

// std::basic_streambuf file:bits/streambuf.tcc line:149
struct PyCallBack_std_streambuf : public std::streambuf {
	using std::streambuf::basic_streambuf;

	void imbue(const class std::locale & a0) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const std::streambuf *>(this), "imbue");
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
	class std::basic_streambuf<char> * setbuf(char * a0, long a1) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const std::streambuf *>(this), "setbuf");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1);
			if (pybind11::detail::cast_is_temporary_value_reference<class std::basic_streambuf<char> *>::value) {
				static pybind11::detail::override_caster_t<class std::basic_streambuf<char> *> caster;
				return pybind11::detail::cast_ref<class std::basic_streambuf<char> *>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<class std::basic_streambuf<char> *>(std::move(o));
		}
		return basic_streambuf::setbuf(a0, a1);
	}
	class std::fpos<__mbstate_t> seekoff(long a0, enum std::_Ios_Seekdir a1, enum std::_Ios_Openmode a2) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const std::streambuf *>(this), "seekoff");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2);
			if (pybind11::detail::cast_is_temporary_value_reference<class std::fpos<__mbstate_t>>::value) {
				static pybind11::detail::override_caster_t<class std::fpos<__mbstate_t>> caster;
				return pybind11::detail::cast_ref<class std::fpos<__mbstate_t>>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<class std::fpos<__mbstate_t>>(std::move(o));
		}
		return basic_streambuf::seekoff(a0, a1, a2);
	}
	class std::fpos<__mbstate_t> seekpos(class std::fpos<__mbstate_t> a0, enum std::_Ios_Openmode a1) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const std::streambuf *>(this), "seekpos");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1);
			if (pybind11::detail::cast_is_temporary_value_reference<class std::fpos<__mbstate_t>>::value) {
				static pybind11::detail::override_caster_t<class std::fpos<__mbstate_t>> caster;
				return pybind11::detail::cast_ref<class std::fpos<__mbstate_t>>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<class std::fpos<__mbstate_t>>(std::move(o));
		}
		return basic_streambuf::seekpos(a0, a1);
	}
	int sync() override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const std::streambuf *>(this), "sync");
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
	long showmanyc() override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const std::streambuf *>(this), "showmanyc");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<long>::value) {
				static pybind11::detail::override_caster_t<long> caster;
				return pybind11::detail::cast_ref<long>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<long>(std::move(o));
		}
		return basic_streambuf::showmanyc();
	}
	long xsgetn(char * a0, long a1) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const std::streambuf *>(this), "xsgetn");
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
	int underflow() override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const std::streambuf *>(this), "underflow");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<int>::value) {
				static pybind11::detail::override_caster_t<int> caster;
				return pybind11::detail::cast_ref<int>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<int>(std::move(o));
		}
		return basic_streambuf::underflow();
	}
	int uflow() override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const std::streambuf *>(this), "uflow");
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
	int pbackfail(int a0) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const std::streambuf *>(this), "pbackfail");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<int>::value) {
				static pybind11::detail::override_caster_t<int> caster;
				return pybind11::detail::cast_ref<int>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<int>(std::move(o));
		}
		return basic_streambuf::pbackfail(a0);
	}
	long xsputn(const char * a0, long a1) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const std::streambuf *>(this), "xsputn");
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
	int overflow(int a0) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const std::streambuf *>(this), "overflow");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<int>::value) {
				static pybind11::detail::override_caster_t<int> caster;
				return pybind11::detail::cast_ref<int>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<int>(std::move(o));
		}
		return basic_streambuf::overflow(a0);
	}
};

void bind_std_locale_classes(std::function< pybind11::module &(std::string const &namespace_) > &M)
{
	{ // std::locale file:bits/locale_classes.h line:62
		pybind11::class_<std::locale, std::shared_ptr<std::locale>> cl(M("std"), "locale", "", pybind11::module_local());
		cl.def( pybind11::init( [](){ return new std::locale(); } ) );
		cl.def( pybind11::init( [](std::locale const &o){ return new std::locale(o); } ) );
		cl.def( pybind11::init<const char *>(), pybind11::arg("__s") );

		cl.def( pybind11::init<const class std::locale &, const char *, int>(), pybind11::arg("__base"), pybind11::arg("__s"), pybind11::arg("__cat") );

		cl.def( pybind11::init<const std::string &>(), pybind11::arg("__s") );

		cl.def( pybind11::init<const class std::locale &, const std::string &, int>(), pybind11::arg("__base"), pybind11::arg("__s"), pybind11::arg("__cat") );

		cl.def( pybind11::init<const class std::locale &, const class std::locale &, int>(), pybind11::arg("__base"), pybind11::arg("__add"), pybind11::arg("__cat") );

		cl.def("assign", (const class std::locale & (std::locale::*)(const class std::locale &)) &std::locale::operator=, "C++: std::locale::operator=(const class std::locale &) --> const class std::locale &", pybind11::return_value_policy::automatic, pybind11::arg("__other"));
		cl.def("name", (std::string (std::locale::*)() const) &std::locale::name, "C++: std::locale::name() const --> std::string");
		cl.def("__eq__", (bool (std::locale::*)(const class std::locale &) const) &std::locale::operator==, "C++: std::locale::operator==(const class std::locale &) const --> bool", pybind11::arg("__other"));
		cl.def("__ne__", (bool (std::locale::*)(const class std::locale &) const) &std::locale::operator!=, "C++: std::locale::operator!=(const class std::locale &) const --> bool", pybind11::arg("__other"));
		cl.def_static("global", (class std::locale (*)(const class std::locale &)) &std::locale::global, "C++: std::locale::global(const class std::locale &) --> class std::locale", pybind11::arg("__loc"));
		cl.def_static("classic", (const class std::locale & (*)()) &std::locale::classic, "C++: std::locale::classic() --> const class std::locale &", pybind11::return_value_policy::automatic);

		{ // std::locale::id file:bits/locale_classes.h line:485
			auto & enclosing_class = cl;
			pybind11::class_<std::locale::id, std::shared_ptr<std::locale::id>> cl(enclosing_class, "id", "", pybind11::module_local());
			cl.def( pybind11::init( [](){ return new std::locale::id(); } ) );
			cl.def("_M_id", (unsigned long (std::locale::id::*)() const) &std::locale::id::_M_id, "C++: std::locale::id::_M_id() const --> unsigned long");
		}

		{ // std::locale::_Impl file:bits/locale_classes.h line:524
			auto & enclosing_class = cl;
			pybind11::class_<std::locale::_Impl, std::locale::_Impl*> cl(enclosing_class, "_Impl", "", pybind11::module_local());
		}

	}
	{ // std::logic_error file:stdexcept line:113
		pybind11::class_<std::logic_error, std::shared_ptr<std::logic_error>, PyCallBack_std_logic_error, std::exception> cl(M("std"), "logic_error", "", pybind11::module_local());
		cl.def( pybind11::init<const std::string &>(), pybind11::arg("__arg") );

		cl.def( pybind11::init<const char *>(), pybind11::arg("") );

		cl.def( pybind11::init( [](PyCallBack_std_logic_error const &o){ return new PyCallBack_std_logic_error(o); } ) );
		cl.def( pybind11::init( [](std::logic_error const &o){ return new std::logic_error(o); } ) );
		cl.def("assign", (class std::logic_error & (std::logic_error::*)(const class std::logic_error &)) &std::logic_error::operator=, "C++: std::logic_error::operator=(const class std::logic_error &) --> class std::logic_error &", pybind11::return_value_policy::automatic, pybind11::arg(""));
		cl.def("what", (const char * (std::logic_error::*)() const) &std::logic_error::what, "C++: std::logic_error::what() const --> const char *", pybind11::return_value_policy::automatic);
	}
	{ // std::invalid_argument file:stdexcept line:168
		pybind11::class_<std::invalid_argument, std::shared_ptr<std::invalid_argument>, PyCallBack_std_invalid_argument, std::logic_error> cl(M("std"), "invalid_argument", "", pybind11::module_local());
		cl.def( pybind11::init<const std::string &>(), pybind11::arg("__arg") );

		cl.def( pybind11::init<const char *>(), pybind11::arg("") );

		cl.def( pybind11::init( [](PyCallBack_std_invalid_argument const &o){ return new PyCallBack_std_invalid_argument(o); } ) );
		cl.def( pybind11::init( [](std::invalid_argument const &o){ return new std::invalid_argument(o); } ) );
		cl.def("assign", (class std::invalid_argument & (std::invalid_argument::*)(const class std::invalid_argument &)) &std::invalid_argument::operator=, "C++: std::invalid_argument::operator=(const class std::invalid_argument &) --> class std::invalid_argument &", pybind11::return_value_policy::automatic, pybind11::arg(""));
	}
	{ // std::runtime_error file:stdexcept line:219
		pybind11::class_<std::runtime_error, std::shared_ptr<std::runtime_error>, PyCallBack_std_runtime_error, std::exception> cl(M("std"), "runtime_error", "", pybind11::module_local());
		cl.def( pybind11::init<const std::string &>(), pybind11::arg("__arg") );

		cl.def( pybind11::init<const char *>(), pybind11::arg("") );

		cl.def( pybind11::init( [](PyCallBack_std_runtime_error const &o){ return new PyCallBack_std_runtime_error(o); } ) );
		cl.def( pybind11::init( [](std::runtime_error const &o){ return new std::runtime_error(o); } ) );
		cl.def("assign", (class std::runtime_error & (std::runtime_error::*)(const class std::runtime_error &)) &std::runtime_error::operator=, "C++: std::runtime_error::operator=(const class std::runtime_error &) --> class std::runtime_error &", pybind11::return_value_policy::automatic, pybind11::arg(""));
		cl.def("what", (const char * (std::runtime_error::*)() const) &std::runtime_error::what, "C++: std::runtime_error::what() const --> const char *", pybind11::return_value_policy::automatic);
	}
	// std::_Ios_Openmode file:bits/ios_base.h line:111
	pybind11::enum_<std::_Ios_Openmode>(M("std"), "_Ios_Openmode", pybind11::arithmetic(), "", pybind11::module_local())
		.value("_S_app", std::_S_app)
		.value("_S_ate", std::_S_ate)
		.value("_S_bin", std::_S_bin)
		.value("_S_in", std::_S_in)
		.value("_S_out", std::_S_out)
		.value("_S_trunc", std::_S_trunc)
		.value("_S_ios_openmode_end", std::_S_ios_openmode_end)
		.value("_S_ios_openmode_max", std::_S_ios_openmode_max)
		.value("_S_ios_openmode_min", std::_S_ios_openmode_min)
		.export_values();

;

	// std::_Ios_Seekdir file:bits/ios_base.h line:193
	pybind11::enum_<std::_Ios_Seekdir>(M("std"), "_Ios_Seekdir", pybind11::arithmetic(), "", pybind11::module_local())
		.value("_S_beg", std::_S_beg)
		.value("_S_cur", std::_S_cur)
		.value("_S_end", std::_S_end)
		.value("_S_ios_seekdir_end", std::_S_ios_seekdir_end)
		.export_values();

;

	{ // std::basic_streambuf file:bits/streambuf.tcc line:149
		pybind11::class_<std::streambuf, std::shared_ptr<std::streambuf>, PyCallBack_std_streambuf> cl(M("std"), "streambuf", "", pybind11::module_local());
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
}
