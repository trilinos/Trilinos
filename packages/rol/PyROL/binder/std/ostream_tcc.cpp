#include <cwchar>
#include <ios>
#include <locale>
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

void bind_std_ostream_tcc(std::function< pybind11::module &(std::string const &namespace_) > &M)
{
	{ // std::basic_ostream file:bits/ostream.tcc line:359
		pybind11::class_<std::ostream, Teuchos::RCP<std::ostream>> cl(M("std"), "ostream", "", pybind11::module_local());
		cl.def( pybind11::init<class std::basic_streambuf<char> *>(), pybind11::arg("__sb") );

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

		{ // std::basic_ostream<char, std::char_traits<char> >::sentry file:ostream line:96
			auto & enclosing_class = cl;
			pybind11::class_<std::basic_ostream<char, std::char_traits<char> >::sentry, Teuchos::RCP<std::basic_ostream<char, std::char_traits<char> >::sentry>> cl(enclosing_class, "sentry", "", pybind11::module_local());
			cl.def( pybind11::init<std::ostream &>(), pybind11::arg("__os") );

		}

	}
}
