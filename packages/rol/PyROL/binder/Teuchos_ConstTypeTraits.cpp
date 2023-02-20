#include <PyROL_Teuchos_Custom.hpp>
#include <Teuchos_ConstTypeTraits.hpp>
#include <Teuchos_FancyOStream.hpp>
#include <Teuchos_RCPDecl.hpp>
#include <Teuchos_TestForException.hpp>
#include <Teuchos_TypeNameTraits.hpp>
#include <Teuchos_any.hpp>
#include <cwchar>
#include <ios>
#include <iterator>
#include <memory>
#include <ostream>
#include <sstream> // __str__
#include <streambuf>
#include <string>
#include <typeinfo>

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

void bind_Teuchos_ConstTypeTraits(std::function< pybind11::module &(std::string const &namespace_) > &M)
{
	{ // Teuchos::ConstTypeTraits file:Teuchos_ConstTypeTraits.hpp line:59
		pybind11::class_<Teuchos::ConstTypeTraits<float>, std::shared_ptr<Teuchos::ConstTypeTraits<float>>> cl(M("Teuchos"), "ConstTypeTraits_float_t", "", pybind11::module_local());
		cl.def( pybind11::init( [](){ return new Teuchos::ConstTypeTraits<float>(); } ) );
	}
	{ // Teuchos::ConstTypeTraits file:Teuchos_ConstTypeTraits.hpp line:59
		pybind11::class_<Teuchos::ConstTypeTraits<double>, std::shared_ptr<Teuchos::ConstTypeTraits<double>>> cl(M("Teuchos"), "ConstTypeTraits_double_t", "", pybind11::module_local());
		cl.def( pybind11::init( [](){ return new Teuchos::ConstTypeTraits<double>(); } ) );
	}
	{ // Teuchos::ConstTypeTraits file:Teuchos_ConstTypeTraits.hpp line:59
		pybind11::class_<Teuchos::ConstTypeTraits<Teuchos::any::placeholder>, std::shared_ptr<Teuchos::ConstTypeTraits<Teuchos::any::placeholder>>> cl(M("Teuchos"), "ConstTypeTraits_Teuchos_any_placeholder_t", "", pybind11::module_local());
		cl.def( pybind11::init( [](){ return new Teuchos::ConstTypeTraits<Teuchos::any::placeholder>(); } ) );
	}
	// Teuchos::demangleName(const std::string &) file:Teuchos_TypeNameTraits.hpp line:77
	M("Teuchos").def("demangleName", (std::string (*)(const std::string &)) &Teuchos::demangleName, "Demangle a C++ name if valid.\n\n The name must have come from typeid(...).name() in order to be\n valid name to pass to this function.\n\n \n\n \n\nC++: Teuchos::demangleName(const std::string &) --> std::string", pybind11::arg("mangledName"));

	{ // Teuchos::TypeNameTraits file:Teuchos_TypeNameTraits.hpp line:85
		pybind11::class_<Teuchos::TypeNameTraits<long long>, std::shared_ptr<Teuchos::TypeNameTraits<long long>>> cl(M("Teuchos"), "TypeNameTraits_long_long_t", "", pybind11::module_local());
		cl.def( pybind11::init( [](){ return new Teuchos::TypeNameTraits<long long>(); } ) );
		cl.def_static("name", (std::string (*)()) &Teuchos::TypeNameTraits<long long>::name, "C++: Teuchos::TypeNameTraits<long long>::name() --> std::string");
		cl.def_static("concreteName", (std::string (*)(const long long &)) &Teuchos::TypeNameTraits<long long>::concreteName, "C++: Teuchos::TypeNameTraits<long long>::concreteName(const long long &) --> std::string", pybind11::arg("t"));
	}
	{ // Teuchos::TypeNameTraits file:Teuchos_TypeNameTraits.hpp line:85
		pybind11::class_<Teuchos::TypeNameTraits<unsigned long long>, std::shared_ptr<Teuchos::TypeNameTraits<unsigned long long>>> cl(M("Teuchos"), "TypeNameTraits_unsigned_long_long_t", "", pybind11::module_local());
		cl.def( pybind11::init( [](){ return new Teuchos::TypeNameTraits<unsigned long long>(); } ) );
		cl.def_static("name", (std::string (*)()) &Teuchos::TypeNameTraits<unsigned long long>::name, "C++: Teuchos::TypeNameTraits<unsigned long long>::name() --> std::string");
		cl.def_static("concreteName", (std::string (*)(const unsigned long long &)) &Teuchos::TypeNameTraits<unsigned long long>::concreteName, "C++: Teuchos::TypeNameTraits<unsigned long long>::concreteName(const unsigned long long &) --> std::string", pybind11::arg("t"));
	}
	{ // Teuchos::TypeNameTraits file:Teuchos_TypeNameTraits.hpp line:85
		pybind11::class_<Teuchos::TypeNameTraits<Teuchos::basic_FancyOStream<char, std::char_traits<char> >>, std::shared_ptr<Teuchos::TypeNameTraits<Teuchos::basic_FancyOStream<char, std::char_traits<char> >>>> cl(M("Teuchos"), "TypeNameTraits_Teuchos_basic_FancyOStream_char_std_char_traits_char_t", "", pybind11::module_local());
		cl.def( pybind11::init( [](){ return new Teuchos::TypeNameTraits<Teuchos::basic_FancyOStream<char, std::char_traits<char> >>(); } ) );
		cl.def_static("name", (std::string (*)()) &Teuchos::TypeNameTraits<Teuchos::basic_FancyOStream<char, std::char_traits<char> > >::name, "C++: Teuchos::TypeNameTraits<Teuchos::basic_FancyOStream<char, std::char_traits<char> > >::name() --> std::string");
		cl.def_static("concreteName", (std::string (*)(const class Teuchos::basic_FancyOStream<char, struct std::char_traits<char> > &)) &Teuchos::TypeNameTraits<Teuchos::basic_FancyOStream<char, std::char_traits<char> > >::concreteName, "C++: Teuchos::TypeNameTraits<Teuchos::basic_FancyOStream<char, std::char_traits<char> > >::concreteName(const class Teuchos::basic_FancyOStream<char, struct std::char_traits<char> > &) --> std::string", pybind11::arg("t"));
	}
	{ // Teuchos::TypeNameTraits file:Teuchos_TypeNameTraits.hpp line:85
		pybind11::class_<Teuchos::TypeNameTraits<Teuchos::any::placeholder>, std::shared_ptr<Teuchos::TypeNameTraits<Teuchos::any::placeholder>>> cl(M("Teuchos"), "TypeNameTraits_Teuchos_any_placeholder_t", "", pybind11::module_local());
		cl.def( pybind11::init( [](){ return new Teuchos::TypeNameTraits<Teuchos::any::placeholder>(); } ) );
		cl.def_static("name", (std::string (*)()) &Teuchos::TypeNameTraits<Teuchos::any::placeholder>::name, "C++: Teuchos::TypeNameTraits<Teuchos::any::placeholder>::name() --> std::string");
		cl.def_static("concreteName", (std::string (*)(const class Teuchos::any::placeholder &)) &Teuchos::TypeNameTraits<Teuchos::any::placeholder>::concreteName, "C++: Teuchos::TypeNameTraits<Teuchos::any::placeholder>::concreteName(const class Teuchos::any::placeholder &) --> std::string", pybind11::arg("t"));
	}
	// Teuchos::typeName(const class Teuchos::any::placeholder &) file:Teuchos_TypeNameTraits.hpp line:115
	M("Teuchos").def("typeName", (std::string (*)(const class Teuchos::any::placeholder &)) &Teuchos::typeName<Teuchos::any::placeholder>, "C++: Teuchos::typeName(const class Teuchos::any::placeholder &) --> std::string", pybind11::arg("t"));

	{ // Teuchos::TypeNameTraits file: line:147
		pybind11::class_<Teuchos::TypeNameTraits<bool>, std::shared_ptr<Teuchos::TypeNameTraits<bool>>> cl(M("Teuchos"), "TypeNameTraits_bool_t", "", pybind11::module_local());
		cl.def( pybind11::init( [](){ return new Teuchos::TypeNameTraits<bool>(); } ) );
		cl.def_static("name", (std::string (*)()) &Teuchos::TypeNameTraits<bool>::name, "C++: Teuchos::TypeNameTraits<bool>::name() --> std::string");
		cl.def_static("concreteName", (std::string (*)(const bool &)) &Teuchos::TypeNameTraits<bool>::concreteName, "C++: Teuchos::TypeNameTraits<bool>::concreteName(const bool &) --> std::string", pybind11::arg(""));
	}
	{ // Teuchos::TypeNameTraits file: line:147
		pybind11::class_<Teuchos::TypeNameTraits<char>, std::shared_ptr<Teuchos::TypeNameTraits<char>>> cl(M("Teuchos"), "TypeNameTraits_char_t", "", pybind11::module_local());
		cl.def( pybind11::init( [](){ return new Teuchos::TypeNameTraits<char>(); } ) );
		cl.def_static("name", (std::string (*)()) &Teuchos::TypeNameTraits<char>::name, "C++: Teuchos::TypeNameTraits<char>::name() --> std::string");
		cl.def_static("concreteName", (std::string (*)(const char &)) &Teuchos::TypeNameTraits<char>::concreteName, "C++: Teuchos::TypeNameTraits<char>::concreteName(const char &) --> std::string", pybind11::arg(""));
	}
	{ // Teuchos::TypeNameTraits file: line:147
		pybind11::class_<Teuchos::TypeNameTraits<signed char>, std::shared_ptr<Teuchos::TypeNameTraits<signed char>>> cl(M("Teuchos"), "TypeNameTraits_signed_char_t", "", pybind11::module_local());
		cl.def( pybind11::init( [](){ return new Teuchos::TypeNameTraits<signed char>(); } ) );
		cl.def_static("name", (std::string (*)()) &Teuchos::TypeNameTraits<signed char>::name, "C++: Teuchos::TypeNameTraits<signed char>::name() --> std::string");
		cl.def_static("concreteName", (std::string (*)(const signed char &)) &Teuchos::TypeNameTraits<signed char>::concreteName, "C++: Teuchos::TypeNameTraits<signed char>::concreteName(const signed char &) --> std::string", pybind11::arg(""));
	}
	{ // Teuchos::TypeNameTraits file: line:147
		pybind11::class_<Teuchos::TypeNameTraits<unsigned char>, std::shared_ptr<Teuchos::TypeNameTraits<unsigned char>>> cl(M("Teuchos"), "TypeNameTraits_unsigned_char_t", "", pybind11::module_local());
		cl.def( pybind11::init( [](){ return new Teuchos::TypeNameTraits<unsigned char>(); } ) );
		cl.def_static("name", (std::string (*)()) &Teuchos::TypeNameTraits<unsigned char>::name, "C++: Teuchos::TypeNameTraits<unsigned char>::name() --> std::string");
		cl.def_static("concreteName", (std::string (*)(const unsigned char &)) &Teuchos::TypeNameTraits<unsigned char>::concreteName, "C++: Teuchos::TypeNameTraits<unsigned char>::concreteName(const unsigned char &) --> std::string", pybind11::arg(""));
	}
	{ // Teuchos::TypeNameTraits file: line:147
		pybind11::class_<Teuchos::TypeNameTraits<short>, std::shared_ptr<Teuchos::TypeNameTraits<short>>> cl(M("Teuchos"), "TypeNameTraits_short_t", "", pybind11::module_local());
		cl.def( pybind11::init( [](){ return new Teuchos::TypeNameTraits<short>(); } ) );
		cl.def_static("name", (std::string (*)()) &Teuchos::TypeNameTraits<short>::name, "C++: Teuchos::TypeNameTraits<short>::name() --> std::string");
		cl.def_static("concreteName", (std::string (*)(const short &)) &Teuchos::TypeNameTraits<short>::concreteName, "C++: Teuchos::TypeNameTraits<short>::concreteName(const short &) --> std::string", pybind11::arg(""));
	}
	{ // Teuchos::TypeNameTraits file: line:147
		pybind11::class_<Teuchos::TypeNameTraits<int>, std::shared_ptr<Teuchos::TypeNameTraits<int>>> cl(M("Teuchos"), "TypeNameTraits_int_t", "", pybind11::module_local());
		cl.def( pybind11::init( [](){ return new Teuchos::TypeNameTraits<int>(); } ) );
		cl.def_static("name", (std::string (*)()) &Teuchos::TypeNameTraits<int>::name, "C++: Teuchos::TypeNameTraits<int>::name() --> std::string");
		cl.def_static("concreteName", (std::string (*)(const int &)) &Teuchos::TypeNameTraits<int>::concreteName, "C++: Teuchos::TypeNameTraits<int>::concreteName(const int &) --> std::string", pybind11::arg(""));
	}
	{ // Teuchos::TypeNameTraits file: line:147
		pybind11::class_<Teuchos::TypeNameTraits<long>, std::shared_ptr<Teuchos::TypeNameTraits<long>>> cl(M("Teuchos"), "TypeNameTraits_long_t", "", pybind11::module_local());
		cl.def( pybind11::init( [](){ return new Teuchos::TypeNameTraits<long>(); } ) );
		cl.def_static("name", (std::string (*)()) &Teuchos::TypeNameTraits<long>::name, "C++: Teuchos::TypeNameTraits<long>::name() --> std::string");
		cl.def_static("concreteName", (std::string (*)(const long &)) &Teuchos::TypeNameTraits<long>::concreteName, "C++: Teuchos::TypeNameTraits<long>::concreteName(const long &) --> std::string", pybind11::arg(""));
	}
	{ // Teuchos::TypeNameTraits file: line:147
		pybind11::class_<Teuchos::TypeNameTraits<unsigned short>, std::shared_ptr<Teuchos::TypeNameTraits<unsigned short>>> cl(M("Teuchos"), "TypeNameTraits_unsigned_short_t", "", pybind11::module_local());
		cl.def( pybind11::init( [](){ return new Teuchos::TypeNameTraits<unsigned short>(); } ) );
		cl.def_static("name", (std::string (*)()) &Teuchos::TypeNameTraits<unsigned short>::name, "C++: Teuchos::TypeNameTraits<unsigned short>::name() --> std::string");
		cl.def_static("concreteName", (std::string (*)(const unsigned short &)) &Teuchos::TypeNameTraits<unsigned short>::concreteName, "C++: Teuchos::TypeNameTraits<unsigned short>::concreteName(const unsigned short &) --> std::string", pybind11::arg(""));
	}
	{ // Teuchos::TypeNameTraits file: line:147
		pybind11::class_<Teuchos::TypeNameTraits<unsigned int>, std::shared_ptr<Teuchos::TypeNameTraits<unsigned int>>> cl(M("Teuchos"), "TypeNameTraits_unsigned_int_t", "", pybind11::module_local());
		cl.def( pybind11::init( [](){ return new Teuchos::TypeNameTraits<unsigned int>(); } ) );
		cl.def_static("name", (std::string (*)()) &Teuchos::TypeNameTraits<unsigned int>::name, "C++: Teuchos::TypeNameTraits<unsigned int>::name() --> std::string");
		cl.def_static("concreteName", (std::string (*)(const unsigned int &)) &Teuchos::TypeNameTraits<unsigned int>::concreteName, "C++: Teuchos::TypeNameTraits<unsigned int>::concreteName(const unsigned int &) --> std::string", pybind11::arg(""));
	}
	{ // Teuchos::TypeNameTraits file: line:147
		pybind11::class_<Teuchos::TypeNameTraits<unsigned long>, std::shared_ptr<Teuchos::TypeNameTraits<unsigned long>>> cl(M("Teuchos"), "TypeNameTraits_unsigned_long_t", "", pybind11::module_local());
		cl.def( pybind11::init( [](){ return new Teuchos::TypeNameTraits<unsigned long>(); } ) );
		cl.def_static("name", (std::string (*)()) &Teuchos::TypeNameTraits<unsigned long>::name, "C++: Teuchos::TypeNameTraits<unsigned long>::name() --> std::string");
		cl.def_static("concreteName", (std::string (*)(const unsigned long &)) &Teuchos::TypeNameTraits<unsigned long>::concreteName, "C++: Teuchos::TypeNameTraits<unsigned long>::concreteName(const unsigned long &) --> std::string", pybind11::arg(""));
	}
	{ // Teuchos::TypeNameTraits file: line:147
		pybind11::class_<Teuchos::TypeNameTraits<float>, std::shared_ptr<Teuchos::TypeNameTraits<float>>> cl(M("Teuchos"), "TypeNameTraits_float_t", "", pybind11::module_local());
		cl.def( pybind11::init( [](){ return new Teuchos::TypeNameTraits<float>(); } ) );
		cl.def_static("name", (std::string (*)()) &Teuchos::TypeNameTraits<float>::name, "C++: Teuchos::TypeNameTraits<float>::name() --> std::string");
		cl.def_static("concreteName", (std::string (*)(const float &)) &Teuchos::TypeNameTraits<float>::concreteName, "C++: Teuchos::TypeNameTraits<float>::concreteName(const float &) --> std::string", pybind11::arg(""));
	}
	{ // Teuchos::TypeNameTraits file: line:147
		pybind11::class_<Teuchos::TypeNameTraits<double>, std::shared_ptr<Teuchos::TypeNameTraits<double>>> cl(M("Teuchos"), "TypeNameTraits_double_t", "", pybind11::module_local());
		cl.def( pybind11::init( [](){ return new Teuchos::TypeNameTraits<double>(); } ) );
		cl.def_static("name", (std::string (*)()) &Teuchos::TypeNameTraits<double>::name, "C++: Teuchos::TypeNameTraits<double>::name() --> std::string");
		cl.def_static("concreteName", (std::string (*)(const double &)) &Teuchos::TypeNameTraits<double>::concreteName, "C++: Teuchos::TypeNameTraits<double>::concreteName(const double &) --> std::string", pybind11::arg(""));
	}
	{ // Teuchos::TypeNameTraits file:Teuchos_TypeNameTraits.hpp line:184
		pybind11::class_<Teuchos::TypeNameTraits<std::string>, std::shared_ptr<Teuchos::TypeNameTraits<std::string>>> cl(M("Teuchos"), "TypeNameTraits_std_string_t", "", pybind11::module_local());
		cl.def( pybind11::init( [](){ return new Teuchos::TypeNameTraits<std::string>(); } ) );
		cl.def_static("name", (std::string (*)()) &Teuchos::TypeNameTraits<std::string >::name, "C++: Teuchos::TypeNameTraits<std::string >::name() --> std::string");
		cl.def_static("concreteName", (std::string (*)(const std::string &)) &Teuchos::TypeNameTraits<std::string >::concreteName, "C++: Teuchos::TypeNameTraits<std::string >::concreteName(const std::string &) --> std::string", pybind11::arg(""));
	}
	{ // Teuchos::TypeNameTraits file:Teuchos_TypeNameTraits.hpp line:193
		pybind11::class_<Teuchos::TypeNameTraits<void *>, std::shared_ptr<Teuchos::TypeNameTraits<void *>>> cl(M("Teuchos"), "TypeNameTraits_void__star__t", "", pybind11::module_local());
		cl.def( pybind11::init( [](){ return new Teuchos::TypeNameTraits<void *>(); } ) );
		cl.def_static("name", (std::string (*)()) &Teuchos::TypeNameTraits<void *>::name, "C++: Teuchos::TypeNameTraits<void *>::name() --> std::string");
		cl.def_static("concreteName", (std::string (*)(const std::string &)) &Teuchos::TypeNameTraits<void *>::concreteName, "C++: Teuchos::TypeNameTraits<void *>::concreteName(const std::string &) --> std::string", pybind11::arg(""));
	}
	{ // Teuchos::TypeNameTraits file:Teuchos_TypeNameTraits.hpp line:206
		pybind11::class_<Teuchos::TypeNameTraits<void>, std::shared_ptr<Teuchos::TypeNameTraits<void>>> cl(M("Teuchos"), "TypeNameTraits_void_t", "", pybind11::module_local());
		cl.def( pybind11::init( [](){ return new Teuchos::TypeNameTraits<void>(); } ) );
		cl.def_static("name", (std::string (*)()) &Teuchos::TypeNameTraits<void>::name, "C++: Teuchos::TypeNameTraits<void>::name() --> std::string");
		cl.def_static("concreteName", (std::string (*)(const std::string &)) &Teuchos::TypeNameTraits<void>::concreteName, "C++: Teuchos::TypeNameTraits<void>::concreteName(const std::string &) --> std::string", pybind11::arg(""));
	}
	// Teuchos::TestForException_incrThrowNumber() file:Teuchos_TestForException.hpp line:61
	M("Teuchos").def("TestForException_incrThrowNumber", (void (*)()) &Teuchos::TestForException_incrThrowNumber, "Increment the throw number.  \n\nC++: Teuchos::TestForException_incrThrowNumber() --> void");

	// Teuchos::TestForException_getThrowNumber() file:Teuchos_TestForException.hpp line:64
	M("Teuchos").def("TestForException_getThrowNumber", (int (*)()) &Teuchos::TestForException_getThrowNumber, "Increment the throw number.  \n\nC++: Teuchos::TestForException_getThrowNumber() --> int");

	// Teuchos::TestForException_break(const std::string &) file:Teuchos_TestForException.hpp line:68
	M("Teuchos").def("TestForException_break", (void (*)(const std::string &)) &Teuchos::TestForException_break, "The only purpose for this function is to set a breakpoint.\n    \n\n\nC++: Teuchos::TestForException_break(const std::string &) --> void", pybind11::arg("msg"));

	// Teuchos::TestForException_setEnableStacktrace(bool) file:Teuchos_TestForException.hpp line:72
	M("Teuchos").def("TestForException_setEnableStacktrace", (void (*)(bool)) &Teuchos::TestForException_setEnableStacktrace, "Set at runtime if stacktracing functionality is enabled when *\n    exceptions are thrown.  \n\n\nC++: Teuchos::TestForException_setEnableStacktrace(bool) --> void", pybind11::arg("enableStrackTrace"));

	// Teuchos::TestForException_getEnableStacktrace() file:Teuchos_TestForException.hpp line:76
	M("Teuchos").def("TestForException_getEnableStacktrace", (bool (*)()) &Teuchos::TestForException_getEnableStacktrace, "Get at runtime if stacktracing functionality is enabled when\n exceptions are thrown. \n\nC++: Teuchos::TestForException_getEnableStacktrace() --> bool");

	// Teuchos::TestForTermination_terminate(const std::string &) file:Teuchos_TestForException.hpp line:79
	M("Teuchos").def("TestForTermination_terminate", (void (*)(const std::string &)) &Teuchos::TestForTermination_terminate, "Prints the message to std::cerr and calls std::terminate. \n\nC++: Teuchos::TestForTermination_terminate(const std::string &) --> void", pybind11::arg("msg"));

}
