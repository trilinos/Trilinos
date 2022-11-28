#include <Teuchos_ArrayRCPDecl.hpp>
#include <Teuchos_ArrayView.hpp>
#include <Teuchos_ArrayViewDecl.hpp>
#include <Teuchos_ENull.hpp>
#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_RCPNode.hpp>
#include <Teuchos_VerbosityLevel.hpp>
#include <Teuchos_as.hpp>
#include <Teuchos_basic_oblackholestream.hpp>
#include <Teuchos_iostream_helpers.hpp>
#include <cwchar>
#include <ios>
#include <istream>
#include <iterator>
#include <locale>
#include <memory>
#include <sstream> // __str__
#include <streambuf>
#include <string>
#include <vector>

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

void bind_Teuchos_iostream_helpers(std::function< pybind11::module &(std::string const &namespace_) > &M)
{
	// Teuchos::enumIstreamExtractor(std::istream &, enum Teuchos::EVerbosityLevel &) file:Teuchos_iostream_helpers.hpp line:58
	M("Teuchos").def("enumIstreamExtractor", (std::istream & (*)(std::istream &, enum Teuchos::EVerbosityLevel &)) &Teuchos::enumIstreamExtractor<Teuchos::EVerbosityLevel>, "C++: Teuchos::enumIstreamExtractor(std::istream &, enum Teuchos::EVerbosityLevel &) --> std::istream &", pybind11::return_value_policy::automatic, pybind11::arg("std_is"), pybind11::arg("enum_value"));

	// Teuchos::EVerbosityLevel file:Teuchos_VerbosityLevel.hpp line:62
	pybind11::enum_<Teuchos::EVerbosityLevel>(M("Teuchos"), "EVerbosityLevel", pybind11::arithmetic(), "Verbosity level.\n\n \n\n ")
		.value("VERB_DEFAULT", Teuchos::VERB_DEFAULT)
		.value("VERB_NONE", Teuchos::VERB_NONE)
		.value("VERB_LOW", Teuchos::VERB_LOW)
		.value("VERB_MEDIUM", Teuchos::VERB_MEDIUM)
		.value("VERB_HIGH", Teuchos::VERB_HIGH)
		.value("VERB_EXTREME", Teuchos::VERB_EXTREME)
		.export_values();

;

	// Teuchos::toString(const enum Teuchos::EVerbosityLevel) file:Teuchos_VerbosityLevel.hpp line:80
	M("Teuchos").def("toString", (std::string (*)(const enum Teuchos::EVerbosityLevel)) &Teuchos::toString, "Return a std::string representation of the verbosity level.\n\n \n\n \n\nC++: Teuchos::toString(const enum Teuchos::EVerbosityLevel) --> std::string", pybind11::arg("verbLevel"));

	// Teuchos::includesVerbLevel(const enum Teuchos::EVerbosityLevel, const enum Teuchos::EVerbosityLevel, const bool) file:Teuchos_VerbosityLevel.hpp line:96
	M("Teuchos").def("includesVerbLevel", [](const enum Teuchos::EVerbosityLevel & a0, const enum Teuchos::EVerbosityLevel & a1) -> bool { return Teuchos::includesVerbLevel(a0, a1); }, "", pybind11::arg("verbLevel"), pybind11::arg("requestedVerbLevel"));
	M("Teuchos").def("includesVerbLevel", (bool (*)(const enum Teuchos::EVerbosityLevel, const enum Teuchos::EVerbosityLevel, const bool)) &Teuchos::includesVerbLevel, "Return true if the verbosity level includes the given level.\n\n \n\n           [in] The verbosity level that is in effect.\n \n\n\n           [in] The verbosity level the client is asking if\n           is included in verbLevel.\n \n\n\n           [in] Set to true if the level in\n           requestedVerbLevel is the default verbosity level.  In\n           this case, if verbLevel==VERB_DEFAULT, then this function\n           will return true.  The default value is false.\n\nC++: Teuchos::includesVerbLevel(const enum Teuchos::EVerbosityLevel, const enum Teuchos::EVerbosityLevel, const bool) --> bool", pybind11::arg("verbLevel"), pybind11::arg("requestedVerbLevel"), pybind11::arg("isDefaultLevel"));

	// Teuchos::incrVerbLevel(const enum Teuchos::EVerbosityLevel, const int) file:Teuchos_VerbosityLevel.hpp line:112
	M("Teuchos").def("incrVerbLevel", (enum Teuchos::EVerbosityLevel (*)(const enum Teuchos::EVerbosityLevel, const int)) &Teuchos::incrVerbLevel, "Return an increased or decreased verbosity level.\n\n \n\n           [in] The base verbosity level.\n \n\n\n           [in] The number of levels to increase (>0) or decrease (<0).\n\n See the function implementation for details on what it does!\n\nC++: Teuchos::incrVerbLevel(const enum Teuchos::EVerbosityLevel, const int) --> enum Teuchos::EVerbosityLevel", pybind11::arg("inputVerbLevel"), pybind11::arg("numLevels"));

	// Teuchos::as(const long &) file:Teuchos_as.hpp line:287
	M("Teuchos").def("as", (int (*)(const long &)) &Teuchos::as<int,long>, "C++: Teuchos::as(const long &) --> int", pybind11::arg("t"));

	// Teuchos::as(const std::string &) file:Teuchos_as.hpp line:287
	M("Teuchos").def("as", (unsigned long (*)(const std::string &)) &Teuchos::as<unsigned long,std::string>, "C++: Teuchos::as(const std::string &) --> unsigned long", pybind11::arg("t"));

	// Teuchos::as(const unsigned long &) file:Teuchos_as.hpp line:287
	M("Teuchos").def("as", (unsigned int (*)(const unsigned long &)) &Teuchos::as<unsigned int,unsigned long>, "C++: Teuchos::as(const unsigned long &) --> unsigned int", pybind11::arg("t"));

	// Teuchos::as(const long &) file:Teuchos_as.hpp line:287
	M("Teuchos").def("as", (short (*)(const long &)) &Teuchos::as<short,long>, "C++: Teuchos::as(const long &) --> short", pybind11::arg("t"));

	// Teuchos::as(const unsigned long &) file:Teuchos_as.hpp line:287
	M("Teuchos").def("as", (unsigned short (*)(const unsigned long &)) &Teuchos::as<unsigned short,unsigned long>, "C++: Teuchos::as(const unsigned long &) --> unsigned short", pybind11::arg("t"));

	// Teuchos::as(const int &) file:Teuchos_as.hpp line:2840
	M("Teuchos").def("as", (long (*)(const int &)) &Teuchos::as<long,int>, "C++: Teuchos::as(const int &) --> long", pybind11::arg("t"));

	// Teuchos::as(const unsigned long &) file:Teuchos_as.hpp line:2840
	M("Teuchos").def("as", (int (*)(const unsigned long &)) &Teuchos::as<int,unsigned long>, "C++: Teuchos::as(const unsigned long &) --> int", pybind11::arg("t"));

	// Teuchos::asSafe(const long &) file:Teuchos_as.hpp line:356
	M("Teuchos").def("asSafe", (int (*)(const long &)) &Teuchos::asSafe<int,long>, "C++: Teuchos::asSafe(const long &) --> int", pybind11::arg("t"));

	// Teuchos::asSafe(const unsigned long &) file:Teuchos_as.hpp line:356
	M("Teuchos").def("asSafe", (unsigned int (*)(const unsigned long &)) &Teuchos::asSafe<unsigned int,unsigned long>, "C++: Teuchos::asSafe(const unsigned long &) --> unsigned int", pybind11::arg("t"));

	// Teuchos::asSafe(const long &) file:Teuchos_as.hpp line:356
	M("Teuchos").def("asSafe", (short (*)(const long &)) &Teuchos::asSafe<short,long>, "C++: Teuchos::asSafe(const long &) --> short", pybind11::arg("t"));

	// Teuchos::asSafe(const unsigned long &) file:Teuchos_as.hpp line:356
	M("Teuchos").def("asSafe", (unsigned short (*)(const unsigned long &)) &Teuchos::asSafe<unsigned short,unsigned long>, "C++: Teuchos::asSafe(const unsigned long &) --> unsigned short", pybind11::arg("t"));

	{ // Teuchos::GlobalMPISession file:Teuchos_GlobalMPISession.hpp line:113
		pybind11::class_<Teuchos::GlobalMPISession, Teuchos::RCP<Teuchos::GlobalMPISession>> cl(M("Teuchos"), "GlobalMPISession", "you would write:\n \n\n\n\n\n\n\n\n\n\n\n This saves you from needing to remember to call MPI_Init() or\n MPI_Finalize().  Also, having the GlobalMPISession object's constructor\n call MPI_Finalize() allows destructors from other objects to call MPI\n functions.  That wold never be possible if you were to directly call\n MPI_Finalize() at the end of main().\n\n This class even works if you have not built Teuchos with MPI support.  In\n that case, it behaves as if MPI_COMM_WORLD had one process, which is\n always the calling process.  Thus, you can use this class to insulate your\n code from needing to know about MPI.  You don't even have to include\n mpi.h, as long as your code doesn't directly use MPI routines or types.\n Teuchos implements wrappers for MPI communicators (see the Teuchos::Comm\n class and its subclasses in the TeuchosComm subpackage) which allow you to\n use a very very small subset of MPI functionality without needing to\n include mpi.h or depend on MPI in any way.\n\n This class also contains the most minimal of other static member functions\n that are needed for only the most simplistic of tasks needed by other\n TeuchosCore software.  For example, you can do a barrier or sum an int\n across processes.  These are needed by the most basic operations involving\n output or determining success or failure across processes for unit tests.\n\n GlobalMPISession's static functions cleverly checks whether MPI has been\n initialized already before calling any MPI functions.  Therefore, you can\n use it in your libraries without requiring that a GlobalMPISession object\n was created in main().");
		cl.def_static("abort", (void (*)()) &Teuchos::GlobalMPISession::abort, "abort the program\n\n Calls MPI_Abort for HAVE_MPI\n Otherwise calls std::abort\n\nC++: Teuchos::GlobalMPISession::abort() --> void");
		cl.def_static("mpiIsInitialized", (bool (*)()) &Teuchos::GlobalMPISession::mpiIsInitialized, "Return whether MPI was initialized.\n\n This is always true if the constructor returned.  If the\n constructor was not called, it may or may not be true, depending\n on whether the user called MPI_Init() themselves.  If the\n constructor was called but threw an exception, then some MPI\n function returned an error code.\n\nC++: Teuchos::GlobalMPISession::mpiIsInitialized() --> bool");
		cl.def_static("mpiIsFinalized", (bool (*)()) &Teuchos::GlobalMPISession::mpiIsFinalized, "Return whether MPI was already finalized.\n\n This is always true if the destructor was called.  If the\n destructor was not called, it may or may not be true, depending\n on whether the user called MPI_Init() themselves.\n\nC++: Teuchos::GlobalMPISession::mpiIsFinalized() --> bool");
		cl.def_static("getRank", (int (*)()) &Teuchos::GlobalMPISession::getRank, "The rank of the calling process in MPI_COMM_WORLD.\n\n \n 0 if MPI has not yet been initialized, else the\n   rank of the calling process in MPI_COMM_WORLD.\n\n You may call this method even if the constructor was never\n called.  Thus, it is safe to use no matter how MPI_Init() was\n called.  However, MPI_Init() must have been called somehow in\n order for this method to return a sensible result.\n\nC++: Teuchos::GlobalMPISession::getRank() --> int");
		cl.def_static("getNProc", (int (*)()) &Teuchos::GlobalMPISession::getNProc, "The number of processes in MPI_COMM_WORLD.\n\n \n 1 if MPI has not yet been initialized, else the\n   number of processes in MPI_COMM_WORLD.\n\n You may call this method even if the constructor was never\n called.  Thus, it is safe to use no matter how MPI_Init() was\n called.  However, MPI_Init() must have been called somehow in\n order for this method to return a sensible result.\n\nC++: Teuchos::GlobalMPISession::getNProc() --> int");
		cl.def_static("barrier", (void (*)()) &Teuchos::GlobalMPISession::barrier, "Call MPI_Barrier() on MPI_COMM_WORLD.\n\n This method must be called collectively on all processes in\n MPI_COMM_WORLD.\n\n \n Users should invoke barrier through the Teuchos::Comm\n   interface.  We only expose this method for Teuchos-internal\n   functionality.\n\nC++: Teuchos::GlobalMPISession::barrier() --> void");
		cl.def_static("sum", (int (*)(int)) &Teuchos::GlobalMPISession::sum, "Sum a set of integers across processes.\n\n This performs an MPI_Allreduce() of localVal over\n MPI_COMM_WORLD, and returns the result (which is the\n same on all processes).\n\n This method must be called collectively on all processes in\n MPI_COMM_WORLD.\n\n \n [in] Value on local process to sum across processes.\n \n\n The global sum (on all processes).\n\n \n Users should invoke reductions through the Teuchos::Comm\n   interface.  We only expose this method for Teuchos-internal\n   functionality.\n\nC++: Teuchos::GlobalMPISession::sum(int) --> int", pybind11::arg("localVal"));
	}
	{ // Teuchos::basic_oblackholestream file:Teuchos_basic_oblackholestream.hpp line:59
		pybind11::class_<Teuchos::basic_oblackholestream<char,std::char_traits<char>>, Teuchos::RCP<Teuchos::basic_oblackholestream<char,std::char_traits<char>>>> cl(M("Teuchos"), "basic_oblackholestream_char_std_char_traits_char_t", "");
		cl.def( pybind11::init( [](){ return new Teuchos::basic_oblackholestream<char,std::char_traits<char>>(); } ) );
	}
}
