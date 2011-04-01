//
//  Suppose we wanted to use Teuchos::ParameterList for the "DEBUG_OUTPUT_STREAM"
//  parameter.  How would we do it?  This test shows how.  The user has to give
//  us a pointer to the stream instead of a reference to the stream but everything
//  else seems to work.
//  
//  We assume the stream is a file, std::cout, std::cerr, std::clog or an ostringstream.
//
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_ParameterEntry.hpp>
#include <sstream>
#include <fstream>
#include <ostream>

using namespace std;

#if 0
// Copy these classes into Teuchos_TypeNameTraits.hpp.
//    We can't template on the ostream, ofstream or ostringstream because 
//    it requires making a copy in the ParameterList.  Assignment is
//    private in these objects.  So we use pointers to the output streams.
namespace Teuchos{
template<> 
class TEUCHOS_LIB_DLL_EXPORT TypeNameTraits<std::ostream*> { 
public: 
  static std::string name() { return ("ostream*"); } 
  static std::string concreteName(const std::ostream *&) { return name(); } 
};

template<> 
class TEUCHOS_LIB_DLL_EXPORT TypeNameTraits<std::ofstream*> { 
public: 
  static std::string name() { return ("ofstream*"); } 
  static std::string concreteName(const std::ofstream *&) { return name(); } 
};

template<> 
class TEUCHOS_LIB_DLL_EXPORT TypeNameTraits<std::ostringstream *> { 
public: 
  static std::string name() { return ("ostringstream*"); } 
  static std::string concreteName(const std::ostringstream *&) { return name(); } 
}; 
}
#endif

void check_set_debug_output_stream(Teuchos::ParameterList &pl, std::ostream **os)
{
  std::ostream **ostreamPtr=NULL;
  std::ofstream **ofstreamPtr=NULL;
  std::ostringstream **ostringstreamPtr=NULL;

  Teuchos::ParameterEntry &entry = pl.getEntry("DEBUG_OSTREAM");

  if (entry.isType<ostream *>())
    *os = entry.getValue(ostreamPtr);

  else if (entry.isType<ofstream *>())
    *os = entry.getValue(ofstreamPtr);

  else if (entry.isType<ostringstream *>())
    *os = entry.getValue(ostringstreamPtr);

  else{
    std::cout << "Invalid debug output stream" << std::endl;
    return;
  }
}

int main(int argc, char *argv[])
{
  std::ostream *debugStream;
  std::ofstream f;
  std::ostringstream oss;

  Teuchos::ParameterList pl("testList");

  // Set debug output to a file.

  f.open("test.txt");
  pl.set("DEBUG_OSTREAM", &f, "Output stream for debugging statements");
  check_set_debug_output_stream(pl, &debugStream);

  *debugStream << "TESTING OUTPUT TO FILE\n";
  f.close();
  std::cout << "Check test.txt for TESTING OUTPUT TO FILE\n";

  // Set debug output to ostringstream.

  pl.set("DEBUG_OSTREAM", &oss, "Output stream for debugging statements");
  check_set_debug_output_stream(pl, &debugStream);

  *debugStream << "TESTING OUTPUT TO STRINGSTREAM\n";
  std::cout << oss.str();

  // Set debug output to std::cout

  pl.set("DEBUG_OSTREAM", &std::cout, "Output stream for debugging statements");
  check_set_debug_output_stream(pl, &debugStream);

  *debugStream << "TESTING OUTPUT TO std::cout\n";

  // Set debug output to std::cerr

  pl.set("DEBUG_OSTREAM", &std::cerr, "Output stream for debugging statements");
  check_set_debug_output_stream(pl, &debugStream);

  *debugStream << "TESTING OUTPUT TO std::cerr\n";

  // Set debug output to std::clog

  pl.set("DEBUG_OSTREAM", &std::clog, "Output stream for debugging statements");
  check_set_debug_output_stream(pl, &debugStream);

  *debugStream << "TESTING OUTPUT TO std::clog\n";
}

