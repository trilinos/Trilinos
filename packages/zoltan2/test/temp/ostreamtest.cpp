//  Suppose we don't use ParameterList for output streams.
//
//  We want caller to give us an ostream for debug output and another for
//  profiling/timing output.  
//
//  We assume the stream is a file, std::cout, std::cerr, std::clog or an ostringstream.
//
#include <Teuchos_FancyOStream.hpp>
#include <Teuchos_RCP.hpp>
#include <sstream>
#include <fstream>
#include <ostream>

using namespace std;

// Our ProblemParameters object would have the information that's in the
// class test_parameters defined below.
//
// Objects in Zoltan2 that want to stream debug or timing data would get the 
// Fancy streamer from the ProblemParameters object. 
//
// An object that wants to indent further than the current indent would:
//     Problem::getDebugFancyOStream().getTabIndentStr()
//     Problem::getDebugFancyOStream().setTabIndentStr() to a bigger indent
//
// The caller's interface would be:
//
// Problem::setDebugOutputStream(std::ostream &os);
// Problem::setTimingOutputStream(std::ostream &os);
//
// instead of setting it in a Teuchos::ParameterList.


class test_parameters{

private:

  Teuchos::RCP<Teuchos::basic_FancyOStream<char> > debugOStream;
  Teuchos::RCP<Teuchos::basic_FancyOStream<char> > timingOStream;

public:

  ~test_parameters(){
    debugOStream = Teuchos::null;
    timingOStream = Teuchos::null;
  }

  void setDebugOutputStream(std::basic_ostream<char> &os){
    Teuchos::RCP<std::basic_ostream<char> > osptr(&os, false);
    debugOStream = Teuchos::null;
    debugOStream = Teuchos::fancyOStream(osptr);
  }

  void setTimingOutputStream(std::ostream &os){
    Teuchos::RCP<std::basic_ostream<char> > osptr(&os, false);
    timingOStream = Teuchos::null;
    timingOStream = Teuchos::fancyOStream(osptr);  
  }
  
  Teuchos::RCP<Teuchos::basic_FancyOStream<char> > &getDebugFancyOStream() {return debugOStream;}
  
  Teuchos::RCP<Teuchos::basic_FancyOStream<char> > &getTimingFancyOStream() {return timingOStream;}
  
  Teuchos::RCP<std::basic_ostream<char> > getDebugOStream() {return debugOStream->getOStream();}
  
  Teuchos::RCP<std::basic_ostream<char> > getTimingOStream() {return timingOStream->getOStream();}

};

int main(int argc, char *argv[])
{
  std::ofstream f;
  std::ostringstream oss;
  Teuchos::RCP<Teuchos::basic_FancyOStream<char> > os;

  test_parameters TestParameters;

  // Set debug output to a file.

  f.open("test.txt");

  TestParameters.setDebugOutputStream(f);

  os = TestParameters.getDebugFancyOStream();
  os->setShowProcRank(true);
  os->setTabIndentStr(std::string("---->"));

  *os << "TESTING OUTPUT TO FILE\n";

  // TODO: set some formatting parameters on the os and then print mroe

  f.close();
  std::cout << "Check test.txt for TESTING OUTPUT TO FILE\n";

  // Set debug output to ostringstream.

  TestParameters.setDebugOutputStream(oss);
  os = TestParameters.getDebugFancyOStream();
  os->setShowLinePrefix(true);

  *os << "TESTING OUTPUT TO STRINGSTREAM\n";
  std::cout << oss.str() << std::endl;

  // Set debug output to std::cout

  TestParameters.setDebugOutputStream(std::cout);
  os = TestParameters.getDebugFancyOStream();

  *os << "TESTING OUTPUT TO std::cout\n" << std::endl;

  os->pushTab();
  *os << "TESTING OUTPUT TO std::cout\n" << std::endl;
  os->pushTab();
  *os << "TESTING OUTPUT TO std::cout\n" << std::endl;
  os->popTab();
  *os << "TESTING OUTPUT TO std::cout\n" << std::endl;
  os->popTab();
  *os << "TESTING OUTPUT TO std::cout\n" << std::endl;

  // Set debug output to std::cerr

  TestParameters.setDebugOutputStream(std::cerr);
  os = TestParameters.getDebugFancyOStream();
  os->setShowAllFrontMatter(true);

  *os << "TESTING OUTPUT TO std::cerr\n" << std::endl;

  // Set debug output to std::clog

  TestParameters.setDebugOutputStream(std::clog);
  os = TestParameters.getDebugFancyOStream();

  *os<< "TESTING OUTPUT TO std::clog\n" << std::endl;
}

