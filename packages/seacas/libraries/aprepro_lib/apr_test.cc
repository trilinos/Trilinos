#include <iostream>
#include <fstream>

#include "aprepro.h"

// This function is used below in the example showing how an
// application can add its own functions to an aprepro instance.
double succ(double i) {
  return ++i;
}

int main(int argc, char *argv[])
{
  bool readfile = false;

  SEAMS::Aprepro aprepro;
  
  // EXAMPLE: Add a function to aprepro...
  SEAMS::symrec *ptr = aprepro.putsym("succ", SEAMS::Aprepro::FUNCTION, 0);
  ptr->value.fnctptr_d = succ;
  ptr->info = "Return the successor to d";
  ptr->syntax = "succ(d)";
  
  for(int ai = 1; ai < argc; ++ai) {
    if (argv[ai] == std::string ("-s")) {
      aprepro.statistics();
    }
    else if (argv[ai] == std::string ("-C")) {
      aprepro.copyright();
    }
    else if (argv[ai] == std::string ("-p")) {
      aprepro.ap_options.trace_parsing = true;
    }
    else if (argv[ai] == std::string ("-d")) {
      aprepro.ap_options.debugging = true;
    }
    else if (argv[ai] == std::string ("-e")) {
      aprepro.ap_options.end_on_exit = true;
    }
    else if (argv[ai] == std::string ("-M")) {
      aprepro.ap_options.info_msg = true;
    }
    else if (argv[ai] == std::string ("-X")) {
      aprepro.ap_options.immutable = true;
      aprepro.stateImmutable = true;
    }
    else if (argv[ai] == std::string ("-W")) {
      aprepro.ap_options.warning_msg = false;
    }
    else if (argv[ai] == std::string ("-v")) {
      std::cerr << "Aprepro version " << aprepro.version() << "\n";
    }
    else if (argv[ai] == std::string ("-i")) {
      // Read from cin and echo each line to cout All results will
      // also be stored in Aprepro::parsing_results() stream if needed
      // at end of file.
      aprepro.ap_options.interactive = true;
      bool result = aprepro.parse_stream(std::cin, "standard input");
      if (result) {
	std::cout << "PARSING RESULTS: " << aprepro.parsing_results().str();
      }
    }
    else {
      // Read and parse a file.  The entire file will be parsed and
      // then the output can be obtained in an std::ostringstream via
      // Aprepro::parsing_results()
      std::fstream infile(argv[ai]);
      if (!infile.good()) {
	std::cerr << "APREPRO: Could not open file: " << argv[ai] << std::endl;
	return 0;
      }

      bool result = aprepro.parse_stream(infile, argv[ai]);
      if (result) {
	std::cout << "PARSING RESULTS: " << aprepro.parsing_results().str();
      }

      readfile = true;
    }
  }

  if (readfile) return 0;
    
  // Read and parse a string's worth of data at a time.
  // Cannot use looping/ifs/... with this method.
  std::string line;
  while( std::cout << "\nexpession: " &&
	 std::getline(std::cin, line) &&
	 !line.empty() ) {
    if (line[0] != '{') 
      line = "{" + line + "}\n";
    else
      line += "\n";
    
    bool result = aprepro.parse_string(line, "input");

    if (result) {
      std::cout << "         : " << aprepro.parsing_results().str();
      aprepro.clear_results();
    }
  }
}
