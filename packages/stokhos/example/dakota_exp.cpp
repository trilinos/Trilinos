#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <cstdlib>

// A simle executable that computes the exponential, reading the input
// value from the input file on the command line, and putting the result
// in the output file on the command line.
//
// The purpose of this is to run Dakota around this executable to experiment
// with polynomial chaos expansions.
int main(int argc, char *argv[]) {

  // Get filenames for Dakota runs
  if (argc != 3) {
    std::cout << "Usage:  dakota_exp.exe [input_file] [output file]" 
	      << std::endl;
    std::exit(-1);
  }

  std::string input_filename = std::string(argv[1]);
  std::string output_filename = std::string(argv[2]);
  std::ifstream input_file(input_filename.c_str());
  int nvar;
  double val;
  std::string name;
  input_file >> nvar >> name;
  double x = 0.0;
  for (int i=0; i<nvar; i++) {
    input_file >> val >> name;
    x += val;
  }
  input_file.close();

  double v = std::exp(x);

  std::ofstream output_file(output_filename.c_str());
  output_file.precision(12);
  output_file.setf(std::ios::scientific);
  output_file << v << " " << v << std::endl;
  output_file.close();

  return 0;
}
