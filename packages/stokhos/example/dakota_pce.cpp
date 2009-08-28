#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>

#include "Stokhos.hpp"

int main(int argc, char *argv[]) {
  Teuchos::Array<double> x;

  // Get filenames for Dakota runs
  if (argc != 3) {
    std::cout << "Usage:  dakota_pce.exe [input_file] [output file]" 
	      << std::endl;
    exit(-1);
  }

  std::string input_filename = std::string(argv[1]);
  std::string output_filename = std::string(argv[2]);
  std::ifstream input_file(input_filename.c_str());
  int nvar;
  std::string name;
  input_file >> nvar >> name;
  x.resize(nvar);
  for (int i=0; i<nvar; i++) {
    input_file >> x[i] >> name;
  }
  input_file.close();

  typedef Stokhos::HermiteBasis<int,double> basis_type;
  int p = 4;
  Teuchos::Array< Teuchos::RCP<const Stokhos::OneDOrthogPolyBasis<int,double> > > bases(nvar); 
  for (int i=0; i<nvar; i++) {
    bases[i] = Teuchos::rcp(new basis_type(p));
  }
  Teuchos::RCP<const Stokhos::OrthogPolyBasis<int,double> > basis = 
    Teuchos::rcp(new Stokhos::CompletePolynomialBasis<int,double>(bases));
  Stokhos::OrthogPolyApprox<int,double> u(basis);

  u[0] = 1.0;
  u[1] = 0.4;
  u[2] = 0.06;
  u[3] = 0.002;

  double uu = u.evaluate(x);
  double v = std::log(uu);
  v = 1.0 / (v*v + 1.0);

  std::ofstream output_file(output_filename.c_str());
  output_file.precision(16);
  output_file.setf(std::ios::scientific);
  output_file << v << " " << v << std::endl;
  output_file.close();

  return 0;
}
