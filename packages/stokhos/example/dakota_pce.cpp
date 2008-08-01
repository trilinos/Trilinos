#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>

#include "Stokhos_OrthogPolyExpansion.hpp"
#include "Stokhos_HermiteEBasis2.hpp"
#include "Stokhos_CompletePolynomialBasis.hpp"

int main(int argc, char *argv[]) {
  std::vector<double> x;

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

  typedef Stokhos::HermiteEBasis2<double> basis_type;
  int p = 1;
  std::vector< Teuchos::RCP<const Stokhos::OrthogPolyBasis<double> > > bases(nvar); 
  std::vector<double> deriv_coeffs(nvar);
  for (int i=0; i<nvar; i++) {
    bases[i] = Teuchos::rcp(new basis_type(p));
    deriv_coeffs[i] = 1.0;
  }
  Teuchos::RCP< Stokhos::CompletePolynomialBasis<double> > basis = 
    Teuchos::rcp(new Stokhos::CompletePolynomialBasis<double>(bases,
							      deriv_coeffs));
  Stokhos::OrthogPolyExpansion<double> he(basis);
  unsigned int sz = basis->size();
  Stokhos::OrthogPolyApprox<double> u(sz);
  u[0] = 0.5;
  if (sz >= 1)
    u[1] = 0.05;
  if (sz >= 2)
    u[2] = 0.05;

  double uu = u.evaluate(*basis,x);
  double v = std::exp(std::exp(uu));

  std::ofstream output_file(output_filename.c_str());
  output_file.precision(15);
  output_file.setf(std::ios::scientific);
  output_file << v << " " << v << std::endl;
  output_file.close();

  return 0;
}
