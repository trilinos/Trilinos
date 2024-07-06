// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER


/** \file
\brief  Unit test (CubatureDirect): correctness of
        integration of monomials for 1D reference cells.
\author Created by P. Bochev and D. Ridzal.
*/

#include "Intrepid2_AdaptiveSparseGrid.hpp"
//#include "Intrepid2_CubatureLineSorted.hpp"
#include "Intrepid2_Utils.hpp"
#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_RefCountPtr.hpp"
#include "Teuchos_GlobalMPISession.hpp"

using namespace Intrepid2;

template<class Scalar>
class StdVector {
private:
  Teuchos::RefCountPtr<std::vector<Scalar> >  std_vec_;

public:

  StdVector( const Teuchos::RefCountPtr<std::vector<Scalar> > & std_vec ) 
  : std_vec_(std_vec) {}

  Teuchos::RefCountPtr<StdVector<Scalar> > Create() const {
    return Teuchos::rcp( new StdVector<Scalar>(
	   Teuchos::rcp(new std::vector<Scalar>(std_vec_->size(),0))));
  }

  void Update( StdVector<Scalar> & s ) {
    int dimension  = (int)(std_vec_->size());
    for (int i=0; i<dimension; i++)
      (*std_vec_)[i] += s[i]; 
  }

  void Update( Scalar alpha, StdVector<Scalar> & s )  {
    int dimension  = (int)(std_vec_->size());
    for (int i=0; i<dimension; i++)
      (*std_vec_)[i] += alpha*s[i];
  }

  Scalar operator[](int i) {
    return (*std_vec_)[i];
  }

  void clear() {
    std_vec_->clear();
  }

  void resize(int n, Scalar p) {
    std_vec_->resize(n,p);
  }

  int size() {
    return (int)std_vec_->size();
  }

  void Set( Scalar alpha )  {
    int dimension  = (int)(std_vec_->size());
    for (int i=0; i<dimension; i++)
      (*std_vec_)[i] = alpha;
  }
};

template<class Scalar, class UserVector>
class ASGdata : 
  public Intrepid2::AdaptiveSparseGridInterface<Scalar,UserVector> {  
public:  
  ~ASGdata() {}

  ASGdata(int dimension,std::vector<EIntrepidBurkardt> rule1D,
	  std::vector<EIntrepidGrowth> growth1D, int maxLevel,
	  bool isNormalized) : AdaptiveSparseGridInterface<Scalar,UserVector>(
	  dimension,rule1D,growth1D,maxLevel,isNormalized) {}

  void eval_integrand(UserVector & output, std::vector<Scalar> & input) {
    output.clear(); output.resize(1,std::exp(-input[0]*input[0])
				  +10.0*std::exp(-input[1]*input[1]));
  }  
  Scalar error_indicator(UserVector & input) {
    int dimension = (int)input.size();
    Scalar norm2  = 0.0;
    for (int i=0; i<dimension; i++)
      norm2 += input[i]*input[i];
    
    Scalar ID = AdaptiveSparseGridInterface<Scalar,UserVector>::
      getInitialDiff();
    norm2 = std::sqrt(norm2)/ID;
    return norm2;
  }
};

long double adaptSG(StdVector<long double> & iv,
	  AdaptiveSparseGridInterface<long double,StdVector<long double> > & 
          problem_data,long double TOL) {

  // Construct a Container for the adapted rule
  int dimension = problem_data.getDimension();
  std::vector<int> index(dimension,1);
  
  // Initialize global error indicator
  long double eta = 1.0;
  
  // Initialize the Active index set
  std::multimap<long double,std::vector<int> > activeIndex;  
  activeIndex.insert(std::pair<long double,std::vector<int> >(eta,index));

  // Initialize the old index set
  std::set<std::vector<int> > oldIndex;

  // Perform Adaptation
  while (eta > TOL) {
    eta = AdaptiveSparseGrid<long double,StdVector<long double> >::refine_grid(
        activeIndex,oldIndex,iv,eta,problem_data);
  }
  return eta;
}

int main(int argc, char *argv[]) {
  
  Teuchos::GlobalMPISession mpiSession(&argc, &argv);
Kokkos::initialize();
  // This little trick lets us print to std::cout only if
  // a (dummy) command-line argument is provided.
  int iprint     = argc - 1;
  Teuchos::RCP<std::ostream> outStream;
  Teuchos::oblackholestream bhs; // outputs nothing
  if (iprint > 0)
    outStream = Teuchos::rcp(&std::cout, false);
  else
    outStream = Teuchos::rcp(&bhs, false);

  // Save the format state of the original std::cout.
  Teuchos::oblackholestream oldFormatState;
  oldFormatState.copyfmt(std::cout);
 
  *outStream \
  << "===============================================================================\n" \
  << "|                                                                             |\n" \
  << "|                         Unit Test (AdaptiveSparseGrid)                      |\n" \
  << "|                                                                             |\n" \
  << "|     1) Integrate a sum of Gaussians in 2D (Gerstner and Griebel).           |\n" \
  << "|                                                                             |\n" \
  << "|  Questions? Contact  Drew Kouri (dpkouri@sandia.gov) or                     |\n" \
  << "|                      Denis Ridzal (dridzal@sandia.gov).                     |\n" \
  << "|                                                                             |\n" \
  << "|  Intrepid's website: http://trilinos.sandia.gov/packages/intrepid           |\n" \
  << "|  Trilinos website:   http://trilinos.sandia.gov                             |\n" \
  << "|                                                                             |\n" \
  << "===============================================================================\n"\
  << "| TEST 20: Integrate an anisotropic sum of Gaussians in 2D                    |\n"\
  << "===============================================================================\n";


  // internal variables:
  int         errorFlag    = 0;
  long double TOL          = INTREPID_TOL;
  int         dimension    = 2;
  int         maxLevel     = 25;
  bool        isNormalized = true;                  

  std::vector<EIntrepidBurkardt> rule1D(dimension,BURK_CLENSHAWCURTIS);
  std::vector<EIntrepidGrowth>   growth1D(dimension,GROWTH_FULLEXP);
 
  ASGdata<long double,StdVector<long double> > problem_data(
    dimension,rule1D,growth1D,maxLevel,isNormalized);  
  Teuchos::RCP<std::vector<long double> > integralValue = 
    Teuchos::rcp(new std::vector<long double>(1,0.0));
  StdVector<long double> sol(integralValue); sol.Set(0.0);
  problem_data.init(sol);

  long double eta = adaptSG(sol,problem_data,TOL); 

  long double analyticInt = (1.0+10.0)*std::sqrt(M_PI)/2.0*erff(1.0);
  long double abstol      = 1.0e1*std::sqrt(INTREPID_TOL);
  long double absdiff     = fabs(analyticInt-sol[0]);
  try { 
    *outStream << "Adaptive Sparse Grid exited with global error " 
	       << std::scientific << std::setprecision(16) << eta << "\n"
	       << "Approx = " << std::scientific << std::setprecision(16) << sol[0] 
	       << ",  Exact = " << std::scientific << std::setprecision(16) << analyticInt << "\n"
	       << "Error = " << std::scientific << std::setprecision(16) << absdiff << "   " 
	       << "<?" << "   " << abstol << "\n"; 
    if (absdiff > abstol) {
      errorFlag++;
      *outStream << std::right << std::setw(104) << "^^^^---FAILURE!\n";
    }
  }
  catch (std::logic_error &err) {    
    *outStream << err.what() << "\n";
    errorFlag = -1;
  };  

  if (errorFlag != 0)
    std::cout << "End Result: TEST FAILED\n";
  else
    std::cout << "End Result: TEST PASSED\n";

  // reset format state of std::cout
  std::cout.copyfmt(oldFormatState);
Kokkos::finalize();
  return errorFlag;
}
