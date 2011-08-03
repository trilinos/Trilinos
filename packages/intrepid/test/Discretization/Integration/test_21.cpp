// @HEADER
// ************************************************************************
//
//                           Intrepid Package
//                 Copyright (2007) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Pavel Bochev (pbboche@sandia.gov) or
//                    Denis Ridzal (dridzal@sandia.gov).
//
// ************************************************************************
// @HEADER


/** \file
\brief  Unit test (CubatureDirect): correctness of
        integration of monomials for 1D reference cells.
\author Created by P. Bochev and D. Ridzal.
*/

#include "Intrepid_AdaptiveSparseGrid.hpp"
//#include "Intrepid_CubatureLineSorted.hpp"
#include "Intrepid_Utils.hpp"
#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_RefCountPtr.hpp"
#include "Teuchos_GlobalMPISession.hpp"

using namespace Intrepid;

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
  public Intrepid::AdaptiveSparseGridInterface<Scalar,UserVector> {  
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

void adaptSG(StdVector<long double> & iv,
	     std::multimap<long double,std::vector<int> > & activeIndex,
	     std::set<std::vector<int> > & oldIndex,
	     AdaptiveSparseGridInterface<long double,StdVector<long double> > & problem_data,
	     long double TOL, int flag) {

  // Construct a Container for the adapted rule
  int dimension = problem_data.getDimension();
  std::vector<int> index(dimension,1);
  
  // Initialize global error indicator
  long double eta = 1.0;
  
  // Initialize the Active index set
  activeIndex.insert(std::pair<long double,std::vector<int> >(eta,index));

  std::vector<EIntrepidBurkardt> rule1D; problem_data.getRule(rule1D);
  std::vector<EIntrepidGrowth> growth1D; problem_data.getGrowth(growth1D);
  bool isNormalized = problem_data.isNormalized();
  CubatureTensorSorted<long double> cubRule(
            dimension,index,rule1D,growth1D,isNormalized);
  std::cout << "BEFORE\n";
  // Perform Adaptation
  while (eta > TOL) {
    if (flag==1) {
      eta = AdaptiveSparseGrid<long double,StdVector<long double> >::refine_grid(activeIndex,oldIndex,iv,eta,problem_data);
     }
    else if (flag==2) {
      eta = AdaptiveSparseGrid<long double,StdVector<long double> >::refine_grid(activeIndex,oldIndex,iv,cubRule,eta,problem_data);
     }
  }
  std::cout << "AFTER\n";
}

int main(int argc, char *argv[]) {
  
  Teuchos::GlobalMPISession mpiSession(&argc, &argv);

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
  << "|     1) Integrate a sum of Gaussians in 2D and compare index sets.           |\n" \
  << "|                                                                             |\n" \
  << "|  Questions? Contact  Drew Kouri (dpkouri@sandia.gov) or                     |\n" \
  << "|                      Denis Ridzal (dridzal@sandia.gov).                     |\n" \
  << "|                                                                             |\n" \
  << "|  Intrepid's website: http://trilinos.sandia.gov/packages/intrepid           |\n" \
  << "|  Trilinos website:   http://trilinos.sandia.gov                             |\n" \
  << "|                                                                             |\n" \
  << "===============================================================================\n"\
  << "| TEST 21: Compare index sets for different instances of refine grid          |\n"\
  << "===============================================================================\n";


  // internal variables:
  int         errorFlag    = 0;
  long double TOL          = INTREPID_TOL;
  int         dimension    = 2;
  int         maxLevel     = 25;
  bool        isNormalized = true;                  

  std::vector<EIntrepidBurkardt> rule1D(dimension,BURK_CLENSHAWCURTIS);
  std::vector<EIntrepidGrowth>   growth1D(dimension,GROWTH_FULLEXP);

  ASGdata<long double,StdVector<long double> > problem_data(dimension,rule1D,growth1D,maxLevel,isNormalized);  
  Teuchos::RCP<std::vector<long double> > integralValue = 
    Teuchos::rcp(new std::vector<long double>(1,0.0));
  StdVector<long double> sol(integralValue); sol.Set(0.0);
  problem_data.init(sol);

  try { 
    
    // Initialize the index sets
    std::multimap<long double,std::vector<int> > activeIndex1;
    std::set<std::vector<int> > oldIndex1;  
    
    adaptSG(sol,activeIndex1,oldIndex1,problem_data,TOL,1);
 
    // Initialize the index sets
    std::multimap<long double,std::vector<int> > activeIndex2;
    std::set<std::vector<int> > oldIndex2;  
    adaptSG(sol,activeIndex2,oldIndex2,problem_data,TOL,2);
    if (activeIndex1.size()!=activeIndex2.size()||oldIndex1.size()!=oldIndex2.size()) {
      errorFlag++;
      *outStream << std::right << std::setw(104) << "^^^^---FAILURE!\n";
    }
    else {
      std::multimap<long double,std::vector<int> > inter1;
      std::set_intersection(activeIndex1.begin(),activeIndex1.end(), 
			    activeIndex2.begin(),activeIndex2.end(),
			    inserter(inter1,inter1.begin()),inter1.value_comp());
      std::set<std::vector<int> > inter2;
      std::set_intersection(oldIndex1.begin(),oldIndex1.end(), 
			    oldIndex2.begin(),oldIndex2.end(),
			    inserter(inter2,inter2.begin()),inter2.value_comp());
      
      if(inter1.size()!=activeIndex1.size()||inter1.size()!=activeIndex2.size()||
	 inter2.size()!=oldIndex1.size()||inter2.size()!=oldIndex2.size()) {
	errorFlag++;
	*outStream << std::right << std::setw(104) << "^^^^---FAILURE!\n";
      }
      else {
	*outStream << "Index sets 1 and " << 2 << " are equal!\n";
      } 
    }
  }
  catch (std::logic_error err) {    
    *outStream << err.what() << "\n";
    errorFlag = -1;
  };  

  if (errorFlag != 0)
    std::cout << "End Result: TEST FAILED\n";
  else
    std::cout << "End Result: TEST PASSED\n";

  // reset format state of std::cout
  std::cout.copyfmt(oldFormatState);

  return errorFlag;
}
