#include "ROL_MultiVectorDefault.hpp"
#include "ROL_StdVector.hpp"
#include "ROL_Types.hpp"

#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_GlobalMPISession.hpp"

#include <iomanip>
#include <iostream>


typedef double RealT;

using namespace ROL;
using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::ArrayRCP;


template<class Real>
Real norm_sum(const MultiVector<Real> &A) {
    int numVectors = A.getNumberOfVectors();
    std::vector<RealT> norms(numVectors);
    A.norms(norms);
    Real sum = 0;
    for(int i=0;i<numVectors;++i) {
        sum += norms[i]; 
    }
    return sum;
}


int main(int argc, char *argv[]) {

    Teuchos::GlobalMPISession mpiSession(&argc,&argv);

    int iprint     = argc - 1;
    Teuchos::RCP<std::ostream> outStream;
    Teuchos::oblackholestream bhs; // outputs nothing
    if (iprint > 0)
        outStream = Teuchos::rcp(&std::cout, false);
    else
        outStream = Teuchos::rcp(&bhs, false);

    int errorFlag = 0;

    try { 

	Teuchos::SerialDenseMatrix<int,RealT> M(2,2,true);
	M(0,0) = 2.0;
	M(0,1) = -1.0;
	M(1,0) = -2.0;
	M(1,1) = 3.0;


	Teuchos::RCP<std::vector<RealT> > w_rcp = Teuchos::rcp(new std::vector<RealT>(2));
	Teuchos::RCP<std::vector<RealT> > x_rcp = Teuchos::rcp(new std::vector<RealT>(2));
	Teuchos::RCP<std::vector<RealT> > y_rcp = Teuchos::rcp(new std::vector<RealT>(2));
	Teuchos::RCP<std::vector<RealT> > z_rcp = Teuchos::rcp(new std::vector<RealT>(2));

	(*w_rcp)[0] = 0.0;
	(*w_rcp)[1] = 1.0;

	(*x_rcp)[0] = 1.0;
	(*x_rcp)[1] = 0.0;

	(*y_rcp)[0] = 3.0;
	(*y_rcp)[1] = 4.0;

	(*z_rcp)[0] = -1.0;
	(*z_rcp)[1] =  1.0;
       
	RCP<Vector<RealT> > w = rcp(new StdVector<RealT>(w_rcp)); 
	RCP<Vector<RealT> > x = rcp(new StdVector<RealT>(x_rcp)); 
	RCP<Vector<RealT> > y = rcp(new StdVector<RealT>(y_rcp)); 
	RCP<Vector<RealT> > z = rcp(new StdVector<RealT>(z_rcp)); 

	ArrayRCP<RCP<Vector<RealT>>> A_rcp(2);
	ArrayRCP<RCP<Vector<RealT>>> B_rcp(2);

	A_rcp[0] = x;     
	A_rcp[1] = y;     

	B_rcp[0] = w;     
	B_rcp[1] = z;     

	RCP<MultiVector<RealT>> A = rcp(new MultiVectorDefault<RealT>(A_rcp));
	RCP<MultiVector<RealT>> B = rcp(new MultiVectorDefault<RealT>(B_rcp));
       
	// Test norm
	if(static_cast<int>(norm_sum(*A)) != 6) {
	    ++errorFlag;
	}

	// Test clone
	RCP<MultiVector<RealT>> C = A->clone();    
	if(norm_sum(*C) != 0) {
	    ++errorFlag;
	}

	// Test deep copy
        RCP<MultiVector<RealT>> D = A->deepCopy();
	if(static_cast<int>(norm_sum(*D)) != 6) {
	    ++errorFlag;
	}
        
        // Test shallow copy
	std::vector<int> index(1);
	index[0] = 0;

        RCP<MultiVector<RealT>> S = A->shallowCopy(index);
	if(static_cast<int>(norm_sum(*S)) != 1) {
	    ++errorFlag;
	}

        // Test scaling
	std::vector<RealT> alpha(2);
	alpha[0] = 4.0;
	alpha[1] = 9.0;
	A->scale(alpha);
	if(static_cast<int>(norm_sum(*A)) != 49) {
	    ++errorFlag;
	}

	// Test matrix multiplication 
	A->gemm(2.0,*B,M,1.0);
	if(static_cast<int>(norm_sum(*A)) != 49) {
	    ++errorFlag;
	}

	// Test set
        A->set(*B);
	if(static_cast<int>(norm_sum(*A)) != 2) {
	    ++errorFlag;
	}

	// Test innerProducts
        Teuchos::SerialDenseMatrix<int,RealT> P(2,2,true);
        B->innerProducts(1.0,*B,P);
        Teuchos::SerialDenseMatrix<int,RealT> Check(2,2,true);
        Check(0,0) = 1.0;   
        Check(0,1) = 1.0;   
        Check(1,0) = 1.0;   
        Check(1,1) = 2.0;   
        if( P != Check ) {
	    ++errorFlag;
        }

	// Test dot products
        std::vector<RealT> dots(2); 
        D->dots(*D,dots);
        if(static_cast<int>(dots[0]) != 1 || 
           static_cast<int>(dots[1]) != 25 ) {
            ++errorFlag;
        }

//    StdVector<RealT> d0 = Teuchos::dyn_cast<StdVector<RealT>>(*D->getVector(0));
//    StdVector<RealT> d1 = Teuchos::dyn_cast<StdVector<RealT>>(*D->getVector(1));

//    std::cout << (*d0.getVector())[0] << std::endl;
//    std::cout << (*d0.getVector())[1] << std::endl;
//    std::cout << (*d1.getVector())[0] << std::endl;
//    std::cout << (*d1.getVector())[1] << std::endl;

    }
    catch(std::logic_error err) {
        *outStream << err.what() << "\n";
        errorFlag = -1000;
    }; // end try
    if (errorFlag != 0)
        std::cout << "End Result: TEST FAILED\n";
    else
        std::cout << "End Result: TEST PASSED\n";

    return 0;
}

