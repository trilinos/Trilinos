// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file  test_03.cpp
    \brief Test MultiVector interface.

*/

#include "ROL_MultiVectorDefault.hpp"
#include "ROL_StdVector.hpp"
#include "ROL_Types.hpp"

#include "ROL_Stream.hpp"
#include "Teuchos_GlobalMPISession.hpp"

#include <iomanip>
#include <iostream>


typedef double RealT;

using namespace ROL;


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
    ROL::Ptr<std::ostream> outStream;
    ROL::nullstream bhs; // outputs nothing
    if (iprint > 0)
        outStream = ROL::makePtrFromRef(std::cout);
    else
        outStream = ROL::makePtrFromRef(bhs);

    int errorFlag = 0;

    try { 

	Teuchos::SerialDenseMatrix<int,RealT> M(2,2,true);
	M(0,0) = 2.0;
	M(0,1) = -1.0;
	M(1,0) = -2.0;
	M(1,1) = 3.0;


	ROL::Ptr<std::vector<RealT> > w_ptr = ROL::makePtr<std::vector<RealT>>(2);
	ROL::Ptr<std::vector<RealT> > x_ptr = ROL::makePtr<std::vector<RealT>>(2);
	ROL::Ptr<std::vector<RealT> > y_ptr = ROL::makePtr<std::vector<RealT>>(2);
	ROL::Ptr<std::vector<RealT> > z_ptr = ROL::makePtr<std::vector<RealT>>(2);

	(*w_ptr)[0] = 0.0;
	(*w_ptr)[1] = 1.0;

	(*x_ptr)[0] = 1.0;
	(*x_ptr)[1] = 0.0;

	(*y_ptr)[0] = 3.0;
	(*y_ptr)[1] = 4.0;

	(*z_ptr)[0] = -1.0;
	(*z_ptr)[1] =  1.0;
       
	ROL::Ptr<Vector<RealT> > w = ROL::makePtr<StdVector<RealT>>(w_ptr); 
	ROL::Ptr<Vector<RealT> > x = ROL::makePtr<StdVector<RealT>>(x_ptr); 
	ROL::Ptr<Vector<RealT> > y = ROL::makePtr<StdVector<RealT>>(y_ptr); 
	ROL::Ptr<Vector<RealT> > z = ROL::makePtr<StdVector<RealT>>(z_ptr); 

	ArrayRCP<ROL::Ptr<Vector<RealT> > > A_ptr(2);
	ArrayRCP<ROL::Ptr<Vector<RealT> > > B_ptr(2);

	A_ptr[0] = x;     
	A_ptr[1] = y;     

	B_ptr[0] = w;     
	B_ptr[1] = z;     

	ROL::Ptr<MultiVector<RealT> > A = ROL::makePtr<MultiVectorDefault<RealT>>(A_ptr);
        ROL::Ptr<MultiVector<RealT> > B = ROL::makePtr<MultiVectorDefault<RealT>>(B_ptr);
       
	// Test norm
	if(static_cast<int>(norm_sum(*A)) != 6) {
            *outStream << "Norm test failed!\n";
	    ++errorFlag;
	}

	// Test clone
	ROL::Ptr<MultiVector<RealT> > C = A->clone();    
	if(norm_sum(*C) != 0) {
            *outStream << "Clone test failed!\n";
	    ++errorFlag;
	}

	// Test deep copy
        ROL::Ptr<MultiVector<RealT> > D = A->deepCopy();
	if(static_cast<int>(norm_sum(*D)) != 6) {
            *outStream << "Deep copy test failed!\n";
	    ++errorFlag;
	}
        
        // Test shallow copy
	std::vector<int> index(1);
	index[0] = 0;

        ROL::Ptr<MultiVector<RealT> > S = A->shallowCopy(index);
	if(static_cast<int>(norm_sum(*S)) != 1) {
            *outStream << "Shallow copy test failed!\n";
	    ++errorFlag;
	}

        // Test scaling
	std::vector<RealT> alpha(2);
	alpha[0] = 4.0;
	alpha[1] = 9.0;
	A->scale(alpha);
	if(static_cast<int>(norm_sum(*A)) != 49) {
            *outStream << "Scaling test failed!\n";
	    ++errorFlag;
	}

	// Test matrix multiplication 
	A->gemm(2.0,*B,M,1.0);
	if(static_cast<int>(norm_sum(*A)) != 53) {
            *outStream << "Matmat multiply test failed!  The norm_sum is " << static_cast<int>(norm_sum(*A)) << ", not equal to 49.\n";
	    ++errorFlag;
	}

	// Test set
        A->set(*B);
	if(static_cast<int>(norm_sum(*A)) != 2) {
            *outStream << "Set test failed!\n";
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
            *outStream << "Inner product test failed!\n";
	    ++errorFlag;
        }

	// Test dot products
        std::vector<RealT> dots(2); 
        D->dots(*D,dots);
        if(static_cast<int>(dots[0]) != 1 || 
           static_cast<int>(dots[1]) != 25 ) {
            *outStream << "Dot product test failed!\n";
            ++errorFlag;
        }

//    StdVector<RealT> d0 = dynamic_cast<StdVector<RealT>>(*D-&>getVector(0));
//    StdVector<RealT> d1 = dynamic_cast<StdVector<RealT>>(*D-&>getVector(1));

//    std::cout << (*d0.getVector())[0] << std::endl;
//    std::cout << (*d0.getVector())[1] << std::endl;
//    std::cout << (*d1.getVector())[0] << std::endl;
//    std::cout << (*d1.getVector())[1] << std::endl;

    }
    catch(std::logic_error& err) {
        *outStream << err.what() << "\n";
        errorFlag = -1000;
    }; // end try

    if (errorFlag != 0) {
        std::cout << "End Result: TEST FAILED\n";
    }
    else {
        std::cout << "End Result: TEST PASSED\n";
    }

    return 0;
}

