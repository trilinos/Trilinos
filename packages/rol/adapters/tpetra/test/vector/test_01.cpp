// @HEADER
// ************************************************************************
//
//               Rapid Optimization Library (ROL) Package
//                 Copyright (2014) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact lead developers:
//              Drew Kouri   (dpkouri@sandia.gov) and
//              Denis Ridzal (dridzal@sandia.gov)
//
// ************************************************************************
// @HEADER

/*! \file  test_01.cpp
    \brief Test Tpetra_MultiVector interface.
*/


#include "ROL_Algorithm.hpp"
#include "ROL_TpetraMultiVector.hpp"
#include "ROL_Types.hpp"
#include "ROL_Zakharov.hpp"

#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_GlobalMPISession.hpp"

#include "Tpetra_DefaultPlatform.hpp"

typedef double RealT;
typedef double ElementT;

typedef Tpetra::Map<>::local_ordinal_type LO;
typedef Tpetra::Map<>::global_ordinal_type GO;
typedef Tpetra::Map<>::node_type Node;
typedef Tpetra::Map<LO, GO, Node> Map;
typedef Tpetra::MultiVector<RealT, LO, GO, Node> MV;

int main(int argc, char *argv[]) {

    
     
    typedef ROL::Ptr<MV> MVP;

    Teuchos::GlobalMPISession mpiSession(&argc, &argv,0);
    typedef Tpetra::DefaultPlatform::DefaultPlatformType Platform;
    Platform &platform = Tpetra::DefaultPlatform::getDefaultPlatform();
    ROL::Ptr<const Teuchos::Comm<int> > comm = platform.getComm();
 

    int iprint = argc - 1;
    Teuchos::oblackholestream bhs; // outputs nothing
    std::ostream& outStream = (iprint > 0) ? std::cout : bhs;

    int errorFlag = 0;

    RealT errtol = ROL::ROL_THRESHOLD<RealT>();

    try {
        // Dimension of the optimization vector

        int dim = 10; 
      
        ROL::Ptr<Map> map = ROL::makePtr<Map>(dim,0,comm);

        // Create Tpetra::MultiVectors (single vectors) 
        MVP x_ptr = ROL::makePtr<MV>(map,1,true); 
        MVP y_ptr = ROL::makePtr<MV>(map,1,true); 

        // Random elements
        x_ptr->randomize(); 

        // Set all values to 2
        y_ptr->putScalar(-2.0);
           
        /*---[ Begin Test of ROL::TpetraMultiVector methods ] ---*/

        // Create ROL vectors
        ROL::TpetraMultiVector<RealT> x(x_ptr); // Testing default parameters here
        ROL::TpetraMultiVector<RealT,LO,GO,Node> y(y_ptr);

        // norm of x
        RealT xnorm = x.norm();
        outStream << "\nNorm of ROL::TpetraMultiVector x: " << xnorm << "\n";

        // norm of y
        RealT ynorm = y.norm();
        outStream << "\nNorm of ROL::TpetraMultiVector y: " << ynorm << "\n"; 

        // scale x
        x.scale(0.5);
        RealT xnorm2 = x.norm();
        outStream << "\nNorm of half of x: " << xnorm2 << "\n";
        if ( std::abs(xnorm/xnorm2 - 2.0) > errtol ) {
            outStream << "---> POSSIBLE ERROR ABOVE!\n";
            errorFlag++;
        }

        // clone z from x, deep copy x into z, norm of z
        ROL::Ptr<ROL::Vector<RealT> > z = x.clone();
        z->set(x);
        RealT znorm = z->norm();
        outStream << "\nNorm of ROL::Vector z (clone of x): " << znorm << "\n";
        if ( std::abs(xnorm2 - znorm) > errtol ) {
            outStream << "---> POSSIBLE ERROR ABOVE!\n";
            errorFlag++;
        }

        // compute norm of x - x - 0
        z->set(x);
        x.scale(-1.0);
        z->plus(x);
        y.zero();
        z->axpy(-1.0, y);
        znorm = z->norm();
        outStream << "\nNorm of (x - x) - 0: " << znorm << "\n";
        if ( std::abs(znorm) > errtol ) {
            outStream << "---> POSSIBLE ERROR ABOVE!\n";
            errorFlag++;
        }

        // set x to first basis vector
        z = x.basis(0);
        znorm = z->norm();
        outStream << "\nNorm of ROL::Vector z (first basis vector): " << znorm << "\n";
        if ( std::abs(znorm-1.0) > errtol ) {
            outStream << "---> POSSIBLE ERROR ABOVE!\n";
            errorFlag++;
        }

        // set x to middle basis vector
        z = x.basis(dim/2);
        znorm = z->norm();
        outStream << "\nNorm of ROL::Vector z ('middle' basis vector): " << znorm << "\n";
        if ( std::abs(znorm-1.0) > errtol ) {
            outStream << "---> POSSIBLE ERROR ABOVE!\n";
            errorFlag++;
        } 

        // set x to last basis vector
        z = x.basis(dim-1);
        znorm = z->norm();
        outStream << "\nNorm of ROL::Vector z (last basis vector): " << znorm << "\n";
        if ( std::abs(znorm-1.0) > errtol ) {
            outStream << "---> POSSIBLE ERROR ABOVE!\n";
            errorFlag++;
        }
        /*---[ End Test of ROL::TpetraMultiVector methods ] ---*/
 
        /*---[ Begin Test of optimization using ROL::TpetraMultiVector ]---*/
        // Set to initial guess
        x_ptr->putScalar(4.0);

        // For constructing Zakharov objective
        MVP k_ptr     = ROL::makePtr<MV>(map,1,true);

        // For gradient and Hessian checks
        MVP xtest_ptr = ROL::makePtr<MV>(map,1,true);
        MVP d_ptr     = ROL::makePtr<MV>(map,1,true);
        MVP v_ptr     = ROL::makePtr<MV>(map,1,true);
        MVP hv_ptr    = ROL::makePtr<MV>(map,1,true);
        MVP ihhv_ptr  = ROL::makePtr<MV>(map,1,true);
 
        int numElem = map->getNodeNumElements();
         
        for( LO lclRow = 0; lclRow < static_cast<LO> (numElem); ++lclRow) {
            const GO gblRow = map->getGlobalElement(lclRow);
            k_ptr->replaceGlobalValue(gblRow,0,gblRow+1.0);
        }

        xtest_ptr->randomize();
        d_ptr->randomize();
        v_ptr->randomize();

        ROL::Ptr<ROL::Vector<RealT> > k = ROL::makePtr<ROL::TpetraMultiVector<RealT,LO,GO,Node>>(k_ptr);

        // Check gradient and Hessian
        ROL::TpetraMultiVector<RealT,LO,GO,Node> xtest(xtest_ptr);
        ROL::TpetraMultiVector<RealT,LO,GO,Node> d(d_ptr);
        ROL::TpetraMultiVector<RealT,LO,GO,Node> v(v_ptr);
        ROL::TpetraMultiVector<RealT,LO,GO,Node> hv(hv_ptr);
        ROL::TpetraMultiVector<RealT,LO,GO,Node> ihhv(ihhv_ptr);

        // Create Objective function 
        ROL::ZOO::Objective_Zakharov<RealT> obj(k);

        
        outStream << "\nChecking Gradient and Hessian: " << znorm << "\n";

        obj.checkGradient(xtest, d, true, outStream);                             outStream << "\n"; 
        obj.checkHessVec(xtest, v, true, outStream);                              outStream << "\n";  
        obj.checkHessSym(xtest, d, v, true, outStream);                           outStream << "\n";
   
        // Check inverse Hessian 
        RealT tol=0;
        obj.hessVec(hv,v,xtest,tol);
        obj.invHessVec(ihhv,hv,xtest,tol);
        ihhv.axpy(-1,v);
        outStream << "Checking inverse Hessian" << std::endl;
        outStream << "||H^{-1}Hv-v|| = " << ihhv.norm() << std::endl;
     
        // Set optimization parameters
        Teuchos::ParameterList parlist;

        // Trust Region Parameters
        parlist.sublist("Step").sublist("Trust Region").set("Subproblem Solver", "Truncated CG");
        parlist.sublist("Step").sublist("Trust Region").set("Initial Radius",    0.5);
        parlist.sublist("Step").sublist("Trust Region").set("Maximum Radius",    2.0);

        // Krylov Parameters
        parlist.sublist("General").sublist("Krylov").set("Absolute Tolerance", 1.e-4);
        parlist.sublist("General").sublist("Krylov").set("Relative Tolerance", 1.e-2);
        parlist.sublist("General").sublist("Krylov").set("Iteration Limit",    10);

        // Define Status Test
        parlist.sublist("Status Test").set("Gradient Tolerance", 1e-12);
        parlist.sublist("Status Test").set("Step Tolerance",     1e-14);
        parlist.sublist("Status Test").set("Iteration Limit",    100);

        // Define Algorithm
        ROL::Algorithm<RealT> algo("Trust Region", parlist);

        // Run Algorithm
        algo.run(x, obj, true, outStream);

        // Compute Error
        RealT abserr = x.norm();
        outStream << std::scientific << "\n   Absolute Error: " << abserr << std::endl;
        if ( abserr > sqrt(ROL::ROL_EPSILON<RealT>()) ) {
            errorFlag += 1;
        }
        /*---[ End Test of optimization using ROL::TpetraMultiVector ]---*/

    }

    catch (std::logic_error err) {
        outStream << err.what() << "\n";
        errorFlag = -1000;
    }; // end try

    if (errorFlag != 0)
        std::cout << "End Result: TEST FAILED\n";
    else
        std::cout << "End Result: TEST PASSED\n";

  return 0;
}
