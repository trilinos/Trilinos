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
    \brief Test TpetraBoundConstraint class
    \author Created by G. von Winckel
*/

#include "ROL_StdBoundConstraint.hpp"
#include "ROL_TpetraBoundConstraint.hpp"

#include "ROL_Types.hpp"

#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_GlobalMPISession.hpp"

#include "Tpetra_DefaultPlatform.hpp"

using Teuchos::RCP;
using Teuchos::rcp; 
using Teuchos::ArrayRCP;

typedef double RealT;

typedef Tpetra::Map<>::local_ordinal_type LO;
typedef Tpetra::Map<>::global_ordinal_type GO;
typedef Tpetra::Map<>::node_type Node;
typedef Tpetra::Map<LO, GO, Node> Map;
typedef Tpetra::MultiVector<RealT, LO, GO, Node> MV;

typedef Tpetra::DefaultPlatform::DefaultPlatformType Platform;

typedef std::vector<RealT> SV;
typedef RCP<MV>            MVP;
typedef RCP<SV>            SVP;


int test(RCP<const Tpetra::Comm<int> > comm, int dim) {

    // Get number of processes
    const int numProc = comm->getSize(); 

    // Total size over all processes
    const int numGblIndices = numProc*dim;

    RCP<Map> map = rcp( new Map(numGblIndices,0,comm) ); 

    int errorFlag = 0;

    // Upper bound is +0.75
    MVP u = rcp( new MV(map,1,true) );
    u->putScalar(0.9);

    // Lower bound is -0.75
    MVP l = rcp( new MV(map,1,true) );
    l->putScalar(-0.9);
 
    // Optimization variable
    MVP x = rcp( new MV(map,1,true) );
    x->randomize();

    ROL::TpetraBoundConstraint<RealT,LO,GO,Node> tcon(l,u);

    ROL::TpetraMultiVector<RealT,LO,GO,Node> X(x);

    // Prune the vector
    tcon.project(X);

    if(!tcon.isFeasible(X)) {
        ++errorFlag; 
    }

    return errorFlag; 
} 




int main(int argc, char *argv[]) {


    int iprint     = argc - 1;
    Teuchos::RCP<std::ostream> outStream;
    Teuchos::oblackholestream bhs; // outputs nothing
    Teuchos::GlobalMPISession mpiSession (&argc, &argv, &bhs);

    if (iprint > 0)
        outStream = Teuchos::rcp(&std::cout, false);
    else
        outStream = Teuchos::rcp(&bhs, false);

    int errorFlag  = 0;

    try {
        
        Platform &platform = Tpetra::DefaultPlatform::getDefaultPlatform();

        RCP<const Tpetra::Comm<int> > comm = platform.getComm();


        // Maximum dimension of test for a given process
        int maxdim = 100; 

        for(int i = 10;i<maxdim;++i) {
            errorFlag += test(comm,10);
        }

        typedef std::vector<int> ivec;  
        Teuchos::RCP<ivec> a_rcp = Teuchos::rcp(new ivec(2,1));
        Teuchos::RCP<ivec> b_rcp = Teuchos::rcp(new ivec(3,2));
        Teuchos::ArrayRCP<Teuchos::RCP<ivec>> v(2); 
        v[0] = a_rcp;
        v[1] = b_rcp;   
        (*b_rcp)[2] = 3;

    }

    catch (std::logic_error err) {
        *outStream << err.what() << "\n";
        errorFlag = -1000;
    }; // end try

    if (errorFlag != 0)
        std::cout << "End Result: TEST FAILED\n";
    else
        std::cout << "End Result: TEST PASSED\n";

    return 0;

}



