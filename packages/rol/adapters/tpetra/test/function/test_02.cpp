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

/*! \file  test_02.cpp
    \brief Test TpetraBoundConstraint project and scaling functions
*/

#include "ROL_TpetraBoundConstraint.hpp"
#include "ROL_Bounds.hpp"

// #include "ROL_Types.hpp"

#include "ROL_StdVector.hpp"
#include "ROL_Stream.hpp"
#include "Teuchos_GlobalMPISession.hpp"

#include "Tpetra_Core.hpp"

using Teuchos::ArrayRCP;

typedef double RealT;

typedef Tpetra::Map<>::local_ordinal_type LO;
typedef Tpetra::Map<>::global_ordinal_type GO;
typedef Tpetra::Map<>::node_type Node;
typedef Tpetra::Map<LO, GO, Node> Map;
typedef Tpetra::MultiVector<RealT, LO, GO, Node> MV;

typedef std::vector<RealT> SV;
typedef ROL::Ptr<MV>            MVP;
typedef ROL::Ptr<SV>            SVP;

RealT calcError(MV &a, const std::vector<RealT> &b, ROL::Ptr<Map> map,
                int numLocalEntries) {
  int j;
  for (int i = 0; i < numLocalEntries; i++) {
    j = map->getGlobalElement(i);
    a.sumIntoLocalValue(i, 0, -b[j]);
  }
  return a.getVector(0)->normInf();
}

void initialize(MV &a, const std::vector<RealT> &b, ROL::Ptr<Map> map,
                int numLocalEntries) {
  int j;
  for (int i = 0; i < numLocalEntries; i++) {
    j = map->getGlobalElement(i);
    a.sumIntoLocalValue(i, 0, b[j]);
  }
}

int testRandomInputs(ROL::Ptr<const Teuchos::Comm<int> > comm, RealT tol,
                     ROL::Ptr<std::ostream> outStream) {

  const int numProcs = comm->getSize();       // number of processes
  const int numLocalEntries  = 1000;          // number of entries per process
  // Total size over all processes
  const int numGlobalEntries = numProcs*numLocalEntries;

  int   errorFlag = 0;
  RealT errorInftyNorm;

  // Generate standard vectors that hold data.
  ROL::Ptr<std::vector<RealT>> vsv
    = ROL::makePtr<std::vector<RealT>>(numGlobalEntries, 0.0);
  ROL::Ptr<std::vector<RealT>> xsv
    = ROL::makePtr<std::vector<RealT>>(numGlobalEntries, 0.0);
  ROL::Ptr<std::vector<RealT>> gsv
    = ROL::makePtr<std::vector<RealT>>(numGlobalEntries, 0.0);
  ROL::Ptr<std::vector<RealT>> lsv
    = ROL::makePtr<std::vector<RealT>>(numGlobalEntries, 0.0);
  ROL::Ptr<std::vector<RealT>> usv
    = ROL::makePtr<std::vector<RealT>>(numGlobalEntries, 0.0);
  // Include space for storing results.
  ROL::Ptr<std::vector<RealT>> outsv
    = ROL::makePtr<std::vector<RealT>>(numGlobalEntries, 0.0);

  // Use these standard vectors to define ROL::StdVectors (or, in the case of l
  // and u, pointers to ROL::Vectors).
  ROL::StdVector<RealT> vs(vsv);
  ROL::StdVector<RealT> xs(xsv);
  ROL::StdVector<RealT> gs(gsv);
  ROL::Ptr<ROL::Vector<RealT>> lsp = ROL::makePtr<ROL::StdVector<RealT>>(lsv);
  ROL::Ptr<ROL::Vector<RealT>> usp = ROL::makePtr<ROL::StdVector<RealT>>(usv);
  ROL::StdVector<RealT> outs(outsv);

  // Initialize.
  lsp->randomize(-10.0, 10.0);
  usp->randomize(  0.0, 10.0);
  usp->plus(*lsp);
  xs.randomize(-20.0, 20.0);
  vs.randomize(- 2.0,  2.0);
  gs.randomize(-20.0, 20.0);

  ROL::Ptr<Map> map = ROL::makePtr<Map>(numGlobalEntries, 0, comm);

  // Build Tpetra Multivectors.
  MVP vmv   = ROL::makePtr<MV>(map, 1, true);
  MVP xmv   = ROL::makePtr<MV>(map, 1, true);
  MVP gmv   = ROL::makePtr<MV>(map, 1, true);
  MVP umv   = ROL::makePtr<MV>(map, 1, true);
  MVP lmv   = ROL::makePtr<MV>(map, 1, true);
  MVP outmv = ROL::makePtr<MV>(map, 1, true);
  ROL::TpetraMultiVector<RealT,LO,GO,Node> vt(vmv);
  ROL::TpetraMultiVector<RealT,LO,GO,Node> xt(xmv);
  ROL::TpetraMultiVector<RealT,LO,GO,Node> gt(gmv);
  ROL::TpetraMultiVector<RealT,LO,GO,Node> outt(outmv);
  initialize(*vmv, *vsv, map, numLocalEntries);
  initialize(*gmv, *gsv, map, numLocalEntries);
  initialize(*umv, *usv, map, numLocalEntries);
  initialize(*lmv, *lsv, map, numLocalEntries);

  ROL::Bounds<RealT>                           boundsBC(lsp, usp);
  ROL::TpetraBoundConstraint<RealT,LO,GO,Node> tpetraBC(lmv, umv);

  // Test 1 - Check that project and projectInterior agree between the
  //          Elementwise and Tpetra implementations.

  outs.set(xs);  // project using a copy of x to preserve the data for later
  initialize(*outmv, *outsv, map, numLocalEntries);
  boundsBC.project(outs);
  tpetraBC.project(outt);
  errorInftyNorm = calcError(*outmv, *outsv, map, numLocalEntries);
  errorFlag += errorInftyNorm > tol;

  *outStream << std::endl;
  *outStream << "Consistency Check at " << numGlobalEntries
             << " Randomly Sampled Points -- " << std::endl
             << " Infinity Norm of | TpetraBoundConstraint - Elementwise |:"
             << std::endl
             << "  Project         = " << errorInftyNorm << std::endl;

  initialize(*xmv, *xsv, map, numLocalEntries);
  boundsBC.projectInterior(xs);
  tpetraBC.projectInterior(xt);
  errorInftyNorm = calcError(*xmv, *xsv, map, numLocalEntries);
  errorFlag += errorInftyNorm > tol;

  *outStream << "  ProjectInterior = " << errorInftyNorm << std::endl;
  *outStream << std::endl;

  initialize(*xmv, *xsv, map, numLocalEntries);

  // Test 2 - Check that applyInverseScalingFunction and
  //          applyScalingFunctionJacobian agree between the Elementwise and
  //          Tpetra implementations.

  boundsBC.applyInverseScalingFunction(outs, vs, xs, gs);
  tpetraBC.applyInverseScalingFunction(outt, vt, xt, gt);
  errorInftyNorm = calcError(*outmv, *outsv, map, numLocalEntries);
  errorFlag += errorInftyNorm > tol;

  *outStream << "Consistency Check at " << numGlobalEntries
             << " Randomly Sampled Points -- " << std::endl
             << " Infinity Norm of | TpetraBoundConstraint - Elementwise |:"
             << std::endl
             << "  Inverse         = " << errorInftyNorm << std::endl;

  boundsBC.applyScalingFunctionJacobian(outs, vs, xs, gs);
  tpetraBC.applyScalingFunctionJacobian(outt, vt, xt, gt);
  errorInftyNorm = calcError(*outmv, *outsv, map, numLocalEntries);
  errorFlag += errorInftyNorm > tol;

  *outStream << "  Jacobian        = " << errorInftyNorm << std::endl;
  *outStream << std::endl;

  return errorFlag;
}

int main(int argc, char *argv[]) {

    int iprint = argc - 1;
    ROL::Ptr<std::ostream> outStream;
    ROL::nullstream bhs; // outputs nothing
    if (iprint > 0)
        outStream = ROL::makePtrFromRef(std::cout);
    else
        outStream = ROL::makePtrFromRef(bhs);

    Teuchos::GlobalMPISession mpiSession(&argc, &argv, &*outStream);

    int errorFlag = 0;

    try {
        ROL::Ptr<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();
        RealT tol = 1e-8;  // tolerance
        errorFlag += testRandomInputs(comm, tol, outStream);
    }
    catch (std::logic_error& err) {
        *outStream << err.what() << "\n";
        errorFlag = -1000;
    }; // end try

    if (errorFlag != 0)
        std::cout << "End Result: TEST FAILED\n";
    else
        std::cout << "End Result: TEST PASSED\n";

    return 0;

}