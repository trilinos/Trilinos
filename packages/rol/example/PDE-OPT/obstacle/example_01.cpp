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

/*! \file  example_01.cpp
    \brief Shows how to solve the obstacle problem.
*/

#include "Teuchos_Comm.hpp"
#include "ROL_Stream.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"

#include "Tpetra_DefaultPlatform.hpp"
#include "Tpetra_Version.hpp"

#include <iostream>
#include <algorithm>

#include "ROL_Algorithm.hpp"
#include "ROL_Bounds.hpp"
#include "ROL_Bundle.hpp"

#include "../TOOLS/meshmanager.hpp"
#include "../TOOLS/pdeobjective.hpp"
#include "../TOOLS/pdevector.hpp"
#include "pde_obstacle.hpp"
#include "energy_objective.hpp"

typedef double RealT;


template <class Real>
Real evaluateUpperBound(const std::vector<Real> & coord) {
  Real distance(0), one(1);
  if (coord[1] < coord[0]) {
    if (coord[1] < one - coord[0]) {
      distance = coord[1];
    }
    else {
      distance = one - coord[0];
    }
  }
  else {
    if (coord[1] < one - coord[0]) {
      distance = coord[0];
    }
    else {
      distance = 1 - coord[1];
    }
  }
  return distance;
}


template <class Real>
void computeUpperBound(const ROL::Ptr<Tpetra::MultiVector<> > & ubVec,
                       const ROL::Ptr<const FE<Real> > & fe,
                       const ROL::Ptr<Intrepid::FieldContainer<Real> > & cellNodes,
                       const ROL::Ptr<Intrepid::FieldContainer<int> > & cellDofs,
                       const Teuchos::Array<int> & cellIds) {
  int c = fe->gradN()->dimension(0);
  int f = fe->gradN()->dimension(1);
  int d = fe->gradN()->dimension(3);
  ROL::Ptr<Intrepid::FieldContainer<Real> > dofPoints =
    ROL::makePtr<Intrepid::FieldContainer<Real>>(c,f,d);
  fe->computeDofCoords(dofPoints, cellNodes);
  
  std::vector<Real> coord(d);
  for (int i = 0; i < c; ++i) {
    for (int j = 0; j < f; ++j) {
      int fidx = (*cellDofs)(cellIds[i],j);
      if (ubVec->getMap()->isNodeGlobalElement(fidx)) {
        for (int k = 0; k < d; ++k) {
          coord[k] = (*dofPoints)(i,j,k);
        }
        ubVec->replaceGlobalValue(fidx,
                                  0,
                                  evaluateUpperBound<Real>(coord));
      }
    }
  }
} 


int main(int argc, char *argv[]) {
  // This little trick lets us print to std::cout only if a (dummy) command-line argument is provided.
  int iprint     = argc - 1;
  ROL::Ptr<std::ostream> outStream;
  ROL::nullstream bhs; // outputs nothing

  /*** Initialize communicator. ***/
  Teuchos::GlobalMPISession mpiSession (&argc, &argv, &bhs);
  ROL::Ptr<const Teuchos::Comm<int> > comm
    = Tpetra::DefaultPlatform::getDefaultPlatform().getComm();
  const int myRank = comm->getRank();
  if ((iprint > 0) && (myRank == 0)) {
    outStream = ROL::makePtrFromRef(std::cout); 
  }
  else {
    outStream = ROL::makePtrFromRef(bhs);
  }
  int errorFlag  = 0;

  // *** Example body.
  try {

    /*** Read in XML input ***/
    std::string filename = "input.xml";
    Teuchos::RCP<Teuchos::ParameterList> parlist = Teuchos::rcp( new Teuchos::ParameterList() );
    Teuchos::updateParametersFromXmlFile( filename, parlist.ptr() );

    /*** Initialize PDE describing the obstacle problem ***/
    ROL::Ptr<MeshManager<RealT> > meshMgr
      = ROL::makePtr<MeshManager_Rectangle<RealT>>(*parlist);
    ROL::Ptr<PDE_Obstacle<RealT> > pde
      = ROL::makePtr<PDE_Obstacle<RealT>>(*parlist);
    ROL::Ptr<EnergyObjective<RealT> > obj
      = ROL::makePtr<EnergyObjective<RealT>>(pde,meshMgr,comm,*parlist,*outStream);
    ROL::Ptr<Assembler<RealT> > assembler = obj->getAssembler();

    // Create state vector and set to zeroes
    ROL::Ptr<Tpetra::MultiVector<> > u_ptr = assembler->createStateVector();
    u_ptr->randomize();
    ROL::Ptr<ROL::Vector<RealT> > up
      = ROL::makePtr<PDE_PrimalSimVector<RealT>>(u_ptr,pde,assembler);
    // Create state direction vector and set to random
    ROL::Ptr<Tpetra::MultiVector<> > du_ptr = assembler->createStateVector();
    du_ptr->randomize();
    ROL::Ptr<ROL::Vector<RealT> > dup
      = ROL::makePtr<PDE_PrimalSimVector<RealT>>(du_ptr,pde,assembler);

    // Build bound constraints
    ROL::Ptr<Tpetra::MultiVector<> > lo_ptr = assembler->createStateVector();
    ROL::Ptr<Tpetra::MultiVector<> > hi_ptr = assembler->createStateVector();
    lo_ptr->putScalar(0.0); hi_ptr->putScalar(1.0);
    computeUpperBound<RealT>(hi_ptr,pde->getFE(),
                             pde->getCellNodes(),
                             assembler->getDofManager()->getCellDofs(),
                             assembler->getCellIds());
    ROL::Ptr<ROL::Vector<RealT> > lop
      = ROL::makePtr<PDE_PrimalSimVector<RealT>>(lo_ptr,pde,assembler);
    ROL::Ptr<ROL::Vector<RealT> > hip
      = ROL::makePtr<PDE_PrimalSimVector<RealT>>(hi_ptr,pde,assembler);
    ROL::Ptr<ROL::BoundConstraint<RealT> > bnd
      = ROL::makePtr<ROL::Bounds<RealT>>(lop,hip);

    // Run derivative checks
    obj->checkGradient(*up,*dup,true,*outStream);
    obj->checkHessVec(*up,*dup,true,*outStream);

    du_ptr->putScalar(0.4);
    up->set(*dup);
    ROL::Algorithm<RealT> algoTR("Trust Region",*parlist,false);
    std::clock_t timerTR = std::clock();
    algoTR.run(*up,*obj,*bnd,true,*outStream);
    *outStream << "Trust Region Time: "
               << static_cast<RealT>(std::clock()-timerTR)/static_cast<RealT>(CLOCKS_PER_SEC)
               << " seconds." << std::endl << std::endl;

    up->set(*dup);
    ROL::Algorithm<RealT> algoPDAS("Primal Dual Active Set",*parlist,false);
    std::clock_t timerPDAS = std::clock();
    algoPDAS.run(*up,*obj,*bnd,true,*outStream);
    *outStream << "PD Active Set Time: "
               << static_cast<RealT>(std::clock()-timerPDAS)/static_cast<RealT>(CLOCKS_PER_SEC)
               << " seconds." << std::endl << std::endl;

    // Output.
    assembler->printMeshData(*outStream);
    assembler->outputTpetraVector(u_ptr,"state.txt");
    assembler->outputTpetraVector(hi_ptr,"upperBound.txt");

    // Get a summary from the time monitor.
    Teuchos::TimeMonitor::summarize();

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
