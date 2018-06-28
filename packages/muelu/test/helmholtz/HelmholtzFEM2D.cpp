// @HEADER
//
// ***********************************************************************
//
//        MueLu: A package for multigrid based preconditioning
//                  Copyright 2012 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Questions? Contact
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Andrey Prokopenko (aprokop@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
#include <iostream>
#include <complex>

// Teuchos
#include <Teuchos_StandardCatchMacros.hpp>

// Xpetra and Galeri
#include <Xpetra_MultiVectorFactory.hpp>
#include <Galeri_XpetraParameters.hpp>
#include <Galeri_XpetraProblemFactory_Helmholtz.hpp>
#include <Galeri_XpetraUtils.hpp>
#include <Galeri_XpetraMaps.hpp>

// MueLu
#include <MueLu_ShiftedLaplacian.hpp>
#include <MueLu_UseDefaultTypesComplex.hpp>
#include <MueLu_MutuallyExclusiveTime.hpp>

int main(int argc, char *argv[]) {
#include <MueLu_UseShortNames.hpp>

  typedef Tpetra::Vector<SC,LO,GO,NO>              TVEC;
  typedef Tpetra::MultiVector<SC,LO,GO,NO>         TMV;
  typedef Tpetra::CrsMatrix<SC,LO,GO,NO>           TCRS;
  typedef Xpetra::CrsMatrix<SC,LO,GO,NO>           XCRS;
  typedef Xpetra::TpetraCrsMatrix<SC,LO,GO,NO>     XTCRS;
  typedef Xpetra::Matrix<SC,LO,GO,NO>              XMAT;
  typedef Xpetra::CrsMatrixWrap<SC,LO,GO,NO>       XWRAP;

  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::TimeMonitor;

  Teuchos::GlobalMPISession mpiSession(&argc, &argv, NULL);

  bool success = false;
  bool verbose = true;
  try {
    RCP< const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();
    Teuchos::CommandLineProcessor clp(false);

    //***********************//
    //   Galeri Parameters   //
    //***********************//

    GO nx, ny, nz;
    GO mx, my, mz;
    double stretchx, stretchy, stretchz, h, delta;
    int PMLXL, PMLXR, PMLYL, PMLYR, PMLZL, PMLZR;
    double omega, shift;
    int model;

    std::ifstream inputfile;
    inputfile.open("helm2D.inp");
    inputfile >> nx       >> ny       >> nz ;
    if(comm->getRank()==0)
      std::cout<<"nx: "<<nx<<"  ny: "<<ny<<"  nz: "<<nz<<std::endl;
    inputfile >> stretchx >> stretchy >> stretchz ;
    if(comm->getRank()==0)
      std::cout<<"stretchx: "<<stretchx<<"  stretchy: "<<stretchy<<"  stretchz: "<<stretchz<<std::endl;
    inputfile >> h        >> delta ;
    if(comm->getRank()==0)
      std::cout<<"h: "<<h<<"  delta: "<<delta<<std::endl;
    inputfile >> PMLXL    >> PMLXR ;
    if(comm->getRank()==0)
      std::cout<<"PMLXL: "<<PMLXL<<"  PMLXR: "<<PMLXR<<std::endl;
    inputfile >> PMLYL    >> PMLYR ;
    if(comm->getRank()==0)
      std::cout<<"PMLYL: "<<PMLYL<<"  PMLYR: "<<PMLYR<<std::endl;
    inputfile >> PMLZL    >> PMLZR ;
    if(comm->getRank()==0)
      std::cout<<"PMLZL: "<<PMLZL<<"  PMLZR: "<<PMLZR<<std::endl;
    inputfile >> omega    >> shift ;
    if(comm->getRank()==0)
      std::cout<<"omega: "<<omega<<"  shift: "<<shift<<std::endl;
    inputfile >> mx       >> my       >> mz ;
    if(comm->getRank()==0)
      std::cout<<"mx: "<<mx<<"  my: "<<my<<"  mz: "<<mz<<std::endl;
    inputfile >> model ;
    if(comm->getRank()==0)
      std::cout<<"velocity model: "<<model<<std::endl;

    Galeri::Xpetra::Parameters<GO> matrixParameters_helmholtz(clp, nx, ny, nz, "HelmholtzFEM2D", 0, stretchx, stretchy, stretchz,
        h, delta, PMLXL, PMLXR, PMLYL, PMLYR, PMLZL, PMLZR, omega, 0.0, mx, my, mz, model);
    Galeri::Xpetra::Parameters<GO> matrixParameters_shifted(clp, nx, ny, nz, "HelmholtzFEM2D", 0, stretchx, stretchy, stretchz,
        h, delta, PMLXL, PMLXR, PMLYL, PMLYR, PMLZL, PMLZR, omega, shift, mx, my, mz, model);
    Xpetra::Parameters             xpetraParameters(clp);

    //****************************************//
    //   Setup Galeri Problems and Matrices   //
    //****************************************//

    RCP<TimeMonitor> globalTimeMonitor = rcp (new TimeMonitor(*TimeMonitor::getNewTimer("ScalingTest: S - Global Time")));
    RCP<TimeMonitor> tm = rcp (new TimeMonitor(*TimeMonitor::getNewTimer("ScalingTest: 1 - Matrix Build")));

    Teuchos::ParameterList pl = matrixParameters_helmholtz.GetParameterList();
    RCP<RealValuedMultiVector> coordinates;
    Teuchos::ParameterList galeriList;
    galeriList.set("nx", pl.get("nx", nx));
    galeriList.set("ny", pl.get("ny", ny));
    galeriList.set("nz", pl.get("nz", nz));
    RCP<const Map> map;

    map = Galeri::Xpetra::CreateMap<LO, GO, Node>(xpetraParameters.GetLib(), "Cartesian2D", comm, galeriList);
    coordinates = Galeri::Xpetra::Utils::CreateCartesianCoordinates<SC, LO, GO, Map, RealValuedMultiVector>("2D", map, matrixParameters_helmholtz.GetParameterList());

    RCP<const Tpetra::Map<LO, GO, NO> > tmap = Xpetra::toTpetra(map);

    Teuchos::ParameterList matrixParams_helmholtz = matrixParameters_helmholtz.GetParameterList();
    Teuchos::ParameterList matrixParams_shifted   = matrixParameters_shifted.GetParameterList();

    RCP<Galeri::Xpetra::Problem_Helmholtz<Map,CrsMatrixWrap,MultiVector> > Pr_helmholtz =
      Galeri::Xpetra::BuildProblem_Helmholtz<SC,LO,GO,Map,CrsMatrixWrap,MultiVector>(matrixParameters_helmholtz.GetMatrixType(), map, matrixParams_helmholtz);
    RCP<Matrix> Kmat, Mmat;
    std::pair< RCP<Matrix>, RCP<Matrix> > system = Pr_helmholtz->BuildMatrices();
    Kmat = system.first;
    Mmat = system.second;
    RCP<Matrix> Amat = Pr_helmholtz->BuildMatrix();

    RCP<Galeri::Xpetra::Problem_Helmholtz<Map,CrsMatrixWrap,MultiVector> > Pr_shifted =
      Galeri::Xpetra::BuildProblem_Helmholtz<SC,LO,GO,Map,CrsMatrixWrap,MultiVector>(matrixParameters_shifted.GetMatrixType(), map, matrixParams_shifted);
    RCP<Matrix> Pmat = Pr_shifted->BuildMatrix();

    comm->barrier();

    tm = Teuchos::null;

    //*************************************************************//
    //   Setup Shifted Laplacian Preconditioner and Belos Solver   //
    //*************************************************************//

    tm = rcp (new TimeMonitor(*TimeMonitor::getNewTimer("ScalingTest: 2 - MueLu Setup")));

    RCP<ShiftedLaplacian> SLSolver = rcp( new ShiftedLaplacian );
    SLSolver -> setstiff(Kmat);
    SLSolver -> setmass(Mmat);
    SLSolver -> setProblemMatrix(Amat);
    SLSolver -> setPreconditioningMatrix(Pmat);
    SLSolver -> initialize();
    SLSolver -> setupFastRAP();

    tm = Teuchos::null;

    //************************************//
    //   Solve linear system with Belos   //
    //************************************//

    tm = rcp (new TimeMonitor(*TimeMonitor::getNewTimer("ScalingTest: 3 - LHS and RHS initialization")));

    RCP<TMV> X = Tpetra::createMultiVector<SC,LO,GO,NO>(tmap,1);
    RCP<TMV> B = Tpetra::createMultiVector<SC,LO,GO,NO>(tmap,1);
    X->putScalar((SC) 0.0);
    B->putScalar((SC) 0.0);
    int pointsourceid=nx*ny/2+nx/2;
    if(map->isNodeGlobalElement(pointsourceid)==true) {
      B->replaceGlobalValue(pointsourceid, 0, (SC) 1.0);
    }

    tm = Teuchos::null;

    tm = rcp (new TimeMonitor(*TimeMonitor::getNewTimer("ScalingTest: 4 - Belos Solve")));

    SLSolver -> solve(B,X);

    tm = Teuchos::null;

    globalTimeMonitor = Teuchos::null;

    TimeMonitor::summarize();

    success = true;
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(verbose, std::cerr, success);
  return ( success ? EXIT_SUCCESS : EXIT_FAILURE );
} //main
