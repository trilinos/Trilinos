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
//                    Jeremie Gaidamour (jngaida@sandia.gov)
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
#include <iostream>
#include <complex>

// Xpetra and Galeri
#include <Xpetra_MultiVectorFactory.hpp>
#include <Galeri_XpetraParameters.hpp>
#include <Galeri_XpetraProblemFactory_Helmholtz.hpp>
#include <Galeri_XpetraUtils.hpp>
#include <Galeri_XpetraMaps.hpp>

// MueLu
#include <MueLu_ShiftedLaplacian.hpp>
#include <MueLu_UseDefaultTypesComplex.hpp>
#include <MueLu_UseShortNames.hpp>
#include <MueLu_MutuallyExclusiveTime.hpp>

typedef Tpetra::Vector<SC,LO,GO,NO>                  TVEC;
typedef Tpetra::MultiVector<SC,LO,GO,NO>             TMV;
typedef Tpetra::CrsMatrix<SC,LO,GO,NO,LMO>           TCRS;
typedef Xpetra::CrsMatrix<SC,LO,GO,NO,LMO>           XCRS;
typedef Xpetra::TpetraCrsMatrix<SC,LO,GO,NO,LMO>     XTCRS; 
typedef Xpetra::Matrix<SC,LO,GO,NO,LMO>              XMAT;
typedef Xpetra::CrsMatrixWrap<SC,LO,GO,NO,LMO>       XWRAP;

int main(int argc, char *argv[]) {

  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::TimeMonitor;

  Teuchos::GlobalMPISession mpiSession(&argc, &argv, NULL);
  RCP< const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();
  Teuchos::CommandLineProcessor clp(false);

  //***********************//
  //   Galeri Parameters   //
  //***********************//

  GO nx, ny, nz;
  int mx, my, mz;
  double stretchx, stretchy, stretchz, h, delta;
  int PMLXL, PMLXR, PMLYL, PMLYR, PMLZL, PMLZR;
  double omega, shift;

  std::ifstream inputfile;
  inputfile.open("helm3D.xml");
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

  Galeri::Xpetra::Parameters<GO> matrixParameters_laplace  (clp, nx, ny, nz, "HelmholtzFEM3D", 0, stretchx, stretchy, stretchz,
							    h, delta, 0,     0,     0,     0,     0,     0,     0.0,   0.0,   mx,   my,   mz  );
  Galeri::Xpetra::Parameters<GO> matrixParameters_helmholtz(clp, nx, ny, nz, "HelmholtzFEM3D", 0, stretchx, stretchy, stretchz,
							    h, delta, PMLXL, PMLXR, PMLYL, PMLYR, PMLZL, PMLZR, omega, 0.0,   mx,   my,   mz  );
  Galeri::Xpetra::Parameters<GO> matrixParameters_shift    (clp, nx, ny, nz, "HelmholtzFEM3D", 0, stretchx, stretchy, stretchz,
							    h, delta, PMLXL, PMLXR, PMLYL, PMLYR, PMLZL, PMLZR, omega, shift, mx,   my,   mz  );
  Xpetra::Parameters             xpetraParameters(clp);

  //****************************************//
  //   Setup Galeri Problems and Matrices   //
  //****************************************//

  RCP<TimeMonitor> globalTimeMonitor = rcp (new TimeMonitor(*TimeMonitor::getNewTimer("ScalingTest: S - Global Time")));
  RCP<TimeMonitor> tm = rcp (new TimeMonitor(*TimeMonitor::getNewTimer("ScalingTest: 1 - Matrix Build")));

  Teuchos::ParameterList pl = matrixParameters_helmholtz.GetParameterList();
  RCP<MultiVector> coordinates;
  Teuchos::ParameterList galeriList;
  galeriList.set("nx", pl.get("nx", nx));
  galeriList.set("ny", pl.get("ny", ny));
  galeriList.set("nz", pl.get("nz", nz));
  RCP<const Map> map;

  if (matrixParameters_helmholtz.GetMatrixType() == "Helmholtz1D") {
    map = MapFactory::Build(xpetraParameters.GetLib(), matrixParameters_helmholtz.GetNumGlobalElements(), 0, comm);
    coordinates = Galeri::Xpetra::Utils::CreateCartesianCoordinates<SC, LO, GO, Map, MultiVector>("1D", map, matrixParameters_helmholtz.GetParameterList());
  }
  else if (matrixParameters_helmholtz.GetMatrixType() == "HelmholtzFEM2D") {
    map = Galeri::Xpetra::CreateMap<LO, GO, Node>(xpetraParameters.GetLib(), "Cartesian2D", comm, galeriList);
    coordinates = Galeri::Xpetra::Utils::CreateCartesianCoordinates<SC, LO, GO, Map, MultiVector>("2D", map, matrixParameters_helmholtz.GetParameterList());
  }
  else if (matrixParameters_helmholtz.GetMatrixType() == "HelmholtzFEM3D") {
    map = Galeri::Xpetra::CreateMap<LO, GO, Node>(xpetraParameters.GetLib(), "Cartesian3D", comm, galeriList);
    coordinates = Galeri::Xpetra::Utils::CreateCartesianCoordinates<SC, LO, GO, Map, MultiVector>("3D", map, matrixParameters_helmholtz.GetParameterList());
  }

  RCP<const Tpetra::Map<LO, GO, NO> > tmap = Xpetra::toTpetra(map);

  Teuchos::ParameterList matrixParams_laplace   = matrixParameters_laplace.GetParameterList();
  Teuchos::ParameterList matrixParams_helmholtz = matrixParameters_helmholtz.GetParameterList();
  Teuchos::ParameterList matrixParams_shift     = matrixParameters_shift.GetParameterList();

  RCP<Galeri::Xpetra::Problem<Map,CrsMatrixWrap,MultiVector> > Pr_laplace =
      Galeri::Xpetra::BuildProblem<SC,LO,GO,Map,CrsMatrixWrap,MultiVector>(matrixParameters_laplace.GetMatrixType(), map, matrixParams_laplace);
  RCP<Matrix> A_laplace = Pr_laplace->BuildMatrix();

  RCP<Galeri::Xpetra::Problem<Map,CrsMatrixWrap,MultiVector> > Pr_helmholtz =
      Galeri::Xpetra::BuildProblem<SC,LO,GO,Map,CrsMatrixWrap,MultiVector>(matrixParameters_helmholtz.GetMatrixType(), map, matrixParams_helmholtz);
  RCP<Matrix> A_helmholtz = Pr_helmholtz->BuildMatrix();

  RCP<Galeri::Xpetra::Problem<Map,CrsMatrixWrap,MultiVector> > Pr_shift =
      Galeri::Xpetra::BuildProblem<SC,LO,GO,Map,CrsMatrixWrap,MultiVector>(matrixParameters_shift.GetMatrixType(), map, matrixParams_shift);
  RCP<Matrix> A_shift = Pr_shift->BuildMatrix();

  comm->barrier();

  tm = Teuchos::null;

  //*************************************************************//
  //   Setup Shifted Laplacian Preconditioner and Belos Solver   //
  //*************************************************************//

  tm = rcp (new TimeMonitor(*TimeMonitor::getNewTimer("ScalingTest: 2 - MueLu Setup")));

  RCP<ShiftedLaplacian> SLSolver = rcp( new ShiftedLaplacian );
  SLSolver -> setLaplacian(A_laplace);
  SLSolver -> setHelmholtz(A_helmholtz);
  SLSolver -> setShiftedLaplacian(A_shift);
  SLSolver -> setup(omega);

  tm = Teuchos::null;
  
  //************************************//
  //   Solve linear system with Belos   //
  //************************************//

  tm = rcp (new TimeMonitor(*TimeMonitor::getNewTimer("ScalingTest: 3 - LHS and RHS initialization")));

  RCP<TMV> X = Tpetra::createMultiVector<SC,LO,GO,NO>(tmap,1);
  RCP<TMV> B = Tpetra::createMultiVector<SC,LO,GO,NO>(tmap,1);  
  X->putScalar((SC) 0.0);
  B->putScalar((SC) 0.0);
  int pointsourceid=nx*ny*nz/2+nx*ny/2+nx/2;
  if(map->isNodeGlobalElement(pointsourceid)==true) {
    B->replaceGlobalValue(pointsourceid, 0, (SC) 1.0);
  }

  tm = Teuchos::null;

  tm = rcp (new TimeMonitor(*TimeMonitor::getNewTimer("ScalingTest: 4 - Belos Solve")));

  SLSolver -> solve(B,X);

  tm = Teuchos::null;

  globalTimeMonitor = Teuchos::null;

  TimeMonitor::summarize();

} //main
