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
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <vector>
#include <iostream>
#include <numeric>

#include "Amesos.h"
#include "Amesos_BaseSolver.h"

#define RegionsSpanProcs  1
#define MultipleRegionsPerProc  2
#ifdef HAVE_MPI
#include "mpi.h"
#endif

#include <Kokkos_DefaultNode.hpp>

#include <MueLu.hpp>
#include <MueLu_EpetraOperator.hpp>
#include <MueLu_DirectSolver.hpp>
#include <MueLu_Hierarchy.hpp>
#include <MueLu_SmootherFactory.hpp>
#include <MueLu_Utilities.hpp>

#include <Xpetra_ConfigDefs.hpp>
#include <Xpetra_EpetraMultiVector.hpp>
#include <Xpetra_Export.hpp>
#include <Xpetra_Import.hpp>
#include <Xpetra_IO.hpp>
#include <Xpetra_Map.hpp>
#include <Xpetra_MultiVector.hpp>
#include <Xpetra_Vector.hpp>

#include <Teuchos_Assert.hpp>
#include <Teuchos_Comm.hpp>
#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_FancyOStream.hpp>
#include <Teuchos_ScalarTraits.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>


#include <SetupRegionHierarchy_def.hpp>
#include <HHG_Utils_def.hpp>

/* This is a driver that is not included in any other file.
 * So, we should be fine to create useful typedefs for Xpetra here and use them in the entire file.
 */
// typedef double Scalar;
// typedef int LocalOrdinal;
// typedef int GlobalOrdinal;
// typedef KokkosClassic::DefaultNode::DefaultNodeType Node;

template<class LocalOrdinal, class GlobalOrdinal, class Node>
using LIDregion_type = LocalOrdinal (*) (void*, const LocalOrdinal, int);

template<class LocalOrdinal, class GlobalOrdinal, class Node>
LocalOrdinal LIDregion(void *ptr, const LocalOrdinal LIDcomp, int whichGrp);

template<class LocalOrdinal, class GlobalOrdinal, class Node>
LocalOrdinal LID2Dregion(void *ptr, const LocalOrdinal LIDcomp, int whichGrp);

extern void edgeGhosts(int ArowPtr[], int Acols[], int &nGhostFound, int ghostCompLIDs[], int edgeLength, int alongX, int ownedEdge, int interiorEdge, int start, int ownedX, int ownedY);

extern void fillCircleSquareData(int ArowPtr[], int Acols[], int ownedX, int ownedY, int Cx, int Cy, int Rx, int Ry, std::vector<int> &appData);

template<class LocalOrdinal, class GlobalOrdinal, class Node>
extern LocalOrdinal LIDregionCircleSquare(void *ptr, const LocalOrdinal compLID,int whichGrp);

template<class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::Array<LocalOrdinal> setLocalNodesPerDim(const std::string& problemType,
                                                 const std::string& caseName, const bool doing1D,
                                                 const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>& rowMap,
                                                 const LocalOrdinal dimX,
                                                 const LocalOrdinal dimY,
                                                 const LocalOrdinal dimZ = 1)
{
  // Number of nodes per x/y/z-direction per processor
  Teuchos::Array<LocalOrdinal> lNodesPerDim(3);

  if (problemType == "structured") {
    if (doing1D) { // One-dimensional problems
      lNodesPerDim[0] = rowMap.getNodeNumElements();
      lNodesPerDim[1] = 1;
      lNodesPerDim[2] = 1;
    }
    else { // Two-dimensional problems

      if (caseName == "caseFifteen")
      {
        lNodesPerDim[0] = 4;
        lNodesPerDim[1] = 4;
        lNodesPerDim[2] = 1;
      }
      else if (caseName == "caseSixteen")
      {
        lNodesPerDim[0] = 7;
        lNodesPerDim[1] = 7;
        lNodesPerDim[2] = 1;
      }
      else if (caseName == "caseSeventeen")
      {
        lNodesPerDim[0] = 31;
        lNodesPerDim[1] = 31;
        lNodesPerDim[2] = 1;
      }
      else if (caseName == "caseEightteen" || caseName == "caseNineteen")
      {
        lNodesPerDim[0] = 16;
        lNodesPerDim[1] = 16;
        lNodesPerDim[2] = 1;
      }
      else if (caseName == "caseTwenty")
      {
        lNodesPerDim[0] = 10;
        lNodesPerDim[1] = 10;
        lNodesPerDim[2] = 1;
      }
      else
      {
        TEUCHOS_TEST_FOR_EXCEPT_MSG(true, "Unknown caseName. Cannot set 'lNodesPerDim' without knowing the spefic case.\n");
      }
    }
  }
  else {

    // Circle-in-a-disk example on 5 processors
    {
      if (rowMap.getComm()->getRank() == 4)
      {
        // This is the unstructured region, so we set dummy values
        lNodesPerDim [0] = -1;
        lNodesPerDim [1] = -1;
        lNodesPerDim [2] = -1;
      }
      else
      {
        lNodesPerDim[0] = dimX;
        lNodesPerDim[1] = dimY;
        lNodesPerDim[2] = dimZ;
      }
    }
  }

  return lNodesPerDim;
}

// Select the type of aggregation in each region
std::string setAggregationTypePerRegion (const std::string& problemType, const int myRank)
{
  std::string aggregationType;

  if (problemType == "structured")
  {
    aggregationType = "structured";
  }
  else
  {
    // Circle-in-a-disk example on 5 processors
    {
      if (myRank == 4)
        aggregationType = "uncoupled";
      else
        aggregationType = "structured";
    }
  }

  return aggregationType;
}

/* To run the region MG solver, first run the Matlab program 'createInput.m'
 * to write a bunch of files with region information to the disk.
 * Then start this executable with the appropriate number of MPI ranks.
 * See unpublished latex document
 * ``Using Trilinos Capabilities to Facilitate Composite-Regional Transitions''
 * for further information.
 */
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
int main_(int argc, char *argv[]) {
#include "Xpetra_UseShortNames.hpp"

#ifdef HAVE_MPI
  MPI_Init(&argc,&argv);
  Epetra_MpiComm CommEpetra(MPI_COMM_WORLD);

  MPI_Group world_group;
  MPI_Comm_group(CommEpetra.Comm(),&world_group);
#else
  Epetra_SerialComm CommEpetra;
#endif

  // wrap communicator into Teuchos::Comm
  Teuchos::RCP<const Teuchos::Comm<int> > Comm = Teuchos::DefaultComm<int>::getComm();

  Xpetra::UnderlyingLib lib = Xpetra::UseEpetra;

  int  myRank;
  FILE *fp;
  char command[40];
  bool doing1D = false;
  int numDimensions = 0;
  int  globalNx, globalNy;
  std::string xmlFileName;
  std::string problemType;
  std::string regionDataDirectory;
  std::string caseName;

  myRank = Comm->getRank();

  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::ParameterList;
  using Teuchos::Array;
  using Teuchos::ArrayRCP;
  using Teuchos::ArrayView;

  typedef typename Teuchos::ScalarTraits<Scalar>::coordinateType real_type;
  typedef Xpetra::MultiVector<real_type,LO,GO,NO> RealValuedMultiVector;

  Comm->barrier();

  // Manually set smoother parameters on level 0
  Array<RCP<Teuchos::ParameterList> > smootherParams(1);
  smootherParams[0] = rcp(new Teuchos::ParameterList);
  smootherParams[0]->set("smoother: type",    "Jacobi");
  smootherParams[0]->set("smoother: sweeps",  20);
  smootherParams[0]->set("smoother: damping", 0.67);

  // read xml filename from command line
  Teuchos::CommandLineProcessor clp;
  {
    // define a help message
    clp.setDocString("Driver for region multigrid\n\nUse Matlab script 'createInput.m' to create necessary input data on the hard dist.\n\nProvide filename of MueLu xml configuration via '--xml=...'.");

    // define command line arguments
    clp.setOption("xml", &xmlFileName, "filename of xml-file with MueLu configuration", true);
    clp.setOption("probType", &problemType, "Problem type [structured, hybrid]", true);
    clp.setOption("regDataDir", &regionDataDirectory, "directory with all region information/files", true);
    clp.setOption("caseName", &caseName, "name of case file", true);

    // force user to specify all options
    clp.recogniseAllOptions(true);

    /* We now parse the command line where argc and argv are passed to
     * the parse method.
     */
    Teuchos::CommandLineProcessor::EParseCommandLineReturn parseReturn = clp.parse(argc, argv);
    if(parseReturn == Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED) {
      return 0;
    }
    if(parseReturn != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL) {
      return 1; // Error!
    }

    // Check for valid command line arguments
    TEUCHOS_TEST_FOR_EXCEPT_MSG(!(problemType == "structured" || problemType == "hybrid"),
        "Unknown problem type. Use either 'structured' or 'hybrid'.\n");
  }

  Comm->barrier();

  // read maxRegPerGID (maximum # of regions shared by any one node) and
  //      maxRegPerProc (maximum # of partial regions owned by any proc)
  //      whichCase (MultipleRegionsPerProc: No regions spans across multiple
  //                                procs but a proc might own multiple regions
  //                 RegionsSpanProcs: processors own only a piece of 1 region
  //                                but regions may span across many procs

  // Provide some feedback to the user
  if (myRank == 0) {

    std::cout << "User input:" << std::endl
        << "  xml-file with MueLu configuration: " << xmlFileName << std::endl
        << "  problem type: " << problemType << std::endl
        << "  Path to directory with region data: " << regionDataDirectory << std::endl;
  }
  Comm->barrier();

  // Define some basic scalar values
  Scalar scalarZero = Teuchos::ScalarTraits<Scalar>::zero();
  Scalar scalarOne = Teuchos::ScalarTraits<Scalar>::one();

  int maxRegPerGID = 0;
  int maxRegPerProc = 0;
  int whichCase = 0;
  int iii = 0;

  // Read region info from file
  {
    std::stringstream fileNameSS;
    fileNameSS << regionDataDirectory << "/myRegionInfo_" << myRank;
    while ((fp = fopen(fileNameSS.str().c_str(),"r") ) == NULL) sleep(1);

    fgets(command,80,fp);
    sscanf(command,"%d",&maxRegPerGID);
    while ( command[iii] != ' ') iii++;
    sscanf(&(command[iii+1]),"%d",&maxRegPerProc);
    while ( command[iii] == ' ') iii++;
    while ( command[iii] != ' ') iii++;
    sscanf(&(command[iii+1]),"%d",&globalNx);
    while ( command[iii] == ' ') iii++;
    while ( command[iii] != ' ') iii++;
    sscanf(&(command[iii+1]),"%d",&globalNy);
    while ( command[iii] == ' ') iii++;
    while ( command[iii] != ' ') iii++;
    if      (command[iii+1] == 'M') whichCase = MultipleRegionsPerProc;
    else if (command[iii+1] == 'R') whichCase = RegionsSpanProcs;
    else {fprintf(stderr,"%d: head messed up %s\n",myRank,command); exit(1);}
  }

  Comm->barrier();

  // check for 1D or 2D problem
  if (globalNy == 1)
    doing1D = true;
  else
    doing1D = false;

  if (doing1D)
    numDimensions = 1;
  else
    numDimensions = 2;

  // ******************************************************************
  // Application Specific Data for LIDregion()
  // ******************************************************************
  struct widget appData;                            // ****************
  std::vector<GlobalOrdinal>  minGIDComp(maxRegPerProc);      // ****************
  std::vector<GlobalOrdinal>  maxGIDComp(maxRegPerProc);      // ****************
  // ******************************************************************
  // ******************************************************************


  std::vector<LocalOrdinal> myRegions; // regions that myRank owns
  Teuchos::RCP<CrsMatrixWrap> AComp = Teuchos::null; // composite form of matrix
  Teuchos::RCP<Map> mapComp = Teuchos::null; // composite map used to build AComp

  // regionsPerGID[i] lists all regions that share the ith composite GID
  // associated with myRank's composite matrix row Map.
  Teuchos::RCP<Xpetra::MultiVector<LocalOrdinal, LocalOrdinal, GlobalOrdinal, Node> > regionsPerGID = Teuchos::null;

  // regionsPerGIDWithGhosts[i] lists all regions that share the ith composite GID
  // associated with myRank's composite matrix col Map.
  Teuchos::RCP<Xpetra::MultiVector<LocalOrdinal, LocalOrdinal, GlobalOrdinal, Node> > regionsPerGIDWithGhosts = Teuchos::null;

  /* rowMapPerGrp[i] and colMapPerGrp[i] are based on the composite maps, i.e.
   * they don't include duplication of interface DOFs.
   *
   * revisedRowMapPerGrp[i] and revisedColMapPerGrp[i] are build manually to
   * account for duplicated, but distinct interface nodes, i.e. they inlucde
   * new GIDs for the duplicated nodes.
   *
   * We make some assumptions:
   * - For the composite as well as the region maps , we assume that
   *   row, column, and range map to be the same. So, we only need a
   *   row and a column map. Range and domain map can be taken as the row map.
   * - For the quasiRegional maps, we only need row and column maps. FillComplete()
   *   will generate some range and domain maps, though we will never use them
   *   since the quasiRegional matrices are only an intermediate step and are
   *   never used for actual calculations.
   */

  std::vector <Teuchos::RCP<Map> > rowMapPerGrp(maxRegPerProc); // row map associated with myRank's ith region in composite layout
  std::vector <Teuchos::RCP<Map> > colMapPerGrp(maxRegPerProc); // column map associated with myRank's ith region in composite layout
  std::vector <Teuchos::RCP<Map> > revisedRowMapPerGrp(maxRegPerProc); // revised row map associated with myRank's ith region for regional layout
  std::vector <Teuchos::RCP<Map> > revisedColMapPerGrp(maxRegPerProc); // revised column map associated with myRank's ith region for regional layout

  std::vector<Teuchos::RCP<Import> > rowImportPerGrp(maxRegPerProc); // row importers per group
  std::vector<Teuchos::RCP<Import> > colImportPerGrp(maxRegPerProc); // column importers per group
  std::vector<Teuchos::RCP<Export> > rowExportPerGrp(maxRegPerProc); // row exporters per group

  std::vector<Teuchos::RCP<Matrix> > quasiRegionGrpMats(maxRegPerProc); // region-wise matrices with quasiRegion maps (= composite GIDs)
  std::vector<Teuchos::RCP<Matrix> > regionGrpMats(maxRegPerProc); // region-wise matrices in true region layout with unique GIDs for replicated interface DOFs

  Teuchos::RCP<Epetra_Map> coarseCompRowMap; // composite row map on the coarse grid
  std::vector<Teuchos::RCP<Map> > coarseRowMapPerGrp(maxRegPerProc); // region-wise row map in true region layout with unique GIDs for replicated interface DOFs
  std::vector<Teuchos::RCP<Map> > coarseQuasiRowMapPerGrp(maxRegPerProc); // region-wise row map in quasiRegion layout with original GIDs from fine level
  std::vector<Teuchos::RCP<Map> > coarseColMapPerGrp(maxRegPerProc); // region-wise columns map in true region layout with unique GIDs for replicated interface DOFs
  std::vector<Teuchos::RCP<Map> > coarseAltColMapPerGrp(maxRegPerProc); // region-wise columns map in true region layout with unique GIDs for replicated interface DOFs
  std::vector<Teuchos::RCP<Import> > coarseRowImportPerGrp(maxRegPerProc); // coarse level row importer per group
  std::vector<Teuchos::RCP<Matrix> > regionGrpProlong(maxRegPerProc); // region-wise prolongator in true region layout with unique GIDs for replicated interface DOFs
  std::vector<Teuchos::RCP<Matrix> > regionAltGrpProlong(maxRegPerProc); // region-wise prolongator in true region layout with unique GIDs for replicated interface DOFs
  std::vector<Teuchos::RCP<Matrix> > regCoarseMatPerGrp(maxRegPerProc); // coarse level operator 'RAP' in region layout

  std::vector<Teuchos::RCP<Vector> > regionInterfaceScaling(maxRegPerProc);
  std::vector<Teuchos::RCP<Vector> > coarseRegionInterfaceScaling(maxRegPerProc);

  Teuchos::RCP<Vector> compX = Teuchos::null; // initial guess for truly composite calculations
  Teuchos::RCP<Vector> compY = Teuchos::null; // result vector for truly composite calculations
  Teuchos::RCP<Vector> regYComp = Teuchos::null; // result vector in composite layout, but computed via regional operations

  Array<Teuchos::RCP<Vector> > quasiRegX(maxRegPerProc); // initial guess associated with myRank's ith region in quasiRegional layout
  Array<Teuchos::RCP<Vector> > quasiRegY(maxRegPerProc); // result vector associated with myRank's ith region in quasiRegional layout
  Array<Teuchos::RCP<Vector> > regX(maxRegPerProc); // initial guess associated with myRank's ith region in regional layout
  Array<Teuchos::RCP<Vector> > regY(maxRegPerProc); // result vector associated with myRank's ith region in regional layout

  std::vector<LocalOrdinal> intIDs; // LIDs of interface DOFs
  std::vector<std::vector<LocalOrdinal> > regIntIDs(maxRegPerProc); // LIDs of interface DOFs

  /* Stuff for multi-level algorithm
   *
   * To allow for multi-level schemes with more than two levels, we need to store
   * maps, matrices, vectors, and stuff like that on each level. Since we call the
   * multi-level scheme recursively, this should be reflected in the design of
   * variables.
   *
   * We use Teuchos::Array<T> to store each quantity on each level.
   */
  int numLevels = 0;
  Array<RCP<Map> > compRowMaps; // composite row maps on each level
  Array<RCP<Map> > compColMaps; // composite columns maps on each level
  Array<std::vector<RCP<Map> > > regRowMaps; // regional row maps on each level
  Array<std::vector<RCP<Map> > > quasiRegRowMaps; // quasiRegional row maps on each level
  Array<std::vector<RCP<Map> > > regColMaps; // regional column maps on each level
  Array<std::vector<RCP<Map> > > quasiRegColMaps; // quasiRegional column maps on each level
  Array<std::vector<RCP<Matrix> > > regMatrices; // regional matrices on each level
  Array<std::vector<RCP<Matrix> > > regProlong; // regional prolongators on each level
  Array<std::vector<RCP<Import> > > regRowImporters; // regional row importers on each level
  Array<Array<RCP<Vector> > > regInterfaceScalings; // regional interface scaling factors on each level

  Teuchos::RCP<Matrix> coarseCompOp = Teuchos::null;

  /* The actual computations start here. It's a sequence of operations to
   * - read region data from files
   * - create composite and region matrices
   * - create a region multigrid hierarchy
   * - run a region-type V-cycle
   */

  Comm->barrier();

  // Load composite map
  {
    std::cout << myRank << " | Loading composite map ..." << std::endl;

    std::stringstream fileNameSS;
    fileNameSS << regionDataDirectory << "/myCompositeMap_" << myRank;
    fp = fopen(fileNameSS.str().c_str(), "r");
    if (fp == NULL)
      std::cout << std::endl << ">>> Check number of MPI ranks!" << std::endl << std::endl;
    TEUCHOS_ASSERT(fp!=NULL);

    Teuchos::Array<GlobalOrdinal> fileData; // composite GIDs
    GlobalOrdinal i;
    while (fscanf(fp, "%d", &i) != EOF)
      fileData.push_back(i);

    Teuchos::ArrayView<GlobalOrdinal> fileDataView = fileData;

    mapComp = MapFactory::Build(lib, Teuchos::OrdinalTraits<GlobalOrdinal>::invalid(), fileDataView,
        Teuchos::OrdinalTraits<GlobalOrdinal>::zero(), Comm);
//    mapComp->describe(*fos, Teuchos::VERB_EXTREME);
  }

  Comm->barrier();

  // Read matrix from file
  {
    std::cout << myRank << " | Reading matrix from file ..." << std::endl;

    std::stringstream fileNameSS;
    fileNameSS << regionDataDirectory << "/Amat.mm";

    AComp = Teuchos::rcp_dynamic_cast<CrsMatrixWrap>(IO::Read(fileNameSS.str(), mapComp));
//    AComp->describe(*fos, Teuchos::VERB_EXTREME);
  }

  Comm->barrier();

  // Load and communicate region assignments
  {
    std::cout << myRank << " | Loading and communicating region assignments ..." << std::endl;

    regionsPerGID = Xpetra::MultiVectorFactory<LocalOrdinal, LocalOrdinal, GlobalOrdinal, Node>::Build(AComp->getRowMap(), maxRegPerGID, true);

    std::stringstream fileNameSS;
    fileNameSS << regionDataDirectory << "/myRegionAssignment_" << myRank;
    fp = fopen(fileNameSS.str().c_str(), "r");
    TEUCHOS_ASSERT(fp!=NULL);

    int k;
    RCP<Xpetra::Vector<LocalOrdinal, LocalOrdinal, GlobalOrdinal, Node> > jthRegions; // pointer to jth column in regionsPerGID
    for (LocalOrdinal i = 0; i < Teuchos::as<LocalOrdinal>(mapComp->getNodeNumElements()); i++) {
      for (int j = 0; j < maxRegPerGID; j++) {
        jthRegions = regionsPerGID->getVectorNonConst(j);

        if (fscanf(fp, "%d", &k) == EOF) {
          fprintf(stderr, "not enough region assignments\n");
          exit(1);
        }
        jthRegions->replaceLocalValue(i, (Scalar) k);

        // identify interface DOFs. An interface DOF is assigned to at least two regions.
        if (j > 0 and k != -1)
          intIDs.push_back(i);
      }
    }

    // make extended Region Assignments
    RCP<Import> Importer = ImportFactory::Build(mapComp, AComp->getColMap());
    regionsPerGIDWithGhosts = Xpetra::MultiVectorFactory<LocalOrdinal, LocalOrdinal, GlobalOrdinal, Node>::Build(AComp->getColMap(), maxRegPerGID, true);
    regionsPerGIDWithGhosts->doImport(*regionsPerGID, *Importer, Xpetra::INSERT);
//    regionsPerGIDWithGhosts->describe(*fos, Teuchos::VERB_EXTREME);
  }

  Comm->barrier();

  // Load Regions
  {
    std::cout << myRank << " | Loading regions ..." << std::endl;

    std::stringstream fileNameSS;
    fileNameSS << regionDataDirectory << "/myRegions_" << myRank;
    fp = fopen(fileNameSS.str().c_str(), "r");
    TEUCHOS_ASSERT(fp!=NULL);

    while (fgets(command, 80, fp) != NULL) {
      int i;
      sscanf(command, "%d", &i);
      myRegions.push_back(i);
    }
  }

  Comm->barrier();

  std::vector<LocalOrdinal> genericVector;
  // Load AppData forLID region
  {
    std::stringstream fileNameSS;
    fileNameSS << regionDataDirectory << "/myAppData_" << myRank;
    fp = fopen(fileNameSS.str().c_str(), "r");
    TEUCHOS_ASSERT(fp!=NULL);

    int retval; // dummy return value for fscanf() to suppress compiler warnings
    appData.maxRegPerProc = maxRegPerProc;

    if (problemType == "structured") {

      // Fill this with dummy entries. Will not be used in fully structured problems.
      genericVector.resize(3);
      genericVector[inpData_isStructured] = 0;
      genericVector[inpData_ownedX] = -1;
      genericVector[inpData_ownedY] = -1;

      if (doing1D) {
        int minGID, maxGID;
        for (int i = 0; i < (int) myRegions.size(); i++) {
          retval = fscanf(fp,"%d%d",&minGID,&maxGID);
          if(retval == 0) {std::cout << "Something probably went wrong while reading minGID and maxGID from file!" << std::endl;}
          minGIDComp[i] = minGID;
          maxGIDComp[i] = maxGID;
        }
        appData.gDim = (int *) malloc(sizeof(int)*3*myRegions.size());
        appData.lDim = (int *) malloc(sizeof(int)*3*myRegions.size());
        appData.lowInd= (int *) malloc(sizeof(int)*3*myRegions.size());
        for (int i = 0; i < (int) myRegions.size(); i++) {
          retval = fscanf(fp,"%d%d%d",&(appData.gDim[3*i]),&(appData.gDim[3*i+1]),&(appData.gDim[3*i+2]));
          retval = fscanf(fp,"%d%d%d",&(appData.lDim[3*i]),&(appData.lDim[3*i+1]),&(appData.lDim[3*i+2]));
          retval = fscanf(fp,"%d%d%d",&(appData.lowInd[3*i]),&(appData.lowInd[3*i+1]),&(appData.lowInd[3*i+2]));
        }
        appData.minGIDComp = minGIDComp.data();
        appData.maxGIDComp = maxGIDComp.data();
      } else {
        appData.nx = globalNx;
        appData.gDim = (int *) malloc(sizeof(int)*3*myRegions.size());
        appData.lDimx= (int *) malloc(sizeof(int)*myRegions.size());
        appData.lDimy= (int *) malloc(sizeof(int)*myRegions.size());
        appData.trueCornerx= (int *) malloc(sizeof(int)*myRegions.size());
        appData.trueCornery= (int *) malloc(sizeof(int)*myRegions.size());
        appData.relcornerx= (int *) malloc(sizeof(int)*myRegions.size());
        appData.relcornery= (int *) malloc(sizeof(int)*myRegions.size());
        int garbage;
        for (int i = 0; i < (int) myRegions.size(); i++) {
          retval = fscanf(fp,"%d%d%d",&(appData.gDim[3*i]),&(appData.gDim[3*i+1]),&(appData.gDim[3*i+2]));
          retval = fscanf(fp,"%d%d%d",&(appData.lDimx[i]),&(appData.lDimy[i]),&garbage);
          retval = fscanf(fp,"%d%d%d",&(appData.relcornerx[i]),&(appData.relcornery[i]),&garbage);
          retval = fscanf(fp,"%d%d%d",&(appData.trueCornerx[i]),&(appData.trueCornery[i]),&garbage);
        }
      }
    } else {
      TEUCHOS_TEST_FOR_EXCEPT_MSG(true, "Hybrid case not migrated to Xpetra, yet.");
//      int isStructured;
//      retval = fscanf(fp,"%d",&isStructured);
//      if (isStructured == 0) {
//         genericVector.resize(3);
//         genericVector[inpData_isStructured] = 0;
//         retval = fscanf(fp,"%d",&(genericVector[inpData_ownedX]));
//         genericVector[inpData_ownedY] = 1;
//      }
//      else {
//         int Cx, Cy, Rx, Ry;
//         int ownedX,ownedY;
//         retval = fscanf(fp,"%d%d",&Rx,&Ry);
//         retval = fscanf(fp,"%d%d",&ownedX,&ownedY);
//         retval = fscanf(fp,"%d%d",&Cx,&Cy);
//         int *rowptr, *cols;  double *values;
//         AComp->ExtractCrsDataPointers(rowptr, cols, values);
//         fillCircleSquareData(rowptr, cols, ownedX,ownedY, Cx, Cy, Rx, Ry,genericVector);
//      }
//      fclose(fp);
    }

    appData.maxRegPerGID = maxRegPerGID;
    appData.myRegions = myRegions.data();
    appData.colMap = (Map*) &*(AComp->getColMap());
    appData.regionsPerGIDWithGhosts = regionsPerGIDWithGhosts;
    appData.myRank = myRank;
  }

  Comm->barrier();

  // Select the appropriate function pointer
  LIDregion_type<LO, GO, Node> myLIDregion;
  if (problemType == "structured") {
    if (doing1D)
      myLIDregion = &LIDregion<LO, GO, NO>;
    else
      myLIDregion = &LID2Dregion<LO, GO, NO>;
  }
  else {
    myLIDregion = &LIDregionCircleSquare<LO, GO, NO>;
  }

  Comm->barrier();

  // Make group region row maps
  MakeGroupRegionRowMaps<SC, LO, GO, NO, widget>(myRank,
                                                 AComp,
                                                 myLIDregion,
                                                 appData,
                                                 myRegions,
                                                 rowMapPerGrp);
  Comm->barrier();

  // Make group region column maps
  MakeGroupRegionColumnMaps<SC, LO, GO, NO, widget>(myRank,
                                                    whichCase,
                                                    AComp,
                                                    myLIDregion,
                                                    appData,
                                                    myRegions,
                                                    rowMapPerGrp,
                                                    colMapPerGrp);
  Comm->barrier();

  // Make extended group region maps
  MakeGroupExtendedMaps<SC, LO, GO, NO, widget>(myRank,
                                                whichCase,
                                                AComp,
                                                mapComp,
                                                myLIDregion,
                                                appData,
                                                myRegions,
                                                rowMapPerGrp,
                                                colMapPerGrp,
                                                revisedRowMapPerGrp,
                                                revisedColMapPerGrp);

  // Setup importers
  for (int j = 0; j < maxRegPerProc; j++) {
    rowImportPerGrp[j] = ImportFactory::Build(mapComp, rowMapPerGrp[j]);
    colImportPerGrp[j] = ImportFactory::Build(mapComp, colMapPerGrp[j]);
  }
  Comm->barrier();

  // Make quasiRegion matrices
  MakeQuasiregionMatrices(AComp,
                          appData.maxRegPerProc,
                          appData.regionsPerGIDWithGhosts,
                          rowMapPerGrp,
                          colMapPerGrp,
                          rowImportPerGrp,
                          quasiRegionGrpMats);
  Comm->barrier();

  // Make region matrices
  MakeRegionMatrices(AComp,
                     AComp->getRowMap(),
                     rowMapPerGrp,
                     revisedRowMapPerGrp,
                     revisedColMapPerGrp,
                     rowImportPerGrp,
                     appData.maxRegPerProc,
                     quasiRegionGrpMats,
                     regionGrpMats);
  Comm->barrier();

  // Setting up parameters before hierarchy construction
  // These need to stay in the driver as they would be provide by an app
  Array<Array<int> > lNodesPerDim(maxRegPerProc);
  Array<std::string> aggregationRegionType(maxRegPerProc);
  Array<RCP<MultiVector> > nullspace(maxRegPerProc);
  Array<RCP<RealValuedMultiVector> > coordinates(maxRegPerProc);
  for(int j = 0; j < maxRegPerProc; ++j) {
    lNodesPerDim[j] = setLocalNodesPerDim(problemType, caseName, doing1D,
                                          *revisedRowMapPerGrp[j],
                                          genericVector[inpData_regionX],
                                          genericVector[inpData_regionY]);

    // Set aggregation type for each region
    aggregationRegionType[j] = setAggregationTypePerRegion(problemType, myRank);

    // create nullspace vector
    nullspace[j] = MultiVectorFactory::Build(revisedRowMapPerGrp[j], 1);
    nullspace[j]->putScalar(scalarOne);

    // create dummy coordinates vector
    coordinates[j] = Xpetra::MultiVectorFactory<real_type, LocalOrdinal, GlobalOrdinal, Node>::Build(revisedRowMapPerGrp[j], 3);
    coordinates[j]->putScalar(scalarOne);
  }

  // Create multigrid hierarchy
  createRegionHierarchy(maxRegPerProc, numDimensions,
                        lNodesPerDim,
                        aggregationRegionType,
                        xmlFileName,
                        nullspace,
                        coordinates,
                        regionGrpMats,
                        mapComp,
                        rowMapPerGrp,
                        colMapPerGrp,
                        revisedRowMapPerGrp,
                        revisedColMapPerGrp,
                        rowImportPerGrp,
                        compRowMaps,
                        compColMaps,
                        regRowMaps,
                        regColMaps,
                        quasiRegRowMaps,
                        quasiRegColMaps,
                        regMatrices,
                        regProlong,
                        regRowImporters,
                        regInterfaceScalings,
                        coarseCompOp,
                        smootherParams);
  Comm->barrier();

  // Run V-cycle
  {
    std::cout << myRank << " | Running V-cycle ..." << std::endl;

    // Extract the number of levels from the prolongator data structure
    numLevels = regProlong.size();

    TEUCHOS_TEST_FOR_EXCEPT_MSG(!(numLevels>0), "We require numLevel > 0. Probably, numLevel has not been set, yet.");

    /* We first use the non-level container variables to setup the fine grid problem.
     * This is ok since the initial setup just mimics the application and the outer
     * Krylov method.
     *
     * We switch to using the level container variables as soon as we enter the
     * recursive part of the algorithm.
     */

    // initial guess for solution
    compX = VectorFactory::Build(mapComp, true);

// -- NOT MIGRATED TO XPETRA YET -- START
// debugging using shadow.m
//double *z;
//compX->ExtractView(&z); for (int kk = 0; kk < compX->MyLength(); kk++) z[kk] = (double) kk*kk;
//compX->Print(std::cout);
// -- NOT MIGRATED TO XPETRA YET -- END
    // forcing vector
    RCP<Vector> compB = VectorFactory::Build(mapComp, true);

    if (doing1D)
    {
      // 1D
      {
        {
          compB->replaceGlobalValue(compB->getGlobalLength() - 1, 0, Teuchos::as<Scalar>(1.0e-3));
        }
//        {
//        compB->putScalar(scalarOne);
//        compB->replaceGlobalValue(Teuchos::OrdinalTraits<GlobalOrdinal>::zero(), scalarOne);
//        }
//        {
//          compB->replaceGlobalValue(Teuchos::as<GlobalOrdinal>(15), scalarOne);
//        }
//        {
//          compB->replaceGlobalValue(Teuchos::as<GlobalOrdinal>(16), scalarOne);
//        }
      }
    }
    else //2D
    {
      std::vector<GlobalOrdinal> dirichletGIDs;

      // 2D
      {
        const int nx = 19; // global number of nodes in x-direction
        for (int i = 0; i < nx; ++i)
          dirichletGIDs.push_back(i);
        for (int i = 0; i < nx; ++i) {
          dirichletGIDs.push_back(i*nx);
          dirichletGIDs.push_back((i+1)*nx - 1);
        }
        for (int i = 0; i < nx; ++i)
          dirichletGIDs.push_back((nx*(nx-1) + i));
      }

//      for (auto gid : dirichletGIDs)
//        std::cout << gid << ", ";
//      std::cout << std::endl;

      {
        compB->putScalar(Teuchos::as<Scalar>(1.0e-3));
        for (size_t i = 0; i < dirichletGIDs.size(); ++i)
          compB->replaceGlobalValue(dirichletGIDs[i], scalarZero);
      }
    }

    // residual vector
    RCP<Vector> compRes = VectorFactory::Build(mapComp, true);
    {
      AComp->apply(*compX, *compRes);
      compRes->update(1.0, *compB, -1.0);
    }

    // transform composite vectors to regional layout
    compositeToRegional(compX, quasiRegX, regX,
        revisedRowMapPerGrp, rowImportPerGrp);

    Array<RCP<Vector> > quasiRegB(maxRegPerProc);
    Array<RCP<Vector> > regB(maxRegPerProc);
    compositeToRegional(compB, quasiRegB, regB,
        revisedRowMapPerGrp, rowImportPerGrp);

//    printRegionalObject<Vector>("regB 0", regB, myRank, *fos);

    Array<RCP<Vector> > regRes(maxRegPerProc);
    for (int j = 0; j < maxRegPerProc; j++) { // step 1
      regRes[j] = VectorFactory::Build(revisedRowMapPerGrp[j], true);
    }

    /////////////////////////////////////////////////////////////////////////
    // SWITCH TO RECURSIVE STYLE --> USE LEVEL CONTAINER VARIABLES
    /////////////////////////////////////////////////////////////////////////

    // define max iteration counts
    const int maxVCycle = 200;

    // Prepare output of residual norm to file
    RCP<std::ofstream> log;
    if (myRank == 0)
    {
      std::string s = "residual_norm.txt";
      log = rcp(new std::ofstream(s.c_str()));
      (*log) << "# num procs = " << Comm->getSize() << "\n"
             << "# iteration | res-norm\n"
             << "#\n";
    }

    // Richardson iterations
    for (int cycle = 0; cycle < maxVCycle; ++cycle) {
      // check for convergence
      {
        ////////////////////////////////////////////////////////////////////////
        // SWITCH BACK TO NON-LEVEL VARIABLES
        ////////////////////////////////////////////////////////////////////////

        computeResidual(regRes, regX, regB, regionGrpMats,
            revisedRowMapPerGrp, rowImportPerGrp);

//        printRegionalObject<Vector>("regB 1", regB, myRank, *fos);

        compRes = VectorFactory::Build(mapComp, true);
        regionalToComposite(regRes, compRes, rowImportPerGrp);
        typename Teuchos::ScalarTraits<Scalar>::magnitudeType normRes = compRes->norm2();

        // Output current residual norm to screen (on proc 0 only)
        if (myRank == 0)
        {
          std::cout << cycle << "\t" << normRes << std::endl;
          (*log) << cycle << "\t" << normRes << "\n";
        }

        if (normRes < 1.0e-12)
          break;
      }

      /////////////////////////////////////////////////////////////////////////
      // SWITCH TO RECURSIVE STYLE --> USE LEVEL CONTAINER VARIABLES
      /////////////////////////////////////////////////////////////////////////

//      printRegionalObject<Vector>("regB 2", regB, myRank, *fos);
      vCycle(0, numLevels,
             regX, regB, regMatrices,
             regProlong, compRowMaps, quasiRegRowMaps, regRowMaps, regRowImporters,
             regInterfaceScalings, smootherParams, coarseCompOp);
    }

    ////////////////////////////////////////////////////////////////////////
    // SWITCH BACK TO NON-LEVEL VARIABLES
    ////////////////////////////////////////////////////////////////////////

    // -----------------------------------------------------------------------
    // Print fine-level solution
    // -----------------------------------------------------------------------
/*
    Comm->barrier();
    sleep(1);

    // ToDo (mayr.mt) Is this the right CombineMode?
    regionalToComposite(regX, compX, rowImportPerGrp);

    std::cout << myRank << " | compX after V-cycle" << std::endl;
    sleep(1);
    compX->describe(*fos, Teuchos::VERB_HIGH);
    sleep(2);

    // Write solution to file for printing
    std::string outFileName = "compX.mm";
    Xpetra::IO<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Write(outFileName, *compX);
*/
  }

  Comm->barrier();

#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  return 0;
}

// Wrapping main_ to have proper template instantiation
int main(int argc, char *argv[]) {

  const int retVal = main_<double, int, int, KokkosClassic::DefaultNode::DefaultNodeType>(argc, argv);

  return retVal;
}



/* Returns local ID (within region curRegion) for the LIDcomp^th composite grid
 * point that myRank owns. If this grid point is not part of curRegion, then
 * -1 is returned.
 *
 * For fully structured 1D problems only.
 */
 template <class LocalOrdinal, class GlobalOrdinal, class Node>
LocalOrdinal LIDregion(void *ptr, const LocalOrdinal LIDcomp, int whichGrp)
{
#include "Xpetra_UseShortNamesOrdinal.hpp"
   struct widget * myWidget = (struct widget *) ptr;

   int        *minGIDComp  = myWidget->minGIDComp;
   int        *maxGIDComp  = myWidget->maxGIDComp;
   int        *myRegions   = myWidget->myRegions;
   Map        *colMap      = myWidget->colMap;
   int        maxRegPerGID = myWidget->maxRegPerGID;
   Teuchos::RCP<Xpetra::MultiVector<LocalOrdinal, LocalOrdinal, GlobalOrdinal, Node> > regionsPerGIDWithGhosts = myWidget->regionsPerGIDWithGhosts;

   int curRegion = myRegions[whichGrp];
   if (LIDcomp >=  Teuchos::as<int>(colMap->getNodeNumElements())) return(-1);

   const Teuchos::ArrayView<const GlobalOrdinal> colGIDsComp = colMap->getNodeElementList();

   if (colGIDsComp[LIDcomp] < minGIDComp[whichGrp]) return(-1);
   if (colGIDsComp[LIDcomp] > maxGIDComp[whichGrp]) return(-1);

   bool found = false;
   for (int j = 0; j <  maxRegPerGID; j++) {
     Teuchos::ArrayRCP<const LocalOrdinal> jthRegions = regionsPerGIDWithGhosts->getData(j);
      if  (jthRegions[LIDcomp] == curRegion ) {
         found = true;
         break;
      }
   }
   if (found == false) return(-1);

   return( colGIDsComp[LIDcomp] - minGIDComp[whichGrp] );
}

/* Returns local ID (within region curRegion) for the LIDcomp^th composite grid
 * point that myRank owns. If this grid point is not part of curRegion, then
 * -1 is returned.
 *
 * For fully structured 2D problems only.
 */
template <class LocalOrdinal, class GlobalOrdinal, class Node>
LocalOrdinal LID2Dregion(void *ptr, LocalOrdinal LIDcomp, int whichGrp)
{
#include "Xpetra_UseShortNamesOrdinal.hpp"
   struct widget * myWidget = (struct widget *) ptr;

//   int       *minGIDComp   = myWidget->minGIDComp;
//   int       *maxGIDComp   = myWidget->maxGIDComp;
   int       *myRegions    = myWidget->myRegions;
   Map       *colMap       = myWidget->colMap;
   int        maxRegPerGID = myWidget->maxRegPerGID;
   int       *trueCornerx  = myWidget->trueCornerx; // global coords of region
   int       *trueCornery  = myWidget->trueCornery; // corner within entire 2D mesh
                                                    // across all regions/procs
   int       *relcornerx   = myWidget->relcornerx;  // coords of corner relative
   int       *relcornery   = myWidget->relcornery;  // to region corner
   int       *lDimx        = myWidget->lDimx;
   int       *lDimy        = myWidget->lDimy;
   Teuchos::RCP<Xpetra::MultiVector<LocalOrdinal, LocalOrdinal, GlobalOrdinal, Node> > regionsPerGIDWithGhosts = myWidget->regionsPerGIDWithGhosts;

   int curRegion = myRegions[whichGrp];
   if (LIDcomp >=  Teuchos::as<int>(colMap->getNodeNumElements())) return(-1);

   const Teuchos::ArrayView<const GlobalOrdinal> colGIDsComp = colMap->getNodeElementList();

   int xGIDComp =  colGIDsComp[LIDcomp]%(myWidget->nx);
   int yGIDComp = (colGIDsComp[LIDcomp] - xGIDComp)/myWidget->nx;

   if (xGIDComp < trueCornerx[whichGrp]+relcornerx[whichGrp]) return(-1);
   if (yGIDComp < trueCornery[whichGrp]+relcornery[whichGrp]) return(-1);
   if (xGIDComp > trueCornerx[whichGrp]+relcornerx[whichGrp]+lDimx[whichGrp]-1) return(-1);
   if (yGIDComp > trueCornery[whichGrp]+relcornery[whichGrp]+lDimy[whichGrp]-1) return(-1);


   bool found = false;
   for (int j = 0; j <  maxRegPerGID; j++) {
     Teuchos::ArrayRCP<const LocalOrdinal> jthRegions = regionsPerGIDWithGhosts->getData(j);
      if  (jthRegions[LIDcomp] == curRegion ) {
         found = true;
         break;
      }
   }
   if (found == false) return(-1);

   return( (yGIDComp - relcornery[whichGrp] - trueCornery[whichGrp])*lDimx[whichGrp]
           + (xGIDComp - relcornerx[whichGrp] - trueCornerx[whichGrp]) );
}



void fillCircleSquareData(int ArowPtr[], int Acols[], int ownedX, int ownedY, int Cx, int Cy, int Rx, int Ry,
std::vector<int> &appData)
{
/* Fills the vector appData so that the function LIDregionCircleSquareWithUnstr2D() can
 * properly map CompositeLIDs to RegionLIDs. Specifically, appData's
 * contents will be given by
 *
 *   appData[inpData_ownedX]  x/y dimensions of rectangle owned by proc, which is
 *   appData[inpData_ownedY]  given as input parameters 'ownedX' and 'ownedY'
 *
 *   appData[inpData_regionX] x/y dimensions of proc's region, which is given
 *   appData[inpData_regionY] as input parameters 'Rx' and 'Ry'
 *
 *   appData[inpData_cornerX] Offset of the lower left corner defined by the
 *   appData[inpData_cornerY] rectangular region piece that is actually owned by
 *                            this processor (in the composite layout). Should be
 *                            either 0 or 1.  So, Cx = Cy=0 means that the processor
 *                            actually owns the lower left corner of the region.
 *                            This is given as input parameters 'Cx' and 'Cy'
 *
 *   appData[k]           Gives the region LID associated with the
 *                        (k-inpData_firstLIDsOfGhosts+ownedX*ownedY)^th
 *                        composite LID  .. for k >= inpData_firstLIDsOfGhosts,
 *                        which is the (k-inpData_firstLIDsOfGhosts)^th ghost
 *                        composite LID.
 *
 *
 * Before filling appData, fillCircleSquareData() first fills
 * ghostCompLIDs with ids of any ghosts associated with a region boundary.
 * The edges are done in the following order: bottom, top, left, right.
 * If a region boundary does not correspond to ghost unknowns in the composite layout,
 * then this boundary is not included in ghostCompLIDs. Once filled, ghostCompLIDs[]
 * should have the following contents:
 *
 *    ghostCompLIDs[startBot:startTop-1]   Ghost composite LIDs for bottom
 *                                         edge of region
 *
 *    ghostCompLIDs[startTop:startLft-1]   Ghost composite LIDs for top
 *                                         edge of region
 *
 *    ghostCompLIDs[startLft:startRgt-1]   Ghost composite LIDs for left edge
 *                                         of region. Does not include corners
 *                                         if they are already included
 *                                         with the bottom or top ghosts
 *
 *    ghostCompLIDs[startRgt:              Ghost composite LIDs for right edge
 *                     nRegionalGhosts ]   of region. Does not include corners
 *                                         if they are already included
 *                                         with the bottom or top ghosts
 *
 * Read comments for edgeGhosts(). The main assumption is that all stencils are full (9
 * point in the interior).
 *
 */

  /* Compute the number of ghosts that lie on the region boundary and where*/
  /* the different region boundaries will start in the vector.             */

  int nRegionalGhosts = 0;
  int startBot = nRegionalGhosts;
  if (Cy==1) {
     nRegionalGhosts= nRegionalGhosts+ownedX;
     if (Cx==1)           nRegionalGhosts=nRegionalGhosts+1;
     if (Rx==ownedX+Cx+1) nRegionalGhosts=nRegionalGhosts+1;
  }
  int startTop = nRegionalGhosts;

  if (Ry == ownedY+Cy+1) {
     nRegionalGhosts= nRegionalGhosts+ownedX;
     if (Cx==1)          nRegionalGhosts=nRegionalGhosts+1;
     if (Rx==ownedX+Cx+1)nRegionalGhosts=nRegionalGhosts+1;
  }
  int startLft = nRegionalGhosts;

  if (Cx==1)             nRegionalGhosts= nRegionalGhosts+ownedY;
  int startRgt = nRegionalGhosts;

  if (Rx == ownedX+Cx+1) nRegionalGhosts= nRegionalGhosts+ownedY;
  std::vector<int> ghostCompLIDs(nRegionalGhosts);

  /* insert ghosts for bottom, top, left, and right edges into ghostCompLIDs */

  int nGhostFound = 0;
  if (Cy==1)            edgeGhosts(ArowPtr, Acols, nGhostFound, ghostCompLIDs.data(), Rx-Cx-1,1,     1,       2,2-Cx,ownedX, ownedY);
  if (Ry == ownedY+Cy+1)edgeGhosts(ArowPtr, Acols, nGhostFound, ghostCompLIDs.data(), Rx-Cx-1,1,ownedY,ownedY-1,2-Cx,ownedX, ownedY);
  if (Cx==1)            edgeGhosts(ArowPtr, Acols, nGhostFound, ghostCompLIDs.data(),ownedY-1,0,     1,       2,   2,ownedX, ownedY);
  if (Rx == ownedX+Cx+1)edgeGhosts(ArowPtr, Acols, nGhostFound, ghostCompLIDs.data(),ownedY-1,0,ownedX,ownedX-1,   2,ownedX, ownedY);

  /* determine the largest ghost LID so that we can allocate enough space */
  int biggest = ownedX*ownedY-1;
  for (int k = 0; k < nRegionalGhosts; k++)
     if (ghostCompLIDs[k] > biggest) biggest = ghostCompLIDs[k];

  // fill appData

  appData.resize(inpData_firstLIDsOfGhosts+biggest-ownedX*ownedY);

  appData[inpData_isStructured] = 1;
  appData[inpData_ownedX] = ownedX;
  appData[inpData_ownedY] = ownedY;
  appData[inpData_regionX] = Rx;
  appData[inpData_regionY] = Ry;
  appData[inpData_cornerX] = Cx;
  appData[inpData_cornerY] = Cy;
  appData[inpData_nGhosts] = biggest-ownedX*ownedY;

  int offset = inpData_firstLIDsOfGhosts-ownedX*ownedY;

  for (int k = 0; k < biggest-ownedX*ownedY; k++) appData[inpData_firstLIDsOfGhosts+k] = -1;
  for (int k=startBot; k < startTop; k++)
     appData[ghostCompLIDs[k]+offset]=k;
  for (int k=startTop; k < startLft; k++)
     appData[ghostCompLIDs[k]+offset] = Rx*(Ry-1)+k-startTop;
  for (int k=startLft; k < startRgt; k++)
     appData[ghostCompLIDs[k]+offset] = Rx*(k+Cy-startLft);
  for (int k=startRgt; k < nRegionalGhosts; k++)
     appData[ghostCompLIDs[k]+offset] = Rx*(k+Cy-startRgt)+Rx-1;

//  for (int i = 0; i < ghostCompLIDs.size(); i++)
//    printf("ghostComp(%d)=%d ", i, ghostCompLIDs[i]);
//  printf("\n");
//  fflush (stdout);

  return;
}

void edgeGhosts(int ArowPtr[], int Acols[], int &nGhostFound, int ghostCompLIDs[], int edgeLength, int alongX, int ownedEdge, int interiorEdge, int start, int ownedX, int ownedY)
{
/*
%
%  Find the local region-oriented ids of a region's shared interface, which is
%  owned by another processor in the composite layout. The situation is basically
%  this
%
%            ?   ?   ?  ?  ?  ?  ?  ?
%             ======================
%            ?|| .   .  .  .  .  .||?
%            ?||                  ||?
%
%  The = and || denote the inter-processor boundary. We know that we have a
%  bunch of externally owned vertices, denoted by ?. The problem is that we
%  are not completely sure which remote local composite id is associated with
%  which vertex on the regular grid layout (as this was done automatically
%  by a fill complete on the composite matrix. To figure this out, this function
%  assumes a 9-pt stencil and that the chunk that each processor owns is at
%  least 3 wide in each dimension.
%
%  The way local region ids are computed is as follows:
%   1) We basically do a series of find(A(row,:)) for matrix rows corresponding
%      to the ownedEdge that is adjacent to the shared interface that we wish
%      to find region ids.
%   2) We assign any non-local columns (those > nOwned) an index, ii, indicating
%      that this column is adjacent to the ii^th point along ownedEdge. Some
%      columns might be adjacent to several points along ownedEdge. These
%      columns end up with the largest ii value, as we over-write the assigned
%      indices. Thus, a typical situation might look like this
%
%            1   2   3  4  5  6  6  6
%             ======================
%            1|| .   .  .  .  .  .||6
%            1||                  ||6
%
%      In this picture, locations corresponding to entries owned by other
%      processors have been assigned a number. The = and || symbols denote
%      the edge of the inter-processor boundary. The .'s lie along ownedEdge
%      So non-local dofs of the kth dot (counting from left to right) are
%      assigned the number k.
%   3) The two lower 1's and the two lower 6's are problematic as we do not
%      wish to consider these vertices as each function invocation is only
%      concerned with the one region edge. So we will stick large numbers
%      (equal to nOwned) in these locations. This is done by probing points
%      that are 1 away (provided externally as interiorEdge) from the corner
%      along orthogonal edges and assigning these a largeNumber. Our example
%      might now look like
%
%            1   2   3  4  5  6  6  6
%             ======================
%            L|| .   .  .  .  .  .||L
%            L|| *               *||L
%            L||                  ||L
%
%      where L denotes newly assigned large numbers and * are the one away
%      point that were just probed to assign large numbers.
%
%   4) We can now determine which remote vertices correspond to the 2D layout.
%      Specifically, we take the 1st vertex along ownedEdge and examine its
%      5 neighbors (e.g., col1,col2,col3,col4,col5). We look at assignedIndices
%      computed in steps 2 and 3. The column with the lowest assignedIndices
%      value (= 1) is the corner point. The column with the next lowest
%      (= 2) is adjacent to this corner point along the desired edge. The
%      column with the next lowest value (= 3) is the next point. We record
%      these 3 column indices in ghostCompLIDs and set assignedIndices for
%      them to a large number (so now 1 2 3 in the above picture would now be
%      replaced by L L L). We now examine the 2nd vertex's remote neighbors
%      which should have the values L L 3 and assign the column associated
%      with the smallest value (= 3) to ghostCompLIDs ... continuing along
%      until we have assigned the entire edge.
%
%  Note: In the above picture we are assuming that all region edges are owned
%  by remote processors. However, this function works as well when either
%  of the orthogonal edges (vertical edges in our example) are owned locally.
%  For example, we might have the case below
%
%                ?   ?  ?  ?  ?  ?  ?
%             ======================
%             || .   .  .  .  .  .||?
%             ||                  ||?
%  Here, the left edge might be a real physical boundary or it might be that
%  the processor owns the shared interface of the region (so it has remotes
%  to the left of the leftmost ||, but there is no need to assign them a
%  local region-oriented id.  In these cases start is not necessarily 1
%  and fullEdgeLength is not necessarily equal to edgeLength
*/

  if (ownedX < 3) { fprintf(stderr,"edges must be longer\n"); exit(1);}
  if (ownedY < 3) { fprintf(stderr,"edges must be longer\n"); exit(1);}

  int  nOwned      = ownedX*ownedY;
  int  largeNumber = nOwned+10;

  std::vector<int> assignedIndices( (ownedX+2)*(ownedY+2) + 1);

  /* Perform steps 1) and 2) described above. */

  int fullEdgeLength, row, *cols, nCols;

  if (alongX)  fullEdgeLength = ownedX;
  else         fullEdgeLength = ownedY;

  for (int ii=1; ii <= fullEdgeLength; ii++) {

    /* non-Ghosts are lexicographically ordered  */

    if   (alongX) row = ownedX*(ownedEdge-1) +     ii   ;
    else          row = ownedX*(   ii    -1) + ownedEdge;
    cols = &(Acols[ArowPtr[row-1]]);
    nCols = ArowPtr[row] - ArowPtr[row-1];
    for (int k = 0; k < nCols; k++) {
       if (cols[k] >= nOwned ) assignedIndices[cols[k]] = ii;
    }
  }

  /* Now assign large numbers (step 3) to 2 closest vertices */
  /* along each orthogonal edge                              */

  /* non-Ghosts are lexicographically ordered  */
  if  (alongX) row = ownedX*(interiorEdge-1)+     1     ;
  else         row =                        interiorEdge;

  cols = &(Acols[ArowPtr[row-1]]);
  nCols = ArowPtr[row] - ArowPtr[row-1];
  bool firstCornerHasGhosts = false;
  for (int k = 0; k < nCols; k++)
    if (cols[k] >= nOwned) {firstCornerHasGhosts = true; assignedIndices[cols[k]]=largeNumber;}

// This is a case not originally considered. I hope it fixes this bug.
// When coding this up, I failed to recognize the case when Cx is 0
// because
if ( (start==2) && (firstCornerHasGhosts == false)) start--;

  /* non-Ghosts are lexicographically ordered  */
  if  (alongX) row = ownedX*(interiorEdge-1) +  edgeLength ;
  else         row = ownedX*( edgeLength -1) + interiorEdge;

  cols = &(Acols[ArowPtr[row-1]]);
  nCols = ArowPtr[row] - ArowPtr[row-1];
  for (int k = 0; k < nCols; k++)
    if (cols[k] >= nOwned) assignedIndices[cols[k]]=largeNumber;

  /* Now perform step 4 by looking at smallest */
  /* assignedIndices along ownedEdge.          */

  int min, kk;

  for (int ii=1; ii <= edgeLength; ii++) {

    /* non-Ghosts are lexicographically ordered  */
    if  (alongX) row = ownedX*(ownedEdge-1) +    ii    ;
    else         row = ownedX*(   ii    -1) + ownedEdge;

    cols = &(Acols[ArowPtr[row-1]]);
    nCols= ArowPtr[row] - ArowPtr[row-1];

    min  = largeNumber-1;  kk = -1;

    for (int k = 0; k < nCols; k++) {
      if (cols[k] >= nOwned ) {
        if (assignedIndices[cols[k]] == min) { fprintf(stderr,"a tie?\n"); exit(1);}
        if (assignedIndices[cols[k]] < min) {
          kk = cols[k];
          min = assignedIndices[cols[k]];
        }
      }
    }
    if ((ii>=start) && (kk != -1)) ghostCompLIDs[nGhostFound++]= kk;
    if (kk != -1) assignedIndices[kk] = largeNumber;

    if (ii==1) {
      for (int kkkk = 1; kkkk <= 2; kkkk++) {
        min  = largeNumber-1;  kk = -1;
        for (int k = 0; k < nCols; k++) {
          if (cols[k] >= nOwned ) {
            if (assignedIndices[cols[k]] == min) { fprintf(stderr,"a Tie?\n"); exit(1);}
            if (assignedIndices[cols[k]] < min) {
              kk = cols[k];
              min = assignedIndices[cols[k]];
            }
          }
        }
        if (kk != -1) {
          ghostCompLIDs[nGhostFound++]= kk;
          assignedIndices[kk] = largeNumber;
        }
      } // for (int kkkk = 1; kkkk <= 2; kkkk++)
    } // if (ii==1)
  } // for (int ii=1; ii <= edgeLength; ii++)
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
LocalOrdinal LIDregionCircleSquare(void *ptr, const LocalOrdinal compLID, int whichGrp)
{
   // Maps composite LIDs to region LIDs for the example
   // corresponding to a circle embedded within a box created
   // by the mkUnstrQuads Matlab functions. Assumes that
   // fillCircleSquareData() has been invoked.

   int* appData      = (int *) ptr;

   if (appData[inpData_isStructured] == 0) return(compLID);

   int* LIDsOfGhosts = (int *) &(appData[inpData_firstLIDsOfGhosts]);
   int  ownedX = appData[inpData_ownedX];
   int  ownedY = appData[inpData_ownedY];
   int  Rx     = appData[inpData_regionX];
   // int  Ry     = appData[inpData_regionY];
   int  Cx     = appData[inpData_cornerX];
   int  Cy     = appData[inpData_cornerY];

   int i,j,ii;
   // local composite ids are assumed to be lexicographical
   // on the owned rectangle. These need to be mapped to
   // the region rectangle.

   if (compLID < ownedX*ownedY) {
      i = (compLID+1)%ownedX;
      if (i==0) i=ownedX;

      j = (compLID+1 - i)/ownedX + 1;

      return(Rx*(j-1+Cy)+i+Cx -1);  // C-style ==> -1
   }
   else {
      ii = compLID - ownedX*ownedY;
      if (ii > appData[inpData_nGhosts] ) return(-1);
      if (LIDsOfGhosts[ii] == -1) return(-1);
      return(LIDsOfGhosts[ii]);
   }
}
