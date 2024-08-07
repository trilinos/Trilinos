// @HEADER
// *****************************************************************************
//   Zoltan2: A package of combinatorial algorithms for scientific computing
//
// Copyright 2012 NTESS and the Zoltan2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file rcbPerformanceZ1.cpp
    \brief rcbPerformance with Zoltan1

    Geometry is a uniform mesh.
*/

#include "Zoltan2_config.h"
#include <zoltan.h>

#include <Zoltan2_Util.hpp>

#include <Teuchos_RCP.hpp>
#include <Teuchos_ArrayView.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_Comm.hpp>
#include <Teuchos_CommHelpers.hpp>

#include <Zoltan2_TestHelpers.hpp>
#include <Tpetra_MultiVector.hpp>
#include <Tpetra_KokkosCompat_DefaultNode.hpp>
#include <GeometricGenerator.hpp>

#include <vector>
#include <string>
#include <ostream>
#include <sstream>
#include <fstream>
using std::string;
using std::vector;
using std::bad_alloc;
using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::Comm;
using Teuchos::ArrayView;
using Teuchos::CommandLineProcessor;

typedef Tpetra::MultiVector<zscalar_t, zlno_t, zgno_t, znode_t> tMVector_t;
typedef Tpetra::Map<zlno_t, zgno_t, znode_t> tMap_t;
typedef tMap_t::node_type znode_t;

//////////////////////////////////////////////////////////////////////////////
// Data structure for data
typedef struct dots {
  vector<vector<float> > weights;
  tMVector_t *coordinates;
} DOTS;

const char param_comment = '#';

string trim_right_copy(
        const string& s,
        const string& delimiters = " \f\n\r\t\v" )
{
    return s.substr( 0, s.find_last_not_of( delimiters ) + 1 );
}

string trim_left_copy(
        const string& s,
        const string& delimiters = " \f\n\r\t\v" )
{
    return s.substr( s.find_first_not_of( delimiters ) );
}

string trim_copy(
        const string& s,
        const string& delimiters = " \f\n\r\t\v" )
{
    return trim_left_copy( trim_right_copy( s, delimiters ), delimiters );
}

void readGeoGenParams(string paramFileName, Teuchos::ParameterList &geoparams, const RCP<const Teuchos::Comm<int> > & comm){
    std::string input = "";
    char inp[25000];
    for(int i = 0; i < 25000; ++i){
        inp[i] = 0;
    }

    bool fail = false;
    if(comm->getRank() == 0){

        std::fstream inParam(paramFileName.c_str());
        if (inParam.fail())
        {
            fail = true;
        }
        if(!fail)
        {
            std::string tmp = "";
            getline (inParam,tmp);
            while (!inParam.eof()){
                if(tmp != ""){
                    tmp = trim_copy(tmp);
                    if(tmp != ""){
                        input += tmp + "\n";
                    }
                }
                getline (inParam,tmp);
            }
            inParam.close();
            for (size_t i = 0; i < input.size(); ++i){
                inp[i] = input[i];
            }
        }
    }



    int size = input.size();
    if(fail){
        size = -1;
    }
    comm->broadcast(0, sizeof(int), (char*) &size);
    if(size == -1){
        throw "File " + paramFileName + " cannot be opened.";
    }
    comm->broadcast(0, size, inp);
    std::istringstream inParam(inp);
    string str;
    getline (inParam,str);
    while (!inParam.eof()){
        if(str[0] != param_comment){
            size_t pos = str.find('=');
            if(pos == string::npos){
                throw  "Invalid Line:" + str  + " in parameter file";
            }
            string paramname = trim_copy(str.substr(0,pos));
            string paramvalue = trim_copy(str.substr(pos + 1));
            geoparams.set(paramname, paramvalue);
        }
        getline (inParam,str);
    }
}


//////////////////////////////////////////////////////////////////////////////
// Zoltan1 query functions

int getNumObj(void *data, int *ierr)
{
  *ierr = 0;
  DOTS *dots = (DOTS *) data;
  return dots->coordinates->getLocalLength();
}

//////////////////////////
int getDim(void *data, int *ierr)
{
  *ierr = 0;
  DOTS *dots = (DOTS *) data;
  int dim =  dots->coordinates->getNumVectors();

  return dim;
}

//////////////////////////
void getObjList(void *data, int numGid, int numLid,
  ZOLTAN_ID_PTR gids, ZOLTAN_ID_PTR lids,
  int num_wgts, float *obj_wgts, int *ierr)
{
  *ierr = 0;
  DOTS *dots = (DOTS *) data;

  size_t localLen = dots->coordinates->getLocalLength();
  const zgno_t *ids =
               dots->coordinates->getMap()->getLocalElementList().getRawPtr();

  if (sizeof(ZOLTAN_ID_TYPE) == sizeof(zgno_t))
    memcpy(gids, ids, sizeof(ZOLTAN_ID_TYPE) * localLen);
  else
    for (size_t i=0; i < localLen; i++)
      gids[i] = static_cast<ZOLTAN_ID_TYPE>(ids[i]);

  if (num_wgts > 0){
    float *wgts = obj_wgts;
    for (size_t i=0; i < localLen; i++)
      for (int w=0; w < num_wgts; w++)
        *wgts++ = dots->weights[w][i];
  }
}

//////////////////////////
void getCoords(void *data, int numGid, int numLid,
  int numObj, ZOLTAN_ID_PTR gids, ZOLTAN_ID_PTR lids,
  int dim, double *coords, int *ierr)
{
  // I know that Zoltan asks for coordinates in gid order.
  if (dim == 3){
  *ierr = 0;
  DOTS *dots = (DOTS *) data;
  double *val = coords;
  const zscalar_t *x = dots->coordinates->getData(0).getRawPtr();
  const zscalar_t *y = dots->coordinates->getData(1).getRawPtr();
  const zscalar_t *z = dots->coordinates->getData(2).getRawPtr();
  for (int i=0; i < numObj; i++){
    *val++ = static_cast<double>(x[i]);
    *val++ = static_cast<double>(y[i]);
    *val++ = static_cast<double>(z[i]);
  }
  }
  else {
      *ierr = 0;
      DOTS *dots = (DOTS *) data;
      double *val = coords;
      const zscalar_t *x = dots->coordinates->getData(0).getRawPtr();
      const zscalar_t *y = dots->coordinates->getData(1).getRawPtr();
      for (int i=0; i < numObj; i++){
        *val++ = static_cast<double>(x[i]);
        *val++ = static_cast<double>(y[i]);
      }


  }
}


//////////////////////////////////////////////////////////////////////////////

enum weightTypes{
  upDown,
  roundRobin,
  increasing,
  numWeightTypes
};

void makeWeights(
  const RCP<const Teuchos::Comm<int> > & comm,
  vector<float> &wgts, weightTypes how, float scale, int rank)
{
  zlno_t len = wgts.size();
  if (how == upDown){
    float val = scale + rank%2;
    for (zlno_t i=0; i < len; i++)
      wgts[i] = val;
  }
  else if (how == roundRobin){
    for (int i=0; i < 10; i++){
      float val = (i + 10)*scale;
      for (zlno_t j=i; j < len; j += 10)
         wgts[j] = val;
    }
  }
  else if (how == increasing){
    float val = scale + rank;
    for (zlno_t i=0; i < len; i++)
      wgts[i] = val;
  }
}

//////////////////////////////////////////////////////////////////////////////
/* Create a mesh of approximately the desired size.
 *
 *  We want 3 dimensions close to equal in length.
 */
tMVector_t* makeMeshCoordinates(
    const RCP<const Teuchos::Comm<int> > & comm,
    zgno_t numGlobalCoords)
{
  int rank = comm->getRank();
  int nprocs = comm->getSize();

  double k = log(numGlobalCoords) / 3;
  double xdimf = exp(k) + 0.5;
  ssize_t xdim = static_cast<ssize_t>(floor(xdimf));
  ssize_t ydim = xdim;
  ssize_t zdim = numGlobalCoords / (xdim*ydim);
  ssize_t num=xdim*ydim*zdim;
  ssize_t diff = numGlobalCoords - num;
  ssize_t newdiff = 0;

  while (diff > 0){
    if (zdim > xdim && zdim > ydim){
      zdim++;
      newdiff = diff - (xdim*ydim);
      if (newdiff < 0)
        if (diff < -newdiff)
          zdim--;
    }
    else if (ydim > xdim && ydim > zdim){
      ydim++;
      newdiff = diff - (xdim*zdim);
      if (newdiff < 0)
        if (diff < -newdiff)
          ydim--;
    }
    else{
      xdim++;
      newdiff = diff - (ydim*zdim);
      if (newdiff < 0)
        if (diff < -newdiff)
          xdim--;
    }

    diff = newdiff;
  }

  num=xdim*ydim*zdim;
  diff = numGlobalCoords - num;
  if (diff < 0)
    diff /= -numGlobalCoords;
  else
    diff /= numGlobalCoords;

  if (rank == 0){
    if (diff > .01)
      std::cout << "Warning: Difference " << diff*100 << " percent" << std::endl;
    std::cout << "Mesh size: " << xdim << "x" << ydim << "x" <<
      zdim << ", " << num << " vertices." << std::endl;
  }

  // Divide coordinates.

  ssize_t numLocalCoords = num / nprocs;
  ssize_t leftOver = num % nprocs;
  ssize_t gid0 = 0;

  if (rank <= leftOver)
    gid0 = rank * (numLocalCoords+1);
  else
    gid0 = (leftOver * (numLocalCoords+1)) +
           ((rank - leftOver) * numLocalCoords);

  if (rank < leftOver)
    numLocalCoords++;

  ssize_t gid1 = gid0 + numLocalCoords;

  zgno_t *ids = new zgno_t[numLocalCoords];
  if (!ids)
    throw bad_alloc();
  ArrayView<zgno_t> idArray(ids, numLocalCoords);
  zgno_t *idptr = ids;

  for (ssize_t i=gid0; i < gid1; i++)
    *idptr++ = zgno_t(i);

  RCP<const tMap_t> idMap = rcp(new tMap_t(num, idArray, 0, comm));

  delete [] ids;

  // Create a Tpetra::MultiVector of coordinates.

  zscalar_t *x = new zscalar_t [numLocalCoords*3];
  if (!x) throw bad_alloc();

  zscalar_t *y = x + numLocalCoords;
  zscalar_t *z = y + numLocalCoords;

  zgno_t xStart = 0;
  zgno_t yStart = 0;
  zgno_t xyPlane = xdim*ydim;
  zgno_t zStart = gid0 / xyPlane;
  zgno_t rem = gid0 % xyPlane;
  if (rem > 0){
    yStart = rem / xdim;
    xStart = rem % xdim;
  }

  zlno_t next = 0;
  for (zscalar_t zval=zStart; next < numLocalCoords && zval < zdim; zval+=1.){
    for (zscalar_t yval=yStart; next < numLocalCoords && yval < ydim; yval+=1.){
      for (zscalar_t xval=xStart; next < numLocalCoords && xval < xdim;xval+=1.){
        x[next] = xval;
        y[next] = yval;
        z[next] = zval;
        next++;
      }
      xStart = 0;
    }
    yStart = 0;
  }

  ArrayView<const zscalar_t> xArray(x, numLocalCoords*3);
  tMVector_t *dots = new tMVector_t(idMap, xArray, numLocalCoords, 3);

  delete [] x;
  return dots;
}


int main(int narg, char *arg[])
{
  // MEMORY_CHECK(true, "Before initializing MPI");

  Tpetra::ScopeGuard tscope(&narg, &arg);
  Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();
  int rank = comm->getRank();
  int nprocs = comm->getSize();
  DOTS dots;

  MEMORY_CHECK(rank==0 || rank==nprocs-1, "After initializing MPI");

  if (rank==0)
    std::cout << "Number of processes: " << nprocs << std::endl;

  // Default values
  zgno_t numGlobalCoords = 1000;
  int nWeights = 0;
  int debugLevel=2;
  string memoryOn("memoryOn");
  string memoryOff("memoryOff");
  bool doMemory=false;
  int numGlobalParts = nprocs;
  int dummyTimer=0;
  bool remap=0;

  string balanceCount("balance_object_count");
  string balanceWeight("balance_object_weight");
  string mcnorm1("multicriteria_minimize_total_weight");
  string mcnorm2("multicriteria_balance_total_maximum");
  string mcnorm3("multicriteria_minimize_maximum_weight");
  string objective(balanceWeight);   // default

  // Process command line input
  CommandLineProcessor commandLine(false, true);
  //commandLine.setOption("size", &numGlobalCoords,
  //  "Approximate number of global coordinates.");
  int input_option = 0;
  commandLine.setOption("input_option", &input_option,
    "whether to use mesh creation, geometric generator, or file input");
  string inputFile = "";

  commandLine.setOption("input_file", &inputFile,
    "the input file for geometric generator or file input");


  commandLine.setOption("size", &numGlobalCoords,
    "Approximate number of global coordinates.");
  commandLine.setOption("numParts", &numGlobalParts,
    "Number of parts (default is one per proc).");
  commandLine.setOption("nWeights", &nWeights,
    "Number of weights per coordinate, zero implies uniform weights.");
  commandLine.setOption("debug", &debugLevel, "Zoltan1 debug level");
  commandLine.setOption("remap", "no-remap", &remap,
    "Zoltan1 REMAP parameter; disabled by default for scalability testing");
  commandLine.setOption("timers", &dummyTimer, "ignored");
  commandLine.setOption(memoryOn.c_str(), memoryOff.c_str(), &doMemory,
    "do memory profiling");

  string doc(balanceCount);
  doc.append(": ignore weights\n");
  doc.append(balanceWeight);
  doc.append(": balance on first weight\n");
  doc.append(mcnorm1);
  doc.append(": given multiple weights, balance their total.\n");
  doc.append(mcnorm3);
  doc.append(": given multiple weights, "
             "balance the maximum for each coordinate.\n");
  doc.append(mcnorm2);
  doc.append(": given multiple weights, balance the L2 norm of the weights.\n");
  commandLine.setOption("objective", &objective,  doc.c_str());

  CommandLineProcessor::EParseCommandLineReturn rc =
    commandLine.parse(narg, arg);



  if (rc != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL) {
    if (rc == Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED) {
      if (rank==0) std::cout << "PASS" << std::endl;
      return 1;
    }
    else {
      if (rank==0) std::cout << "FAIL" << std::endl;
      return 0;
    }
  }

  //MEMORY_CHECK(doMemory && rank==0, "After processing parameters");

  // Create the data structure
  size_t numLocalCoords = 0;
  if (input_option == 0){
      dots.coordinates = makeMeshCoordinates(comm, numGlobalCoords);
      numLocalCoords = dots.coordinates->getLocalLength();

#if 0
      comm->barrier();
      for (int p=0; p < nprocs; p++){
          if (p==rank){
              std::cout << "Rank " << rank << ", " << numLocalCoords << "coords" << std::endl;
              const zscalar_t *x = coordinates->getData(0).getRawPtr();
              const zscalar_t *y = coordinates->getData(1).getRawPtr();
              const zscalar_t *z = coordinates->getData(2).getRawPtr();
              for (zlno_t i=0; i < numLocalCoords; i++)
                  std::cout << " " << x[i] << " " << y[i] << " " << z[i] << std::endl;
          }
          std::cout.flush();
          comm->barrier();
      }
#endif

      if (nWeights > 0){

          dots.weights.resize(nWeights);

          int wt = 0;
          float scale = 1.0;
          for (int i=0; i < nWeights; i++){
              dots.weights[i].resize(numLocalCoords);
              makeWeights(comm, dots.weights[i], weightTypes(wt++), scale, rank);

              if (wt == numWeightTypes){
                  wt = 0;
                  scale++;
              }
          }
      }
  }
  else if(input_option == 1){
      Teuchos::ParameterList geoparams("geo params");
      readGeoGenParams(inputFile, geoparams, comm);
      GeometricGen::GeometricGenerator<zscalar_t, zlno_t, zgno_t, znode_t> *gg = new GeometricGen::GeometricGenerator<zscalar_t, zlno_t, zgno_t, znode_t>(geoparams,comm);

      int coord_dim = gg->getCoordinateDimension();
      nWeights = gg->getNumWeights();
      numLocalCoords = gg->getNumLocalCoords();
      numGlobalCoords = gg->getNumGlobalCoords();
      zscalar_t **coords = new zscalar_t * [coord_dim];
      for(int i = 0; i < coord_dim; ++i){
          coords[i] = new zscalar_t[numLocalCoords];
      }
      gg->getLocalCoordinatesCopy(coords);
      zscalar_t **weight = NULL;
      if(nWeights){
          weight= new zscalar_t * [nWeights];
          for(int i = 0; i < nWeights; ++i){
              weight[i] = new zscalar_t[numLocalCoords];
          }
          gg->getLocalWeightsCopy(weight);
      }

      delete gg;

      RCP<Tpetra::Map<zlno_t, zgno_t, znode_t> > mp = rcp(
              new Tpetra::Map<zlno_t, zgno_t, znode_t> (numGlobalCoords, numLocalCoords, 0, comm));

      Teuchos::Array<Teuchos::ArrayView<const zscalar_t> > coordView(coord_dim);
      for (int i=0; i < coord_dim; i++){
          if(numLocalCoords > 0){
              Teuchos::ArrayView<const zscalar_t> a(coords[i], numLocalCoords);
              coordView[i] = a;
          } else{
              Teuchos::ArrayView<const zscalar_t> a;
              coordView[i] = a;
          }
      }

      tMVector_t *tmVector = new tMVector_t( mp, coordView.view(0, coord_dim), coord_dim);

      dots.coordinates = tmVector;
      dots.weights.resize(nWeights);

      if(nWeights){
          for (int i = 0; i < nWeights;++i){
              for (zlno_t j = 0; j < zlno_t(numLocalCoords); ++j){
                  dots.weights[i].push_back(weight[i][j]);
              }
          }
      }
      if(nWeights){
          for(int i = 0; i < nWeights; ++i)
              delete [] weight[i];
          delete [] weight;
      }
  }
  else {

      UserInputForTests uinput(testDataFilePath, inputFile, comm, true);
      RCP<tMVector_t> coords = uinput.getUICoordinates();
      tMVector_t *newMulti = new tMVector_t(*coords);
      dots.coordinates = newMulti;
      numLocalCoords = coords->getLocalLength();
      numGlobalCoords = coords->getGlobalLength();
  }

  MEMORY_CHECK(doMemory && rank==0, "After creating input");

  // Now call Zoltan to partition the problem.

  float ver;
  int aok = Zoltan_Initialize(narg, arg, &ver);

  if (aok != 0){
    printf("Zoltan_Initialize failed\n");
    exit(0);
  }

  struct Zoltan_Struct *zz;
  zz = Zoltan_Create(MPI_COMM_WORLD);

  Zoltan_Set_Param(zz, "LB_METHOD", "RCB");
  Zoltan_Set_Param(zz, "LB_APPROACH", "PARTITION");
  Zoltan_Set_Param(zz, "CHECK_GEOM", "0");
  Zoltan_Set_Param(zz, "NUM_GID_ENTRIES", "1");
  Zoltan_Set_Param(zz, "NUM_LID_ENTRIES", "0");
  Zoltan_Set_Param(zz, "RETURN_LISTS", "PART");
  std::ostringstream oss;
  oss << numGlobalParts;
  Zoltan_Set_Param(zz, "NUM_GLOBAL_PARTS", oss.str().c_str());
  oss.str("");
  oss << debugLevel;
  Zoltan_Set_Param(zz, "DEBUG_LEVEL", oss.str().c_str());

  if (remap)
    Zoltan_Set_Param(zz, "REMAP", "1");
  else
    Zoltan_Set_Param(zz, "REMAP", "0");

  if (objective != balanceCount){
    oss.str("");
    oss << nWeights;
    Zoltan_Set_Param(zz, "OBJ_WEIGHT_DIM", oss.str().c_str());

    if (objective == mcnorm1)
      Zoltan_Set_Param(zz, "RCB_MULTICRITERIA_NORM", "1");
    else if (objective == mcnorm2)
      Zoltan_Set_Param(zz, "RCB_MULTICRITERIA_NORM", "2");
    else if (objective == mcnorm3)
      Zoltan_Set_Param(zz, "RCB_MULTICRITERIA_NORM", "3");
  }
  else{
    Zoltan_Set_Param(zz, "OBJ_WEIGHT_DIM", "0");
  }

  Zoltan_Set_Num_Obj_Fn(zz, getNumObj, &dots);
  Zoltan_Set_Obj_List_Fn(zz, getObjList, &dots);
  Zoltan_Set_Num_Geom_Fn(zz, getDim, &dots);
  Zoltan_Set_Geom_Multi_Fn(zz, getCoords, &dots);

  int changes, numGidEntries, numLidEntries, numImport, numExport;
  ZOLTAN_ID_PTR importGlobalGids, importLocalGids;
  ZOLTAN_ID_PTR exportGlobalGids, exportLocalGids;
  int *importProcs, *importToPart, *exportProcs, *exportToPart;

  MEMORY_CHECK(doMemory && rank==0, "Before Zoltan_LB_Partition");

  if (rank == 0) std::cout << "Calling Zoltan_LB_Partition" << std::endl;
  aok = Zoltan_LB_Partition(zz, &changes, &numGidEntries, &numLidEntries,
                            &numImport, &importGlobalGids, &importLocalGids,
                            &importProcs, &importToPart,
                            &numExport, &exportGlobalGids, &exportLocalGids,
                            &exportProcs, &exportToPart);
  if (rank == 0) std::cout << "Returned from Zoltan_LB_Partition" << std::endl;

  MEMORY_CHECK(doMemory && rank==0, "After Zoltan_LB_Partition");

  /* Print the load-balance stats here */

  float *sumWgtPerPart = new float[numGlobalParts];
  float *gsumWgtPerPart = new float[numGlobalParts];
  for (int i = 0; i < numGlobalParts; i++) sumWgtPerPart[i] = 0.;

  for (size_t i = 0; i < numLocalCoords; i++)
    sumWgtPerPart[exportToPart[i]] += (nWeights ? dots.weights[0][i]: 1.);

  Teuchos::reduceAll<int, float>(*comm, Teuchos::REDUCE_SUM, numGlobalParts,
                                  sumWgtPerPart, gsumWgtPerPart);

  float maxSumWgtPerPart = 0.;
  float minSumWgtPerPart = std::numeric_limits<float>::max();
  float totWgt = 0.;
  int maxSumWgtPart=0, minSumWgtPart=0;
  for (int i = 0; i < numGlobalParts; i++) {
    if (gsumWgtPerPart[i] > maxSumWgtPerPart) {
      maxSumWgtPerPart = gsumWgtPerPart[i];
      maxSumWgtPart = i;
    }
    if (gsumWgtPerPart[i] < minSumWgtPerPart) {
      minSumWgtPerPart = gsumWgtPerPart[i];
      minSumWgtPart = i;
    }
    totWgt += gsumWgtPerPart[i];
  }

  if (rank == 0)
    std::cout << std::endl << std::endl
              << "Part loads (per part for " << numGlobalParts << " parts):"
              << std::endl
              << "   min = " << minSumWgtPerPart
                             << " in part " << minSumWgtPart << std::endl
              << "   max = " << maxSumWgtPerPart
                             << " in part " << maxSumWgtPart << std::endl
              << "   tot = " << totWgt << std::endl
              << "   avg = " << totWgt / numGlobalParts
              << std::endl << std::endl << std::endl;

  delete [] sumWgtPerPart;
  delete [] gsumWgtPerPart;

  Zoltan_Destroy(&zz);
  MEMORY_CHECK(doMemory && rank==0, "After Zoltan_Destroy");

  delete dots.coordinates;
  for (int i = 0; i < nWeights; i++)
    dots.weights[i].clear();
  dots.weights.clear();

  MEMORY_CHECK(doMemory && rank==0, "After destroying input");

  if (rank==0){
    if (aok != 0)
      std::cout << "FAIL" << std::endl;
    else
      std::cout << "PASS" << std::endl;
  }

  return 0;
}
