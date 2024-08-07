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
    \todo  get the imbalance when done
*/

#include <Zoltan2_TestHelpers.hpp>

#include <zoltan.h>
#include <Teuchos_CommandLineProcessor.hpp>

#include <vector>
#include <ostream>
#include <sstream>
#include <string>

#include <Zoltan2_TestHelpers.hpp>
#include <Zoltan2_XpetraMultiVectorAdapter.hpp>
#include <Zoltan2_PartitioningSolution.hpp>
#include <Zoltan2_PartitioningProblem.hpp>
#include <GeometricGenerator.hpp>
#include <Zoltan2_EvaluatePartition.hpp>

using std::vector;
using std::bad_alloc;
using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::Comm;
using Teuchos::ArrayView;
using Teuchos::ArrayRCP;
using Teuchos::CommandLineProcessor;

typedef Tpetra::MultiVector<zscalar_t, zlno_t, zgno_t, znode_t> tMVector_t;
typedef Tpetra::Map<zlno_t, zgno_t, znode_t> tMap_t;

static ArrayRCP<ArrayRCP<zscalar_t> > weights;
static RCP<tMVector_t> coordinates;


const char param_comment = '#';

std::string trim_right_copy(
    const std::string& s,
    const std::string& delimiters = " \f\n\r\t\v" )
{
  return s.substr( 0, s.find_last_not_of( delimiters ) + 1 );
}

std::string trim_left_copy(
    const std::string& s,
    const std::string& delimiters = " \f\n\r\t\v" )
{
  return s.substr( s.find_first_not_of( delimiters ) );
}

std::string trim_copy(
    const std::string& s,
    const std::string& delimiters = " \f\n\r\t\v" )
{
  return trim_left_copy( trim_right_copy( s, delimiters ), delimiters );
}

// Zoltan1 query functions
bool getArgumentValue(std::string &argumentid, double &argumentValue, std::string argumentline){
  std::stringstream stream(std::stringstream::in | std::stringstream::out);
  stream << argumentline;
  getline(stream, argumentid, '=');
  if (stream.eof()){
    return false;
  }
  stream >> argumentValue;
  return true;
}

std::string convert_to_string(char *args){
  std::string tmp = "";
  for(int i = 0; args[i] != 0; i++)
    tmp += args[i];
  return tmp;
}
int getNumObj(void *data, int *ierr)
{
  *ierr = 0;

  return coordinates->getLocalLength();
}

int getDim(void *data, int *ierr)
{
  *ierr = 0;
  return 3;
}

void getObjList(void *data, int numGid, int numLid,
  zgno_t * gids, zgno_t * lids, 
  int num_wgts, float *obj_wgts, int *ierr)
{
  *ierr = 0;
  size_t localLen = coordinates->getLocalLength();
  const zgno_t *ids = coordinates->getMap()->getLocalElementList().getRawPtr();
  zgno_t *idsNonConst = const_cast<zgno_t *>(ids);

  if (sizeof(zgno_t) == sizeof(zgno_t)){
    memcpy(gids, idsNonConst, sizeof(zgno_t) * localLen);
  }
  else{
    for (size_t i=0; i < localLen; i++)
      gids[i] = static_cast<zgno_t>(idsNonConst[i]);
  }

  if (num_wgts > 0){
    float *wgts = obj_wgts;
    for (size_t i=0; i < localLen; i++)
      for (int w=0; w < num_wgts; w++)
        *wgts++ = static_cast<float>(weights[w][i]);
  }
}

void getCoords(void *data, int numGid, int numLid,
  int numObj, zgno_t * gids, zgno_t * lids,
  int dim, double *coords, int *ierr)
{
  // I know that Zoltan asks for coordinates in gid order.
  *ierr = 0;
  double *val = coords;
  const zscalar_t *x = coordinates->getData(0).getRawPtr();
  const zscalar_t *y = coordinates->getData(1).getRawPtr();
  const zscalar_t *z = coordinates->getData(2).getRawPtr();
  for (zlno_t i=0; i < numObj; i++){
    *val++ = static_cast<double>(x[i]);
    *val++ = static_cast<double>(y[i]);
    *val++ = static_cast<double>(z[i]);
  }
}



enum weightTypes{
  upDown,
  roundRobin,
  increasing,
  numWeightTypes
};

ArrayRCP<zscalar_t> makeWeights(
  const RCP<const Teuchos::Comm<int> > & comm,
  zlno_t len, weightTypes how, zscalar_t scale, int rank)
{
  zscalar_t *wgts = new zscalar_t [len];
  if (!wgts)
    throw bad_alloc();

  ArrayRCP<zscalar_t> wgtArray(wgts, 0, len, true);

  if (how == upDown){
    zscalar_t val = scale + rank%2;
    for (zlno_t i=0; i < len; i++)
      wgts[i] = val;
  }
  else if (how == roundRobin){
    for (int i=0; i < 10; i++){
      zscalar_t val = (i + 10)*scale;
      for (int j=i; j < len; j += 10)
         wgts[j] = val;
    }
  }
  else if (how == increasing){
    zscalar_t val = scale + rank;
    for (zlno_t i=0; i < len; i++)
      wgts[i] = val;
  }

  return wgtArray;
}

/*! \brief Create a mesh of approximately the desired size.
 *
 *  We want 3 dimensions close to equal in length.
 */
const RCP<tMVector_t> getMeshCoordinates(
    const RCP<const Teuchos::Comm<int> > & comm,
    zgno_t numGlobalCoords)
{
  int rank = comm->getRank();
  int nprocs = comm->getSize();

  double k = log(numGlobalCoords) / 3;
  double xdimf = exp(k) + 0.5;
  zgno_t xdim = static_cast<int>(floor(xdimf));
  zgno_t ydim = xdim;
  zgno_t zdim = numGlobalCoords / (xdim*ydim);
  zgno_t num=xdim*ydim*zdim;
  zgno_t diff = numGlobalCoords - num;
  zgno_t newdiff = 0;

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

  zgno_t numLocalCoords = num / nprocs;
  zgno_t leftOver = num % nprocs;
  zgno_t gid0 = 0;

  if (rank <= leftOver)
    gid0 = zgno_t(rank) * (numLocalCoords+1);
  else
    gid0 = (leftOver * (numLocalCoords+1)) + 
           ((zgno_t(rank) - leftOver) * numLocalCoords);

  if (rank < leftOver)
    numLocalCoords++;

  zgno_t gid1 = gid0 + numLocalCoords;

  zgno_t *ids = new zgno_t [numLocalCoords];
  if (!ids)
    throw bad_alloc();
  ArrayRCP<zgno_t> idArray(ids, 0, numLocalCoords, true);

  for (zgno_t i=gid0; i < gid1; i++)
    *ids++ = i;   

  RCP<const tMap_t> idMap = rcp(
    new tMap_t(num, idArray.view(0, numLocalCoords), 0, comm));

  // Create a Tpetra::MultiVector of coordinates.

  zscalar_t *x = new zscalar_t [numLocalCoords*3]; 
  if (!x)
    throw bad_alloc();
  ArrayRCP<zscalar_t> coordArray(x, 0, numLocalCoords*3, true);

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
  for (zscalar_t zval=zStart; next < numLocalCoords && zval < zdim; zval++){
    for (zscalar_t yval=yStart; next < numLocalCoords && yval < ydim; yval++){
      for (zscalar_t xval=xStart; next < numLocalCoords && xval < xdim; xval++){
        x[next] = xval;
        y[next] = yval;
        z[next] = zval;
        next++;
      }
      xStart = 0;
    }
    yStart = 0;
  }

  ArrayView<const zscalar_t> xArray(x, numLocalCoords);
  ArrayView<const zscalar_t> yArray(y, numLocalCoords);
  ArrayView<const zscalar_t> zArray(z, numLocalCoords);
  ArrayRCP<ArrayView<const zscalar_t> > coordinates =
    arcp(new ArrayView<const zscalar_t> [3], 0, 3);
  coordinates[0] = xArray;
  coordinates[1] = yArray;
  coordinates[2] = zArray;

  ArrayRCP<const ArrayView<const zscalar_t> > constCoords =
   coordinates.getConst();

  RCP<tMVector_t> meshCoords = rcp(new tMVector_t(
    idMap, constCoords.view(0,3), 3));

  return meshCoords;
}


void getArgVals(int narg, char **arg,   int &numParts,
    std::string &paramFile){

  for(int i = 0; i < narg; ++i){
    std::string tmp = convert_to_string(arg[i]);
    std::string identifier = "";
    long long int value = -1; double fval = -1;
    if(!getArgumentValue(identifier, fval, tmp)) continue;
    value = (long long int) (fval);
    if(identifier == "C"){
      if(value > 0){
        numParts=value;
      } else {
        throw  "Invalid argument at " + tmp;
      }
    }
    else if(identifier == "PF"){
      std::stringstream stream(std::stringstream::in | std::stringstream::out);
      stream << tmp;
      getline(stream, paramFile, '=');
      stream >> paramFile;
    }
    else {
      throw "Invalid argument at " + tmp;
    }

  }

}

void readGeoGenParams(std::string paramFileName, Teuchos::ParameterList &geoparams, const RCP<const Teuchos::Comm<int> > & comm){
  std::string input = "";
  char inp[25000];
  for(int i = 0; i < 25000; ++i){
    inp[i] = 0;
  }
  if(comm->getRank() == 0){
    std::fstream inParam(paramFileName.c_str());

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


  int size = input.size();

  //MPI_Bcast(&size,1,MPI_INT, 0, MPI_COMM_WORLD);
  //MPI_Bcast(inp,size,MPI_CHAR, 0, MPI_COMM_WORLD);
  //Teuchos::broadcast<int, char>(comm, 0,inp);

  comm->broadcast(0, sizeof(int), (char*) &size);
  comm->broadcast(0, size, inp);
  //Teuchos::broadcast<int,std::string>(comm,0, &input);
  std::istringstream inParam(inp);
  std::string str;
  getline (inParam,str);
  while (!inParam.eof()){
    if(str[0] != param_comment){
      size_t pos = str.find('=');
      if(pos == std::string::npos){
        throw  "Invalid Line:" + str  + " in parameter file";
      }
      std::string paramname = trim_copy(str.substr(0,pos));
      std::string paramvalue = trim_copy(str.substr(pos + 1));
      geoparams.set(paramname, paramvalue);
    }
    getline (inParam,str);
  }
}

int main(int narg, char *arg[])
{
  // MEMORY_CHECK(true, "Before initializing MPI");
  Tpetra::ScopeGuard tscope(&narg, &arg);
  Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();

  int rank = comm->getRank();
  int nprocs = comm->getSize();

  MEMORY_CHECK(rank==0 || rank==nprocs-1, "After initializing MPI");

  if (rank==0)
    std::cout << "Number of processes: " << nprocs << std::endl;

  // Default values
  //double numGlobalCoords = 1000;
  //int nWeights = 0;
  int debugLevel=2;          // for timing
  std::string memoryOn("memoryOn");
  std::string memoryOff("memoryOff");
  bool doMemory=false;

  int numGlobalParts = 512 * 512;
  std::string paramFile = "p2.txt";
  //getArgVals(narg, arg,   numGlobalParts, paramFile);

  int dummyTimer=0;



  Teuchos::ParameterList geoparams("geo params");

  readGeoGenParams(paramFile, geoparams, comm);
#ifdef HAVE_ZOLTAN2_OMP
  double begin = omp_get_wtime();
#endif
  GeometricGenerator<zscalar_t, zlno_t, zgno_t, znode_t> *gg = new GeometricGenerator<zscalar_t, zlno_t, zgno_t, znode_t>(geoparams,comm);
#ifdef HAVE_ZOLTAN2_OMP
  double end = omp_get_wtime();
  std::cout << "GeometricGen Time:" << end - begin << std::endl;
#endif
  int coord_dim = gg->getCoordinateDimension();
  int nWeights = gg->getNumWeights();
  zlno_t numLocalPoints = gg->getNumLocalCoords();
  zgno_t numGlobalPoints = gg->getNumGlobalCoords();
  zscalar_t **coords = new zscalar_t * [coord_dim];
  for(int i = 0; i < coord_dim; ++i){
    coords[i] = new zscalar_t[numLocalPoints];
  }
  gg->getLocalCoordinatesCopy(coords);
  zscalar_t **weight = NULL;
  if(nWeights){
    weight= new zscalar_t * [nWeights];
    for(int i = 0; i < nWeights; ++i){
      weight[i] = new zscalar_t[numLocalPoints];
    }
    gg->getLocalWeightsCopy(weight);
  }

  zgno_t globalSize = static_cast<zgno_t>(numGlobalPoints);
  delete gg;

  std::cout << "coord_dim:" << coord_dim << " nWeights:" << nWeights << " numLocalPoints:" << numLocalPoints << " numGlobalPoints:" << numGlobalPoints << std::endl;

  CommandLineProcessor commandLine(false, true);
  commandLine.setOption("size", &numGlobalPoints,
    "Approximate number of global coordinates.");
  commandLine.setOption("numParts", &numGlobalParts, 
    "Number of parts (default is one per proc).");
  commandLine.setOption("numWeights", &nWeights,
    "Number of weights per coordinate, zero implies uniform weights.");
  commandLine.setOption("debug", &debugLevel, "Zoltan1 debug level");
  commandLine.setOption("timers", &dummyTimer, "ignored");
  commandLine.setOption(memoryOn.c_str(), memoryOff.c_str(), &doMemory,
    "do memory profiling");

  std::string balanceCount("balance_object_count");
  std::string balanceWeight("balance_object_weight");
  std::string mcnorm1("multicriteria_minimize_total_weight");
  std::string mcnorm2("multicriteria_balance_total_maximum");
  std::string mcnorm3("multicriteria_minimize_maximum_weight");

  std::string objective(balanceWeight);   // default

  std::string doc(balanceCount);
  doc.append(": ignore weights\n");

  doc.append(balanceWeight);
  doc.append(": balance on first weight\n");

  doc.append(mcnorm1);
  doc.append(": given multiple weights, balance their total.\n");

  doc.append(mcnorm3);
  doc.append(": given multiple weights, balance the maximum for each coordinate.\n");

  doc.append(mcnorm2);
  doc.append(": given multiple weights, balance the L2 norm of the weights.\n");

  commandLine.setOption("objective", &objective,  doc.c_str());

  CommandLineProcessor::EParseCommandLineReturn rc = 
    commandLine.parse(narg, arg);

  if (rc != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL){
    if (rc == Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED){
      if (rank==0)
        std::cout << "PASS" << std::endl;
      return 1;
    }
    else{
      if (rank==0)
        std::cout << "FAIL" << std::endl;
      return 0;
    }
  }

  //MEMORY_CHECK(doMemory && rank==0, "After processing parameters");




  RCP<Tpetra::Map<zlno_t, zgno_t, znode_t> > mp = rcp(
      new Tpetra::Map<zlno_t, zgno_t, znode_t> (numGlobalPoints, numLocalPoints, 0, comm));

  Teuchos::Array<Teuchos::ArrayView<const zscalar_t> > coordView(coord_dim);
  for (int i=0; i < coord_dim; i++){
    if(numLocalPoints > 0){
      Teuchos::ArrayView<const zscalar_t> a(coords[i], numLocalPoints);
      coordView[i] = a;
    } else{
      Teuchos::ArrayView<const zscalar_t> a;
      coordView[i] = a;
    }
  }

  coordinates = RCP< Tpetra::MultiVector<zscalar_t, zlno_t, zgno_t, znode_t> >(
      new Tpetra::MultiVector<zscalar_t, zlno_t, zgno_t, znode_t>( mp, coordView.view(0, coord_dim), coord_dim));

  //typedef Tpetra::MultiVector<zscalar_t, zlno_t, zgno_t, znode_t> tMVector_t;

  RCP<const tMVector_t> coordsConst = Teuchos::rcp_const_cast<const tMVector_t>(coordinates);



  //coordinates = getMeshCoordinates(comm, globalSize);
  size_t numLocalCoords = coordinates->getLocalLength();

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

    weights = arcp(new ArrayRCP<zscalar_t> [nWeights],
      0, nWeights, true);
    for(int i = 0; i < nWeights; ++i){
      //weights[i] = ArrayRCP<zscalar_t>(weight[i]);
    }
  }

  MEMORY_CHECK(doMemory && rank==0, "After creating input");

  // Now call Zoltan to partition the problem.

  float ver;
  int aok = Zoltan_Initialize(narg, arg, &ver);

  if (aok != 0){
    printf("sorry...\n");
    exit(0);
  }

  struct Zoltan_Struct *zz;
  zz = Zoltan_Create(MPI_COMM_WORLD);
  
  Zoltan_Set_Param(zz, "LB_METHOD", "RCB");
  Zoltan_Set_Param(zz, "LB_APPROACH", "PARTITION");
  Zoltan_Set_Param(zz, "CHECK_GEOM", "0");
  Zoltan_Set_Param(zz, "NUM_GID_ENTRIES", "1"); // compiled with ULONG option
  Zoltan_Set_Param(zz, "NUM_LID_ENTRIES", "0");
  Zoltan_Set_Param(zz, "RETURN_LISTS", "ALL");

  //Zoltan_Set_Param(zz, "FINAL_OUTPUT", "1");
  Zoltan_Set_Param(zz, "IMBALANCE_TOL", "1.0");
  std::ostringstream oss;
  oss << numGlobalParts;
  Zoltan_Set_Param(zz, "NUM_GLOBAL_PARTS", oss.str().c_str());
  oss.str("");
  oss << debugLevel;
  Zoltan_Set_Param(zz, "DEBUG_LEVEL", oss.str().c_str());

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

  Zoltan_Set_Num_Obj_Fn(zz, getNumObj, NULL);
  Zoltan_Set_Obj_List_Fn(zz, getObjList,NULL);
  Zoltan_Set_Num_Geom_Fn(zz, getDim, NULL);
  Zoltan_Set_Geom_Multi_Fn(zz, getCoords, NULL);

  int changes, numGidEntries, numLidEntries, numImport, numExport;
  zgno_t * importGlobalGids, * importLocalGids;
  zgno_t * exportGlobalGids, * exportLocalGids;
  int *importProcs, *importToPart, *exportProcs, *exportToPart;

  MEMORY_CHECK(doMemory && rank==0, "Before Zoltan_LB_Partition");


  aok = Zoltan_LB_Partition(zz, /* input (all remaining fields are output) */
        &changes,        /* 1 if partitioning was changed, 0 otherwise */
        &numGidEntries,  /* Number of integers used for a global ID */
        &numLidEntries,  /* Number of integers used for a local ID */
        &numImport,      /* Number of vertices to be sent to me */
        &importGlobalGids,  /* Global IDs of vertices to be sent to me */
        &importLocalGids,   /* Local IDs of vertices to be sent to me */
        &importProcs,    /* Process rank for source of each incoming vertex */
        &importToPart,   /* New partition for each incoming vertex */
        &numExport,      /* Number of vertices I must send to other processes*/
        &exportGlobalGids,  /* Global IDs of the vertices I must send */
        &exportLocalGids,   /* Local IDs of the vertices I must send */
        &exportProcs,    /* Process to which I send each of the vertices */
        &exportToPart);  /* Partition to which each vertex will belong */

  MEMORY_CHECK(doMemory && rank==0, "After Zoltan_LB_Partition");
  Zoltan_LB_Free_Part(importGlobalGids, importLocalGids, importProcs, importToPart);
  Zoltan_LB_Free_Part(exportGlobalGids, exportLocalGids, exportProcs, exportToPart);
  Zoltan_Destroy(&zz);
  MEMORY_CHECK(doMemory && rank==0, "After Zoltan_Destroy");

  if (rank==0){
    if (aok != 0)
      std::cout << "FAIL" << std::endl;
    else
      std::cout << "PASS" << std::endl;
  }

  return 0;
}
