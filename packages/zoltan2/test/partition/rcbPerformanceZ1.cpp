// @HEADER
//
// ***********************************************************************
//
//   Zoltan2: A package of combinatorial algorithms for scientific computing
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
// Questions? Contact Karen Devine      (kddevin@sandia.gov)
//                    Erik Boman        (egboman@sandia.gov)
//                    Siva Rajamanickam (srajama@sandia.gov)
//
// ***********************************************************************
//
// @HEADER

/*! \file rcbPerformanceZ1.cpp
    \brief rcbPerformance with Zoltan1

    Geometry is a uniform mesh.
    \todo  get the imbalance when done
*/

#include <Zoltan2_TestHelpers.hpp>

#define HAVE_ZOLTAN2_ZOLTAN
#ifdef HAVE_ZOLTAN2_ZOLTAN
#include <zoltan.h>
#include <Teuchos_CommandLineProcessor.hpp>

#include <vector>
#include <ostream>
#include <sstream>
#include <string>

#include <Zoltan2_TestHelpers.hpp>
#include <Zoltan2_BasicCoordinateInput.hpp>
#include <Zoltan2_XpetraMultiVectorInput.hpp>
#include <Zoltan2_PartitioningSolution.hpp>
#include <Zoltan2_PartitioningProblem.hpp>
#include <GeometricGenerator.hpp>
#include <Zoltan2_PartitioningSolutionQuality.hpp>
using namespace std;
using std::vector;
using std::cout;
using std::cerr;
using std::endl;
using std::bad_alloc;
using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::Comm;
using Teuchos::ArrayView;
using Teuchos::ArrayRCP;
using Teuchos::CommandLineProcessor;

typedef Tpetra::MultiVector<scalar_t, lno_t, gno_t, node_t> tMVector_t;
typedef Tpetra::Map<lno_t, gno_t, node_t> tMap_t;

static ArrayRCP<ArrayRCP<scalar_t> > weights;
static RCP<tMVector_t> coordinates;


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

// Zoltan1 query functions
bool getArgumentValue(string &argumentid, double &argumentValue, string argumentline){
  stringstream stream(stringstream::in | stringstream::out);
  stream << argumentline;
  getline(stream, argumentid, '=');
  if (stream.eof()){
    return false;
  }
  stream >> argumentValue;
  return true;
}

string convert_to_string(char *args){
  string tmp = "";
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
  gid_t * gids, gid_t * lids, 
  int wgt_dim, float *obj_wgts, int *ierr)
{
  *ierr = 0;
  size_t localLen = coordinates->getLocalLength();
  const gno_t *ids = coordinates->getMap()->getNodeElementList().getRawPtr();
  gno_t *idsNonConst = const_cast<gno_t *>(ids);

  if (sizeof(gid_t) == sizeof(gno_t)){
    memcpy(gids, idsNonConst, sizeof(gid_t) * localLen);
  }
  else{
    for (size_t i=0; i < localLen; i++)
      gids[i] = static_cast<gid_t>(idsNonConst[i]);
  }

  if (wgt_dim > 0){
    float *wgts = obj_wgts;
    for (size_t i=0; i < localLen; i++)
      for (int w=0; w < wgt_dim; w++)
        *wgts++ = static_cast<float>(weights[w][i]);
  }
}

void getCoordinates(void *data, int numGid, int numLid,
  int numObj, gid_t * gids, gid_t * lids,
  int dim, double *coords, int *ierr)
{
  // I know that Zoltan asks for coordinates in gid order.
  *ierr = 0;
  double *val = coords;
  const scalar_t *x = coordinates->getData(0).getRawPtr();
  const scalar_t *y = coordinates->getData(1).getRawPtr();
  const scalar_t *z = coordinates->getData(2).getRawPtr();
  for (lno_t i=0; i < numObj; i++){
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

ArrayRCP<scalar_t> makeWeights(
  const RCP<const Teuchos::Comm<int> > & comm,
  lno_t len, weightTypes how, scalar_t scale, int rank)
{
  scalar_t *wgts = new scalar_t [len];
  if (!wgts)
    throw bad_alloc();

  ArrayRCP<scalar_t> wgtArray(wgts, 0, len, true);

  if (how == upDown){
    scalar_t val = scale + rank%2;
    for (lno_t i=0; i < len; i++)
      wgts[i] = val;
  }
  else if (how == roundRobin){
    for (int i=0; i < 10; i++){
      scalar_t val = (i + 10)*scale;
      for (int j=i; j < len; j += 10)
         wgts[j] = val;
    }
  }
  else if (how == increasing){
    scalar_t val = scale + rank;
    for (lno_t i=0; i < len; i++)
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
    gno_t numGlobalCoords)
{
  int rank = comm->getRank();
  int nprocs = comm->getSize();

  double k = log(numGlobalCoords) / 3;
  double xdimf = exp(k) + 0.5;
  gno_t xdim = static_cast<int>(floor(xdimf));
  gno_t ydim = xdim;
  gno_t zdim = numGlobalCoords / (xdim*ydim);
  gno_t num=xdim*ydim*zdim;
  gno_t diff = numGlobalCoords - num;
  gno_t newdiff = 0;

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
      cout << "Warning: Difference " << diff*100 << " percent" << endl;
    cout << "Mesh size: " << xdim << "x" << ydim << "x" <<
      zdim << ", " << num << " vertices." << endl;
  }

  // Divide coordinates.

  gno_t numLocalCoords = num / nprocs;
  gno_t leftOver = num % nprocs;
  gno_t gid0 = 0;

  if (rank <= leftOver)
    gid0 = gno_t(rank) * (numLocalCoords+1);
  else
    gid0 = (leftOver * (numLocalCoords+1)) + 
           ((gno_t(rank) - leftOver) * numLocalCoords);

  if (rank < leftOver)
    numLocalCoords++;

  gno_t gid1 = gid0 + numLocalCoords;

  gno_t *ids = new gno_t [numLocalCoords];
  if (!ids)
    throw bad_alloc();
  ArrayRCP<gno_t> idArray(ids, 0, numLocalCoords, true);

  for (gno_t i=gid0; i < gid1; i++)
    *ids++ = i;   

  RCP<const tMap_t> idMap = rcp(
    new tMap_t(num, idArray.view(0, numLocalCoords), 0, comm));

  // Create a Tpetra::MultiVector of coordinates.

  scalar_t *x = new scalar_t [numLocalCoords*3]; 
  if (!x)
    throw bad_alloc();
  ArrayRCP<scalar_t> coordArray(x, 0, numLocalCoords*3, true);

  scalar_t *y = x + numLocalCoords;
  scalar_t *z = y + numLocalCoords;

  gno_t xStart = 0;
  gno_t yStart = 0;
  gno_t xyPlane = xdim*ydim;
  gno_t zStart = gid0 / xyPlane;
  gno_t rem = gid0 % xyPlane;
  if (rem > 0){
    yStart = rem / xdim;
    xStart = rem % xdim;
  }

  lno_t next = 0;
  for (scalar_t zval=zStart; next < numLocalCoords && zval < zdim; zval++){
    for (scalar_t yval=yStart; next < numLocalCoords && yval < ydim; yval++){
      for (scalar_t xval=xStart; next < numLocalCoords && xval < xdim; xval++){
        x[next] = xval;
        y[next] = yval;
        z[next] = zval;
        next++;
      }
      xStart = 0;
    }
    yStart = 0;
  }

  ArrayView<const scalar_t> xArray(x, numLocalCoords);
  ArrayView<const scalar_t> yArray(y, numLocalCoords);
  ArrayView<const scalar_t> zArray(z, numLocalCoords);
  ArrayRCP<ArrayView<const scalar_t> > coordinates =
    arcp(new ArrayView<const scalar_t> [3], 0, 3);
  coordinates[0] = xArray;
  coordinates[1] = yArray;
  coordinates[2] = zArray;

  ArrayRCP<const ArrayView<const scalar_t> > constCoords =
   coordinates.getConst();

  RCP<tMVector_t> meshCoords = rcp(new tMVector_t(
    idMap, constCoords.view(0,3), 3));

  return meshCoords;
}


void getArgVals(int argc, char **argv,   partId_t &numParts,
    std::string &paramFile){

  for(int i = 0; i < argc; ++i){
    string tmp = convert_to_string(argv[i]);
    string identifier = "";
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
      stringstream stream(stringstream::in | stringstream::out);
      stream << tmp;
      getline(stream, paramFile, '=');
      stream >> paramFile;
    }
    else {
      throw "Invalid argument at " + tmp;
    }

  }

}

void readGeoGenParams(string paramFileName, Teuchos::ParameterList &geoparams, const RCP<const Teuchos::Comm<int> > & comm){
  std::string input = "";
  char inp[25000];
  for(int i = 0; i < 25000; ++i){
    inp[i] = 0;
  }
  if(comm->getRank() == 0){
    fstream inParam(paramFileName.c_str());

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
  //Teuchos::broadcast<int,string>(comm,0, &input);
  istringstream inParam(inp);
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

int main(int argc, char *argv[])
{
  // MEMORY_CHECK(true, "Before initializing MPI");

  Teuchos::GlobalMPISession session(&argc, &argv, NULL);
  RCP<const Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();
  int rank = comm->getRank();
  int nprocs = comm->getSize();

  MEMORY_CHECK(rank==0 || rank==nprocs-1, "After initializing MPI");

  if (rank==0)
    cout << "Number of processes: " << nprocs << endl;

  // Default values
  //double numGlobalCoords = 1000;
  //int weightDim = 0;
  int debugLevel=2;          // for timing
  string memoryOn("memoryOn");
  string memoryOff("memoryOff");
  bool doMemory=false;

  int numGlobalParts = 512 * 512;
  string paramFile = "p2.txt";
  //getArgVals(argc, argv,   numGlobalParts, paramFile);

  int dummyTimer=0;



  Teuchos::ParameterList geoparams("geo params");

  readGeoGenParams(paramFile, geoparams, comm);
#ifdef HAVE_ZOLTAN2_OMP
  double begin = omp_get_wtime();
#endif
  GeometricGenerator<scalar_t, lno_t, gno_t, node_t> *gg = new GeometricGenerator<scalar_t, lno_t, gno_t, node_t>(geoparams,comm);
#ifdef HAVE_ZOLTAN2_OMP
  double end = omp_get_wtime();
  cout << "GeometricGen Time:" << end - begin << endl;
#endif
  int coord_dim = gg->getCoordinateDimension();
  int weight_dim = gg->getWeightDimension();
  lno_t numLocalPoints = gg->getNumLocalCoords();
  gno_t numGlobalPoints = gg->getNumGlobalCoords();
  scalar_t **coords = new scalar_t * [coord_dim];
  for(int i = 0; i < coord_dim; ++i){
    coords[i] = new scalar_t[numLocalPoints];
  }
  gg->getLocalCoordinatesCopy(coords);
  scalar_t **weight = NULL;
  if(weight_dim){
    weight= new scalar_t * [weight_dim];
    for(int i = 0; i < weight_dim; ++i){
      weight[i] = new scalar_t[numLocalPoints];
    }
    gg->getLocalWeightsCopy(weight);
  }

  gno_t globalSize = static_cast<gno_t>(numGlobalPoints);
  delete gg;

  cout << "coord_dim:" << coord_dim << " weight_dim:" << weight_dim << " numLocalPoints:" << numLocalPoints << " numGlobalPoints:" << numGlobalPoints << endl;

  CommandLineProcessor commandLine(false, true);
  commandLine.setOption("size", &numGlobalPoints,
    "Approximate number of global coordinates.");
  commandLine.setOption("numParts", &numGlobalParts, 
    "Number of parts (default is one per proc).");
  commandLine.setOption("weightDim", &weight_dim,
    "Number of weights per coordinate, zero implies uniform weights.");
  commandLine.setOption("debug", &debugLevel, "Zoltan1 debug level");
  commandLine.setOption("timers", &dummyTimer, "ignored");
  commandLine.setOption(memoryOn.c_str(), memoryOff.c_str(), &doMemory,
    "do memory profiling");

  string balanceCount("balance_object_count");
  string balanceWeight("balance_object_weight");
  string mcnorm1("multicriteria_minimize_total_weight");
  string mcnorm2("multicriteria_balance_total_maximum");
  string mcnorm3("multicriteria_minimize_maximum_weight");

  string objective(balanceWeight);   // default

  string doc(balanceCount);
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
    commandLine.parse(argc, argv);

  if (rc != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL){
    if (rc == Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED){
      if (rank==0)
        cout << "PASS" << endl;
      return 1;
    }
    else{
      if (rank==0)
        cout << "FAIL" << endl;
      return 0;
    }
  }

  //MEMORY_CHECK(doMemory && rank==0, "After processing parameters");




  RCP<Tpetra::Map<lno_t, gno_t, node_t> > mp = rcp(
      new Tpetra::Map<lno_t, gno_t, node_t> (numGlobalPoints, numLocalPoints, 0, comm));

  Teuchos::Array<Teuchos::ArrayView<const scalar_t> > coordView(coord_dim);
  for (int i=0; i < coord_dim; i++){
    if(numLocalPoints > 0){
      Teuchos::ArrayView<const scalar_t> a(coords[i], numLocalPoints);
      coordView[i] = a;
    } else{
      Teuchos::ArrayView<const scalar_t> a;
      coordView[i] = a;
    }
  }

  coordinates = RCP< Tpetra::MultiVector<scalar_t, lno_t, gno_t, node_t> >(
      new Tpetra::MultiVector<scalar_t, lno_t, gno_t, node_t>( mp, coordView.view(0, coord_dim), coord_dim));

  //typedef Tpetra::MultiVector<scalar_t, lno_t, gno_t, node_t> tMVector_t;

  RCP<const tMVector_t> coordsConst = Teuchos::rcp_const_cast<const tMVector_t>(coordinates);



  //coordinates = getMeshCoordinates(comm, globalSize);
  size_t numLocalCoords = coordinates->getLocalLength();

#if 0
  comm->barrier();
  for (int p=0; p < nprocs; p++){
    if (p==rank){
      cout << "Rank " << rank << ", " << numLocalCoords << "coords" << endl;
      const scalar_t *x = coordinates->getData(0).getRawPtr();
      const scalar_t *y = coordinates->getData(1).getRawPtr();
      const scalar_t *z = coordinates->getData(2).getRawPtr();
      for (lno_t i=0; i < numLocalCoords; i++)
        cout << " " << x[i] << " " << y[i] << " " << z[i] << endl;
    }
    cout.flush();
    comm->barrier();
  }
#endif

  if (weight_dim > 0){

    weights = arcp(new ArrayRCP<scalar_t> [weight_dim],
      0, weight_dim, true);
    for(int i = 0; i < weight_dim; ++i){
      //weights[i] = ArrayRCP<scalar_t>(weight[i]);
    }
  }

  MEMORY_CHECK(doMemory && rank==0, "After creating input");

  // Now call Zoltan to partition the problem.

  float ver;
  int aok = Zoltan_Initialize(argc, argv, &ver);

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
    oss << weight_dim;
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
  Zoltan_Set_Geom_Multi_Fn(zz, getCoordinates,NULL);

  int changes, numGidEntries, numLidEntries, numImport, numExport;
  gid_t * importGlobalGids, * importLocalGids;
  gid_t * exportGlobalGids, * exportLocalGids;
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

#else
#include <iostream>
int main(int argc, char *argv[])
{
  std::cout << "Test did not run due to faulty configuration." << std::endl;
  std::cout << "PASS" << std::endl;
}
#endif
