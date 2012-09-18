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

/*! \file PQJagged.cpp
    \brief An example of partitioning coordinates with PQJagged.
    \todo add more cases to this test.
 */

#include <Zoltan2_TestHelpers.hpp>
#include <Zoltan2_BasicCoordinateInput.hpp>
#include <Zoltan2_XpetraMultiVectorInput.hpp>
#include <Zoltan2_PartitioningSolution.hpp>
#include <Zoltan2_PartitioningProblem.hpp>
#include <GeometricGenerator.hpp>
#include <vector>

#include <Zoltan2_PartitioningSolutionQuality.hpp>


#include <Teuchos_LAPACK.hpp>
#include <fstream>
#include <string>
using namespace std;
using Teuchos::RCP;
using Teuchos::rcp;



typedef Tpetra::MultiVector<scalar_t, lno_t, gno_t, node_t> tMVector_t;
typedef Zoltan2::BasicUserTypes<scalar_t, gno_t, lno_t, gno_t> myTypes_t;


/*! \test PQJaggedTest.cpp
    An example of the use of the PQJagged algorithm to partition coordinate data.
 */


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

    fstream inParam(paramFileName.c_str());
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

  //MPI_Bcast(&size,1,MPI_INT, 0, MPI_COMM_WORLD);
  //MPI_Bcast(inp,size,MPI_CHAR, 0, MPI_COMM_WORLD);
  //Teuchos::broadcast<int, char>(comm, 0,inp);

  comm->broadcast(0, sizeof(int), (char*) &size);
  if(size == -1){
    throw "File " + paramFileName + " cannot be opened.";
  }
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

void GeometricGen(const RCP<const Teuchos::Comm<int> > & comm, partId_t numParts, float imbalance, std::string paramFile, std::string pqParts/*, std::string paramFile*/, partId_t k,
    bool force_binary, bool force_linear){

  Teuchos::ParameterList geoparams("geo params");

  readGeoGenParams(paramFile, geoparams, comm);

#ifdef HAVE_ZOLTAN2_OMP
  double begin = omp_get_wtime();
#endif
  GeometricGenerator<scalar_t, lno_t, gno_t, node_t> *gg = new GeometricGenerator<scalar_t, lno_t, gno_t, node_t>(geoparams,comm);
#ifdef HAVE_ZOLTAN2_OMP
  double end = omp_get_wtime();
  //cout << "GeometricGen Time:" << end - begin << endl;
#endif
  int coord_dim = gg->getCoordinateDimension();
  int weight_dim = gg->getWeightDimension();
  lno_t numLocalPoints = gg->getNumLocalCoords(); gno_t numGlobalPoints = gg->getNumGlobalCoords();
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

  delete gg;

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

  RCP< Tpetra::MultiVector<scalar_t, lno_t, gno_t, node_t> >tmVector = RCP< Tpetra::MultiVector<scalar_t, lno_t, gno_t, node_t> >(
      new Tpetra::MultiVector<scalar_t, lno_t, gno_t, node_t>( mp, coordView.view(0, coord_dim), coord_dim));


  RCP<const tMVector_t> coordsConst = Teuchos::rcp_const_cast<const tMVector_t>(tmVector);
  vector<const scalar_t *> weights;
  if(weight_dim){
    for (int i = 0; i < weight_dim;++i){
      weights.push_back(weight[i]);
    }
  }
  vector <int> stride;

#if 0
  typedef Zoltan2::BasicCoordinateInput<tMVector_t> inputAdapter_t;
  inputAdapter_t ia(localCount, globalIds, x, y, z, 1, 1, 1);
#else
  typedef Zoltan2::XpetraMultiVectorInput<tMVector_t> inputAdapter_t;
  //inputAdapter_t ia(coordsConst);
  inputAdapter_t ia(coordsConst,weights, stride);
#endif


  Teuchos::ParameterList params("test params");

  params.set("pqParts", pqParts);
  params.set("timer_output_stream" , "std::cout");
  params.set("num_global_parts", numParts);
  params.set("algorithm", "pqjagged");
  params.set("compute_metrics", "true");
  params.set("imbalance_tolerance", double(imbalance));
  params.set("parallel_part_calculation_count", k);
  if(force_binary){
    params.set("force_binary_search", "yes");
  }
  else {
    params.set("force_binary_search", "no");
  }
  if(force_linear){
    params.set("force_linear_search", "yes");
  }
  else {
    params.set("force_linear_search", "no");
  }
  params.set("bisection_num_test_cuts", 7);

#ifdef HAVE_ZOLTAN2_MPI


  Zoltan2::PartitioningProblem<inputAdapter_t> problem(&ia, &params,
      MPI_COMM_WORLD);


#else
  Zoltan2::PartitioningProblem<inputAdapter_t> problem(&ia, &params);
#endif

  problem.solve();

  //const Zoltan2::PartitioningSolution<inputAdapter_t> &solution =
  //    problem.getSolution();


  if (comm->getRank() == 0){

    problem.printMetrics(cout);

    cout << "testFromDataFile is done " << endl;
  }
  problem.printTimers();
  if(weight_dim){
    for(int i = 0; i < weight_dim; ++i)
    delete [] weight[i];
    delete [] weight;
  }
  if(coord_dim){
    for(int i = 0; i < coord_dim; ++i)
    delete [] coords[i];
    delete [] coords;
  }
}

void testFromDataFile(const RCP<const Teuchos::Comm<int> > & comm, partId_t numParts, float imbalance, std::string fname, std::string pqParts, partId_t k,
    bool force_binary, bool force_linear)
{
  //std::string fname("simple");
  cout << "running " << fname << endl;

  UserInputForTests uinput(testDataFilePath, fname, comm, true);

  RCP<tMVector_t> coords = uinput.getCoordinates();

#if 0
  size_t localCount = coords->getLocalLength();
  int dim = coords->getNumVectors();

  scalar_t *x=NULL, *y=NULL, *z=NULL;
  x = coords->getDataNonConst(0).getRawPtr();

  if (dim > 1){
    y = coords->getDataNonConst(1).getRawPtr();
    if (dim > 2)
      z = coords->getDataNonConst(2).getRawPtr();
  }

  const gno_t *globalIds = coords->getMap()->getNodeElementList().getRawPtr();

  typedef Zoltan2::BasicCoordinateInput<tMVector_t> inputAdapter_t;
  inputAdapter_t ia(localCount, globalIds, x, y, z, 1, 1, 1);
#else
  RCP<const tMVector_t> coordsConst = rcp_const_cast<const tMVector_t>(coords);

  typedef Zoltan2::XpetraMultiVectorInput<tMVector_t> inputAdapter_t;
  inputAdapter_t ia(coordsConst);
#endif


  Teuchos::ParameterList params("test params");


  params.set("pqParts", pqParts);
  params.set("timer_output_stream" , "std::cout");
  params.set("parallel_part_calculation_count", k);
  if(force_binary){
    params.set("force_binary_search", "yes");
  }
  else {
    params.set("force_binary_search", "no");
  }
  if(force_linear){
    params.set("force_linear_search", "yes");
  }
  else {
    params.set("force_linear_search", "no");
  }

  params.set("num_global_parts", numParts);
  params.set("algorithm", "pqjagged");
  params.set("compute_metrics", "true");
  params.set("imbalance_tolerance", double(imbalance));
  params.set("bisection_num_test_cuts", 7);

#ifdef HAVE_ZOLTAN2_MPI
  Zoltan2::PartitioningProblem<inputAdapter_t> problem(&ia, &params,
      MPI_COMM_WORLD);
#else
  Zoltan2::PartitioningProblem<inputAdapter_t> problem(&ia, &params);
#endif

  problem.solve();

  //const Zoltan2::PartitioningSolution<inputAdapter_t> &solution =
      //problem.getSolution();

  if (comm->getRank() == 0){
    problem.printMetrics(cout);

    cout << "testFromDataFile is done " << endl;
  }

  problem.printTimers();


}






string convert_to_string(char *args){
  string tmp = "";
  for(int i = 0; args[i] != 0; i++)
    tmp += args[i];
  return tmp;
}
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

void getArgVals(int argc, char **argv,   partId_t &numParts, float &imbalance ,
     string &pqParts, int &opt,std::string &fname/*, std::string &paramFile*/, partId_t &k,
     bool &force_binary , bool &force_linear){

  bool isCset = false;
  bool isPset = false;
  bool isFset = false;

  for(int i = 0; i < argc; ++i){
    string tmp = convert_to_string(argv[i]);
    string identifier = "";
    long long int value = -1; double fval = -1;
    if(!getArgumentValue(identifier, fval, tmp)) continue;
    value = (long long int) (fval);

    if(identifier == "C"){

      if(value > 0){


        numParts=value;
        isCset = true;
      } else {
        throw  "Invalid argument at " + tmp;
      }
    } else
      if(identifier == "P"){

        stringstream stream(stringstream::in | stringstream::out);

        stream << tmp;
        string ttmp;
        getline(stream, ttmp, '=');
        stream >> pqParts;
        isPset = true;
      }else if(identifier == "I"){
        if(fval > 0){
          imbalance=fval;
        } else {
          throw "Invalid argument at " + tmp;
        }
      } else if(identifier == "F"){
        stringstream stream(stringstream::in | stringstream::out);
        stream << tmp;
        getline(stream, fname, '=');

        stream >> fname;
        isFset = true;
      }
/*
      else if(identifier == "PF"){
        stringstream stream(stringstream::in | stringstream::out);
        stream << tmp;
        getline(stream, paramFile, '=');
        stream >> paramFile;

      }
      */
      else if(identifier == "EL"){
        if(value == 0){
          force_linear = false;
        } else if(value == 1){
          force_linear = true;
        } else {
          throw "Invalid argument at " + tmp;
        }
      }
      else if(identifier == "EB"){
        if(value == 0){
          force_binary = false;
        } else if(value == 1){
          force_binary = true;
        } else {
          throw "Invalid argument at " + tmp;
        }
      }
      else if(identifier == "O"){
        if(value >= 0 && value <= 3){
          opt = value;
        } else {
          throw "Invalid argument at " + tmp;
        }
      }else if(identifier == "K"){
        if(value >=0 ){
          k = value;
        } else {
          throw "Invalid argument at " + tmp;
        }
      }
      else {
        throw "Invalid argument at " + tmp;
      }

  }
  if(!(isCset && isPset && isFset)){
    throw "P, C and F are mandatory arguments.";
  }

}

void print_usage(char *executable){
  cout << "\nUsage:" << endl;
  cout << executable << " arglist" << endl;
  cout << "arglist:" << endl;
  cout << "\tC=numParts: numParts > 0" << endl;
  cout << "\tP=pqJaggedPart: Example: P=512,512" << endl;
  cout << "\tI=imbalance: Example I=1.03 (ignored for now.)" << endl;
  cout << "\tF=filePath: When O=0 the path of the coordinate input file, for O>1 the path to the geometric generator parameter file." << endl;
  cout << "\tEL=force linear search: EL=0 as default, EL=1 will force to do linear search regardless of the part number" << endl;
  cout << "\tEB=force binary search: EB=0 as default, EB=1 will force to do binary search regardless of the part number" << endl;
  cout << "\tO=input option: O=0 for reading coordinate from file, O>0 for generating coordinate from coordinate generator file. Default will run geometric generator." << endl;
  cout << "\tK=concurrent part calculation input: K>0." << endl;
  cout << "Example:\n" << executable << " P=2,2,2 C=8 F=simple O=0" << endl;
}

int main(int argc, char *argv[])
{
  Teuchos::GlobalMPISession session(&argc, &argv);

  RCP<const Teuchos::Comm<int> > tcomm = Teuchos::DefaultComm<int>::getComm();
  int rank = tcomm->getRank();


  partId_t numParts = 10; float imbalance = 1.03;
  partId_t k = 1;

  string pqParts = "";
  int opt = 1;
  std::string fname = "";
  //std::string paramFile = "";
  bool force_binary = false, force_linear = false;

  try{
    try {
    getArgVals(argc, argv,   numParts, imbalance ,
        pqParts, opt,fname/*, paramFile*/, k, force_binary, force_linear);
    }
    catch(std::string s){
      if(tcomm->getRank() == 0){
        print_usage(argv[0]);
      }
      throw s;
    }

    switch (opt){
    case 0:
      testFromDataFile(tcomm,numParts, imbalance,fname,pqParts, k, force_binary, force_linear);
      break;
    default:
      GeometricGen(tcomm, numParts, imbalance, fname, pqParts/*, paramFile*/, k, force_binary, force_linear);
      break;
    }

    if (rank == 0)
      std::cout << "PASS" << std::endl;
  }


  catch(std::string s){
    if (rank == 0)
    cerr << s << endl;
  }

  catch(char * s){
    if (rank == 0)
    cerr << s << endl;
  }

  catch(char const* s){
    if (rank == 0)
    cerr << s << endl;
  }
  //if (rank == 0)
  //  serialTest(numParts, numCoords, imbalance);

}
