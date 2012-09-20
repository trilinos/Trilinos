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

/*! \file UserInputForTests.hpp
 *  \brief Generate input for testing purposes.
 */

#ifndef USERINPUTFORTESTS
#define USERINPUTFORTESTS

#include <Zoltan2_XpetraTraits.hpp>
#include <Zoltan2_ChacoReader.hpp>

#include <Tpetra_MultiVector.hpp>
#include <Xpetra_Vector.hpp>
#include <Xpetra_CrsMatrix.hpp>
#include <Xpetra_CrsGraph.hpp>

#include <MatrixMarket_Tpetra.hpp>
#include <Galeri_XpetraMatrixFactory.hpp>
#include <Galeri_XpetraParameters.hpp>

#include <Kokkos_DefaultNode.hpp>

//#include <Xpetra_EpetraUtils.hpp>
#ifdef HAVE_ZOLTAN2_MPI
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif

using Teuchos::RCP;
using Teuchos::ArrayRCP;
using Teuchos::ArrayView;
using Teuchos::Array;
using Teuchos::Comm;
using Teuchos::rcp;
using Teuchos::arcp;
using Teuchos::rcp_const_cast;
using std::string;

/*! \brief A class that generates typical user input for testing.
 *
 *  Given:
 *  \li The name of matrix market file in the data directory, or
 *  \li The name of a chaco file in the Zoltan1 test directory, or
 *  \li x, y and z dimensions and a problem type
 *
 *  Retrieve the data in any of the following forms:
 *   \li a Tpetra::CrsMatrix
 *   \li a Tpetra::CrsGraph
 *   \li a Xpetra::CrsMatrix
 *   \li a Xpetra::CrsGraph
 *   \li a Epetra_CrsMatrix (if built with double, int, int)
 *   \li a Epetra_CrsGraph  (if built with double, int, int)
 *
 *  Retrieve any of the these with the same map but with random data:
 *   \li a Tpetra::Vector
 *   \li a Tpetra::MultiVector
 *   \li a Xpetra::Vector
 *   \li a Xpetra::MultiVector
 *
 *  Coordinates and/or weights will be retrieved if they are found.
 *  They are provided in a Tpetra::MultiVector.
 *
 *  Random data can also be generated.
 *
 *  \todo for very large files, each process reads in part of the file
 *  \todo Zoltan1 mtx and mtxp files
 */

class UserInputForTests
{
public:

  typedef Tpetra::CrsMatrix<scalar_t, lno_t, gno_t, node_t> tcrsMatrix_t;
  typedef Tpetra::CrsGraph<lno_t, gno_t, node_t> tcrsGraph_t;
  typedef Tpetra::Vector<scalar_t, lno_t, gno_t, node_t> tVector_t;
  typedef Tpetra::MultiVector<scalar_t, lno_t, gno_t, node_t> tMVector_t;

  typedef Xpetra::CrsMatrix<scalar_t, lno_t, gno_t, node_t> xcrsMatrix_t;
  typedef Xpetra::CrsGraph<lno_t, gno_t, node_t> xcrsGraph_t;
  typedef Xpetra::Vector<scalar_t, lno_t, gno_t, node_t> xVector_t;
  typedef Xpetra::MultiVector<scalar_t, lno_t, gno_t, node_t> xMVector_t;

  typedef Tpetra::Map<lno_t, gno_t, node_t> map_t;
  typedef Tpetra::Export<lno_t, gno_t, node_t> export_t;
  typedef Tpetra::Import<lno_t, gno_t, node_t> import_t;
  typedef Kokkos::DefaultNode::DefaultNodeType default_node_t;

  /*! \brief Constructor that reads in a matrix/graph from disk.
   *   \param path is the path to the test data.  In the case of
   *       Zoltan2 test data it is the path to the "data" directory.
   *       In the case of Zoltan1 test data, it is the path to the
   *       "test" directory.
   *   \param  testData is the root name of the data file or files
   *       of interest.
   *   \param  c  is the communicator for the processes that
   *         share the data.
   *   \param  debugInfo if true process zero will print out status.
   *
   *  For example, if \c path is the path to the Zoltan1 test
   *  directory and \c testData is \c brack2_3, then we'll read
   *  in ch_brack2_3/brack2_3.graph and ch_brack2_3/brack2_3.coords.
   */

  UserInputForTests(string path, string testData, 
    const RCP<const Comm<int> > &c, bool debugInfo=false);

  /*! \brief Constructor that generates an arbitrary sized sparse matrix.
   *   \param x the x dimension of the mesh that generates the matrix.
   *   \param y the y dimension of the mesh that generates the matrix.
   *   \param z the z dimension of the mesh that generates the matrix.
   *   \param problemType the type of problem that will generate a
   *             sparse matrix from the mesh.  If the problemType is
   *             empty we'll pick a default.
   *   \param  c  is the communicator for the processes that
   *         share the data.
   *   \param  debugInfo if true process zero will print out status.
   *
   * Problems can be "Laplace1D", "Laplace2D", "Star2D", "BigStar2D", 
   * "Laplace3D", "Brick3D" and "Identity".
   * See Galeri::Xpetra::CreateCrsMatrix() for more information
   * about problem types.
   */
  UserInputForTests(int x, int y, int z, string matrixType,
    const RCP<const Comm<int> > &c, bool debugInfo=false);

  /*! \brief Generate lists of random scalars.
   */
  static void getRandomData(unsigned int seed, lno_t length,
    scalar_t min, scalar_t max, ArrayView<ArrayRCP<scalar_t > > data);

  RCP<tMVector_t> getCoordinates();

  RCP<tMVector_t> getWeights();

  RCP<tMVector_t> getEdgeWeights();

  RCP<tcrsMatrix_t> getTpetraCrsMatrix();

  RCP<tcrsGraph_t> getTpetraCrsGraph();

  RCP<tVector_t> getTpetraVector();

  RCP<tMVector_t> getTpetraMultiVector(int nvec);

  RCP<xcrsMatrix_t> getXpetraCrsMatrix();

  RCP<xcrsGraph_t> getXpetraCrsGraph();

  RCP<xVector_t> getXpetraVector();

  RCP<xMVector_t> getXpetraMultiVector(int nvec);

#ifdef HAVE_EPETRA_DATA_TYPES
  RCP<Epetra_CrsGraph> getEpetraCrsGraph();

  RCP<Epetra_CrsMatrix> getEpetraCrsMatrix();

  RCP<Epetra_Vector> getEpetraVector();

  RCP<Epetra_MultiVector> getEpetraMultiVector(int nvec);
#endif

private:

  bool verbose_;

  const RCP<const Comm<int> > tcomm_;

  RCP<tcrsMatrix_t> M_; 
  RCP<xcrsMatrix_t> xM_; 

  RCP<tMVector_t> xyz_;
  RCP<tMVector_t> vtxWeights_;
  RCP<tMVector_t> edgWeights_;

#ifdef HAVE_EPETRA_DATA_TYPES
  RCP<const Epetra_Comm> ecomm_;
  RCP<Epetra_CrsMatrix> eM_; 
  RCP<Epetra_CrsGraph> eG_; 
#endif

  // Read a Matrix Market file into M_
  // using Tpetra::MatrixMarket::Reader.
  // If there are "Tim Davis" style coordinates
  // that go with the file,  read those into xyz_.

  void readMatrixMarketFile(string path, string testData);

  // Build matrix M_ from a mesh and a problem type
  // with Galeri::Xpetra.

  void buildCrsMatrix(int xdim, int ydim, int zdim, string type);

  // Read a Zoltan1 Chaco or Matrix Market file
  // into M_.  If it has geometric coordinates,
  // read them into xyz_.  If it has weights,
  // read those into vtxWeights_ and edgWeights_.

  void readZoltanTestData(string path, string testData);

  // Read Zoltan data that is in a .graph file.
  void getChacoGraph(FILE *fptr, string name);

  // Read Zoltan data that is in a .coords file.
  void getChacoCoords(FILE *fptr, string name);

  static string shortPathName(std::ostringstream &oss, string dirName);
};

UserInputForTests::UserInputForTests(string path, string testData, 
  const RCP<const Comm<int> > &c, bool debugInfo):
    verbose_(debugInfo),
    tcomm_(c), M_(), xM_(), xyz_(), vtxWeights_(), edgWeights_()
#ifdef HAVE_EPETRA_DATA_TYPES
    ,ecomm_(), eM_(), eG_()
#endif
{
  bool zoltan1 = false;
  string::size_type loc = path.find("/data/");  // Zoltan2 data
  if (loc == string::npos)
    zoltan1 = true;

  if (zoltan1)
    readZoltanTestData(path, testData);
  else
    readMatrixMarketFile(path, testData);
  
#ifdef HAVE_EPETRA_DATA_TYPES
  ecomm_ = Xpetra::toEpetra(c);
#endif
}

UserInputForTests::UserInputForTests(int x, int y, int z, 
  string matrixType, const RCP<const Comm<int> > &c, bool debugInfo):
    verbose_(debugInfo),
    tcomm_(c), M_(), xM_(), xyz_(), vtxWeights_(), edgWeights_()
#ifdef HAVE_EPETRA_DATA_TYPES
    ,ecomm_(), eM_(), eG_() 
#endif
{
  if (matrixType.size() == 0){
    int dim = 0;
    if (x > 0) dim++;
    if (y > 0) dim++;
    if (z > 0) dim++;
    if (dim == 1)
      matrixType = string("Laplace1D");
    else if (dim == 2)
      matrixType = string("Laplace2D");
    else if (dim == 3)
      matrixType = string("Laplace3D");
    else
      throw std::runtime_error("input");

    if (verbose_ && tcomm_->getRank() == 0)
      std::cout << "UserInputForTests, Matrix type : " << matrixType << std::endl;
  }

  buildCrsMatrix(x, y, z, matrixType);
  
#ifdef HAVE_EPETRA_DATA_TYPES
  ecomm_ = Xpetra::toEpetra(c);
#endif
}

RCP<UserInputForTests::tMVector_t> UserInputForTests::getCoordinates()
{ 
  if (xyz_.is_null())
    throw std::runtime_error("could not read coord file");
  return xyz_;
}

RCP<UserInputForTests::tMVector_t> UserInputForTests::getWeights()
{ 
  return vtxWeights_;
}

RCP<UserInputForTests::tMVector_t> UserInputForTests::getEdgeWeights()
{ 
  return edgWeights_;
}

RCP<UserInputForTests::tcrsMatrix_t> UserInputForTests::getTpetraCrsMatrix() 
{ 
  if (M_.is_null())
    throw std::runtime_error("could not read mtx file");
  return M_;
}

RCP<UserInputForTests::tcrsGraph_t> UserInputForTests::getTpetraCrsGraph() 
{ 
  if (M_.is_null())
    throw std::runtime_error("could not read mtx file");
  return rcp_const_cast<tcrsGraph_t>(M_->getCrsGraph());
}

RCP<UserInputForTests::tVector_t> UserInputForTests::getTpetraVector() 
{ 
  RCP<tVector_t> V = rcp(new tVector_t(M_->getRowMap(),  1));
  V->randomize();
  
  return V;
}

RCP<UserInputForTests::tMVector_t> UserInputForTests::getTpetraMultiVector(int nvec) 
{ 
  RCP<tMVector_t> mV = rcp(new tMVector_t(M_->getRowMap(), nvec));
  mV->randomize();
  
  return mV;
}

RCP<UserInputForTests::xcrsMatrix_t> UserInputForTests::getXpetraCrsMatrix() 
{ 
  if (M_.is_null())
    throw std::runtime_error("could not read mtx file");
  return xM_;
}

RCP<UserInputForTests::xcrsGraph_t> UserInputForTests::getXpetraCrsGraph() 
{ 
  if (M_.is_null())
    throw std::runtime_error("could not read mtx file");
  return rcp_const_cast<xcrsGraph_t>(xM_->getCrsGraph());
}

RCP<UserInputForTests::xVector_t> UserInputForTests::getXpetraVector() 
{ 
  RCP<const tVector_t> tV = getTpetraVector();
  RCP<const xVector_t> xV =
    Zoltan2::XpetraTraits<tVector_t>::convertToXpetra(tV);
  return rcp_const_cast<xVector_t>(xV);
}

RCP<UserInputForTests::xMVector_t> UserInputForTests::getXpetraMultiVector(int nvec) 
{ 
  RCP<const tMVector_t> tMV = getTpetraMultiVector(nvec);
  RCP<const xMVector_t> xMV =
    Zoltan2::XpetraTraits<tMVector_t>::convertToXpetra(tMV);
  return rcp_const_cast<xMVector_t>(xMV);
}

#ifdef HAVE_EPETRA_DATA_TYPES
RCP<Epetra_CrsGraph> UserInputForTests::getEpetraCrsGraph()
{
  if (M_.is_null())
    throw std::runtime_error("could not read mtx file");
  RCP<const tcrsGraph_t> tgraph = M_->getCrsGraph();
  RCP<const Tpetra::Map<lno_t, gno_t> > trowMap = tgraph->getRowMap();
  RCP<const Tpetra::Map<lno_t, gno_t> > tcolMap = tgraph->getColMap();

  int nElts = static_cast<int>(trowMap->getGlobalNumElements());
  int nMyElts = static_cast<int>(trowMap->getNodeNumElements());
  int base = trowMap->getIndexBase();
  ArrayView<const int> gids = trowMap->getNodeElementList();

  Epetra_BlockMap erowMap(nElts, nMyElts,
    gids.getRawPtr(), 1, base, *ecomm_);

  Array<int> rowSize(nMyElts);
  for (int i=0; i < nMyElts; i++){
    rowSize[i] = static_cast<int>(M_->getNumEntriesInLocalRow(i+base));
  }

  size_t maxRow = M_->getNodeMaxNumRowEntries();
  Array<int> colGids(maxRow);
  ArrayView<const int> colLid;

  eG_ = rcp(new Epetra_CrsGraph(Copy, erowMap, 
    rowSize.getRawPtr(), true));

  for (int i=0; i < nMyElts; i++){
    tgraph->getLocalRowView(i+base, colLid);
    for (int j=0; j < colLid.size(); j++)
      colGids[j] = tcolMap->getGlobalElement(colLid[j]);
    eG_->InsertGlobalIndices(gids[i], rowSize[i], colGids.getRawPtr());
  }
  eG_->FillComplete();
  return eG_;
}

RCP<Epetra_CrsMatrix> UserInputForTests::getEpetraCrsMatrix()
{
  if (M_.is_null())
    throw std::runtime_error("could not read mtx file");
  RCP<Epetra_CrsGraph> egraph = getEpetraCrsGraph();
  eM_ = rcp(new Epetra_CrsMatrix(Copy, *egraph));

  size_t maxRow = M_->getNodeMaxNumRowEntries();
  int nrows = egraph->NumMyRows();
  int base = egraph->IndexBase();
  const Epetra_BlockMap &rowMap = egraph->RowMap();
  const Epetra_BlockMap &colMap = egraph->ColMap();
  Array<int> colGid(maxRow);

  for (int i=0; i < nrows; i++){
    ArrayView<const int> colLid;
    ArrayView<const scalar_t> nz;
    M_->getLocalRowView(i+base, colLid, nz);
    size_t rowSize = colLid.size();
    int rowGid = rowMap.GID(i+base);
    for (size_t j=0; j < rowSize; j++){
      colGid[j] = colMap.GID(colLid[j]);
    }
    eM_->InsertGlobalValues(
      rowGid, rowSize, nz.getRawPtr(), colGid.getRawPtr());
  }
  eM_->FillComplete();
  return eM_;
}

RCP<Epetra_Vector> UserInputForTests::getEpetraVector() 
{ 
  RCP<Epetra_CrsGraph> egraph = getEpetraCrsGraph();
  RCP<Epetra_Vector> V = rcp(new Epetra_Vector(egraph->RowMap()));
  V->Random();
  return V;
}

RCP<Epetra_MultiVector> UserInputForTests::getEpetraMultiVector(int nvec) 
{ 
  RCP<Epetra_CrsGraph> egraph = getEpetraCrsGraph();
  RCP<Epetra_MultiVector> mV = 
    rcp(new Epetra_MultiVector(egraph->RowMap(), nvec));
  mV->Random();
  return mV;
}
#endif

void UserInputForTests::getRandomData(unsigned int seed, lno_t length, 
    scalar_t min, scalar_t max,
    ArrayView<ArrayRCP<scalar_t > > data)
{
  if (length < 1)
    return;

  size_t dim = data.size();
  for (size_t i=0; i < dim; i++){
    scalar_t *tmp = new scalar_t [length];
    if (!tmp)
       throw (std::bad_alloc());
    data[i] = Teuchos::arcp(tmp, 0, length, true);
  }

  scalar_t scalingFactor = (max-min) / RAND_MAX;
  srand(seed);
  for (size_t i=0; i < dim; i++){
    scalar_t *x = data[i].getRawPtr();
    for (lno_t j=0; j < length; j++)
      *x++ = min + (scalar_t(rand()) * scalingFactor);
  }
}

void UserInputForTests::readMatrixMarketFile(string path, string testData)
{
  std::ostringstream fname;
  fname << path << "/" << testData << ".mtx";

  RCP<Kokkos::DefaultNode::DefaultNodeType> dnode 
    = Kokkos::DefaultNode::getDefaultNode();

  if (verbose_ && tcomm_->getRank() == 0)
    std::cout << "UserInputForTests, Read: " << 
      shortPathName(fname, "test") << std::endl;

  // This reader has some problems.  "Pattern" matrices
  // cannot be read.  Until the
  // reader is fixed, we'll have to get inputs that are consistent with
  // the reader. (Tpetra bug 5611 and 5624)
 
  bool aok = true;
  try{
    M_ = Tpetra::MatrixMarket::Reader<tcrsMatrix_t>::readSparseFile(
      fname.str(), tcomm_, dnode, true, true, false);
  }
  catch (std::exception &e) {
    //TEST_FAIL_AND_THROW(*tcomm_, 1, e.what());
    aok = false;
  }

  if (aok){
    RCP<const xcrsMatrix_t> xm = 
      Zoltan2::XpetraTraits<tcrsMatrix_t>::convertToXpetra(M_);
    xM_ = rcp_const_cast<xcrsMatrix_t>(xm);
  }
  else{
    if (tcomm_->getRank() == 0)
      std::cout << "UserInputForTests unable to read matrix." << std::endl;
  }

  // Open the coordinate file.

  fname.str("");
  fname << path << "/" << testData << "_coord.mtx";

  size_t coordDim = 0, numGlobalCoords = 0;
  size_t msg[2]={0,0};
  ArrayRCP<ArrayRCP<scalar_t> > xyz;
  std::ifstream coordFile;

  if (tcomm_->getRank() == 0){

    if (verbose_)
      std::cout << "UserInputForTests, Read: " << 
         shortPathName(fname, "test") << std::endl;
  
    int fail = 0;
    try{
      coordFile.open(fname.str().c_str());
    }
    catch (std::exception &e){ // there is no coordinate file
      fail = 1;
    }
  
    if (!fail){
  
      // Read past banner to number and dimension of coordinates.
    
      char c[256];
      bool done=false;
    
      while (!done && !fail && coordFile.good()){
        coordFile.getline(c, 256);
        if (!c[0])
          fail = 1;
        else if (c[0] == '%')
          continue;
        else {
          done=true;
          std::istringstream s(c);
          s >> numGlobalCoords >> coordDim;
          if (!s.eof() || numGlobalCoords < 1 || coordDim < 1)
            fail=1;
        }
      }

      if (done){
    
        // Read in the coordinates.
      
        xyz = Teuchos::arcp(new ArrayRCP<scalar_t> [coordDim], 0, coordDim);
      
        for (size_t dim=0; !fail && dim < coordDim; dim++){
          size_t idx;
          scalar_t *tmp = new scalar_t [numGlobalCoords];
          if (!tmp)
            fail = 1;
          else{
            xyz[dim] = Teuchos::arcp(tmp, 0, numGlobalCoords);
        
            for (idx=0; !coordFile.eof() && idx < numGlobalCoords; idx++){
              coordFile.getline(c, 256);
              std::istringstream s(c);
              s >> tmp[idx];
            }
      
            if (idx < numGlobalCoords)
              fail = 1;
          }
        }

        if (fail){
          ArrayRCP<scalar_t> emptyArray;
          for (size_t dim=0; dim < coordDim; dim++)
            xyz[dim] = emptyArray;   // free the memory

          coordDim = 0;
        }
      }
      else{
        fail = 1;
      }
    
      coordFile.close();
    }

    msg[0] = coordDim;
    msg[1] = numGlobalCoords;
  }

  // Broadcast coordinate dimension
  Teuchos::broadcast<int, size_t>(*tcomm_, 0, 2, msg);

  coordDim = msg[0];
  numGlobalCoords= msg[1];

  if (coordDim == 0)
    return;

  gno_t base;
  RCP<const map_t> toMap;

  if (!M_.is_null()){
    base = M_->getIndexBase();
    const RCP<const map_t> &mapM = M_->getRowMap();
    toMap = mapM;
  }
  else{
    base = 0;
    toMap = rcp(new map_t(numGlobalCoords, base, tcomm_));
  }

  // Export coordinates to their owners

  xyz_ = rcp(new tMVector_t(toMap, coordDim));

  ArrayRCP<ArrayView<const scalar_t> > coordLists(coordDim);

  if (tcomm_->getRank() == 0){

    for (size_t dim=0; dim < coordDim; dim++)
      coordLists[dim] = xyz[dim].view(0, numGlobalCoords);

    gno_t *tmp = new gno_t [numGlobalCoords];
    if (!tmp)
      throw std::bad_alloc();

    ArrayRCP<const gno_t> rowIds = Teuchos::arcp(tmp, 0, numGlobalCoords);

    gno_t basePlusNumGlobalCoords = base+numGlobalCoords;
    for (gno_t id=base; id < basePlusNumGlobalCoords; id++)
      *tmp++ = id;

    RCP<const map_t> fromMap = rcp(new map_t(numGlobalCoords, 
      rowIds.view(0, numGlobalCoords), base, tcomm_));

    tMVector_t allCoords(fromMap, coordLists.view(0, coordDim), coordDim);

    export_t exporter(fromMap, toMap);

    xyz_->doExport(allCoords, exporter, Tpetra::INSERT);
  }
  else{

    RCP<const map_t> fromMap = rcp(new map_t(numGlobalCoords, 
      ArrayView<gno_t>(), base, tcomm_));

    tMVector_t allCoords(fromMap, coordLists.view(0, coordDim), coordDim);

    export_t exporter(fromMap, toMap);

    xyz_->doExport(allCoords, exporter, Tpetra::INSERT);
  }
}

void UserInputForTests::buildCrsMatrix(int xdim, int ydim, int zdim, 
  string problemType)
{
  Teuchos::CommandLineProcessor tclp;
  Galeri::Xpetra::Parameters<gno_t> params(tclp,
     xdim, ydim, zdim, problemType);

  RCP<const Tpetra::Map<lno_t, gno_t> > map =
    rcp(new Tpetra::Map<lno_t, gno_t>(
      params.GetNumGlobalElements(), 0, tcomm_));

  if (verbose_ && tcomm_->getRank() == 0){
    std::cout << "UserInputForTests, Create matrix with " << problemType;
    std::cout << " (and" << xdim;
    if (zdim > 0)
      std::cout << " x " << ydim << " x " << zdim << std::endl;
    else if (ydim > 0)
      std::cout << " x"  << ydim << " x 1" << std::endl;
    else
      std::cout << "x 1 x 1" << std::endl;

    std::cout << " mesh)" << std::endl;
  }

  try{
    M_ = Galeri::Xpetra::CreateCrsMatrix<scalar_t, lno_t, gno_t, 
      Tpetra::Map<lno_t, gno_t>, 
      Tpetra::CrsMatrix<scalar_t, lno_t, gno_t> >(params.GetMatrixType(),
         map, params.GetParameterList()); 
  }
  catch (std::exception &e) {    // Probably not enough memory
    TEST_FAIL_AND_THROW(*tcomm_, 1, e.what());
  }

  RCP<const xcrsMatrix_t> xm = 
    Zoltan2::XpetraTraits<tcrsMatrix_t>::convertToXpetra(M_);
  xM_ = rcp_const_cast<xcrsMatrix_t>(xm);

  // Compute the coordinates for the matrix rows.

  if (verbose_ && tcomm_->getRank() == 0)
    std::cout << 
     "UserInputForTests, Implied matrix row coordinates computed" << 
       std::endl;
  
  ArrayView<const gno_t> gids = map->getNodeElementList();
  lno_t count = gids.size();
  int dim = 3;
  size_t pos = problemType.find("2D");
  if (pos != string::npos)
    dim = 2;
  else if (problemType == string("Laplace1D") || 
           problemType == string("Identity"))
    dim = 1;

  Array<ArrayRCP<scalar_t> > coordinates(dim);

  if (count > 0){
    for (int i=0; i < dim; i++){
      scalar_t *c = new scalar_t [count];
      if (!c)
        throw(std::bad_alloc());
      coordinates[i] = Teuchos::arcp(c, 0, count, true);
    }
  
    if (dim==3){
      scalar_t *x = coordinates[0].getRawPtr();
      scalar_t *y = coordinates[1].getRawPtr();
      scalar_t *z = coordinates[2].getRawPtr();
      gno_t xySize = xdim * ydim;
      for (lno_t i=0; i < count; i++){
        gno_t iz = gids[i] / xySize;
        gno_t xy = gids[i] - iz*xySize;
        z[i] = scalar_t(iz);
        y[i] = scalar_t(xy / xdim);
        x[i] = scalar_t(xy % xdim);
      }
    }
    else if (dim==2){
      scalar_t *x = coordinates[0].getRawPtr();
      scalar_t *y = coordinates[1].getRawPtr();
      for (lno_t i=0; i < count; i++){
        y[i] = scalar_t(gids[i] / xdim);
        x[i] = scalar_t(gids[i] % xdim);
      }
    }
    else{
      scalar_t *x = coordinates[0].getRawPtr();
      for (lno_t i=0; i < count; i++)
        x[i] = scalar_t(gids[i]);
    }
  }

  Array<ArrayView<const scalar_t> > coordView(dim);
  if (count > 0)
    for (int i=0; i < dim; i++)
      coordView[i] = coordinates[i].view(0,count);

  xyz_ = rcp(new tMVector_t(map, coordView.view(0, dim), dim));
}

void UserInputForTests::readZoltanTestData(string path, string testData)
{
  int rank = tcomm_->getRank();
  FILE *graphFile = NULL;
  FILE *coordFile = NULL;
  int fileInfo[2];

  if (rank == 0){
    std::ostringstream chGraphFileName;
    chGraphFileName << path << "/ch_" << testData << "/" << testData << ".graph";
    std::ostringstream chCoordFileName;
    chCoordFileName << path << "/ch_" << testData << "/" << testData << ".coords";
    memset(fileInfo, 0, sizeof(int) * 2);

    graphFile = fopen(chGraphFileName.str().c_str(), "r");
    if (graphFile){
      fileInfo[0] = 1;
      if (verbose_ && tcomm_->getRank() == 0)
        std::cout << "UserInputForTests, open " << 
          shortPathName(chGraphFileName, "zoltan") << std::endl;
      
      coordFile = fopen(chCoordFileName.str().c_str(), "r");
      if (coordFile){
        fileInfo[1] = 1;
        if (verbose_ && tcomm_->getRank() == 0)
          std::cout << "UserInputForTests, open " << 
            shortPathName(chCoordFileName, "zoltan") << std::endl;
      }
    }
  }

  Teuchos::broadcast<int, int>(*tcomm_, 0, 2, fileInfo);

  bool haveGraph = (fileInfo[0] == 1);
  bool haveCoords = (fileInfo[1] == 1);

  if (haveGraph){
    // Writes M_, vtxWeights_, and edgWeights_ and closes file.
    try{
      getChacoGraph(graphFile, testData);
    }
    Z2_FORWARD_EXCEPTIONS

    if (haveCoords){
      // Writes xyz_ and closes the file.
      try{
        getChacoCoords(coordFile, testData);
      }
      Z2_FORWARD_EXCEPTIONS
    }
  }

  RCP<const xcrsMatrix_t> xm =
    Zoltan2::XpetraTraits<tcrsMatrix_t>::convertToXpetra(M_);
  xM_ = rcp_const_cast<xcrsMatrix_t>(xm);
}

void UserInputForTests::getChacoGraph(FILE *fptr, string fname)
{
  int rank = tcomm_->getRank();
  int graphCounts[5];
  int nvtxs=0, nedges=0;
  int vwgt_dim=0, ewgt_dim=0;
  int *start = NULL, *adj = NULL;
  float *ewgts = NULL, *vwgts = NULL;
  size_t *nzPerRow = NULL;
  size_t maxRowLen = 0;
  gno_t base = 0;
  ArrayRCP<const size_t> rowSizes;
  int fail = 0;
  bool haveEdges = true;

  if (rank == 0){

    memset(graphCounts, 0, 5*sizeof(int));
  
    // This function is in the Zoltan C-library.
  
    // Reads in the file and closes it when done.
    char *nonConstName = new char [fname.size() + 1];
    strcpy(nonConstName, fname.c_str());

    fail = Zoltan2::chaco_input_graph(fptr, nonConstName,
      &start, &adj, &nvtxs, &vwgt_dim, &vwgts, &ewgt_dim, &ewgts);
    delete [] nonConstName;

    // There are Zoltan2 test graphs that have no edges.

    if (start == NULL)
      haveEdges = false;

    if (verbose_){
      std::cout << "UserInputForTests, " << nvtxs << " vertices,";
      if (haveEdges)
        std::cout << start[nvtxs] << " edges,";
      else
        std::cout << "no edges,";
      std::cout << vwgt_dim << " vertex weights, ";
      std::cout << ewgt_dim << " edge weights" << std::endl;
    }

    if (nvtxs==0)
      fail = true;

    if (fail){
      Teuchos::broadcast<int, int>(*tcomm_, 0, 5, graphCounts);
      throw std::runtime_error("Unable to read chaco file");
    }

    if (haveEdges)
      nedges = start[nvtxs];

    nzPerRow = new size_t [nvtxs];
    if (!nzPerRow)
      throw std::bad_alloc();
    rowSizes = arcp(nzPerRow, 0, nvtxs, true);

    if (haveEdges){
      for (int i=0; i < nvtxs; i++){
        nzPerRow[i] = start[i+1] - start[i];
        if (nzPerRow[i] > maxRowLen)
          maxRowLen = nzPerRow[i];
      }
    }
    else{
      memset(nzPerRow, 0, sizeof(size_t) * nvtxs);
    }

    if (haveEdges){
      free(start);
      start = NULL;
    }
  
    // Make sure base gid is zero.

    if (nedges){
      int chbase = adj[0];
      for (int i=1; i < nedges; i++)
        if (adj[i] < chbase)
          chbase = adj[i];
  
      if (chbase > 0){
        for (int i=0; i < nedges; i++)
          adj[i] -= chbase;
      }
    }

    graphCounts[0] = nvtxs;
    graphCounts[1] = nedges;
    graphCounts[2] = vwgt_dim;
    graphCounts[3] = ewgt_dim;
    graphCounts[4] = maxRowLen;  // size_t maxRowLen will fit; it is <= (int-int)
  }
  
  Teuchos::broadcast<int, int>(*tcomm_, 0, 5, graphCounts);

  if (graphCounts[0] == 0)
    throw std::runtime_error("Unable to read chaco file");

  haveEdges = (graphCounts[1] > 0);

  RCP<tcrsMatrix_t> fromMatrix;
  RCP<const map_t> fromMap;

  // Create a Tpetra::CrsMatrix where rank 0 has entire matrix.

  if (rank == 0){
    fromMap = rcp(new map_t(nvtxs, nvtxs, base, tcomm_));

    fromMatrix = 
      rcp(new tcrsMatrix_t(fromMap, rowSizes, Tpetra::StaticProfile));

    if (haveEdges){

      gno_t *edgeIds = new gno_t [nedges];
      if (nedges && !edgeIds)
        throw std::bad_alloc();
      for (int i=0; i < nedges; i++)
         edgeIds[i] = adj[i];
  
      free(adj);
      adj = NULL; 

      gno_t *nextId = edgeIds;
      Array<scalar_t> values(maxRowLen, 1.0);
  
      for (gno_t i=0; i < nvtxs; i++){
        if (nzPerRow[i] > 0){
          ArrayView<const gno_t> rowNz(nextId, nzPerRow[i]);
          fromMatrix->insertGlobalValues(i, rowNz, values.view(0,nzPerRow[i]));
          nextId += nzPerRow[i];
        }
      }
  
      delete [] edgeIds;
      edgeIds = NULL;
    }

    fromMatrix->fillComplete();
  }
  else{
    nvtxs = graphCounts[0];
    nedges = graphCounts[1];
    vwgt_dim = graphCounts[2];
    ewgt_dim  = graphCounts[3];
    maxRowLen = graphCounts[4];

    // Create a Tpetra::CrsMatrix where rank 0 has entire matrix.

    fromMap = rcp(new map_t(nvtxs, 0, base, tcomm_));

    fromMatrix = 
      rcp(new tcrsMatrix_t(fromMap, rowSizes, Tpetra::StaticProfile));

    fromMatrix->fillComplete();
  }

  // Create a Tpetra::CrsMatrix with default row distribution.

  RCP<const map_t> toMap = rcp(new map_t(nvtxs, base, tcomm_));
  RCP<tcrsMatrix_t> toMatrix = rcp(new tcrsMatrix_t(toMap, maxRowLen));

  // Import the data.

  RCP<import_t> importer = rcp(new import_t(fromMap, toMap));
  toMatrix->doImport(*fromMatrix, *importer, Tpetra::INSERT);
  toMatrix->fillComplete();

  M_ = toMatrix;

  // Vertex weights, if any

  typedef ArrayRCP<const ArrayView<const scalar_t> > arrayArray_t;

  if (vwgt_dim > 0){

    ArrayRCP<scalar_t> weightBuf;
    ArrayView<const scalar_t> *wgts = new ArrayView<const scalar_t> [vwgt_dim];

    if (rank == 0){
      size_t len = vwgt_dim * nvtxs;
      scalar_t *buf = new scalar_t [len];
      if (!buf) throw std::bad_alloc();
      weightBuf = arcp(buf, 0, len, true);

      for (int wdim=0; wdim < vwgt_dim; wdim++){
        wgts[wdim] = ArrayView<const scalar_t>(buf, nvtxs);
        float *vw = vwgts + wdim;
        for (int i=0; i < nvtxs; i++, vw += vwgt_dim)
          buf[i] = *vw;
        buf += nvtxs;
      }

      free(vwgts);
      vwgts = NULL;
    }

    arrayArray_t vweights = arcp(wgts, 0, vwgt_dim, true);

    RCP<tMVector_t> fromVertexWeights = 
      rcp(new tMVector_t(fromMap, vweights.view(0, vwgt_dim), vwgt_dim));

    RCP<tMVector_t> toVertexWeights = rcp(new tMVector_t(toMap, vwgt_dim));

    toVertexWeights->doImport(*fromVertexWeights, *importer, Tpetra::INSERT);

    vtxWeights_ = toVertexWeights;
  }

  // Edge weights, if any

  if (haveEdges && ewgt_dim > 0){

    ArrayRCP<scalar_t> weightBuf;
    ArrayView<const scalar_t> *wgts = new ArrayView<const scalar_t> [ewgt_dim];

    toMap = rcp(new map_t(nedges, M_->getNodeNumEntries(), base, tcomm_));

    if (rank == 0){
      size_t len = ewgt_dim * nedges;
      scalar_t *buf = new scalar_t [len];
      if (!buf) throw std::bad_alloc();
      weightBuf = arcp(buf, 0, len, true);

      for (int wdim=0; wdim < ewgt_dim; wdim++){
        wgts[wdim] = ArrayView<const scalar_t>(buf, nedges);
        float *ew = ewgts + wdim;
        for (int i=0; i < nedges; i++, ew += ewgt_dim)
          buf[i] = *ew;
        buf += nedges;
      }

      free(ewgts);
      ewgts = NULL;
      fromMap = rcp(new map_t(nedges, nedges, base, tcomm_));
    }
    else{
      fromMap = rcp(new map_t(nedges, 0, base, tcomm_));
    }

    arrayArray_t eweights = arcp(wgts, 0, ewgt_dim, true);

    RCP<tMVector_t> fromEdgeWeights = 
      rcp(new tMVector_t(fromMap, eweights.view(0, ewgt_dim), ewgt_dim));

    RCP<tMVector_t> toEdgeWeights = rcp(new tMVector_t(toMap, ewgt_dim));

    RCP<import_t> edgeImporter = rcp(new import_t(fromMap, toMap));

    toEdgeWeights->doImport(*fromEdgeWeights, *edgeImporter, Tpetra::INSERT);

    edgWeights_ = toEdgeWeights;
  }
}

void UserInputForTests::getChacoCoords(FILE *fptr, string fname)
{
  int rank = tcomm_->getRank();
  int ndim=0;
  float *x=NULL, *y=NULL, *z=NULL;
  int fail = 0;

  size_t localNumVtx = M_->getNodeNumRows();
  size_t globalNumVtx = M_->getGlobalNumRows();

  if (rank == 0){
  
    // This function is in the Zoltan C-library.
  
    // Reads in the file and closes it when done.
    char *nonConstName = new char [fname.size() + 1];
    strcpy(nonConstName, fname.c_str());
    fail = Zoltan2::chaco_input_geom(fptr, nonConstName, globalNumVtx,
      &ndim, &x, &y, &z);
    delete [] nonConstName;

    if (fail)
      ndim = 0;

    if (verbose_){
      std::cout << "UserInputForTests, read " << globalNumVtx;
      std::cout << " " << ndim << "-dimensional coordinates." << std::endl;
    }
  }
  
  Teuchos::broadcast<int, int>(*tcomm_, 0, 1, &ndim);

  if (ndim == 0)
    throw std::runtime_error("Can't read coordinate file");

  ArrayRCP<ArrayRCP<const scalar_t> > coords(ndim);
  lno_t len = 0;

  if (rank == 0){

    for (int dim=0; dim < ndim; dim++){
      scalar_t *v = new scalar_t [globalNumVtx];
      if (!v)
        throw std::bad_alloc();
      coords[dim] = arcp<const scalar_t>(v, 0, globalNumVtx, true);
      float *val = (dim==0 ? x : (dim==1 ? y : z));
      for (size_t i=0; i < globalNumVtx; i++)
        v[i] = scalar_t(val[i]);

      free(val);
    }

    len = globalNumVtx;;
  }

  RCP<const map_t> fromMap = rcp(new map_t(globalNumVtx, len, 0, tcomm_));
  RCP<const map_t> toMap = rcp(new map_t(globalNumVtx, localNumVtx, 0, tcomm_));
  RCP<import_t> importer = rcp(new import_t(fromMap, toMap));

  Array<ArrayView<const scalar_t> > coordData;
  for (int dim=0; dim < ndim; dim++)
    coordData.push_back(coords[dim].view(0, len));

  RCP<tMVector_t> fromCoords = 
    rcp(new tMVector_t(fromMap, coordData.view(0, ndim), ndim));

  RCP<tMVector_t> toCoords = rcp(new tMVector_t(toMap, ndim));

  toCoords->doImport(*fromCoords, *importer, Tpetra::INSERT);

  xyz_ = toCoords;

}

string UserInputForTests::shortPathName(std::ostringstream &oss, string dirName)
{
  string shorterName(oss.str());
  string slash("/");
  string token = slash + dirName + slash;
  
  size_t pos = shorterName.find(token);
  if (pos != string::npos)
    shorterName = shorterName.substr(pos+1);
  return shorterName;
}

#endif
