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

#include "Zoltan2_TestHelpers.hpp"
#include <Zoltan2_XpetraTraits.hpp>
#include <Zoltan2_Typedefs.hpp>

#include <Tpetra_MultiVector.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_Map.hpp>
#include <Xpetra_Vector.hpp>
#include <Xpetra_CrsMatrix.hpp>
#include <Xpetra_CrsGraph.hpp>

#include <MatrixMarket_Tpetra.hpp>

#ifdef HAVE_ZOLTAN2_GALERI
#include <Galeri_XpetraProblemFactory.hpp>
#include <Galeri_XpetraParameters.hpp>
#endif

#include <Kokkos_DefaultNode.hpp>

#include "GeometricGenerator.hpp"
#include <fstream>
#include <string>

#include <TpetraExt_MatrixMatrix_def.hpp>

// pamgen required includes
#include "Zoltan2_PamgenMeshStructure.hpp"


using Teuchos::RCP;
using Teuchos::ArrayRCP;
using Teuchos::ArrayView;
using Teuchos::Array;
using Teuchos::Comm;
using Teuchos::rcp;
using Teuchos::arcp;
using Teuchos::rcp_const_cast;
using Teuchos::ParameterList;
using namespace Zoltan2_TestingFramework;

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

typedef enum USERINPUT_FILE_FORMATS{MATRIX_MARKET, CHACO, GEOMGEN, PAMGEN} USERINPUT_FILE_FORMATS;

class UserInputForTests
{
public:

  typedef Tpetra::Map<zlno_t, zgno_t, znode_t> map_t;
  typedef Tpetra::Export<zlno_t, zgno_t, znode_t> export_t;
  typedef Tpetra::Import<zlno_t, zgno_t, znode_t> import_t;
  typedef map_t::node_type default_znode_t;


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
   *   \param  distributeInput if true, initial data will be distributed
   *           across processes.  Currently distributeInput=false is
   *           supported only for Zoltan input, not MatrixMarket files.
   *
   *  For example, if \c path is the path to the Zoltan1 test
   *  directory and \c testData is \c brack2_3, then we'll read
   *  in ch_brack2_3/brack2_3.graph and ch_brack2_3/brack2_3.coords.
   */

  UserInputForTests(string path, string testData,
                    const RCP<const Comm<int> > &c, bool debugInfo=false,
                    bool distributeInput=true);

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
   * See Galeri::Xpetra::BuildProblem() for more information
   * about problem types.
   */
  UserInputForTests(int x, int y, int z, string matrixType,
                    const RCP<const Comm<int> > &c, bool debugInfo=false,
                    bool distributeInput=true);

  /*! \brief Constructor that can generate a sp matrix or read a
   *      matrix/graph file defined in an input parameter list.
   *   \param pList the parameter list containing input options
   *   \param problemType the type of problem that will generate a
   *             sparse matrix from the mesh.  If the problemType is
   *             empty we'll pick a default.
   *   \param  c  is the communicator for the processes that
   *         share the data.
   *   \param  debugInfo if true process zero will print out status.
   *
   * Problems can be "Laplace1D", "Laplace2D", "Star2D", "BigStar2D",
   * "Laplace3D", "Brick3D" and "Identity".
   * See Galeri::Xpetra::BuildProblem() for more information
   * about problem types.
   */
  UserInputForTests(const ParameterList &pList,
                    const RCP<const Comm<int> > &c);

  /*! \brief Generate lists of random scalars.
   */
  static void getUIRandomData(unsigned int seed, zlno_t length,
                              zscalar_t min, zscalar_t max, ArrayView<ArrayRCP<zscalar_t > > data);

  RCP<tMVector_t> getUICoordinates();

  RCP<tMVector_t> getUIWeights();

  RCP<tMVector_t> getUIEdgeWeights();

  RCP<tcrsMatrix_t> getUITpetraCrsMatrix();

  RCP<tcrsGraph_t> getUITpetraCrsGraph();

  RCP<tVector_t> getUITpetraVector();

  RCP<tMVector_t> getUITpetraMultiVector(int nvec);

  RCP<xcrsMatrix_t> getUIXpetraCrsMatrix();

  RCP<xcrsGraph_t> getUIXpetraCrsGraph();

  RCP<xVector_t> getUIXpetraVector();

  RCP<xMVector_t> getUIXpetraMultiVector(int nvec);

#ifdef HAVE_ZOLTAN2_PAMGEN
  PamgenMesh * getPamGenMesh(){return this->pamgen_mesh.operator->();}
#endif

#ifdef HAVE_EPETRA_DATA_TYPES
  RCP<Epetra_CrsGraph> getUIEpetraCrsGraph();

  RCP<Epetra_CrsMatrix> getUIEpetraCrsMatrix();

  RCP<Epetra_Vector> getUIEpetraVector();

  RCP<Epetra_MultiVector> getUIEpetraMultiVector(int nvec);
#endif
  bool hasInput();

  bool hasInputDataType(const string &input_type);

  bool hasUICoordinates();

  bool hasUIWeights();

  bool hasUIEdgeWeights();

  bool hasUITpetraCrsMatrix();

  bool hasUITpetraCrsGraph();

  bool hasUITpetraVector();

  bool hasUITpetraMultiVector();

  bool hasUIXpetraCrsMatrix();

  bool hasUIXpetraCrsGraph();

  bool hasUIXpetraVector();

  bool hasUIXpetraMultiVector();

  bool hasPamgenMesh();
#ifdef HAVE_EPETRA_DATA_TYPES
  bool hasUIEpetraCrsGraph();

  bool hasUIEpetraCrsMatrix();

  bool hasUIEpetraVector();

  bool hasUIEpetraMultiVector();

#endif

private:

  bool verbose_;

  const RCP<const Comm<int> > tcomm_;

  bool havePamgenMesh;
#ifdef HAVE_ZOLTAN2_PAMGEN
  RCP<PamgenMesh> pamgen_mesh;
#endif

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

  void readMatrixMarketFile(string path, string testData,bool distributeInput = true);

  // Build matrix M_ from a mesh and a problem type
  // with Galeri::Xpetra.

  void buildCrsMatrix(int xdim, int ydim, int zdim, string type,
                      bool distributeInput);

  // Read a Zoltan1 Chaco or Matrix Market file
  // into M_.  If it has geometric coordinates,
  // read them into xyz_.  If it has weights,
  // read those into vtxWeights_ and edgWeights_.
  void readZoltanTestData(string path, string testData,
                          bool distributeInput);

  // Modify the Maps of an input matrix to make them non-contiguous
  RCP<tcrsMatrix_t> modifyMatrixGIDs(RCP<tcrsMatrix_t> &in);
  inline zgno_t newID(const zgno_t id) { return id * 2 + 10001; }

  // Read Zoltan data that is in a .graph file.
  void getUIChacoGraph(FILE *fptr, bool haveAssign, FILE *assignFile,
                       string name, bool distributeInput);

  // Read Zoltan data that is in a .coords file.
  void getUIChacoCoords(FILE *fptr, string name);

  //  Chaco reader code: This code is copied from zoltan/ch.
  //  It might benefit from a rewrite and simplification.

  //  Chaco reader helper functions:  copied from zoltan/ch
  static const int CHACO_LINE_LENGTH=200;
  char chaco_line[CHACO_LINE_LENGTH];    /* space to hold values */
  int chaco_offset;        /* offset into line for next data */
  int chaco_break_pnt;    /* place in sequence to pause */
  int chaco_save_pnt;        /* place in sequence to save */

  double chaco_read_val(FILE* infile, int *end_flag);
  int chaco_read_int(FILE* infile, int *end_flag);
  void chaco_flush_line(FILE*);

  // Chaco graph reader:  copied from zoltan/ch
  int chaco_input_graph(FILE *fin, const char *inname, int **start,
                        int **adjacency, int  *nvtxs, int  *nVwgts,
                        float **vweights, int  *nEwgts, float **eweights);

  // Chaco coordinate reader:  copied from zoltan/ch
  int chaco_input_geom(FILE *fingeom, const char *geomname, int nvtxs,
                       int  *igeom, double **x, double **y, double **z);

  // Chaco coordinate reader:  copied from zoltan/ch
  int chaco_input_assign(FILE *finassign, const char *assignname, int nvtxs,
                         short *assignments);


  // Read a GeomGen.txt file into M_
  // Read coordinates into xyz_.
  // If iti has weights read those to vtxWeights_
  // and edgeWeights_
  void readGeometricGenTestData(string path, string testData);

  // Geometry Gnearatory helper function
  void readGeoGenParams(string paramFileName,
                        ParameterList &geoparams);

  // utility methods used when reading geom gen files

  static string trim_right_copy(const string& s,
                                const string& delimiters = " \f\n\r\t\v" );

  static string trim_left_copy(const string& s,
                               const string& delimiters = " \f\n\r\t\v" );

  static string trim_copy(const string& s,
                          const string& delimiters = " \f\n\r\t\v" );


  // Read a pamgen mesh
  void readPamgenMeshFile(string path, string testData);
#ifdef HAVE_ZOLTAN2_PAMGEN
  void setPamgenAdjacencyGraph();
  void setPamgenCoordinateMV();
#endif
};

UserInputForTests::UserInputForTests(string path, string testData,
                                     const RCP<const Comm<int> > &c,
                                     bool debugInfo, bool distributeInput):
verbose_(debugInfo), tcomm_(c), havePamgenMesh(false),
M_(), xM_(), xyz_(), vtxWeights_(), edgWeights_(),
#ifdef HAVE_EPETRA_DATA_TYPES
ecomm_(), eM_(), eG_(),
#endif
chaco_offset(0), chaco_break_pnt(CHACO_LINE_LENGTH)
{
  bool zoltan1 = false;
  string::size_type loc = path.find("/zoltan/test/");  // Zoltan1 data
  if (loc != string::npos)
    zoltan1 = true;

  if (zoltan1)
    readZoltanTestData(path, testData, distributeInput);
  else
    readMatrixMarketFile(path, testData);

#ifdef HAVE_EPETRA_DATA_TYPES
  ecomm_ = Xpetra::toEpetra(c);
#endif
}

UserInputForTests::UserInputForTests(int x, int y, int z,
                                     string matrixType,
                                     const RCP<const Comm<int> > &c,
                                     bool debugInfo,
                                     bool distributeInput):
verbose_(debugInfo), tcomm_(c), havePamgenMesh(false),
M_(), xM_(), xyz_(), vtxWeights_(), edgWeights_(),
#ifdef HAVE_EPETRA_DATA_TYPES
ecomm_(), eM_(), eG_(),
#endif
chaco_offset(0), chaco_break_pnt(CHACO_LINE_LENGTH)
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

  buildCrsMatrix(x, y, z, matrixType, distributeInput);

#ifdef HAVE_EPETRA_DATA_TYPES
  ecomm_ = Xpetra::toEpetra(c);
#endif
}

UserInputForTests::UserInputForTests(const ParameterList &pList,
                                     const RCP<const Comm<int> > &c):
tcomm_(c), havePamgenMesh(false),
M_(), xM_(), xyz_(), vtxWeights_(), edgWeights_(),
#ifdef HAVE_EPETRA_DATA_TYPES
ecomm_(), eM_(), eG_(),
#endif
chaco_offset(0), chaco_break_pnt(CHACO_LINE_LENGTH)
{

  // get options
  bool distributeInput = true, debugInfo = true;

  if(pList.isParameter("distribute input"))
    distributeInput = pList.get<bool>("distribute input");

  if(pList.isParameter("debug"))
    debugInfo = pList.get<bool>("debug");
  this->verbose_ = debugInfo;

  if(pList.isParameter("input file"))
  {

    // get input path
    string path(".");
    if(pList.isParameter("input path"))
      path = pList.get<string>("input path");

    string testData = pList.get<string>("input file");

    // find out if we are working from the zoltan1 test diretory
    USERINPUT_FILE_FORMATS file_format = MATRIX_MARKET;

    // find out if we are using the geometric generator
    if(pList.isParameter("file type") && pList.get<string>("file type") == "Geometric Generator")
      file_format = GEOMGEN;
    else if(pList.isParameter("file type") && pList.get<string>("file type") == "Pamgen")
    {
      file_format = PAMGEN;
    }
    else if(pList.isParameter("file type") && pList.get<string>("file type") == "Chaco")
      file_format = CHACO; // this flag calls read ZoltanTestData, which calls the chaco readers...

    // read the input file
    switch (file_format) {
      case GEOMGEN: readGeometricGenTestData(path,testData); break;
      case PAMGEN: readPamgenMeshFile(path,testData); break;
      case CHACO: readZoltanTestData(path, testData, distributeInput); break;
      default: readMatrixMarketFile(path, testData, distributeInput); break;
    }

  }else if(pList.isParameter("x") || pList.isParameter("y") || pList.isParameter("z")){

    int x,y,z;
    x = y = z = 0;
    if(pList.isParameter("x")) x = pList.get<int>("x");
    if(pList.isParameter("y")) y = pList.get<int>("y");
    if(pList.isParameter("z")) z = pList.get<int>("z");

    string problemType = "";
    if(pList.isParameter("equation type")) problemType = pList.get<string>("equation type");

    if (problemType.size() == 0){ /** default behavior */
      int dim = 0;
      if (x > 0) dim++;
      if (y > 0) dim++;
      if (z > 0) dim++;
      if (dim == 1)
        problemType = string("Laplace1D");
      else if (dim == 2)
        problemType = string("Laplace2D");
      else if (dim == 3)
        problemType = string("Laplace3D");
      else
        throw std::runtime_error("input");

      if (verbose_ && tcomm_->getRank() == 0)
        std::cout << "UserInputForTests, Matrix type : " << problemType << std::endl;
    }


    buildCrsMatrix(x, y, z, problemType, distributeInput);

  }else{
    std::cerr << "Input file block undefined!" << std::endl;
  }

#ifdef HAVE_EPETRA_DATA_TYPES
  ecomm_ = Xpetra::toEpetra(c);
#endif

}


RCP<Zoltan2_TestingFramework::tMVector_t> UserInputForTests::getUICoordinates()
{
  if (xyz_.is_null())
    throw std::runtime_error("could not read coord file");
  return xyz_;
}

RCP<Zoltan2_TestingFramework::tMVector_t> UserInputForTests::getUIWeights()
{
  return vtxWeights_;
}

RCP<Zoltan2_TestingFramework::tMVector_t> UserInputForTests::getUIEdgeWeights()
{
  return edgWeights_;
}

RCP<Zoltan2_TestingFramework::tcrsMatrix_t> UserInputForTests::getUITpetraCrsMatrix()
{
  if (M_.is_null())
    throw std::runtime_error("could not read mtx file");
  return M_;
}

RCP<Zoltan2_TestingFramework::tcrsGraph_t> UserInputForTests::getUITpetraCrsGraph()
{
  if (M_.is_null())
    throw std::runtime_error("could not read mtx file");
  return rcp_const_cast<tcrsGraph_t>(M_->getCrsGraph());
}

RCP<Zoltan2_TestingFramework::tVector_t> UserInputForTests::getUITpetraVector()
{
  RCP<tVector_t> V = rcp(new tVector_t(M_->getRowMap(),  1));
  V->randomize();

  return V;
}

RCP<Zoltan2_TestingFramework::tMVector_t> UserInputForTests::getUITpetraMultiVector(int nvec)
{
  RCP<tMVector_t> mV = rcp(new tMVector_t(M_->getRowMap(), nvec));
  mV->randomize();

  return mV;
}

RCP<Zoltan2_TestingFramework::xcrsMatrix_t> UserInputForTests::getUIXpetraCrsMatrix()
{
  if (M_.is_null())
    throw std::runtime_error("could not read mtx file");
  return xM_;
}

RCP<Zoltan2_TestingFramework::xcrsGraph_t> UserInputForTests::getUIXpetraCrsGraph()
{
  if (M_.is_null())
    throw std::runtime_error("could not read mtx file");
  return rcp_const_cast<xcrsGraph_t>(xM_->getCrsGraph());
}

RCP<Zoltan2_TestingFramework::xVector_t> UserInputForTests::getUIXpetraVector()
{
  return Zoltan2::XpetraTraits<tVector_t>::convertToXpetra(getUITpetraVector());
}

RCP<Zoltan2_TestingFramework::xMVector_t> UserInputForTests::getUIXpetraMultiVector(int nvec)
{
  RCP<tMVector_t> tMV = getUITpetraMultiVector(nvec);
  return Zoltan2::XpetraTraits<tMVector_t>::convertToXpetra(tMV);
}

#ifdef HAVE_EPETRA_DATA_TYPES
RCP<Epetra_CrsGraph> UserInputForTests::getUIEpetraCrsGraph()
{
  if (M_.is_null())
    throw std::runtime_error("could not read mtx file");
  RCP<const tcrsGraph_t> tgraph = M_->getCrsGraph();
  RCP<const Tpetra::Map<zlno_t, zgno_t> > trowMap = tgraph->getRowMap();
  RCP<const Tpetra::Map<zlno_t, zgno_t> > tcolMap = tgraph->getColMap();

  int nElts = static_cast<int>(trowMap->getGlobalNumElements());
  int nMyElts = static_cast<int>(trowMap->getNodeNumElements());
  int base = 0;
  ArrayView<const int> gids = trowMap->getNodeElementList();

  Epetra_BlockMap erowMap(nElts, nMyElts,
                          gids.getRawPtr(), 1, base, *ecomm_);

  Array<int> rowSize(nMyElts);
  for (int i=0; i < nMyElts; i++){
    rowSize[i] = static_cast<int>(M_->getNumEntriesInLocalRow(i));
  }

  size_t maxRow = M_->getNodeMaxNumRowEntries();
  Array<int> colGids(maxRow);
  ArrayView<const int> colLid;

  eG_ = rcp(new Epetra_CrsGraph(Copy, erowMap,
                                rowSize.getRawPtr(), true));

  for (int i=0; i < nMyElts; i++){
    tgraph->getLocalRowView(i, colLid);
    for (int j=0; j < colLid.size(); j++)
      colGids[j] = tcolMap->getGlobalElement(colLid[j]);
    eG_->InsertGlobalIndices(gids[i], rowSize[i], colGids.getRawPtr());
  }
  eG_->FillComplete();
  return eG_;
}

RCP<Epetra_CrsMatrix> UserInputForTests::getUIEpetraCrsMatrix()
{
  if (M_.is_null())
    throw std::runtime_error("could not read mtx file");
  RCP<Epetra_CrsGraph> egraph = getUIEpetraCrsGraph();
  eM_ = rcp(new Epetra_CrsMatrix(Copy, *egraph));

  size_t maxRow = M_->getNodeMaxNumRowEntries();
  int nrows = egraph->NumMyRows();
  const Epetra_BlockMap &rowMap = egraph->RowMap();
  const Epetra_BlockMap &colMap = egraph->ColMap();
  Array<int> colGid(maxRow);

  for (int i=0; i < nrows; i++){
    ArrayView<const int> colLid;
    ArrayView<const zscalar_t> nz;
    M_->getLocalRowView(i, colLid, nz);
    size_t rowSize = colLid.size();
    int rowGid = rowMap.GID(i);
    for (size_t j=0; j < rowSize; j++){
      colGid[j] = colMap.GID(colLid[j]);
    }
    eM_->InsertGlobalValues(rowGid, (int)rowSize, nz.getRawPtr(), colGid.getRawPtr());
  }
  eM_->FillComplete();
  return eM_;
}

RCP<Epetra_Vector> UserInputForTests::getUIEpetraVector()
{
  RCP<Epetra_CrsGraph> egraph = getUIEpetraCrsGraph();
  RCP<Epetra_Vector> V = rcp(new Epetra_Vector(egraph->RowMap()));
  V->Random();
  return V;
}

RCP<Epetra_MultiVector> UserInputForTests::getUIEpetraMultiVector(int nvec)
{
  RCP<Epetra_CrsGraph> egraph = getUIEpetraCrsGraph();
  RCP<Epetra_MultiVector> mV =
  rcp(new Epetra_MultiVector(egraph->RowMap(), nvec));
  mV->Random();
  return mV;
}
#endif

bool UserInputForTests::hasInput()
{
  // find out if an input source has been loaded
  return  this->hasUICoordinates() || \
          this->hasUITpetraCrsMatrix() || \
          this->hasUITpetraCrsGraph() || \
          this->hasPamgenMesh();
}

bool UserInputForTests::hasInputDataType(const string &input_type)
{
  if(input_type == "coordinates")
    return this->hasUICoordinates();
  else if(input_type == "tpetra_vector")
    return this->hasUITpetraVector();
  else if(input_type == "tpetra_multivector")
    return this->hasUITpetraMultiVector();
  else if(input_type == "tpetra_crs_graph")
    return this->hasUITpetraCrsGraph();
  else if(input_type == "tpetra_crs_matrix")
    return this->hasUITpetraCrsMatrix();
  else if(input_type == "xpetra_vector")
    return this->hasUIXpetraVector();
  else if(input_type == "xpetra_multivector")
    return this->hasUIXpetraMultiVector();
  else if(input_type == "xpetra_crs_graph")
    return this->hasUIXpetraCrsGraph();
  else if(input_type == "xpetra_crs_matrix")
    return this->hasUIXpetraCrsMatrix();
#ifdef HAVE_EPETRA_DATA_TYPES
  else if(input_type == "epetra_vector")
    return this->hasUIEpetraVector();
  else if(input_type == "epetra_multivector")
    return this->hasUIEpetraMultiVector();
  else if(input_type == "epetra_crs_graph")
    return this->hasUIEpetraCrsGraph();
  else if(input_type == "epetra_crs_matrix")
    return this->hasUIEpetraCrsMatrix();
#endif

  return false;
}

bool UserInputForTests::hasUICoordinates()
{
  return xyz_.is_null() ? false : true;
}

bool UserInputForTests::hasUIWeights()
{
  return vtxWeights_.is_null() ? false : true;
}

bool UserInputForTests::hasUIEdgeWeights()
{
  return edgWeights_.is_null() ? false : true;
}

bool UserInputForTests::hasUITpetraCrsMatrix()
{
  return M_.is_null() ? false : true;
}

bool UserInputForTests::hasUITpetraCrsGraph()
{
  return M_.is_null() ? false : true;
}

bool UserInputForTests::hasUITpetraVector()
{
  return true;
}

bool UserInputForTests::hasUITpetraMultiVector()
{
  return true;
}

bool UserInputForTests::hasUIXpetraCrsMatrix()
{
  return M_.is_null() ? false : true;
}

bool UserInputForTests::hasUIXpetraCrsGraph()
{
  return M_.is_null() ? false : true;
}

bool UserInputForTests::hasUIXpetraVector()
{
  return true;
}

bool UserInputForTests::hasUIXpetraMultiVector()
{
  return true;
}

bool UserInputForTests::hasPamgenMesh()
{
  return this->havePamgenMesh;
}

#ifdef HAVE_EPETRA_DATA_TYPES
bool UserInputForTests::hasUIEpetraCrsGraph()
{
  return M_.is_null() ? false : true;
}

bool UserInputForTests::hasUIEpetraCrsMatrix()
{
  return hasUIEpetraCrsGraph();
}

bool UserInputForTests::hasUIEpetraVector()
{
  return hasUIEpetraCrsGraph();
}

bool UserInputForTests::hasUIEpetraMultiVector()
{
  return hasUIEpetraCrsGraph();
}
#endif

void UserInputForTests::getUIRandomData(unsigned int seed, zlno_t length,
                                        zscalar_t min, zscalar_t max,
                                        ArrayView<ArrayRCP<zscalar_t > > data)
{
  if (length < 1)
    return;

  size_t dim = data.size();
  for (size_t i=0; i < dim; i++){
    zscalar_t *tmp = new zscalar_t [length];
    if (!tmp)
      throw (std::bad_alloc());
    data[i] = Teuchos::arcp(tmp, 0, length, true);
  }

  zscalar_t scalingFactor = (max-min) / RAND_MAX;
  srand(seed);
  for (size_t i=0; i < dim; i++){
    zscalar_t *x = data[i].getRawPtr();
    for (zlno_t j=0; j < length; j++)
      *x++ = min + (zscalar_t(rand()) * scalingFactor);
  }
}

// utility methods used when reading geom gen files

string UserInputForTests::trim_right_copy(
                                          const string& s,
                                          const string& delimiters)
{
  return s.substr( 0, s.find_last_not_of( delimiters ) + 1 );
}

string UserInputForTests::trim_left_copy(
                                         const string& s,
                                         const string& delimiters)
{
  return s.substr( s.find_first_not_of( delimiters ) );
}

string UserInputForTests::trim_copy(
                                    const string& s,
                                    const string& delimiters)
{
  return trim_left_copy( trim_right_copy( s, delimiters ), delimiters );
}

void UserInputForTests::readGeometricGenTestData(string path,
                                                 string testData)
{

  std::ostringstream fname;
  fname << path << "/" << testData << ".txt";

  if (verbose_ && tcomm_->getRank() == 0)
    std::cout << "UserInputForTests, Read: " << fname.str() << std::endl;

  Teuchos::ParameterList geoparams("geo params");
  readGeoGenParams(fname.str(),geoparams);

  geometricgen_t * gg = new geometricgen_t(geoparams, this->tcomm_);

  // get coordinate and point info
  int coord_dim = gg->getCoordinateDimension();
  int numWeightsPerCoord = gg->getNumWeights();
  zlno_t numLocalPoints = gg->getNumLocalCoords();
  zgno_t numGlobalPoints = gg->getNumGlobalCoords();

  // allocate an array of coordinate arrays
  zscalar_t **coords = new zscalar_t * [coord_dim];
  for(int i = 0; i < coord_dim; ++i){
    coords[i] = new zscalar_t[numLocalPoints];
  }

  // get a copy of the data
  gg->getLocalCoordinatesCopy(coords);

  // get an array of arrays of weight data (if any)
  zscalar_t **weight = NULL;
  if (numWeightsPerCoord) {
    // memory allocation
    weight = new zscalar_t * [numWeightsPerCoord];
    for(int i = 0; i < numWeightsPerCoord; ++i){
      weight[i] = new zscalar_t[numLocalPoints];
    }

    // get a copy of the weight data
    gg->getLocalWeightsCopy(weight);
  }

  delete gg; // free up memory from geom gen


  // make a Tpetra map
  RCP<Tpetra::Map<zlno_t, zgno_t, znode_t> > mp =
  rcp(new Tpetra::Map<zlno_t, zgno_t, znode_t>(numGlobalPoints, numLocalPoints, 0, this->tcomm_));

  // make an array of array views containing the coordinate data
  Teuchos::Array<Teuchos::ArrayView<const zscalar_t> > coordView(coord_dim);
  for (int i=0; i < coord_dim; i++){
    if(numLocalPoints > 0){
      Teuchos::ArrayView<const zscalar_t> a(coords[i], numLocalPoints);
      coordView[i] = a;
    }
    else {
      Teuchos::ArrayView<const zscalar_t> a;
      coordView[i] = a;
    }
  }

  // set the xyz_ multivector
  xyz_ = RCP<tMVector_t>(new
                         tMVector_t(mp, coordView.view(0, coord_dim),
                                    coord_dim));

  // set the vtx weights
  if (numWeightsPerCoord) {
    // make an array of array views containing the weight data
    Teuchos::Array<Teuchos::ArrayView<const zscalar_t> > weightView(numWeightsPerCoord);
    for (int i=0; i < numWeightsPerCoord; i++){
      if(numLocalPoints > 0){
        Teuchos::ArrayView<const zscalar_t> a(weight[i], numLocalPoints);
        weightView[i] = a;
      }
      else {
        Teuchos::ArrayView<const zscalar_t> a;
        weightView[i] = a;
      }
    }

    vtxWeights_ = RCP<tMVector_t>(new tMVector_t(mp, weightView.view(0, numWeightsPerCoord),
                                                 numWeightsPerCoord));
  }
}

void UserInputForTests::readGeoGenParams(string paramFileName,
                                         ParameterList &geoparams){

  const char param_comment = '#';

  std::string input = "";
  char inp[25000];
  for(int i = 0; i < 25000; ++i){
    inp[i] = 0;
  }

  bool fail = false;
  if(this->tcomm_->getRank() == 0){

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



  int size = (int)input.size();
  if(fail){
    size = -1;
  }
  this->tcomm_->broadcast(0, sizeof(int), (char*) &size);
  if(size == -1){
    throw "File " + paramFileName + " cannot be opened.";
  }
  this->tcomm_->broadcast(0, size, inp);
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

/////////////////////////////////////////////////////////////////////////////
RCP<tcrsMatrix_t> UserInputForTests::modifyMatrixGIDs(
  RCP<tcrsMatrix_t> &inMatrix
)
{
  // Produce a new matrix with the same structure as inMatrix,
  // but whose row/column GIDs are non-contiguous values.
  // In this case GID g in inMatrix becomes g*2+1 in outMatrix.

  // Create the map for the new matrix: same structure as inMap but with
  // the GIDs modified.
  RCP<const map_t> inMap = inMatrix->getRowMap();

  size_t nRows = inMap->getNodeNumElements();
  auto inRows = inMap->getMyGlobalIndices();
  Teuchos::Array<zgno_t> outRows(nRows);
  for (size_t i = 0; i < nRows; i++) {
    outRows[i] = newID(inRows[i]);
  }

  Tpetra::global_size_t nGlobalRows = inMap->getGlobalNumElements();
  RCP<map_t> outMap = rcp(new map_t(nGlobalRows, outRows(), 0,
                                    inMap->getComm()));

#ifdef INCLUDE_LENGTHY_OUTPUT
  // Sanity check output
  {
    std::cout << inMap->getComm()->getRank() << " KDDKDD "
              << "nGlobal " << inMap->getGlobalNumElements() << " "
                            << outMap->getGlobalNumElements() << "; "
              << "nLocal  " << inMap->getNodeNumElements() << " "
                            << outMap->getNodeNumElements() << "; "
              << std::endl;
    std::cout << inMap->getComm()->getRank() << " KDDKDD ";
    for (size_t i = 0; i < nRows; i++)
      std::cout << "(" << inMap->getMyGlobalIndices()[i] << ", "
                << outMap->getMyGlobalIndices()[i] << ") ";
    std::cout << std::endl;
  }
#endif // INCLUDE_LENGTHY_OUTPUT

  // Create a new matrix using the new map
  // Get the length of the longest row; allocate memory.
  size_t rowLen = inMatrix->getNodeMaxNumRowEntries();
  RCP<tcrsMatrix_t> outMatrix = rcp(new tcrsMatrix_t(outMap, rowLen));

  Teuchos::Array<zgno_t> indices(rowLen);
  Teuchos::Array<zscalar_t> values(rowLen);

  for (size_t i = 0; i < nRows; i++) {
    size_t nEntries;
    zgno_t inGid = inMap->getGlobalElement(i);
    inMatrix->getGlobalRowCopy(inGid, indices, values, nEntries);
    for (size_t j = 0; j < nEntries; j++)
      indices[j] = newID(indices[j]);

    zgno_t outGid = outMap->getGlobalElement(i);
    outMatrix->insertGlobalValues(outGid, indices(0, nEntries),
                                          values(0, nEntries));
  }
  outMatrix->fillComplete();

#ifdef INCLUDE_LENGTHY_OUTPUT
  // Sanity check output
  {
    std::cout << inMap->getComm()->getRank() << " KDDKDD Rows "
              << "nGlobal " << inMatrix->getGlobalNumRows() << " "
                            << outMatrix->getGlobalNumRows() << "; "
              << "nLocal  " << inMatrix->getNodeNumRows() << " "
                            << outMatrix->getNodeNumRows() << std::endl;
    std::cout << inMap->getComm()->getRank() << " KDDKDD NNZS "
              << "nGlobal " << inMatrix->getGlobalNumEntries() << " "
                            << outMatrix->getGlobalNumEntries() << "; "
              << "nLocal  " << inMatrix->getNodeNumEntries() << " "
                            << outMatrix->getNodeNumEntries() << std::endl;

    size_t nIn, nOut;
    Teuchos::Array<zgno_t> in(rowLen), out(rowLen);
    Teuchos::Array<zscalar_t> inval(rowLen), outval(rowLen);

    for (size_t i = 0; i < nRows; i++) {
      std::cout << inMap->getComm()->getRank() << " KDDKDD " << i << " nnz(";
      inMatrix->getGlobalRowCopy(inMap->getGlobalElement(i), in, inval, nIn);
      outMatrix->getGlobalRowCopy(outMap->getGlobalElement(i), out, outval,
                                  nOut);

      std::cout << nIn << ", " << nOut << "): ";
      for (size_t j = 0; j < nIn; j++) {
        std::cout << "(" << in[j] << " " << inval[j] << ", "
                  << out[j] << " " << outval[j] << ") ";
      }
      std::cout << std::endl;
    }
  }
#endif // INCLUDE_LENGTHY_OUTPUT

  return outMatrix;
}

/////////////////////////////////////////////////////////////////////////////

void UserInputForTests::readMatrixMarketFile(
  string path,
  string testData,
  bool distributeInput
)
{
  std::ostringstream fname;
  fname << path << "/" << testData << ".mtx";

  if (verbose_ && tcomm_->getRank() == 0)
    std::cout << "UserInputForTests, Read: " << fname.str() << std::endl;

  // FIXME (mfh 01 Aug 2016) Tpetra::MatrixMarket::Reader has a graph
  // ("pattern" matrix) reader.  Call its readSparseGraphFile method.

  RCP<tcrsMatrix_t> toMatrix;
  RCP<tcrsMatrix_t> fromMatrix;
  bool aok = true;
  try{
    typedef Tpetra::MatrixMarket::Reader<tcrsMatrix_t> reader_type;
    fromMatrix = reader_type::readSparseFile(fname.str(), tcomm_,
                                             true, false, false);
#ifdef KDD_NOT_READY_YET
    // See note below about modifying coordinate IDs as well.
    //if (makeNonContiguous)
      fromMatrix = modifyMatrixGIDs(fromMatrix);
#endif

    if(!distributeInput)
    {
      if (verbose_ && tcomm_->getRank() == 0)
        std::cout << "Constructing serial distribution of matrix" <<  std::endl;
      // need to make a serial map and then import the data to redistribute it
      RCP<const map_t> fromMap = fromMatrix->getRowMap();

      size_t numGlobalCoords = fromMap->getGlobalNumElements();
      size_t numLocalCoords = this->tcomm_->getRank() == 0 ? numGlobalCoords : 0;
      RCP<const map_t> toMap = rcp(new map_t(numGlobalCoords,numLocalCoords, 0, tcomm_));

      RCP<import_t> importer = rcp(new import_t(fromMap, toMap));
      toMatrix = rcp(new tcrsMatrix_t(toMap,0));
      toMatrix->doImport(*fromMatrix, *importer, Tpetra::INSERT);
      toMatrix->fillComplete();

    }else{
      toMatrix = fromMatrix;
    }
  }catch (std::exception &e) {
    if (tcomm_->getRank() == 0) {
      std::cout << "UserInputForTests unable to read matrix market file:"
                << fname.str() << std::endl;
      std::cout << e.what() << std::endl;
    }
    aok = false;
  }
  TEST_FAIL_AND_THROW(*tcomm_, aok,
                      "UserInputForTests unable to read matrix market file");

  M_ = toMatrix;
#ifdef INCLUDE_LENGTHY_OUTPUT
  std::cout << tcomm_->getRank() << " KDDKDD " << M_->getNodeNumRows()
            << " " << M_->getGlobalNumRows()
            << " " << M_->getNodeNumEntries()
            << " " << M_->getGlobalNumEntries() << std::endl;
#endif // INCLUDE_LENGTHY_OUTPUT

  xM_ = Zoltan2::XpetraTraits<tcrsMatrix_t>::convertToXpetra(M_);

  // Open the coordinate file.

  fname.str("");
  fname << path << "/" << testData << "_coord.mtx";

  size_t coordDim = 0, numGlobalCoords = 0;
  size_t msg[2]={0,0};
  ArrayRCP<ArrayRCP<zscalar_t> > xyz;
  std::ifstream coordFile;

  if (tcomm_->getRank() == 0){

    if (verbose_)
      std::cout << "UserInputForTests, Read: " <<
      fname.str() << std::endl;

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

        xyz = Teuchos::arcp(new ArrayRCP<zscalar_t> [coordDim], 0, coordDim);

        for (size_t dim=0; !fail && dim < coordDim; dim++){
          size_t idx;
          zscalar_t *tmp = new zscalar_t [numGlobalCoords];
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
          ArrayRCP<zscalar_t> emptyArray;
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
  numGlobalCoords = msg[1];

  if (coordDim == 0)
    return;

  zgno_t base = 0;
  RCP<const map_t> toMap;

  if (!M_.is_null()){
    const RCP<const map_t> &mapM = M_->getRowMap();
    toMap = mapM;
  }
  else{
    if (verbose_ && tcomm_->getRank() == 0)
    {
      std::cout << "Matrix was null. ";
      std::cout << "Constructing distribution map for coordinate vector."
                <<  std::endl;
    }

    if(!distributeInput)
    {
      if (verbose_ && tcomm_->getRank() == 0)
        std::cout << "Constructing serial distribution map for coordinates."
                  <<  std::endl;

      size_t numLocalCoords = this->tcomm_->getRank()==0 ? numGlobalCoords : 0;
      toMap = rcp(new map_t(numGlobalCoords,numLocalCoords, base, tcomm_));
    }else{
      toMap = rcp(new map_t(numGlobalCoords, base, tcomm_));
    }
  }

  // Export coordinates to their owners

  xyz_ = rcp(new tMVector_t(toMap, coordDim));

  ArrayRCP<ArrayView<const zscalar_t> > coordLists(coordDim);

  if (tcomm_->getRank() == 0){

    for (size_t dim=0; dim < coordDim; dim++)
      coordLists[dim] = xyz[dim].view(0, numGlobalCoords);

    zgno_t *tmp = new zgno_t [numGlobalCoords];
    if (!tmp)
      throw std::bad_alloc();

    ArrayRCP<const zgno_t> rowIds = Teuchos::arcp(tmp, 0, numGlobalCoords);

#ifdef KDD_NOT_READY_YET
    // TODO if modifyMatrixGIDs, we need to modify ids here as well
    for (zgno_t id=0; id < zgno_t(numGlobalCoords); id++)
      *tmp++ = newID(id);
#else
    for (zgno_t id=0; id < zgno_t(numGlobalCoords); id++)
      *tmp++ = id;
#endif

    RCP<const map_t> fromMap = rcp(new map_t(numGlobalCoords,
                                             rowIds.view(0, numGlobalCoords),
                                             base, tcomm_));

    tMVector_t allCoords(fromMap, coordLists.view(0, coordDim), coordDim);

    export_t exporter(fromMap, toMap);

    xyz_->doExport(allCoords, exporter, Tpetra::INSERT);
  }
  else{

    RCP<const map_t> fromMap = rcp(new map_t(numGlobalCoords,
                                             ArrayView<zgno_t>(), base, tcomm_));

    tMVector_t allCoords(fromMap, coordLists.view(0, coordDim), coordDim);

    export_t exporter(fromMap, toMap);

    xyz_->doExport(allCoords, exporter, Tpetra::INSERT);
  }
}

void UserInputForTests::buildCrsMatrix(int xdim, int ydim, int zdim,
                                       string problemType, bool distributeInput)
{
#ifdef HAVE_ZOLTAN2_GALERI
  Teuchos::CommandLineProcessor tclp;
  Galeri::Xpetra::Parameters<zgno_t> params(tclp,
                                            xdim, ydim, zdim, problemType);

  RCP<const Tpetra::Map<zlno_t, zgno_t> > map;
  if (distributeInput)
    map = rcp(new Tpetra::Map<zlno_t, zgno_t>(params.GetNumGlobalElements(),
                                              0, tcomm_));
  else {
    // All data initially on rank 0
    size_t nGlobalElements = params.GetNumGlobalElements();
    size_t nLocalElements = ((tcomm_->getRank() == 0) ? nGlobalElements : 0);
    map = rcp(new Tpetra::Map<zlno_t, zgno_t>(nGlobalElements, nLocalElements, 0,
                                              tcomm_));
  }

  if (verbose_ && tcomm_->getRank() == 0){

    std::cout << "Matrix is " << (distributeInput ? "" : "not");
    std::cout << "distributed." << std::endl;

    std::cout << "UserInputForTests, Create matrix with " << problemType;
    std::cout << " (and " << xdim;
    if (zdim > 0)
      std::cout << " x " << ydim << " x " << zdim;
    else if (ydim > 0)
      std::cout << " x"  << ydim << " x 1";
    else
      std::cout << "x 1 x 1";

    std::cout << " mesh)" << std::endl;

  }

  bool aok = true;
  try{
    RCP<Galeri::Xpetra::Problem<Tpetra::Map<zlno_t, zgno_t>, Tpetra::CrsMatrix<zscalar_t, zlno_t, zgno_t>, Tpetra::MultiVector<zscalar_t, zlno_t, zgno_t> > > Pr =
    Galeri::Xpetra::BuildProblem<zscalar_t, zlno_t, zgno_t, Tpetra::Map<zlno_t, zgno_t>, Tpetra::CrsMatrix<zscalar_t, zlno_t, zgno_t>, Tpetra::MultiVector<zscalar_t, zlno_t, zgno_t> >
    (params.GetMatrixType(), map, params.GetParameterList());
    M_ = Pr->BuildMatrix();
  }
  catch (std::exception &e) {    // Probably not enough memory
    aok = false;
  }
  TEST_FAIL_AND_THROW(*tcomm_, aok,
                      "UserInputForTests Galeri::Xpetra::BuildProblem failed");

  xM_ = Zoltan2::XpetraTraits<tcrsMatrix_t>::convertToXpetra(M_);

  // Compute the coordinates for the matrix rows.

  if (verbose_ && tcomm_->getRank() == 0)
    std::cout <<
    "UserInputForTests, Implied matrix row coordinates computed" <<
    std::endl;

  ArrayView<const zgno_t> gids = map->getNodeElementList();
  zlno_t count = static_cast<zlno_t>(gids.size());
  int dim = 3;
  size_t pos = problemType.find("2D");
  if (pos != string::npos)
    dim = 2;
  else if (problemType == string("Laplace1D") ||
           problemType == string("Identity"))
    dim = 1;

  Array<ArrayRCP<zscalar_t> > coordinates(dim);

  if (count > 0){
    for (int i=0; i < dim; i++){
      zscalar_t *c = new zscalar_t [count];
      if (!c)
        throw(std::bad_alloc());
      coordinates[i] = Teuchos::arcp(c, 0, count, true);
    }

    if (dim==3){
      zscalar_t *x = coordinates[0].getRawPtr();
      zscalar_t *y = coordinates[1].getRawPtr();
      zscalar_t *z = coordinates[2].getRawPtr();
      zgno_t xySize = xdim * ydim;
      for (zlno_t i=0; i < count; i++){
        zgno_t iz = gids[i] / xySize;
        zgno_t xy = gids[i] - iz*xySize;
        z[i] = zscalar_t(iz);
        y[i] = zscalar_t(xy / xdim);
        x[i] = zscalar_t(xy % xdim);
      }
    }
    else if (dim==2){
      zscalar_t *x = coordinates[0].getRawPtr();
      zscalar_t *y = coordinates[1].getRawPtr();
      for (zlno_t i=0; i < count; i++){
        y[i] = zscalar_t(gids[i] / xdim);
        x[i] = zscalar_t(gids[i] % xdim);
      }
    }
    else{
      zscalar_t *x = coordinates[0].getRawPtr();
      for (zlno_t i=0; i < count; i++)
        x[i] = zscalar_t(gids[i]);
    }
  }

  Array<ArrayView<const zscalar_t> > coordView(dim);
  if (count > 0)
    for (int i=0; i < dim; i++)
      coordView[i] = coordinates[i].view(0,count);

  xyz_ = rcp(new tMVector_t(map, coordView.view(0, dim), dim));
#else
  throw std::runtime_error("Galeri input requested but Trilinos is "
                           "not built with Galeri.");
#endif
}

void UserInputForTests::readZoltanTestData(string path, string testData,
                                           bool distributeInput)
{
  int rank = tcomm_->getRank();
  FILE *graphFile = NULL;
  FILE *coordFile = NULL;
  FILE *assignFile = NULL;
  int fileInfo[3];

  for (int i = 0; i < CHACO_LINE_LENGTH; i++) chaco_line[i] = '\0';

  if (rank == 0){
    // set chacho graph file name
    std::ostringstream chGraphFileName;
    chGraphFileName << path << "/" << testData << ".graph";

    // set chaco graph
    std::ostringstream chCoordFileName;
    chCoordFileName << path << "/" << testData << ".coords";

    // set chaco graph
    std::ostringstream chAssignFileName;
    chAssignFileName << path << "/" << testData << ".assign";

    // open file
    graphFile = fopen(chGraphFileName.str().c_str(), "r");

    if(!graphFile) // maybe the user is using the default zoltan1 path convention
    {
      chGraphFileName.str("");
      chCoordFileName.str("");
      // try constructing zoltan1 paths
      chGraphFileName << path << "/ch_" << testData << "/" << testData << ".graph";
      chCoordFileName << path << "/ch_" << testData << "/" << testData << ".coords";
      chAssignFileName << path << "/ch_" << testData << "/" << testData << ".assign";
      // try to open the graph file again, if this doesn't open
      // the user has not provided a valid path to the file
      graphFile = fopen(chGraphFileName.str().c_str(), "r");
    }

    memset(fileInfo, 0, sizeof(int) * 3); // set fileinfo to 0's
    if (graphFile){
      fileInfo[0] = 1;
      if (verbose_ && tcomm_->getRank() == 0)
        std::cout << "UserInputForTests, open " <<
        chGraphFileName.str () << std::endl;

      coordFile = fopen(chCoordFileName.str().c_str(), "r");
      if (coordFile){
        fileInfo[1] = 1;
        if (verbose_ && tcomm_->getRank() == 0)
          std::cout << "UserInputForTests, open " <<
          chCoordFileName.str () << std::endl;
      }

      assignFile = fopen(chAssignFileName.str().c_str(), "r");
      if (assignFile){
        fileInfo[2] = 1;
        if (verbose_ && tcomm_->getRank() == 0)
          std::cout << "UserInputForTests, open " <<
          chAssignFileName.str () << std::endl;
      }
    }else{
      if (verbose_ && tcomm_->getRank() == 0){
        std::cout << "UserInputForTests, unable to open file: ";
        std::cout << chGraphFileName.str() << std::endl;
      }
    }
  }

  // broadcast whether we have graphs and coords to all processes
  Teuchos::broadcast<int, int>(*tcomm_, 0, 3, fileInfo);

  bool haveGraph = (fileInfo[0] == 1);
  bool haveCoords = (fileInfo[1] == 1);
  bool haveAssign = (fileInfo[2] == 1);

  if (haveGraph){
    // builds M_, vtxWeights_, and edgWeights_ and closes file.
    try{
      getUIChacoGraph(graphFile, haveAssign, assignFile,
                      testData, distributeInput);
    }
    Z2_FORWARD_EXCEPTIONS

    if (haveCoords){
      // builds xyz_ and closes the file.
      try{
        getUIChacoCoords(coordFile, testData);
      }
      Z2_FORWARD_EXCEPTIONS
    }
  }

  xM_ = Zoltan2::XpetraTraits<tcrsMatrix_t>::convertToXpetra(M_);
}

void UserInputForTests::getUIChacoGraph(FILE *fptr, bool haveAssign,
                                        FILE *assignFile, string fname,
                                        bool distributeInput)
{
  int rank = tcomm_->getRank();
  int graphCounts[5];
  int nvtxs=0, nedges=0;
  int nVwgts=0, nEwgts=0;
  int *start = NULL, *adj = NULL;
  float *ewgts = NULL, *vwgts = NULL;
  size_t *nzPerRow = NULL;
  size_t maxRowLen = 0;
  zgno_t base = 0;
  ArrayRCP<const size_t> rowSizes;
  int fail = 0;
  bool haveEdges = true;

  if (rank == 0){

    memset(graphCounts, 0, 5*sizeof(int));

    // Reads in the file and closes it when done.
    fail = chaco_input_graph(fptr, fname.c_str(), &start, &adj,
                             &nvtxs, &nVwgts, &vwgts, &nEwgts, &ewgts);

    // There are Zoltan2 test graphs that have no edges.

    // nEwgts must be 1 or 0 - add error

    if (start == NULL)
      haveEdges = false;

    if (verbose_)
    {
      std::cout << "UserInputForTests, " << nvtxs << " vertices,";
      if (haveEdges)
        std::cout << start[nvtxs] << " edges,";
      else
        std::cout << "no edges,";
      std::cout << nVwgts << " vertex weights, ";
      std::cout << nEwgts << " edge weights" << std::endl;
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
    graphCounts[2] = nVwgts;
    graphCounts[3] = nEwgts;
    graphCounts[4] = (int)maxRowLen; // size_t maxRowLen will fit; it is <= (int-int)
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
      rcp(new tcrsMatrix_t(fromMap, rowSizes(), Tpetra::StaticProfile));

    if (haveEdges){

      zgno_t *edgeIds = new zgno_t [nedges];
      if (nedges && !edgeIds)
        throw std::bad_alloc();
      for (int i=0; i < nedges; i++)
        edgeIds[i] = adj[i];

      free(adj);
      adj = NULL;

      zgno_t *nextId = edgeIds;
      Array<zscalar_t> values(nedges, 1.0);
      if (nedges > 0 && nEwgts > 0) {
        for (int i=0; i < nedges; i++)
          values[i] = ewgts[i];
        free(ewgts);
        ewgts = NULL;
      }

      for (int i=0; i < nvtxs; i++){
        if (nzPerRow[i] > 0){
          ArrayView<const zgno_t> rowNz(nextId, nzPerRow[i]);
          fromMatrix->insertGlobalValues(i, rowNz, values.view(start[i], start[i+1] - start[i]));
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
    nVwgts = graphCounts[2];
    nEwgts  = graphCounts[3];
    maxRowLen = graphCounts[4];

    // Create a Tpetra::CrsMatrix where rank 0 has entire matrix.

    fromMap = rcp(new map_t(nvtxs, 0, base, tcomm_));

    fromMatrix =
      rcp(new tcrsMatrix_t(fromMap, rowSizes(), Tpetra::StaticProfile));

    fromMatrix->fillComplete();
  }

#ifdef KDDKDDPRINT
  if (rank == 0) {
    size_t sz = fromMatrix->getNodeMaxNumRowEntries();
    Teuchos::Array<zgno_t> indices(sz);
    Teuchos::Array<zscalar_t> values(sz);
    for (size_t i = 0; i < fromMatrix->getNodeNumRows(); i++) {
      zgno_t gid = fromMatrix->getRowMap()->getGlobalElement(i);
      size_t num;
      fromMatrix->getGlobalRowCopy(gid, indices(), values(), num);
      std::cout << "ROW " << gid << ": ";
      for (size_t j = 0; j < num; j++)
        std::cout << indices[j] << " ";
      std::cout << std::endl;
    }
  }
#endif

  RCP<const map_t> toMap;
  RCP<tcrsMatrix_t> toMatrix;
  RCP<import_t> importer;

  if (distributeInput) {
    if (haveAssign) {
      // Read assignments from Chaco assignment file
      short *assignments = new short[nvtxs];
      if (rank == 0) {
        fail = chaco_input_assign(assignFile, fname.c_str(), nvtxs, assignments);
      }
      // Broadcast coordinate dimension
      Teuchos::broadcast<int, short>(*tcomm_, 0, nvtxs, assignments);

      // Build map with my vertices
      Teuchos::Array<zgno_t> mine;
      for (int i = 0; i < nvtxs; i++) {
        if (assignments[i] == rank)
          mine.push_back(i);
      }

      Tpetra::global_size_t dummy =
              Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid();
      toMap = rcp(new map_t(dummy, mine(), base, tcomm_));
      delete [] assignments;
    }
    else {
      // Create a Tpetra::Map with default row distribution.
      toMap = rcp(new map_t(nvtxs, base, tcomm_));
    }
    toMatrix = rcp(new tcrsMatrix_t(toMap, maxRowLen));

    // Import the data.
    importer = rcp(new import_t(fromMap, toMap));
    toMatrix->doImport(*fromMatrix, *importer, Tpetra::INSERT);
    toMatrix->fillComplete();
  }
  else {
    toMap = fromMap;
    toMatrix = fromMatrix;
  }

  M_ = toMatrix;

  // Vertex weights, if any

  typedef ArrayRCP<const ArrayView<const zscalar_t> > arrayArray_t;

  if (nVwgts > 0){

    ArrayRCP<zscalar_t> weightBuf;
    ArrayView<const zscalar_t> *wgts = new ArrayView<const zscalar_t> [nVwgts];

    if (rank == 0){
      size_t len = nVwgts * nvtxs;
      zscalar_t *buf = new zscalar_t [len];
      if (!buf) throw std::bad_alloc();
      weightBuf = arcp(buf, 0, len, true);

      for (int widx=0; widx < nVwgts; widx++){
        wgts[widx] = ArrayView<const zscalar_t>(buf, nvtxs);
        float *vw = vwgts + widx;
        for (int i=0; i < nvtxs; i++, vw += nVwgts)
          buf[i] = *vw;
        buf += nvtxs;
      }

      free(vwgts);
      vwgts = NULL;
    }

    arrayArray_t vweights = arcp(wgts, 0, nVwgts, true);

    RCP<tMVector_t> fromVertexWeights =
    rcp(new tMVector_t(fromMap, vweights.view(0, nVwgts), nVwgts));

    RCP<tMVector_t> toVertexWeights;
    if (distributeInput) {
      toVertexWeights = rcp(new tMVector_t(toMap, nVwgts));
      toVertexWeights->doImport(*fromVertexWeights, *importer, Tpetra::INSERT);
    }
    else
      toVertexWeights = fromVertexWeights;

    vtxWeights_ = toVertexWeights;
  }

  // Edge weights, if any

  if (haveEdges && nEwgts > 0){

    // No longer distributing edge weights; they are the matrix values
    /*
    ArrayRCP<zscalar_t> weightBuf;
    ArrayView<const zscalar_t> *wgts = new ArrayView<const zscalar_t> [nEwgts];

    toMap = rcp(new map_t(nedges, M_->getNodeNumEntries(), base, tcomm_));

    if (rank == 0){
      size_t len = nEwgts * nedges;
      zscalar_t *buf = new zscalar_t [len];
      if (!buf) throw std::bad_alloc();
      weightBuf = arcp(buf, 0, len, true);

      for (int widx=0; widx < nEwgts; widx++){
        wgts[widx] = ArrayView<const zscalar_t>(buf, nedges);
        float *ew = ewgts + widx;
        for (int i=0; i < nedges; i++, ew += nEwgts)
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

    arrayArray_t eweights = arcp(wgts, 0, nEwgts, true);

    RCP<tMVector_t> fromEdgeWeights;
    RCP<tMVector_t> toEdgeWeights;
    RCP<import_t> edgeImporter;

    if (distributeInput) {
      fromEdgeWeights =
      rcp(new tMVector_t(fromMap, eweights.view(0, nEwgts), nEwgts));
      toEdgeWeights = rcp(new tMVector_t(toMap, nEwgts));
      edgeImporter = rcp(new import_t(fromMap, toMap));
      toEdgeWeights->doImport(*fromEdgeWeights, *edgeImporter, Tpetra::INSERT);
    }
    else
      toEdgeWeights = fromEdgeWeights;

    edgWeights_ = toEdgeWeights;
    */

    toMap = rcp(new map_t(nedges, M_->getNodeNumEntries(), base, tcomm_));
    edgWeights_ = rcp(new tMVector_t(toMap, nEwgts));

    size_t maxSize = M_->getNodeMaxNumRowEntries();
    Array<zlno_t> colind(maxSize);
    Array<zscalar_t> vals(maxSize);
    size_t nEntries;

    for (size_t i = 0, idx = 0; i < M_->getNodeNumRows(); i++) {
      M_->getLocalRowCopy(i, colind, vals, nEntries);
      for (size_t j = 0; j < nEntries; j++) {
        edgWeights_->replaceLocalValue(idx, 0, vals[j]); // Assuming nEwgts==1
        idx++;
      }
    }
  }

  if (start) {
    free(start);
    start = NULL;
  }
}

void UserInputForTests::getUIChacoCoords(FILE *fptr, string fname)
{
  int rank = tcomm_->getRank();
  int ndim=0;
  double *x=NULL, *y=NULL, *z=NULL;
  int fail = 0;

  size_t globalNumVtx = M_->getGlobalNumRows();

  if (rank == 0){

    // This function is from the Zoltan C-library.

    // Reads in the file and closes it when done.
    fail = chaco_input_geom(fptr, fname.c_str(), (int)globalNumVtx,
                            &ndim, &x, &y, &z);

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

  ArrayRCP<ArrayRCP<const zscalar_t> > coords(ndim);
  zlno_t len = 0;

  if (rank == 0){

    for (int dim=0; dim < ndim; dim++){
      zscalar_t *v = new zscalar_t [globalNumVtx];
      if (!v)
        throw std::bad_alloc();
      coords[dim] = arcp<const zscalar_t>(v, 0, globalNumVtx, true);
      double *val = (dim==0 ? x : (dim==1 ? y : z));
      for (size_t i=0; i < globalNumVtx; i++)
        v[i] = zscalar_t(val[i]);

      free(val);
    }

    len = static_cast<zlno_t>(globalNumVtx);
  }

  RCP<const map_t> fromMap = rcp(new map_t(globalNumVtx, len, 0, tcomm_));
  RCP<const map_t> toMap = M_->getRowMap();
  RCP<import_t> importer = rcp(new import_t(fromMap, toMap));

  Array<ArrayView<const zscalar_t> > coordData;
  for (int dim=0; dim < ndim; dim++)
    coordData.push_back(coords[dim].view(0, len));

  RCP<tMVector_t> fromCoords =
  rcp(new tMVector_t(fromMap, coordData.view(0, ndim), ndim));

  RCP<tMVector_t> toCoords = rcp(new tMVector_t(toMap, ndim));

  toCoords->doImport(*fromCoords, *importer, Tpetra::INSERT);

  xyz_ = toCoords;

}

///////////  Chaco reader -- copied from zoltan/ch //////////////////////////
///////////  This code might benefit from rewrite and simplification. ///////

double  UserInputForTests::chaco_read_val(
                                          FILE* infile,        /* file to read value from */
                                          int  *end_flag       /* 0 => OK, 1 => EOL, -1 => EOF */
)
{
  double  val;        /* return value */
  char   *ptr;        /* ptr to next string to read */
  char   *ptr2;       /* ptr to next string to read */
  int     length;     /* length of line to read */
  int     length_left;/* length of line still around */
  int     white_seen; /* have I detected white space yet? */
  int     done;       /* checking for end of scan */
  int     i;          /* loop counter */

  *end_flag = 0;

  if (chaco_offset == 0 || chaco_offset >= chaco_break_pnt) {
    if (chaco_offset >= chaco_break_pnt) { /* Copy rest of line back to beginning. */
      length_left = CHACO_LINE_LENGTH - chaco_save_pnt - 1;
      ptr2 = chaco_line;
      ptr = &chaco_line[chaco_save_pnt];
      for (i=length_left; i; i--) *ptr2++ = *ptr++;
      length = chaco_save_pnt + 1;
    }
    else {
      length = CHACO_LINE_LENGTH;
      length_left = 0;
    }

    /* Now read next line, or next segment of current one. */
    ptr2 = fgets(&chaco_line[length_left], length, infile);

    if (ptr2 == (char *) NULL) {    /* We've hit end of file. */
      *end_flag = -1;
      return((double) 0.0);
    }

    if ((chaco_line[CHACO_LINE_LENGTH - 2] != '\n') && (chaco_line[CHACO_LINE_LENGTH - 2] != '\f')
        && (strlen(chaco_line) == CHACO_LINE_LENGTH - 1)){
      /* Line too long.  Find last safe place in chaco_line. */
      chaco_break_pnt = CHACO_LINE_LENGTH - 1;
      chaco_save_pnt = chaco_break_pnt;
      white_seen = 0;
      done = 0;
      while (!done) {
        --chaco_break_pnt;
        if (chaco_line[chaco_break_pnt] != '\0') {
          if (isspace((int)(chaco_line[chaco_break_pnt]))) {
            if (!white_seen) {
              chaco_save_pnt = chaco_break_pnt + 1;
              white_seen = 1;
            }
          }
          else if (white_seen) {
            done= 1;
          }
        }
      }
    }
    else {
      chaco_break_pnt = CHACO_LINE_LENGTH;
    }

    chaco_offset = 0;
  }

  while (isspace((int)(chaco_line[chaco_offset])) && chaco_offset < CHACO_LINE_LENGTH) chaco_offset++;
  if (chaco_line[chaco_offset] == '%' || chaco_line[chaco_offset] == '#') {
    *end_flag = 1;
    if (chaco_break_pnt < CHACO_LINE_LENGTH) {
      chaco_flush_line(infile);
    }
    return((double) 0.0);
  }

  ptr = &(chaco_line[chaco_offset]);
  val = strtod(ptr, &ptr2);

  if (ptr2 == ptr) {    /* End of input line. */
    chaco_offset = 0;
    *end_flag = 1;
    return((double) 0.0);
  }
  else {
    chaco_offset = (int) (ptr2 - chaco_line) / sizeof(char);
  }

  return(val);
}


int UserInputForTests::chaco_read_int(
                                      FILE *infile,        /* file to read value from */
                                      int    *end_flag         /* 0 => OK, 1 => EOL, -1 => EOF */
)
{
  int     val;        /* return value */
  char   *ptr;        /* ptr to next string to read */
  char   *ptr2;        /* ptr to next string to read */
  int     length;        /* length of line to read */
  int     length_left;    /* length of line still around */
  int     white_seen;    /* have I detected white space yet? */
  int     done;        /* checking for end of scan */
  int     i;        /* loop counter */

  *end_flag = 0;

  if (chaco_offset == 0 || chaco_offset >= chaco_break_pnt) {
    if (chaco_offset >= chaco_break_pnt) { /* Copy rest of line back to beginning. */
      length_left = CHACO_LINE_LENGTH - chaco_save_pnt - 1;
      ptr2 = chaco_line;
      ptr = &chaco_line[chaco_save_pnt];
      for (i=length_left; i; i--) *ptr2++ = *ptr++;
      length = chaco_save_pnt + 1;
    }
    else {
      length = CHACO_LINE_LENGTH;
      length_left = 0;
    }

    /* Now read next line, or next segment of current one. */
    ptr2 = fgets(&chaco_line[length_left], length, infile);

    if (ptr2 == (char *) NULL) {    /* We've hit end of file. */
      *end_flag = -1;
      return(0);
    }

    if ((chaco_line[CHACO_LINE_LENGTH - 2] != '\n') && (chaco_line[CHACO_LINE_LENGTH - 2] != '\f')
        && (strlen(chaco_line) == CHACO_LINE_LENGTH - 1)){
      /* Line too long.  Find last safe place in line. */
      chaco_break_pnt = CHACO_LINE_LENGTH - 1;
      chaco_save_pnt = chaco_break_pnt;
      white_seen = 0;
      done = 0;
      while (!done) {
        --chaco_break_pnt;
        if (chaco_line[chaco_break_pnt] != '\0') {
          if (isspace((int)(chaco_line[chaco_break_pnt]))) {
            if (!white_seen) {
              chaco_save_pnt = chaco_break_pnt + 1;
              white_seen = 1;
            }
          }
          else if (white_seen) {
            done= 1;
          }
        }
      }
    }
    else {
      chaco_break_pnt = CHACO_LINE_LENGTH;
    }

    chaco_offset = 0;
  }

  while (isspace((int)(chaco_line[chaco_offset])) && chaco_offset < CHACO_LINE_LENGTH) chaco_offset++;
  if (chaco_line[chaco_offset] == '%' || chaco_line[chaco_offset] == '#') {
    *end_flag = 1;
    if (chaco_break_pnt < CHACO_LINE_LENGTH) {
      chaco_flush_line(infile);
    }
    return(0);
  }

  ptr = &(chaco_line[chaco_offset]);
  val = (int) strtol(ptr, &ptr2, 10);

  if (ptr2 == ptr) {    /* End of input chaco_line. */
    chaco_offset = 0;
    *end_flag = 1;
    return(0);
  }
  else {
    chaco_offset = (int) (ptr2 - chaco_line) / sizeof(char);
  }

  return(val);
}

void UserInputForTests::chaco_flush_line(
                                         FILE *infile         /* file to read value from */
)
{
  char    c;        /* character being read */

  c = fgetc(infile);
  while (c != '\n' && c != '\f')
    c = fgetc(infile);
}

int UserInputForTests::chaco_input_graph(
  FILE *fin,            /* input file */
  const char *inname,        /* name of input file */
  int **start,        /* start of edge list for each vertex */
  int **adjacency,        /* edge list data */
  int *nvtxs,        /* number of vertices in graph */
  int *nVwgts,        /* # of vertex weights per node */
  float **vweights,        /* vertex weight list data */
  int *nEwgts,        /* # of edge weights per edge */
  float **eweights         /* edge weight list data */
)
{
  int    *adjptr;        /* loops through adjacency data */
  float  *ewptr;        /* loops through edge weight data */
  int     narcs;        /* number of edges expected in graph */
  int     nedges;        /* twice number of edges really in graph */
  int     nedge;        /* loops through edges for each vertex */
  int     flag;        /* condition indicator */
  int     skip_flag;    /* should this edge be ignored? */
  int     end_flag;        /* indicates end of line or file */
  int     vtx;        /* vertex in graph */
  int     line_num;        /* line number in input file */
  int     sum_edges;    /* total number of edges read so far */
  int     option = 0;    /* input option */
  int     using_ewgts;    /* are edge weights in input file? */
  int     using_vwgts;    /* are vertex weights in input file? */
  int     vtxnums;        /* are vertex numbers in input file? */
  int     vertex;        /* current vertex being read */
  int     new_vertex;    /* new vertex being read */
  float   weight;        /* weight being read */
  float   eweight;        /* edge weight being read */
  int     neighbor;        /* neighbor of current vertex */
  int     error_flag;    /* error reading input? */
  int     j;        /* loop counters */

  /* Read first line  of input (= nvtxs, narcs, option). */
  /* The (decimal) digits of the option variable mean: 1's digit not zero => input
   edge weights 10's digit not zero => input vertex weights 100's digit not zero
   => include vertex numbers */

  *start = NULL;
  *adjacency = NULL;
  *vweights = NULL;
  *eweights = NULL;

  error_flag = 0;
  line_num = 0;

  /* Read any leading comment lines */
  end_flag = 1;
  while (end_flag == 1) {
    *nvtxs = chaco_read_int(fin, &end_flag);
    ++line_num;
  }
  if (*nvtxs <= 0) {
    printf("ERROR in graph file `%s':", inname);
    printf(" Invalid number of vertices (%d).\n", *nvtxs);
    fclose(fin);
    return(1);
  }

  narcs = chaco_read_int(fin, &end_flag);
  if (narcs < 0) {
    printf("ERROR in graph file `%s':", inname);
    printf(" Invalid number of expected edges (%d).\n", narcs);
    fclose(fin);
    return(1);
  }

  /*  Check if vertex or edge weights are used */
  if (!end_flag) {
    option = chaco_read_int(fin, &end_flag);
  }
  using_ewgts = option - 10 * (option / 10);
  option /= 10;
  using_vwgts = option - 10 * (option / 10);
  option /= 10;
  vtxnums = option - 10 * (option / 10);

  /* Get weight info from Chaco option */
  (*nVwgts) = using_vwgts;
  (*nEwgts) = using_ewgts;

  /* Read numbers of weights if they are specified separately */
  if (!end_flag && using_vwgts==1){
    j = chaco_read_int(fin, &end_flag);
    if (!end_flag) (*nVwgts) = j;
  }
  if (!end_flag && using_ewgts==1){
    j = chaco_read_int(fin, &end_flag);
    if (!end_flag) (*nEwgts) = j;
  }

  /* Discard rest of line */
  while (!end_flag)
    j = chaco_read_int(fin, &end_flag);

  /* Allocate space for rows and columns. */
  *start = (int *) malloc((unsigned) (*nvtxs + 1) * sizeof(int));
  if (narcs != 0)
    *adjacency = (int *) malloc((unsigned) (2 * narcs + 1) * sizeof(int));
  else
    *adjacency = NULL;

  if (using_vwgts)
    *vweights = (float *) malloc((unsigned) (*nvtxs) * (*nVwgts) * sizeof(float));
  else
    *vweights = NULL;

  if (using_ewgts)
    *eweights = (float *)
    malloc((unsigned) (2 * narcs + 1) * (*nEwgts) * sizeof(float));
  else
    *eweights = NULL;

  adjptr = *adjacency;
  ewptr = *eweights;

  sum_edges = 0;
  nedges = 0;
  (*start)[0] = 0;
  vertex = 0;
  vtx = 0;
  new_vertex = 1;
  while ((using_vwgts || vtxnums || narcs) && end_flag != -1) {
    ++line_num;

    /* If multiple input lines per vertex, read vertex number. */
    if (vtxnums) {
      j = chaco_read_int(fin, &end_flag);
      if (end_flag) {
        if (vertex == *nvtxs)
          break;
        printf("ERROR in graph file `%s':", inname);
        printf(" no vertex number in line %d.\n", line_num);
        fclose(fin);
        return (1);
      }
      if (j != vertex && j != vertex + 1) {
        printf("ERROR in graph file `%s':", inname);
        printf(" out-of-order vertex number in line %d.\n", line_num);
        fclose(fin);
        return (1);
      }
      if (j != vertex) {
        new_vertex = 1;
        vertex = j;
      }
      else
        new_vertex = 0;
    }
    else
      vertex = ++vtx;

    if (vertex > *nvtxs)
      break;

    /* If vertices are weighted, read vertex weight. */
    if (using_vwgts && new_vertex) {
      for (j=0; j<(*nVwgts); j++){
        weight = chaco_read_val(fin, &end_flag);
        if (end_flag) {
          printf("ERROR in graph file `%s':", inname);
          printf(" not enough weights for vertex %d.\n", vertex);
          fclose(fin);
          return (1);
        }
        (*vweights)[(vertex-1)*(*nVwgts)+j] = weight;
      }
    }

    nedge = 0;

    /* Read number of adjacent vertex. */
    neighbor = chaco_read_int(fin, &end_flag);

    while (!end_flag) {
      skip_flag = 0;

      if (using_ewgts) {    /* Read edge weight if it's being input. */
        for (j=0; j<(*nEwgts); j++){
          eweight = chaco_read_val(fin, &end_flag);

          if (end_flag) {
            printf("ERROR in graph file `%s':", inname);
            printf(" not enough weights for edge (%d,%d).\n", vertex, neighbor);
            fclose(fin);
            return (1);
          }

          else {
            *ewptr++ = eweight;
          }
        }
      }

      /* Add edge to data structure. */
      if (!skip_flag) {
        if (++nedges > 2*narcs) {
          printf("ERROR in graph file `%s':", inname);
          printf(" at least %d adjacencies entered, but nedges = %d\n",
                 nedges, narcs);
          fclose(fin);
          return (1);
        }
        *adjptr++ = neighbor;
        nedge++;
      }

      /* Read number of next adjacent vertex. */
      neighbor = chaco_read_int(fin, &end_flag);
    }

    sum_edges += nedge;
    (*start)[vertex] = sum_edges;
  }

  /* Make sure there's nothing else in file. */
  flag = 0;
  while (!flag && end_flag != -1) {
    chaco_read_int(fin, &end_flag);
    if (!end_flag)
      flag = 1;
  }

  (*start)[*nvtxs] = sum_edges;

  if (vertex != 0) {        /* Normal file was read. */
    if (narcs) {
    }
    else { /* no edges, but did have vertex weights or vertex numbers */
      free(*start);
      *start = NULL;
      if (*adjacency != NULL)
        free(*adjacency);
      *adjacency = NULL;
      if (*eweights != NULL)
        free(*eweights);
      *eweights = NULL;
    }
  }

  else {
    /* Graph was empty */
    free(*start);
    if (*adjacency != NULL)
      free(*adjacency);
    if (*vweights != NULL)
      free(*vweights);
    if (*eweights != NULL)
      free(*eweights);
    *start = NULL;
    *adjacency = NULL;
  }

  fclose(fin);

  return (error_flag);
}


int UserInputForTests::chaco_input_geom(
  FILE *fingeom,        /* geometry input file */
  const char *geomname,        /* name of geometry file */
  int     nvtxs,        /* number of coordinates to read */
  int    *igeom,        /* dimensionality of geometry */
  double   **x,             /* coordinates of vertices */
  double   **y,
  double   **z
)
{
  double  xc, yc, zc =0;    /* first x, y, z coordinate */
  int     nread;        /* number of lines of coordinates read */
  int     flag;        /* any bad data at end of file? */
  int     line_num;        /* counts input lines in file */
  int     end_flag;        /* return conditional */
  int     ndims;        /* number of values in an input line */
  int     i=0;        /* loop counter */

  *x = *y = *z = NULL;
  line_num = 0;
  end_flag = 1;
  while (end_flag == 1) {
    xc = chaco_read_val(fingeom, &end_flag);
    ++line_num;
  }

  if (end_flag == -1) {
    printf("No values found in geometry file `%s'\n", geomname);
    fclose(fingeom);
    return (1);
  }

  ndims = 1;
  yc = chaco_read_val(fingeom, &end_flag);
  if (end_flag == 0) {
    ndims = 2;
    zc = chaco_read_val(fingeom, &end_flag);
    if (end_flag == 0) {
      ndims = 3;
      chaco_read_val(fingeom, &end_flag);
      if (!end_flag) {
        printf("Too many values on input line of geometry file `%s'\n",
               geomname);

        printf(" Maximum dimensionality is 3\n");
        fclose(fingeom);
        return (1);
      }
    }
  }

  *igeom = ndims;

  *x = (double *) malloc((unsigned) nvtxs * sizeof(double));
  (*x)[0] = xc;
  if (ndims > 1) {
    *y = (double *) malloc((unsigned) nvtxs * sizeof(double));
    (*y)[0] = yc;
  }
  if (ndims > 2) {
    *z = (double *) malloc((unsigned) nvtxs * sizeof(double));
    (*z)[0] = zc;
  }

  for (nread = 1; nread < nvtxs; nread++) {
    ++line_num;
    if (ndims == 1) {
      i = fscanf(fingeom, "%lf", &((*x)[nread]));
    }
    else if (ndims == 2) {
      i = fscanf(fingeom, "%lf%lf", &((*x)[nread]), &((*y)[nread]));
    }
    else if (ndims == 3) {
      i = fscanf(fingeom, "%lf%lf%lf", &((*x)[nread]), &((*y)[nread]),
                 &((*z)[nread]));
    }

    if (i == EOF) {
      printf("Too few lines of values in geometry file; nvtxs=%d, but only %d read\n",
             nvtxs, nread);
      fclose(fingeom);
      return (1);
    }
    else if (i != ndims) {
      printf("Wrong number of values in line %d of geometry file `%s'\n",
             line_num, geomname);
      fclose(fingeom);
      return (1);
    }
  }

  /* Check for spurious extra stuff in file. */
  flag = 0;
  end_flag = 0;
  while (!flag && end_flag != -1) {
    chaco_read_val(fingeom, &end_flag);
    if (!end_flag)
      flag = 1;
  }

  fclose(fingeom);

  return (0);
}

// Chaco input assignments from filename.assign; copied from Zoltan

int UserInputForTests::chaco_input_assign(
  FILE *finassign,              /* input assignment file */
  const char *inassignname,             /* name of input assignment file */
  int nvtxs,            /* number of vertices to output */
  short *assignment)            /* values to be printed */
{
    int       flag;             /* logical conditional */
    int       end_flag;         /* return flag from read_int() */
    int       i, j;             /* loop counter */

    /* Get the assignment vector one line at a time, checking as you go. */
    /* First read past any comments at top. */
    end_flag = 1;
    while (end_flag == 1) {
        assignment[0] = chaco_read_int(finassign, &end_flag);
    }

    if (assignment[0] < 0) {
        printf("ERROR: Entry %d in assignment file `%s' less than zero (%d)\n",
               1, inassignname, assignment[0]);
        fclose(finassign);
        return (1);
    }

    if (end_flag == -1) {
        printf("ERROR: No values found in assignment file `%s'\n", inassignname);
        fclose(finassign);
        return (1);
    }

    flag = 0;
    int np = this->tcomm_->getSize();
    if (assignment[0] >= np) flag = assignment[0];
    srand(this->tcomm_->getRank());
    for (i = 1; i < nvtxs; i++) {
        j = fscanf(finassign, "%hd", &(assignment[i]));
        if (j != 1) {
            printf("ERROR: Too few values in assignment file `%s'.\n", inassignname);
            fclose(finassign);
            return (1);
        }
        if (assignment[i] < 0) {
            printf("ERROR: Entry %d in assignment file `%s' less than zero (%d)\n",
                   i+1, inassignname, assignment[i]);
            fclose(finassign);
            return (1);
        }
        if (assignment[i] >= np) {    // warn since perhaps an error -- initial part 
                                      // assignment is greater than number of processors
            if (assignment[i] > flag)
                flag = assignment[i];
            assignment[i] = rand() % np;  // randomly assign vtx to a proc in this case
        }
    }
    srand(rand());

    if (flag) {
        printf("WARNING: Possible error in assignment file `%s'\n",
                         inassignname);
        printf("         Max assignment set (%d) greater than "
                         "max processor rank (%d)\n", flag, np-1);
        printf("         Some vertices given random initial assignments\n");
    }

    /* Check for spurious extra stuff in file. */
    flag = 0;
    end_flag = 0;
    while (!flag && end_flag != -1) {
        chaco_read_int(finassign, &end_flag);
        if (!end_flag)
            flag = 1;
    }
    if (flag) {
        printf("WARNING: Possible error in assignment file `%s'\n", inassignname);
        printf("         Numerical data found after expected end of file\n");
    }

    fclose(finassign);
    return (0);
}

// Pamgen Reader
void UserInputForTests::readPamgenMeshFile(string path, string testData)
{
#ifdef HAVE_ZOLTAN2_PAMGEN
  int rank = this->tcomm_->getRank();
  if (verbose_ && tcomm_->getRank() == 0)
    std::cout << "UserInputForTestsBD::readPamgenFile, Read: " << testData << std::endl;

  size_t len;
  std::fstream file;
  int dimension;
  if (rank == 0){
    // set file name
    std::ostringstream meshFileName;
    meshFileName << path << "/" << testData << ".pmgen";
    // open file

    file.open(meshFileName.str(), std::ios::in);

    if(!file.is_open()) // may be a problem with path or filename
    {
      if(verbose_ && tcomm_->getRank() == 0)
      {
        std::cout << "Unable to open pamgen mesh: ";
        std::cout << meshFileName.str();
        std::cout <<"\nPlease check file path and name." << std::endl;
      }
      len = 0; // broadcaset 0 length ->will cause exit
    }else{
      // write to character array
      // get size of file
      file.seekg (0,file.end);
      len = file.tellg();
      file.seekg (0);

      // get dimension
      dimension = 2;
      std::string line;
      while(std::getline(file,line))
      {
        if( line.find("nz") != std::string::npos ||
           line.find("nphi") != std::string::npos)
        {
          dimension = 3;
          break;
        }
      }

      file.clear();
      file.seekg(0, std::ios::beg);
    }
  }

  // broadcast the file size
  this->tcomm_->broadcast(0,sizeof(int), (char *)&dimension);
  this->tcomm_->broadcast(0,sizeof(size_t),(char *)&len);
  this->tcomm_->barrier();

  if(len == 0){
    if(verbose_ && tcomm_->getRank() == 0)
      std::cout << "Pamgen Mesh file size == 0, exiting UserInputForTests early." << std::endl;
    return;
  }

  char * file_data = new char[len+1];
  file_data[len] = '\0'; // critical to null terminate buffer
  if(rank == 0){
    file.read(file_data,len); // if proc 0 then read file
  }

  // broadcast the file to the world
  this->tcomm_->broadcast(0,(int)len+1,file_data);
  this->tcomm_->barrier();

  // Create the PamgenMesh

  this->pamgen_mesh = rcp(new PamgenMesh);
  this->havePamgenMesh = true;
  pamgen_mesh->createMesh(file_data,dimension,this->tcomm_);

  // save mesh info
  pamgen_mesh->storeMesh();
  this->tcomm_->barrier();

  // set coordinates
  this->setPamgenCoordinateMV();

  // set adjacency graph
  this->setPamgenAdjacencyGraph();

  this->tcomm_->barrier();
  if(rank == 0) file.close();
  delete [] file_data;
#else
  throw std::runtime_error("Pamgen requested but Trilinos "
                           "not built with Pamgen");
#endif
}

#ifdef HAVE_ZOLTAN2_PAMGEN
void UserInputForTests::setPamgenCoordinateMV()
{
  int dimension = pamgen_mesh->num_dim;
  // get coordinate and point info;
//  zlno_t numLocalPoints = pamgen_mesh->num_nodes;
//  zgno_t numGlobalPoints = pamgen_mesh->num_nodes_global;
  zgno_t numelements = pamgen_mesh->num_elem;
  zgno_t numGlobalElements = pamgen_mesh->num_elems_global;
  // allocate and set an array of coordinate arrays
  zscalar_t **elem_coords = new zscalar_t * [dimension];
  for(int i = 0; i < dimension; ++i){
    elem_coords[i] = new zscalar_t[numelements];
    double *tmp = &pamgen_mesh->element_coord[i*numelements];
    for (int j = 0; j < numelements; j++) elem_coords[i][j] = tmp[j];
  }

  // make a Tpetra map
  RCP<Tpetra::Map<zlno_t, zgno_t, znode_t> > mp;
  //   mp = rcp(new map_t(numGlobalElements, numelements, 0, this->tcomm_)); // constructo 1

//  Array<zgno_t>::size_type numEltsPerProc = numelements;
  Array<zgno_t> elementList(numelements);
  for (Array<zgno_t>::size_type k = 0; k < numelements; ++k) {
    elementList[k] = pamgen_mesh->element_order_map[k];
  }

  mp = rcp (new map_t (numGlobalElements, elementList, 0, this->tcomm_)); // constructor 2


  // make an array of array views containing the coordinate data
  Teuchos::Array<Teuchos::ArrayView<const zscalar_t> > coordView(dimension);
  for (int i = 0; i < dimension; i++){
    if(numelements > 0){
      Teuchos::ArrayView<const zscalar_t> a(elem_coords[i], numelements);
      coordView[i] = a;
    }
    else {
      Teuchos::ArrayView<const zscalar_t> a;
      coordView[i] = a;
    }
  }

  // set the xyz_ multivector
  xyz_ = RCP<tMVector_t>(new
                         tMVector_t(mp, coordView.view(0, dimension),
                                    dimension));
}

void UserInputForTests::setPamgenAdjacencyGraph()
{
//  int rank = this->tcomm_->getRank();
//  if(rank == 0) std::cout << "Making a graph from our pamgen mesh...." << std::endl;

  // Define Types
//  typedef zlno_t lno_t;
//  typedef zgno_t gno_t;

  // get info for setting up map
  size_t local_nodes = (size_t)this->pamgen_mesh->num_nodes;
  size_t local_els = (size_t)this->pamgen_mesh->num_elem;

  size_t global_els = (size_t)this->pamgen_mesh->num_elems_global; // global rows
  size_t global_nodes = (size_t)this->pamgen_mesh->num_nodes_global; //global columns
  // make map with global elements assigned to this mesh
  // make range map
//  if(rank == 0) std::cout << "Building Rowmap: " << std::endl;
  RCP<const map_t> rowMap = rcp(new map_t(global_els,0,this->tcomm_));
  RCP<const map_t> rangeMap = rowMap;

  // make domain map
  RCP<const map_t> domainMap = rcp(new map_t(global_nodes,0,this->tcomm_));

  // Get max number of nodes per element
  int blks = this->pamgen_mesh->num_elem_blk;
  int max_nodes_per_el = 0;
  for(int i = 0; i < blks; i++)
    if (this->pamgen_mesh->nodes_per_element[i] > max_nodes_per_el)
      max_nodes_per_el = this->pamgen_mesh->nodes_per_element[i];

  // make the element-node adjacency matrix
  Teuchos::RCP<tcrsMatrix_t> C = rcp(new tcrsMatrix_t(rowMap,max_nodes_per_el));

  Array<zgno_t> g_el_ids(local_els);
  for (size_t k = 0; k < local_els; ++k) {
    g_el_ids[k] = pamgen_mesh->global_element_numbers[k]-1;
  }

  Array<zgno_t> g_node_ids(local_nodes);
  for (size_t k = 0; k < local_nodes; ++k) {
    g_node_ids[k] = pamgen_mesh->global_node_numbers[k]-1;
  }

  zlno_t el_no = 0;
  zscalar_t one = static_cast<zscalar_t>(1);

//  if(rank == 0) std::cout << "Writing C... " << std::endl;
  for(int i = 0; i < blks; i++)
  {
    int el_per_block = this->pamgen_mesh->elements[i];
    int nodes_per_el = this->pamgen_mesh->nodes_per_element[i];
    int * connect = this->pamgen_mesh->elmt_node_linkage[i];

    for(int j = 0; j < el_per_block; j++)
    {
      const zgno_t gid = static_cast<zgno_t>(g_el_ids[el_no]);
      for(int k = 0; k < nodes_per_el; k++)
      {
        int g_node_i = g_node_ids[connect[j*nodes_per_el+k]-1];
        C->insertGlobalValues(gid,
                              Teuchos::tuple<zgno_t>(g_node_i),
                              Teuchos::tuple<zscalar_t>(one));
      }
      el_no++;
    }
  }

  C->fillComplete(domainMap, rangeMap);


  // Compute product C*C'
//  if(rank == 0) std::cout << "Compute Multiplication C... " << std::endl;
  RCP<tcrsMatrix_t> A = rcp(new tcrsMatrix_t(rowMap,0));
  Tpetra::MatrixMatrix::Multiply(*C, false, *C, true, *A);

  // remove entries not adjacent
  // make graph
//  if(rank == 0) std::cout << "Writing M_... " << std::endl;
  this->M_ = rcp(new tcrsMatrix_t(rowMap, A->getGlobalMaxNumRowEntries()));

//  if(rank == 0) std::cout << "\nSetting graph of connectivity..." << std::endl;
  Teuchos::ArrayView<const zgno_t> rowMapElementList =
                                        rowMap->getNodeElementList();
  for (Teuchos_Ordinal ii = 0; ii < rowMapElementList.size(); ii++)
  {
    zgno_t gid = rowMapElementList[ii];
    size_t numEntriesInRow = A->getNumEntriesInGlobalRow (gid);
    Array<zscalar_t> rowvals (numEntriesInRow);
    Array<zgno_t> rowinds (numEntriesInRow);

    // modified
    Array<zscalar_t> mod_rowvals;
    Array<zgno_t> mod_rowinds;
    A->getGlobalRowCopy (gid, rowinds (), rowvals (), numEntriesInRow);
    for (size_t i = 0; i < numEntriesInRow; i++) {
//      if (rowvals[i] == 2*(this->pamgen_mesh->num_dim-1))
//      {
      if (rowvals[i] >= 1)
      {
        mod_rowvals.push_back(one);
        mod_rowinds.push_back(rowinds[i]);
      }
    }
     this->M_->insertGlobalValues(gid, mod_rowinds, mod_rowvals);
  }

  this->M_->fillComplete();
  this->xM_ = Zoltan2::XpetraTraits<tcrsMatrix_t>::convertToXpetra(M_);
  //  if(rank == 0) std::cout << "Completed M... " << std::endl;

}

#endif

#endif
