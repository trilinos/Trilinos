// @HEADER
// ***********************************************************************
//                Copyright message goes here.
// ***********************************************************************

/*! \file UerInputForTests.hpp
 *  \brief Generate input for testing purposes.
 *
 *  Given:
 *  \li The name of matrix market file, or
 *  \li x, y and z dimensions and a problem type
 *
 *  Retrieve any of the following:
 *   \li a Tpetra::CrsMatrix
 *   \li a Tpetra::CrsGraph
 *   \li a Tpetra::Vector
 *   \li a Tpetra::MultiVector
 *   \li a Xpetra::CrsMatrix
 *   \li a Xpetra::CrsGraph
 *   \li a Xpetra::Vector
 *   \li a Xpetra::MultiVector
 *   \li a Epetra_CrsMatrix (if built with double, int, int)
 *   \li a Epetra_CrsGraph  (if built with double, int, int)
 *   \li the coordinates  (matrix market input only)
 *
 *  \todo compute coordinates for the x,y,z problems.
 *  \todo for very large files, each process reads in part of the file
 */

#include <Zoltan2_XpetraTraits.hpp>

#include <Teuchos_RCP.hpp>
#include <Teuchos_ArrayView.hpp>
#include <Teuchos_Array.hpp>
#include <Teuchos_Comm.hpp>

#include <Epetra_Vector.h>
#include <Epetra_CrsMatrix.h>

#include <Tpetra_Vector.hpp>
#include <Tpetra_CrsMatrix.hpp>

#include <Xpetra_Vector.hpp>
#include <Xpetra_CrsMatrix.hpp>

#include <MatrixMarket_Tpetra.hpp>
#include <MueLu_MatrixFactory.hpp>
#include <MueLu_GalleryParameters.hpp>

#include <Xpetra_EpetraUtils.hpp>
#ifdef HAVE_ZOLTAN2_MPI
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif

#include <bitset>

using Teuchos::RCP;
using Teuchos::ArrayRCP;
using Teuchos::ArrayView;
using Teuchos::Array;
using Teuchos::Comm;
using Teuchos::rcp;
using Teuchos::rcp_const_cast;

enum testOutputType {
  OBJECT_DATA,
  OBJECT_COORDINATES,
  NUM_TEST_OUTPUT_TYPE
};

typedef std::bitset<NUM_TEST_OUTPUT_TYPE> outputFlag_t;
outputFlag_t defaultFlags(OBJECT_DATA);

static int z2Test_read_mtx_coords(
  std::string &fname, ArrayRCP<ArrayRCP<scalar_t> > &xyz);

/*! \brief A class that generates typical user input for testing.
 */

class UserInputForTests
{
private:
  typedef Tpetra::CrsMatrix<scalar_t, lno_t, gno_t, node_t> tcrsMatrix_t;
  typedef Tpetra::CrsGraph<lno_t, gno_t, node_t> tcrsGraph_t;
  typedef Tpetra::Vector<scalar_t, lno_t, gno_t, node_t> tVector_t;
  typedef Tpetra::MultiVector<scalar_t, lno_t, gno_t, node_t> tMVector_t;

  typedef Xpetra::CrsMatrix<scalar_t, lno_t, gno_t, node_t> xcrsMatrix_t;
  typedef Xpetra::CrsGraph<lno_t, gno_t, node_t> xcrsGraph_t;
  typedef Xpetra::Vector<scalar_t, lno_t, gno_t, node_t> xVector_t;
  typedef Xpetra::MultiVector<scalar_t, lno_t, gno_t, node_t> xMVector_t;

  gno_t xdim_, ydim_, zdim_;

  std::string mtype_;

  std::string fname_;
  RCP<Zoltan2::default_node_t> node_;

  const RCP<const Comm<int> > tcomm_;

#ifdef HAVE_EPETRA_DATA_TYPES
  RCP<const Epetra_Comm> ecomm_;
#endif

  RCP<tcrsMatrix_t> M_; 
  RCP<xcrsMatrix_t> xM_; 

#ifdef HAVE_EPETRA_DATA_TYPES
  RCP<Epetra_CrsMatrix> eM_; 
  RCP<Epetra_CrsGraph> eG_; 
#endif

  outputFlag_t flags_;
  RCP<tMVector_t> xyz_;

  void readMatrixMarketFile()
  {
    using Tpetra::Map;
    using Tpetra::MultiVector;
    using Tpetra::Export;
 
    if (flags_.test(OBJECT_DATA)){
      try{
        M_ = Tpetra::MatrixMarket::Reader<tcrsMatrix_t>::readSparseFile(fname_, 
          tcomm_, node_);
      }
      catch (std::exception &e) {
        TEST_FAIL_AND_THROW(*tcomm_, 1, e.what());
      }
      RCP<const xcrsMatrix_t> xm = 
        Zoltan2::XpetraTraits<tcrsMatrix_t>::convertToXpetra(M_);
      xM_ = rcp_const_cast<xcrsMatrix_t>(xm);
    }

    if (!flags_.test(OBJECT_COORDINATES))
      return;

    // Rank 0 reads coordinate file and exports it.

    typedef Map<lno_t, gno_t, node_t> map_t;
    typedef Export<lno_t, gno_t, node_t> export_t;

    gno_t globalNrows = M_->getGlobalNumRows();
    gno_t base = M_->getIndexBase();
    const RCP<const map_t> &toMap = M_->getRowMap();
    int coordDim = 0;
    ArrayRCP<ArrayRCP<scalar_t> > xyz;

    if (tcomm_->getRank() == 0){
      bool fail = z2Test_read_mtx_coords(fname_, xyz); // Get coordinates

      if (fail)
        coordDim = 0;
      else
        coordDim = xyz.size();
    }

    // Broadcast coordinate dimension

    Teuchos::broadcast<int, int>(*tcomm_, 0, 1, &coordDim);

    if (coordDim == 0)
      throw std::runtime_error("No coordinates or not enough memory.");

    // Export coordinates to their owners

    xyz_ = rcp(new tMVector_t(toMap, coordDim));

    ArrayRCP<ArrayView<const scalar_t> > coordLists(coordDim);

    if (tcomm_->getRank() == 0){

      for (int dim=0; dim < coordDim; dim++)
        coordLists[dim] = xyz[dim].view(0, globalNrows);

      gno_t *tmp = new gno_t [globalNrows];
      if (!tmp)
        throw std::bad_alloc();

      ArrayRCP<const gno_t> rowIds = Teuchos::arcp(tmp, 0, globalNrows);

      for (gno_t id=base; id < base+globalNrows; id++)
        *tmp++ = id;

      RCP<const map_t> fromMap = rcp(new map_t(globalNrows, 
        rowIds.view(0, globalNrows), base, tcomm_));

      tMVector_t allCoords(fromMap, coordLists.view(0, coordDim), coordDim);

      export_t exporter(fromMap, toMap);

      xyz_->doExport(allCoords, exporter, Tpetra::REPLACE);
    }
    else{

      RCP<const map_t> fromMap = rcp(new map_t(globalNrows, 
        ArrayView<gno_t>(), base, tcomm_));

      tMVector_t allCoords(fromMap, coordLists.view(0, coordDim), coordDim);

      export_t exporter(fromMap, toMap);

      xyz_->doExport(allCoords, exporter, Tpetra::REPLACE);
    }
  }

  void buildCrsMatrix()
  {
    Teuchos::CommandLineProcessor tclp;
    MueLu::Gallery::Parameters<gno_t> params(tclp,
       xdim_, ydim_, zdim_, mtype_);

    RCP<const Tpetra::Map<lno_t, gno_t> > map =
      rcp(new Tpetra::Map<lno_t, gno_t>(
        params.GetNumGlobalElements(), 0, tcomm_));

    try{
      M_ = MueLu::Gallery::CreateCrsMatrix<scalar_t, lno_t, gno_t, 
        Tpetra::Map<lno_t, gno_t>, Tpetra::CrsMatrix<scalar_t, lno_t, gno_t> >
          (params.GetMatrixType(), map, params.GetParameterList()); 
    }
    catch (std::exception &e) {    // Probably not enough memory
      TEST_FAIL_AND_THROW(*tcomm_, 1, e.what());
    }

    RCP<const xcrsMatrix_t> xm = 
      Zoltan2::XpetraTraits<tcrsMatrix_t>::convertToXpetra(M_);
    xM_ = rcp_const_cast<xcrsMatrix_t>(xm);
  }

  void createMatrix()
  {
    if (M_.is_null()){
      if (xdim_ > 0){
        buildCrsMatrix();
      }
      else if (fname_.size() > 0){
        readMatrixMarketFile();
      }
      else{
        throw std::logic_error("programming error");
      }
    }
  }

public:

#ifdef HAVE_EPETRA_DATA_TYPES
  // Constructor for a user object created from a Matrix
  // Market file.

  UserInputForTests(std::string s, const RCP<const Comm<int> > &c,
    outputFlag_t flags=defaultFlags): 
    xdim_(0), ydim_(0), zdim_(0), mtype_(), fname_(s),
     node_(Kokkos::DefaultNode::getDefaultNode()), 
     tcomm_(c), ecomm_(),
     M_(), xM_(), eM_(), eG_(), flags_(flags), xyz_()
  {
    ecomm_ = Xpetra::toEpetra(c);
  }

  // Constructor for a user object created in memory using
  // a MueLue::Gallery factory.

  UserInputForTests(int x, int y, int z, const RCP<const Comm<int> > &c,
    std::string matrixType=std::string("Laplace3D"),
    outputFlag_t flags=defaultFlags): 
     xdim_(x), ydim_(y), zdim_(z), mtype_(matrixType), fname_(),
     node_(Kokkos::DefaultNode::getDefaultNode()),
     tcomm_(c), ecomm_(),
     M_(), xM_(), eM_(), eG_(), flags_(flags), xyz_()
  {
    if (flags_.test(OBJECT_COORDINATES))
      std::cout << "Coordinates for meshes not supported yet" << std::endl;
    ecomm_ = Xpetra::toEpetra(c);
  }
#else
  // Constructor for a user object created from a Matrix
  // Market file.

  UserInputForTests(std::string s, const RCP<const Comm<int> > &c,
    outputFlag_t flags=defaultFlags):
    xdim_(0), ydim_(0), zdim_(0), mtype_(), fname_(s), 
     node_(Kokkos::DefaultNode::getDefaultNode()),
     tcomm_(c), M_(), xM_(), flags_(flags), xyz_()
  {
  }

  // Constructor for a user object created in memory using
  // a MueLue::Gallery factory.

  UserInputForTests(gno_t x, gno_t y, gno_t z, const RCP<const Comm<int> > &c,
    std::string matrixType=std::string("Laplace3D"),
    outputFlag_t flags=defaultFlags): 
     xdim_(x), ydim_(y), zdim_(z), mtype_(matrixType), fname_(), 
     node_(Kokkos::DefaultNode::getDefaultNode()),
     tcomm_(c), M_(), xM_(), flags_(flags), xyz_()
  {
    if (flags_.test(OBJECT_COORDINATES))
      std::cout << "Coordinates for meshes not supported yet" << std::endl;
  }
#endif
  
  RCP<tMVector_t> getCoordinates()
  { 
    if (flags_.test(OBJECT_COORDINATES)){
      if (M_.is_null())
        createMatrix();
    }
    return xyz_;
  }

  RCP<tcrsMatrix_t> getTpetraCrsMatrix() 
  { 
    if (!flags_.test(OBJECT_DATA)) return M_;
    if (M_.is_null())
     createMatrix();
    return M_;
  }

  RCP<tcrsGraph_t> getTpetraCrsGraph() 
  { 
    if (!flags_.test(OBJECT_DATA)) 
      throw std::runtime_error("object data wasn't requested.");
    if (M_.is_null())
     createMatrix();
    return rcp_const_cast<tcrsGraph_t>(M_->getCrsGraph());
  }

  RCP<tVector_t> getTpetraVector() 
  { 
    if (!flags_.test(OBJECT_DATA)) 
      throw std::runtime_error("object data wasn't requested.");
    if (M_.is_null())
     createMatrix();
    RCP<tVector_t> V = rcp(new tVector_t(M_->getRowMap(),  1));
    V->randomize();
    
    return V;
  }

  RCP<tMVector_t> getTpetraMultiVector(int nvec) 
  { 
    if (!flags_.test(OBJECT_DATA)) 
      throw std::runtime_error("object data wasn't requested.");
    if (M_.is_null())
     createMatrix();
    RCP<tMVector_t> mV = rcp(new tMVector_t(M_->getRowMap(), nvec));
    mV->randomize();
    
    return mV;
  }

  RCP<xcrsMatrix_t> getXpetraCrsMatrix() 
  { 
    if (!flags_.test(OBJECT_DATA)) 
      throw std::runtime_error("object data wasn't requested.");
    if (xM_.is_null())
     createMatrix();
    return xM_;
  }

  RCP<xcrsGraph_t> getXpetraCrsGraph() 
  { 
    if (!flags_.test(OBJECT_DATA)) 
      throw std::runtime_error("object data wasn't requested.");
    if (xM_.is_null())
     createMatrix();
    return rcp_const_cast<xcrsGraph_t>(xM_->getCrsGraph());
  }

  RCP<xVector_t> getXpetraVector() 
  { 
    if (!flags_.test(OBJECT_DATA)) 
      throw std::runtime_error("object data wasn't requested.");
    RCP<const tVector_t> tV = getTpetraVector();
    RCP<const xVector_t> xV =
      Zoltan2::XpetraTraits<tVector_t>::convertToXpetra(tV);
    return rcp_const_cast<xVector_t>(xV);
  }

  RCP<xMVector_t> getXpetraMultiVector(int nvec) 
  { 
    if (!flags_.test(OBJECT_DATA)) 
      throw std::runtime_error("object data wasn't requested.");
    RCP<const tMVector_t> tMV = getTpetraMultiVector(nvec);
    RCP<const xMVector_t> xMV =
      Zoltan2::XpetraTraits<tMVector_t>::convertToXpetra(tMV);
    return rcp_const_cast<xMVector_t>(xMV);
  }

#ifdef HAVE_EPETRA_DATA_TYPES
  RCP<Epetra_CrsGraph> getEpetraCrsGraph()
  {
    if (!flags_.test(OBJECT_DATA)) 
      throw std::runtime_error("object data wasn't requested.");
    if (eG_.is_null()){
      if (M_.is_null())
        createMatrix();

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
    }
    return eG_;
  }

  RCP<Epetra_CrsMatrix> getEpetraCrsMatrix()
  {
    if (!flags_.test(OBJECT_DATA)) 
      throw std::runtime_error("object data wasn't requested.");
    if (eM_.is_null()){
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
    }
    return eM_;
  }

  RCP<Epetra_Vector> getEpetraVector() 
  { 
    if (!flags_.test(OBJECT_DATA)) 
      throw std::runtime_error("object data wasn't requested.");
    RCP<Epetra_CrsGraph> egraph = getEpetraCrsGraph();
    RCP<Epetra_Vector> V = 
      rcp(new Epetra_Vector(egraph->RowMap()));
    V->Random();
    return V;
  }

  RCP<Epetra_MultiVector> getEpetraMultiVector(int nvec) 
  { 
    if (!flags_.test(OBJECT_DATA)) 
      throw std::runtime_error("object data wasn't requested.");
    RCP<Epetra_CrsGraph> egraph = getEpetraCrsGraph();
    RCP<Epetra_MultiVector> mV = 
      rcp(new Epetra_MultiVector(egraph->RowMap(), nvec));
    mV->Random();
    return mV;
  }
#endif
};

static int z2Test_read_mtx_coords(std::string &fname, 
  ArrayRCP<ArrayRCP<scalar_t> > &xyz)

{
  // Open the coordinate file.
  // If fname is "old_name.mtx", the new name is "old_name_coords.mtx".

  std::ifstream coordFile;
  std::ostringstream newName;
  std::string::size_type loc = fname.find('.');
  int fail=1;

  if (loc != std::string::npos){
    newName << fname.substr(0, loc) << "_coord" << fname.substr(loc);

    try{
      coordFile.open(newName.str().c_str());
      fail = 0;
    }
    catch (std::exception &e){ // there is no coordinate file
    }
  }

  if (fail)
    return 1;

  // Read past banner to number and dimension of coordinates.

  char c[256];
  gno_t numCoords;
  int coordDim;
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
      s >> numCoords >> coordDim;
      if (!s.eof() || numCoords < 1 || coordDim < 1)
        fail=1;
    }
  }

  if (fail || !done){
    coordFile.close();
    return 1;
  }

  // Read in the coordinates.

  xyz = Teuchos::arcp(new ArrayRCP<scalar_t> [coordDim], 0, coordDim);

  for (int dim=0; !fail && dim < coordDim; dim++){
    lno_t idx;
    scalar_t *tmp = new scalar_t [numCoords];
    if (!tmp)
      fail = 1;
    else{
      xyz[dim] = Teuchos::arcp(tmp, 0, numCoords);
  
      for (idx=0; !coordFile.eof() && idx < numCoords; idx++){
        coordFile.getline(c, 256);
        std::istringstream s(c);
        s >> tmp[idx];
      }

      if (idx < numCoords)
        fail = 1;
    }
  }

  coordFile.close();

  if (fail){
    ArrayRCP<scalar_t> emptyArray;
    for (int dim=0; dim < coordDim; dim++)
      xyz[dim] = emptyArray;   // free the memory
    return 1;
  }
  
  return 0;
}
