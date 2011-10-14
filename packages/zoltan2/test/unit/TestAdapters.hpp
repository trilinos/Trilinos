// @HEADER
// ***********************************************************************
//
//         Zoltan2: Sandia Partitioning Ordering & Coloring Library
//
//                Copyright message goes here.   TODO
//
// ***********************************************************************
//
// Create xpetra, tpetra, or epetra graph or matrix input adapters for
// testing purposes.  Two choices:
//
//   1. Read the matrix from a MatrixMarket file.
//   2. Build the matrix in-core using MueLu::Gallery.
//
// Note: These adapters use the default Kokkos::Node.

#include <iostream>
#include <vector>
#include <string>
#include <Zoltan2_config.h>

#include <Teuchos_RCP.hpp>
#include <Teuchos_ArrayRCP.hpp>
#include <Teuchos_ArrayView.hpp>
#include <Teuchos_Array.hpp>
#include <Teuchos_CommHelpers.hpp>
#include <MatrixMarket_Tpetra.hpp>
#include <Xpetra_EpetraUtils.hpp>

#include <Zoltan2_EpetraCrsGraphInput.hpp>
#include <Zoltan2_TpetraCrsGraphInput.hpp>
#include <Zoltan2_EpetraCrsMatrixInput.hpp>
#include <Zoltan2_TpetraCrsMatrixInput.hpp>
#include <Zoltan2_XpetraCrsMatrixInput.hpp>
#include <Zoltan2_XpetraCrsGraphInput.hpp>

#include <Epetra_SerialComm.h>
#include <Teuchos_DefaultComm.hpp>
#ifdef HAVE_MPI
#include <Zoltan2_Util.hpp>
#include <Epetra_MpiComm.h>
#endif

#include <MueLu_MatrixFactory.hpp>
#include <MueLu_GalleryParameters.hpp>


#define TEST_FAIL_AND_THROW(comm, ok, s){ \
int gval, lval=( (ok) ? 0 : 1);       \
Teuchos::reduceAll<int,int>(comm, Teuchos::REDUCE_SUM, 1, &lval, &gval);\
if (gval){ \
  throw std::runtime_error(std::string(s)); \
} \
}

#define TEST_FAIL_AND_EXIT(comm, ok, s, code){ \
int gval, lval=( (ok) ? 0 : 1);       \
Teuchos::reduceAll<int,int>(comm, Teuchos::REDUCE_SUM, 1, &lval, &gval);\
if (gval){ \
  if ((comm).getRank() == 0){\
    std::cerr << "Error: " << s << std::endl;\
    std::cout << "FAIL" << std::endl;\
  } \
  exit(code);\
} \
}

using Teuchos::RCP;
using Teuchos::ArrayRCP;
using Teuchos::ArrayView;
using Teuchos::Array;
using Teuchos::Comm;
using Teuchos::rcp_implicit_cast;
using Teuchos::rcp;

template <Z2CLASS_TEMPLATE>
class TestAdapters{

private:
    typedef Tpetra::CrsMatrix<Scalar, LNO, GNO> tcrsMatrix_t;
    typedef Xpetra::CrsMatrix<Scalar, LNO, GNO> xcrsMatrix_t;
    typedef Xpetra::TpetraCrsMatrix<Scalar, LNO, GNO> xtcrsMatrix_t;
    typedef Tpetra::CrsGraph<LNO, GNO> tcrsGraph_t;
    typedef Xpetra::CrsGraph<LNO, GNO> xcrsGraph_t;
    typedef Tpetra::Map<LNO, GNO> map_t;

    typedef Zoltan2::EpetraCrsGraphInput<Epetra_CrsGraph> EpetraCrsGraphInput;
    typedef Zoltan2::EpetraCrsMatrixInput<Epetra_CrsMatrix> EpetraCrsMatrixInput;
    typedef Zoltan2::TpetraCrsGraphInput<tcrsGraph_t> TpetraCrsGraphInput;
    typedef Zoltan2::TpetraCrsMatrixInput<tcrsMatrix_t> TpetraCrsMatrixInput;
    typedef Zoltan2::XpetraCrsGraphInput<xcrsGraph_t> XpetraCrsGraphInput;
    typedef Zoltan2::XpetraCrsMatrixInput<xcrsMatrix_t> XpetraCrsMatrixInput;

    GNO xdim_, ydim_, zdim_;

    std::string fname_;
    RCP<const Comm<int> > tcomm_; 
    RCP<Zoltan2::default_node_t> node_;

    RCP<tcrsMatrix_t> M_; 
    RCP<xcrsMatrix_t> xM_; 

    RCP<TpetraCrsGraphInput > tgi_;
    RCP<TpetraCrsMatrixInput > tmi_;
    RCP<XpetraCrsGraphInput > xgi_;
    RCP<XpetraCrsMatrixInput > xmi_;

    // not implemented, waiting for Tpetra::Map bug fix
    RCP<TpetraCrsMatrixInput > tmi_64_;

#ifdef HAVE_MALLINFO
    size_t mBytes_;
#endif

    void readMatrixMarketFile()
    {
#ifdef HAVE_MALLINFO
        mBytes_ = Zoltan2::getAllocatedMemory();
#endif
      try{
        M_ = Tpetra::MatrixMarket::Reader<tcrsMatrix_t>::readSparseFile(
                 fname_, tcomm_, node_);
      }
      catch (std::exception &e) {
        TEST_FAIL_AND_THROW(*tcomm_, 1, e.what());
      }
#ifdef HAVE_MALLINFO
        mBytes_ = Zoltan2::getAllocatedMemory() - mBytes_;
#endif
      RCP<xtcrsMatrix_t> xtcrs = rcp(new xtcrsMatrix_t(M_));
      xM_ = rcp_implicit_cast<xcrsMatrix_t>(xtcrs);
    }

    void buildCrsMatrix()
    {
#ifdef HAVE_MALLINFO
        mBytes_ = Zoltan2::getAllocatedMemory();
#endif
      Teuchos::CommandLineProcessor tclp;
      MueLu::Gallery::Parameters<GNO> params(tclp,
         xdim_, ydim_, zdim_, std::string("Laplace3D"));
 
      RCP<const Tpetra::Map<LNO, GNO> > map =
        rcp(new Tpetra::Map<LNO, GNO>(
          params.GetNumGlobalElements(), 0, tcomm_));

      try{
        M_ = MueLu::Gallery::CreateCrsMatrix<Scalar, LNO, GNO, 
          Tpetra::Map<LNO, GNO>, Tpetra::CrsMatrix<Scalar, LNO, GNO> >(
            params.GetMatrixType(), map, params.GetParameterList()); 
      }
      catch (std::exception &e) {    // Probably not enough memory
        TEST_FAIL_AND_THROW(*tcomm_, 1, e.what());
      }
#ifdef HAVE_MALLINFO
        mBytes_ = Zoltan2::getAllocatedMemory() - mBytes_;
#endif
      RCP<xtcrsMatrix_t> xtcrs = rcp(new xtcrsMatrix_t(M_));
      xM_ = rcp_implicit_cast<xcrsMatrix_t>(xtcrs);
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
    // Constructor for an InputAdapter created from a Matrix
    // Market file.
  
    TestAdapters(std::string s, const RCP<const Comm<int> > &c): 
      xdim_(0), ydim_(0), zdim_(0), fname_(s), tcomm_(c),
       node_(Kokkos::DefaultNode::getDefaultNode()), M_(), xM_(),
        tgi_(), tmi_(), xgi_(), xmi_(), tmi_64_() {}

    // Constructor for an InputAdapter created in memory using
    // a MueLue::Gallery factory.

    TestAdapters(GNO x, GNO y, GNO z, const RCP<const Comm<int> > &c): 
       xdim_(x), ydim_(y), zdim_(z), fname_(), tcomm_(c),
       node_(Kokkos::DefaultNode::getDefaultNode()), M_(), xM_(),
        tgi_(), tmi_(), xgi_(), xmi_(), tmi_64_() {}
    

    RCP<tcrsMatrix_t> getTpetraMatrix() 
    { 
      if (M_.is_null())
       createMatrix();
      return M_;
    }

    RCP<xcrsMatrix_t> getXpetraMatrix() 
    { 
      if (xM_.is_null())
       createMatrix();
      return xM_;
    }

#ifdef HAVE_MALLINFO
    // A count of memory bytes should be a size_t, but
    // mallinfo() predates size_t.
    int getMatrixSize()
    {
      if (M_.is_null())
       createMatrix();
      return mBytes_;
    }
#endif

    RCP<EpetraCrsGraphInput > getEpetraCrsGraphInputAdapter()
    {
      // Defined in specialization below.
      return Teuchos::null;
    }

    RCP<EpetraCrsMatrixInput > getEpetraCrsMatrixInputAdapter()
    {
      // Defined in specialization below.
      return Teuchos::null;
    }

  
    RCP<TpetraCrsGraphInput > getTpetraCrsGraphInputAdapter()
    {
      if (tgi_.is_null()){
        if (M_.is_null())
          createMatrix();
        tgi_ = rcp(new TpetraCrsGraphInput(M_->getCrsGraph()));
      }
      return tgi_;
    }
    
    RCP<TpetraCrsMatrixInput > getTpetraCrsMatrixInputAdapter()
    {
      if (tmi_.is_null()){
        if (M_.is_null())
          createMatrix();
        tmi_ = rcp(new TpetraCrsMatrixInput(M_));
      }
      return tmi_;
    }
    
    RCP<TpetraCrsMatrixInput > getTpetraCrsMatrixInputAdapter64(
      bool DestroyM=false)
    {
      if (sizeof(GNO) < 8){
        throw std::runtime_error("global IDs are less than 8 bytes");
      }

      throw std::runtime_error("not done yet");

      if (tmi_64_.is_null()){
        if (M_.is_null())
          createMatrix();

        // We're creating a new matrix which is the original
        // matrix with global IDs that use the high order bytes 
        // of an 8 byte ID.

        GNO base = M_->getIndexBase();
        GNO idOffset = 0x70f000000000; 
        global_size_t nrows = M_->getNodeNumRows();
        global_size_t maxnnz = M_->getNodeMaxNumRowEntries();
        global_size_t ngrows = M_->getGlobalNumRows();
        RCP<const map_t > rowMap = M_->getRowMap();
        RCP<const map_t > colMap = M_->getColMap();

        Array<GNO> newRowGNOs(nrows);
        ArrayView<const GNO> oldRowGNOs = rowMap->getNodeElementList();
        for (size_t i=0; i < nrows; i++){
          newRowGNOs[i] = oldRowGNOs[i]+ idOffset;
        }

        ArrayRCP<size_t> entriesPerRow(nrows);
        for (size_t i=0; i < nrows; i++){
          entriesPerRow[i] = M_->getNumEntriesInLocalRow(i);
        }

        RCP<map_t > newMap = 
          rcp(new map_t(ngrows, newRowGNOs, base, tcomm_));

        RCP<tcrsMatrix_t > newM64 =
          rcp(new tcrsMatrix_t(newMap, entriesPerRow, 
            Tpetra::StaticProfile));
 
        ArrayView<const LNO> lids;
        ArrayView<const Scalar> nzvals;
        Array<GNO> gids(maxnnz);

        for (size_t i=0; i < nrows; i++){
          M_->getLocalRowView(LNO(i), lids, nzvals);
          size_t numIndices = lids.size();
          for (size_t j=0; j < numIndices; j++){
            gids[j] = colMap->getGlobalElement(lids[j]) + idOffset;
          }
          newM64->insertGlobalValues( newMap->getGlobalElement(i),
            gids.view(0,numIndices), nzvals);
        }

        newM64->fillComplete();

        if (DestroyM)
          M_.release();
        
        tmi_64_ = rcp(new TpetraCrsMatrixInput(newM64));
      }
      return tmi_64_;
    }

    RCP<XpetraCrsMatrixInput > getXpetraCrsMatrixInputAdapter()
    {
      if (xmi_.is_null()){
        if (xM_.is_null())
          createMatrix();

        xmi_ = rcp(new XpetraCrsMatrixInput(xM_));
      }
      return xmi_;
    }
    
    RCP<XpetraCrsGraphInput > getXpetraCrsGraphInputAdapter()
    {
      if (xgi_.is_null()){
        if (xM_.is_null())
          createMatrix();

        xgi_ = rcp(new XpetraCrsGraphInput(xM_->getCrsGraph()));
      }
      return xgi_;
    }
};

//
// Specialization for Epetra_CrsGraph and Epetra_CrsMatrix
//

template <>
class TestAdapters<double,int,int,int,int,Kokkos::DefaultNode::DefaultNodeType>
{
private:
    typedef Tpetra::CrsMatrix<double, int, int> tcrsMatrix_t;
    typedef Xpetra::CrsMatrix<double, int, int> xcrsMatrix_t;
    typedef Xpetra::TpetraCrsMatrix<double, int, int> xtcrsMatrix_t;
    typedef Tpetra::CrsGraph<int, int> tcrsGraph_t;
    typedef Xpetra::CrsGraph<int, int> xcrsGraph_t;
    typedef Tpetra::Map<int, int> map_t;

    typedef Zoltan2::EpetraCrsGraphInput<Epetra_CrsGraph> EpetraCrsGraphInput;
    typedef Zoltan2::EpetraCrsMatrixInput<Epetra_CrsMatrix> EpetraCrsMatrixInput;
    typedef Zoltan2::TpetraCrsGraphInput<tcrsGraph_t> TpetraCrsGraphInput;
    typedef Zoltan2::TpetraCrsMatrixInput<tcrsMatrix_t> TpetraCrsMatrixInput;
    typedef Zoltan2::XpetraCrsGraphInput<xcrsGraph_t> XpetraCrsGraphInput;
    typedef Zoltan2::XpetraCrsMatrixInput<xcrsMatrix_t> XpetraCrsMatrixInput;

    int xdim_, ydim_, zdim_;

    std::string fname_;
    const RCP<const Comm<int> > tcomm_;
    RCP<const Epetra_Comm> ecomm_;

    RCP<Zoltan2::default_node_t> node_;

    RCP<tcrsMatrix_t> M_; 
    RCP<xcrsMatrix_t> xM_; 

    RCP<EpetraCrsGraphInput > egi_;
    RCP<EpetraCrsMatrixInput > emi_;
    RCP<TpetraCrsGraphInput > tgi_;
    RCP<TpetraCrsMatrixInput > tmi_;
    RCP<XpetraCrsGraphInput > xgi_;
    RCP<XpetraCrsMatrixInput > xmi_;

    void readMatrixMarketFile()
    {
      try{
        M_ = Tpetra::MatrixMarket::Reader<tcrsMatrix_t>::readSparseFile(
                 fname_, tcomm_, node_);
      }
      catch (std::exception &e) {
        TEST_FAIL_AND_THROW(*tcomm_, 1, e.what());
      }
      RCP<xtcrsMatrix_t> xtcrs = rcp(new xtcrsMatrix_t(M_));
      xM_ = rcp_implicit_cast<xcrsMatrix_t>(xtcrs);
    }

    void buildCrsMatrix()
    {
      Teuchos::CommandLineProcessor tclp;
      MueLu::Gallery::Parameters<int> params(tclp,
         xdim_, ydim_, zdim_, std::string("Laplace3D"));

      RCP<const Tpetra::Map<int, int> > map =
        rcp(new Tpetra::Map<int, int>(
          params.GetNumGlobalElements(), 0, tcomm_));

      try{
        // Note: MueLu::Gallery creats a matrix using the default node.
        M_ = MueLu::Gallery::CreateCrsMatrix<double, int, int,
          Tpetra::Map<int, int>, Tpetra::CrsMatrix<double, int, int> >(
            params.GetMatrixType(), map, params.GetParameterList());
      }
      catch (std::exception &e) {    // Probably not enough memory
        TEST_FAIL_AND_THROW(*tcomm_, 1, e.what());
      }
      RCP<xtcrsMatrix_t> xtcrs = rcp(new xtcrsMatrix_t(M_));
      xM_ = rcp_implicit_cast<xcrsMatrix_t>(xtcrs);
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
    ~TestAdapters() { }

    // Constructor for an InputAdapter created from a Matrix
    // Market file.

    TestAdapters(std::string s, const RCP<const Comm<int> > &c): 
       xdim_(0), ydim_(0), zdim_(0), fname_(s), tcomm_(c), ecomm_(), 
       node_(Kokkos::DefaultNode::getDefaultNode()), M_(), xM_(),
        egi_(), emi_(), tgi_(), tmi_(), xgi_(), xmi_() 
    { 
      ecomm_ = Xpetra::toEpetra(c); 
    }

    // Constructor for an InputAdapter created in memory using
    // a MueLue::Gallery factory.

    TestAdapters(int x, int y, int z, const RCP<const Comm<int> > &c): 
       xdim_(x), ydim_(y), zdim_(z), fname_(), tcomm_(c), ecomm_(), 
       node_(Kokkos::DefaultNode::getDefaultNode()), M_(), xM_(),
        egi_(), emi_(), tgi_(), tmi_(), xgi_(), xmi_() 
    { 
      ecomm_ = Xpetra::toEpetra(c); 
    }

    RCP<tcrsMatrix_t> getTpetraMatrix() 
    { 
      if (M_.is_null())
       createMatrix();
      return M_;
    }

    RCP<xcrsMatrix_t> getXpetraMatrix() 
    { 
      if (xM_.is_null())
       createMatrix();
      return xM_;
    }

    RCP<EpetraCrsGraphInput > getEpetraCrsGraphInputAdapter()
    {
      if (egi_.is_null()){
        if (M_.is_null())
          createMatrix();
        RCP<const tcrsGraph_t> tgraph = M_->getCrsGraph();
        RCP<const map_t > trowMap = tgraph->getRowMap();
        RCP<const map_t > tcolMap = tgraph->getColMap();

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
      
        RCP<Epetra_CrsGraph> egraph = rcp(
          new Epetra_CrsGraph(Copy, erowMap, rowSize.getRawPtr(), true));

        for (int i=0; i < nMyElts; i++){
          tgraph->getLocalRowView(i+base, colLid);
          for (int j=0; j < colLid.size(); j++)
            colGids[j] = tcolMap->getGlobalElement(colLid[j]);
          egraph->InsertGlobalIndices(gids[i], rowSize[i], colGids.getRawPtr());
        }
        egraph->FillComplete();

        egi_ = rcp(new EpetraCrsGraphInput(egraph));
      }
      return egi_;
    }

    RCP<EpetraCrsMatrixInput > getEpetraCrsMatrixInputAdapter()
    {
      if (emi_.is_null()){
        if (M_.is_null())
          createMatrix();
        RCP<EpetraCrsGraphInput> graphAdapter = 
          getEpetraCrsGraphInputAdapter();

        RCP<const Epetra_CrsGraph> egraph = graphAdapter->getGraph();

        RCP<Epetra_CrsMatrix> matrix = rcp(
          new Epetra_CrsMatrix(Copy, *egraph));

        size_t maxRow = M_->getNodeMaxNumRowEntries();
        int nrows = egraph->NumMyRows();
        int base = egraph->IndexBase();
        const Epetra_BlockMap &rowMap = egraph->RowMap();
        const Epetra_BlockMap &colMap = egraph->ColMap();
        Array<int> colGid(maxRow);

        for (int i=0; i < nrows; i++){
          ArrayView<const int> colLid;
          ArrayView<const double> nz;
          M_->getLocalRowView(i+base, colLid, nz);
          size_t rowSize = colLid.size();
          int rowGid = rowMap.GID(i+base);
          for (size_t j=0; j < rowSize; j++){
            colGid[j] = colMap.GID(colLid[j]);
          }
          matrix->InsertGlobalValues(
            rowGid, rowSize, nz.getRawPtr(), colGid.getRawPtr());
        }
        matrix->FillComplete();

        emi_ = rcp(new EpetraCrsMatrixInput(matrix));
      }
      return emi_;
    }
    
    RCP<TpetraCrsGraphInput > getTpetraCrsGraphInputAdapter()
    {
      if (tgi_.is_null()){
        if (M_.is_null())
          createMatrix();
        RCP<const tcrsGraph_t> graph = M_->getCrsGraph();
        tgi_ = rcp(new TpetraCrsGraphInput(graph));
      }
      return tgi_;
    }
    
    RCP<TpetraCrsMatrixInput > getTpetraCrsMatrixInputAdapter()
    {
      if (tmi_.is_null()){
        if (M_.is_null())
          createMatrix();
        tmi_ = rcp(new TpetraCrsMatrixInput(M_));
      }
      return tmi_;
    }

    RCP<XpetraCrsMatrixInput > getXpetraCrsMatrixInputAdapter()
    {
      if (xmi_.is_null()){
        if (xM_.is_null())
          createMatrix();

        xmi_ = rcp(new XpetraCrsMatrixInput(xM_));
      }
      return xmi_;
    }
    
    RCP<XpetraCrsGraphInput > getXpetraCrsGraphInputAdapter()
    {
      if (xgi_.is_null()){
        if (xM_.is_null())
          createMatrix();

        xgi_ = rcp(new XpetraCrsGraphInput(xM_->getCrsGraph()));
      }
      return xgi_;
    }
};
