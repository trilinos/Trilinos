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
  if ((comm).getRank() == 0)\
    std::cerr << "Error: " << s << std::endl;\
  exit(code);\
} \
}

template <Z2CLASS_TEMPLATE>
class TestAdapters{

private:
    typedef Tpetra::CrsMatrix<Scalar, LNO, GNO> tcrsMatrix_t;
    typedef Tpetra::CrsGraph<LNO, GNO> tcrsGraph_t;
    typedef Tpetra::Map<LNO, GNO> map_t;

    typedef Zoltan2::EpetraCrsGraphInput<Epetra_CrsGraph> EpetraCrsGraphInput;
    typedef Zoltan2::EpetraCrsMatrixInput<Epetra_CrsMatrix> EpetraCrsMatrixInput;
    typedef Zoltan2::TpetraCrsGraphInput<tcrsGraph_t> TpetraCrsGraphInput;
    typedef Zoltan2::TpetraCrsMatrixInput<tcrsMatrix_t> TpetraCrsMatrixInput;
    typedef Zoltan2::XpetraCrsGraphInput<tcrsGraph_t> XpetraCrsGraphInput;
    typedef Zoltan2::XpetraCrsMatrixInput<tcrsMatrix_t> XpetraCrsMatrixInput;

    GNO xdim_, ydim_, zdim_;

    std::string fname_;
    Teuchos::RCP<Teuchos::Comm<int> > tcomm_; 
    Teuchos::RCP<Zoltan2::default_node_t> node_;

    Teuchos::RCP<tcrsMatrix_t> M_; 

    Teuchos::RCP<TpetraCrsGraphInput > tgi_;
    Teuchos::RCP<TpetraCrsMatrixInput > tmi_;
    Teuchos::RCP<TpetraCrsMatrixInput > tmi_64_;
    Teuchos::RCP<XpetraCrsGraphInput > xgi_;
    Teuchos::RCP<XpetraCrsMatrixInput > xmi_;

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
    }

    void buildCrsMatrix()
    {
#ifdef HAVE_MALLINFO
        mBytes_ = Zoltan2::getAllocatedMemory();
#endif
      Teuchos::CommandLineProcessor tclp;
      MueLu::Gallery::Parameters<GNO> params(tclp,
         xdim_, ydim_, zdim_, std::string("Laplace3D"));
 
      Teuchos::RCP<const Tpetra::Map<LNO, GNO> > map =
        Teuchos::rcp(new Tpetra::Map<LNO, GNO>(
          params.GetNumGlobalElements(), 0, tcomm_));

      try{
        // Note: MueLu::Gallery creats a matrix using the default node.
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
  
    TestAdapters(std::string s): xdim_(0), ydim_(0), zdim_(0),
       fname_(s), tcomm_(Teuchos::rcp(new Teuchos::SerialComm<int>)), 
       node_(Kokkos::DefaultNode::getDefaultNode()), M_(), 
        tgi_(), tmi_(), tmi_64_(), xgi_(), xmi_() {}

    // Constructor for an InputAdapter created in memory using
    // a MueLue::Gallery factory.

    TestAdapters(GNO x, GNO y, GNO z): xdim_(x), ydim_(y), zdim_(z),
       fname_(), tcomm_(Teuchos::rcp(new Teuchos::SerialComm<int>)), 
       node_(Kokkos::DefaultNode::getDefaultNode()), M_(), 
        tgi_(), tmi_(), tmi_64_(), xgi_(), xmi_() {}
    

#ifdef HAVE_MPI
    // Must have communicator before creating adapters.
    void setMpiCommunicator(MPI_Comm comm)
    {
      tcomm_ = Zoltan2::getTeuchosMpiComm<int>(comm);
    }
#endif

    void setNode(Teuchos::RCP<Node> n) { node_ = n; }

    Teuchos::RCP<tcrsMatrix_t> getMatrix() 
    { 
      if (M_.is_null())
       createMatrix();
      return M_;
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

    Teuchos::RCP<EpetraCrsGraphInput > getEpetraCrsGraphInputAdapter()
    {
      // Defined in specialized version below.
      return Teuchos::null;
    }

    Teuchos::RCP<EpetraCrsMatrixInput > getEpetraCrsMatrixInputAdapter()
    {
      // Defined in specialized version below.
      return Teuchos::null;
    }

  
    Teuchos::RCP<TpetraCrsGraphInput > getTpetraCrsGraphInputAdapter()
    {
      if (tgi_.is_null()){
        if (M_.is_null())
          createMatrix();
        tgi_ = Teuchos::rcp(new TpetraCrsGraphInput(M_->getCrsGraph()));
      }
      return tgi_;
    }
    
    Teuchos::RCP<TpetraCrsMatrixInput > getTpetraCrsMatrixInputAdapter()
    {
      if (tmi_.is_null()){
        if (M_.is_null())
          createMatrix();
        tmi_ = Teuchos::rcp(new TpetraCrsMatrixInput(M_));
      }
      return tmi_;
    }
    
    Teuchos::RCP<TpetraCrsMatrixInput > getTpetraCrsMatrixInputAdapter64(
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
        Teuchos::RCP<const map_t > rowMap = M_->getRowMap();
        Teuchos::RCP<const map_t > colMap = M_->getColMap();

        Teuchos::Array<GNO> newRowGNOs(nrows);
        Teuchos::ArrayView<const GNO> oldRowGNOs = rowMap->getNodeElementList();
        for (size_t i=0; i < nrows; i++){
          newRowGNOs[i] = oldRowGNOs[i]+ idOffset;
        }

        Teuchos::ArrayRCP<size_t> entriesPerRow(nrows);
        for (size_t i=0; i < nrows; i++){
          entriesPerRow[i] = M_->getNumEntriesInLocalRow(i);
        }

        Teuchos::RCP<map_t > newMap = 
          Teuchos::rcp(new map_t(ngrows, newRowGNOs, base, tcomm_));

        Teuchos::RCP<tcrsMatrix_t > newM64 =
          Teuchos::rcp(new tcrsMatrix_t(newMap, entriesPerRow, 
            Tpetra::StaticProfile));
 
        Teuchos::ArrayView<const LNO> lids;
        Teuchos::ArrayView<const Scalar> nzvals;
        Teuchos::Array<GNO> gids(maxnnz);

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
        
        tmi_64_ = Teuchos::rcp(new TpetraCrsMatrixInput(newM64));
      }
      return tmi_64_;
    }

    Teuchos::RCP<XpetraCrsMatrixInput > getXpetraCrsMatrixInputAdapter()
    {
      if (xmi_.is_null()){
        Teuchos::RCP<TpetraCrsMatrixInput> tmatrix = 
          getTpetraCrsMatrixInputAdapter();

        xmi_ = Teuchos::rcp_implicit_cast<XpetraCrsMatrixInput>(tmatrix);
      }
      return xmi_;
    }
    
    Teuchos::RCP<XpetraCrsGraphInput > getXpetraCrsGraphInputAdapter()
    {
      if (xgi_.is_null()){
        Teuchos::RCP<TpetraCrsGraphInput> tgraph = 
          getTpetraCrsGraphInputAdapter();

        xgi_ = Teuchos::rcp_implicit_cast<XpetraCrsGraphInput>(tgraph);
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
    typedef Kokkos::DefaultNode::DefaultNodeType nodeType;
    typedef Epetra_CrsMatrix ecrsMatrix_t;
    typedef Tpetra::CrsMatrix<double,int,int,nodeType> tcrsMatrix_t;
    typedef Tpetra::CrsGraph<int,int,nodeType> tcrsGraph_t;

    typedef Zoltan2::EpetraCrsGraphInput<Epetra_CrsGraph> EpetraCrsGraphInput;   
    typedef Zoltan2::EpetraCrsMatrixInput<ecrsMatrix_t> EpetraCrsMatrixInput;  
    typedef Zoltan2::TpetraCrsGraphInput<tcrsGraph_t> TpetraCrsGraphInput;
    typedef Zoltan2::TpetraCrsMatrixInput<tcrsMatrix_t> TpetraCrsMatrixInput;
    typedef Zoltan2::XpetraCrsGraphInput<Tpetra::CrsGraph<int,int,nodeType> > XpetraCrsGraphInput;
    typedef Zoltan2::XpetraCrsMatrixInput<tcrsMatrix_t> XpetraCrsMatrixInput;
    typedef Tpetra::Map<int,int,nodeType> map_t;

    int xdim_, ydim_, zdim_;

    std::string fname_;
    Teuchos::RCP<Teuchos::Comm<int> > tcomm_;
#ifdef HAVE_MPI
    Epetra_MpiComm *ecomm_;
#else
    Epetra_SerialComm *ecomm_;
#endif

    Teuchos::RCP<Zoltan2::default_node_t> node_;

    Teuchos::RCP<tcrsMatrix_t> M_; 

    Teuchos::RCP<EpetraCrsGraphInput > egi_;
    Teuchos::RCP<EpetraCrsMatrixInput > emi_;
    Teuchos::RCP<TpetraCrsGraphInput > tgi_;
    Teuchos::RCP<TpetraCrsMatrixInput > tmi_;
    Teuchos::RCP<XpetraCrsGraphInput > xgi_;
    Teuchos::RCP<XpetraCrsMatrixInput > xmi_;

    void readMatrixMarketFile()
    {
      try{
        M_ = Tpetra::MatrixMarket::Reader<tcrsMatrix_t>::readSparseFile(
                 fname_, tcomm_, node_);
      }
      catch (std::exception &e) {
        TEST_FAIL_AND_THROW(*tcomm_, 1, e.what());
      }
    }

    void buildCrsMatrix()
    {
      Teuchos::CommandLineProcessor tclp;
      MueLu::Gallery::Parameters<int> params(tclp,
         xdim_, ydim_, zdim_, std::string("Laplace3D"));

      Teuchos::RCP<const Tpetra::Map<int, int> > map =
        Teuchos::rcp(new Tpetra::Map<int, int>(
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
    ~TestAdapters() { if (ecomm_) delete ecomm_; }

    // Constructor for an InputAdapter created from a Matrix
    // Market file.

    TestAdapters(std::string s): xdim_(0), ydim_(0), zdim_(0),
       fname_(s), tcomm_(), ecomm_(NULL), 
       node_(Kokkos::DefaultNode::getDefaultNode()), M_(),
        egi_(), emi_(), tgi_(), tmi_(), xgi_(), xmi_() 
    {
#ifdef HAVE_MPI
      ecomm_ = new Epetra_MpiComm(MPI_COMM_WORLD);
      tcomm_ = Zoltan2::getTeuchosMpiComm<int>(MPI_COMM_WORLD);
#else
      ecomm_ = new Epetra_SerialComm;
      tcomm_ = Teuchos::rcp(new Teuchos::SerialComm<int>):
#endif
    }

    // Constructor for an InputAdapter created in memory using
    // a MueLue::Gallery factory.

    TestAdapters(int x, int y, int z): xdim_(x), ydim_(y), zdim_(z),
       fname_(), tcomm_(), ecomm_(NULL), 
       node_(Kokkos::DefaultNode::getDefaultNode()), M_(),
        egi_(), emi_(), tgi_(), tmi_(), xgi_(), xmi_() 
    {
#ifdef HAVE_MPI
      ecomm_ = new Epetra_MpiComm(MPI_COMM_WORLD);
      tcomm_ = Zoltan2::getTeuchosMpiComm<int>(MPI_COMM_WORLD);
#else
      ecomm_ = new Epetra_SerialComm;
      tcomm_ = Teuchos::rcp(new Teuchos::SerialComm<int>):
#endif
    }

#ifdef HAVE_MPI
    // Must have communicator before creating adapters.
    void setMpiCommunicator(MPI_Comm comm)
    {
      tcomm_ = Zoltan2::getTeuchosMpiComm<int>(comm);
      delete ecomm_;
      ecomm_ = new Epetra_MpiComm(comm);
    }
#endif

    Teuchos::RCP<tcrsMatrix_t> getMatrix() 
    { 
      if (M_.is_null())
       createMatrix();
      return M_;
    }
  
    Teuchos::RCP<EpetraCrsGraphInput > getEpetraCrsGraphInputAdapter()
    {
      if (egi_.is_null()){
        if (M_.is_null())
          createMatrix();
        Teuchos::RCP<const tcrsGraph_t> tgraph = M_->getCrsGraph();
        Teuchos::RCP<const map_t > trowMap = tgraph->getRowMap();
        Teuchos::RCP<const map_t > tcolMap = tgraph->getColMap();

        int nElts = static_cast<int>(trowMap->getGlobalNumElements());
        int nMyElts = static_cast<int>(trowMap->getNodeNumElements());
        int base = trowMap->getIndexBase();
        Teuchos::ArrayView<const int> gids = trowMap->getNodeElementList();
        
        Epetra_BlockMap erowMap(nElts, nMyElts, 
          gids.getRawPtr(), 1, base, *ecomm_);

        Teuchos::Array<int> rowSize(nMyElts);
        for (int i=0; i < nMyElts; i++){
          rowSize[i] = static_cast<int>(M_->getNumEntriesInLocalRow(i+base));
        }

        size_t maxRow = M_->getNodeMaxNumRowEntries();
        Teuchos::Array<int> colGids(maxRow);
        Teuchos::ArrayView<const int> colLid;
      
        Teuchos::RCP<Epetra_CrsGraph> egraph = Teuchos::rcp(
          new Epetra_CrsGraph(Copy, erowMap, rowSize.getRawPtr(), true));

        for (int i=0; i < nMyElts; i++){
          tgraph->getLocalRowView(i+base, colLid);
          for (int j=0; j < colLid.size(); j++)
            colGids[j] = tcolMap->getGlobalElement(colLid[j]);
          egraph->InsertGlobalIndices(gids[i], rowSize[i], colGids.getRawPtr());
        }
        egraph->FillComplete();

        egi_ = Teuchos::rcp(new EpetraCrsGraphInput(egraph));
      }
      return egi_;
    }

    Teuchos::RCP<EpetraCrsMatrixInput > getEpetraCrsMatrixInputAdapter()
    {
      if (emi_.is_null()){
        Teuchos::RCP<EpetraCrsGraphInput> graphAdapter = 
          getEpetraCrsGraphInputAdapter();

        Teuchos::RCP<const Epetra_CrsGraph> egraph = graphAdapter->getGraph();

        Teuchos::RCP<Epetra_CrsMatrix> matrix = Teuchos::rcp(
          new Epetra_CrsMatrix(Copy, *egraph));

        size_t maxRow = M_->getNodeMaxNumRowEntries();
        int nrows = egraph->NumMyRows();
        int base = egraph->IndexBase();
        const Epetra_BlockMap &rowMap = egraph->RowMap();
        const Epetra_BlockMap &colMap = egraph->ColMap();
        Teuchos::Array<int> colGid(maxRow);

        for (int i=0; i < nrows; i++){
          Teuchos::ArrayView<const int> colLid;
          Teuchos::ArrayView<const double> nz;
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

        emi_ = Teuchos::rcp(new EpetraCrsMatrixInput(matrix));
      }
      return emi_;
    }
    
    Teuchos::RCP<TpetraCrsGraphInput > getTpetraCrsGraphInputAdapter()
    {
      if (tgi_.is_null()){
        if (M_.is_null())
          createMatrix();
        Teuchos::RCP<const tcrsGraph_t> graph = M_->getCrsGraph();
        tgi_ = Teuchos::rcp(new TpetraCrsGraphInput(graph));
      }
      return tgi_;
    }
    
    Teuchos::RCP<TpetraCrsMatrixInput > getTpetraCrsMatrixInputAdapter()
    {
      if (tmi_.is_null()){
        if (M_.is_null())
          createMatrix();
        tmi_ = Teuchos::rcp(new TpetraCrsMatrixInput(M_));
      }
      return tmi_;
    }
    
    Teuchos::RCP<XpetraCrsMatrixInput > getXpetraCrsMatrixInputAdapter()
    {
      if (xmi_.is_null()){
        if (M_.is_null())
          createMatrix();
        Teuchos::RCP<Xpetra::TpetraCrsMatrix<double,int,int> > xM = 
          Teuchos::rcp(new Xpetra::TpetraCrsMatrix<double,int,int>(
                                   Teuchos::rcp_const_cast<tcrsMatrix_t>(M_)));
        xmi_ = Teuchos::rcp(new XpetraCrsMatrixInput(xM));
      }
      return xmi_;
    }
    
    Teuchos::RCP<XpetraCrsGraphInput > getXpetraCrsGraphInputAdapter()
    {
      if (xgi_.is_null()){
        if (M_.is_null())
          createMatrix();
        Teuchos::RCP<Xpetra::TpetraCrsMatrix<double,int,int> > xM = 
          Teuchos::rcp(new Xpetra::TpetraCrsMatrix<double,int,int>(
                                   Teuchos::rcp_const_cast<tcrsMatrix_t>(M_)));
        Teuchos::RCP<const Xpetra::CrsGraph<int,int> > graph = 
          xM->getCrsGraph();

        xgi_ = Teuchos::rcp(new XpetraCrsGraphInput(graph));
      }
      return xgi_;
    }
};
