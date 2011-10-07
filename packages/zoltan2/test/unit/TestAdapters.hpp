// @HEADER
// ***********************************************************************
//
//         Zoltan2: Sandia Partitioning Ordering & Coloring Library
//
//                Copyright message goes here.   TODO
//
// ***********************************************************************
//
// Create xpetra, tpetra, or epetra graph or matrix input adapters  
//   from matrix market files for testing.

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

#include <Epetra_CrsGraph.h>
#include <Zoltan2_EpetraCrsGraphInput.hpp>

#include <Tpetra_Map.hpp>
#include <Tpetra_CrsGraph.hpp>
#include <Zoltan2_TpetraCrsGraphInput.hpp>

#include <Epetra_CrsMatrix.h>
#include <Zoltan2_EpetraCrsMatrixInput.hpp>

#include <Tpetra_CrsMatrix.hpp>
#include <Zoltan2_TpetraCrsMatrixInput.hpp>

#include <Xpetra_CrsMatrix.hpp>
#include <Zoltan2_XpetraCrsMatrixInput.hpp>

#include <Xpetra_CrsGraph.hpp>
#include <Zoltan2_XpetraCrsGraphInput.hpp>

#include <Epetra_SerialComm.h>
#include <Teuchos_DefaultComm.hpp>
#ifdef HAVE_MPI
#include <Zoltan2_Util.hpp>
#include <Epetra_MpiComm.h>
#endif


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
    typedef Zoltan2::EpetraCrsGraphInput EpetraCrsGraphInput;
    typedef Zoltan2::EpetraCrsMatrixInput EpetraCrsMatrixInput;
    typedef Zoltan2::TpetraCrsGraphInput<LNO, GNO, LID, GID, Node> 
      TpetraCrsGraphInput;
    typedef Zoltan2::TpetraCrsMatrixInput<Z2PARAM_TEMPLATE> 
      TpetraCrsMatrixInput;
    typedef Zoltan2::XpetraCrsGraphInput<LNO, GNO, LID, GID, Node> 
      XpetraCrsGraphInput;
    typedef Zoltan2::XpetraCrsMatrixInput<Z2PARAM_TEMPLATE> 
      XpetraCrsMatrixInput;
    typedef Tpetra::CrsMatrix<Scalar, LNO, GNO, Node> crsMatrix_t;
    typedef Tpetra::Map<LNO, GNO, Node> map_t;

    std::string _fname;
    Teuchos::RCP<Teuchos::Comm<int> > _tcomm; 
    Teuchos::RCP<Node > _node;
    
    Teuchos::RCP<crsMatrix_t> _M; 

    Teuchos::RCP<TpetraCrsGraphInput > _tgi;
    Teuchos::RCP<TpetraCrsMatrixInput > _tmi;
    Teuchos::RCP<TpetraCrsMatrixInput > _tmi64;
    Teuchos::RCP<XpetraCrsGraphInput > _xgi;
    Teuchos::RCP<XpetraCrsMatrixInput > _xmi;

    void createMatrix()
    {
      try{
        _M = Tpetra::MatrixMarket::Reader<crsMatrix_t>::readSparseFile(
                 _fname, _tcomm,  _node);
      }
      catch (std::exception &e) {
        TEST_FAIL_AND_THROW(*_tcomm, 1, e.what());
      }
    }

public:
    // Constructor  
    // TODO not sure that Kokkos::DefaultNode::getDefaultNode() will
    //             compile for all Node types. 
    TestAdapters(std::string s): _fname(s), 
       _tcomm(Teuchos::rcp(new Teuchos::SerialComm<int>)), 
       _node(Kokkos::DefaultNode::getDefaultNode()), _M(), 
        _tgi(), _tmi(), _xgi(), _xmi() {}

#ifdef HAVE_MPI
    // Must have communicator before creating adapters.
    void setMpiCommunicator(MPI_Comm comm)
    {
      _tcomm = Zoltan2::getTeuchosMpiComm<int>(comm);
    }
#endif

    void setNode(Teuchos::RCP<Node> n) { _node = n; }

    Teuchos::RCP<crsMatrix_t> getMatrix() 
    { 
      if (_M.is_null())
       createMatrix();
      return _M;
    }

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
      if (_tgi.is_null()){
        if (_M.is_null())
          createMatrix();
        _tgi = Teuchos::rcp(new TpetraCrsGraphInput(_M->getCrsGraph()));
      }
      return _tgi;
    }
    
    Teuchos::RCP<TpetraCrsMatrixInput > getTpetraCrsMatrixInputAdapter()
    {
      if (_tmi.is_null()){
        if (_M.is_null())
          createMatrix();
        _tmi = Teuchos::rcp(new TpetraCrsMatrixInput(_M));
      }
      return _tmi;
    }
    
    Teuchos::RCP<TpetraCrsMatrixInput > getTpetraCrsMatrixInputAdapter64(
      bool DestroyM=false)
    {
      if (sizeof(GNO) < 8){
        throw std::runtime_error("global IDs are less than 8 bytes");
      }

      throw std::runtime_error("not done yet");

      if (_tmi64.is_null()){
        if (_M.is_null())
          createMatrix();

        // We're creating a new matrix which is the original
        // matrix with global IDs that use the high order bytes 
        // of an 8 byte ID.

        GNO base = _M->getIndexBase();
        GNO idOffset = 0x70f000000000; 
        global_size_t nrows = _M->getNodeNumRows();
        global_size_t maxnnz = _M->getNodeMaxNumRowEntries();
        global_size_t ngrows = _M->getGlobalNumRows();
        Teuchos::RCP<const map_t > rowMap = _M->getRowMap();
        Teuchos::RCP<const map_t > colMap = _M->getColMap();

        Teuchos::Array<GNO> newRowGNOs(nrows);
        Teuchos::ArrayView<const GNO> oldRowGNOs = rowMap->getNodeElementList();
        for (size_t i=0; i < nrows; i++){
          newRowGNOs[i] = oldRowGNOs[i]+ idOffset;
        }

        Teuchos::ArrayRCP<size_t> entriesPerRow(nrows);
        for (size_t i=0; i < nrows; i++){
          entriesPerRow[i] = _M->getNumEntriesInLocalRow(i);
        }

        Teuchos::RCP<map_t > newMap = 
          Teuchos::rcp(new map_t(ngrows, newRowGNOs, base, _tcomm));

        Teuchos::RCP<crsMatrix_t > newM64 =
          Teuchos::rcp(new crsMatrix_t(newMap, entriesPerRow, 
            Tpetra::StaticProfile));
 
        Teuchos::ArrayView<const LNO> lids;
        Teuchos::ArrayView<const Scalar> nzvals;
        Teuchos::Array<GNO> gids(maxnnz);

        for (size_t i=0; i < nrows; i++){
          _M->getLocalRowView(LNO(i), lids, nzvals);
          size_t numIndices = lids.size();
          for (size_t j=0; j < numIndices; j++){
            gids[j] = colMap->getGlobalElement(lids[j]) + idOffset;
          }
          newM64->insertGlobalValues( newMap->getGlobalElement(i),
            gids.view(0,numIndices), nzvals);
        }

        newM64->fillComplete();

        if (DestroyM)
          _M.release();
        
        _tmi64 = Teuchos::rcp(new TpetraCrsMatrixInput(newM64));
      }
      return _tmi64;
    }

    Teuchos::RCP<XpetraCrsMatrixInput > getXpetraCrsMatrixInputAdapter()
    {
      if (_xmi.is_null()){
        Teuchos::RCP<TpetraCrsMatrixInput> tmatrix = 
          getTpetraCrsMatrixInputAdapter();

        _xmi = Teuchos::rcp_implicit_cast<XpetraCrsMatrixInput>(tmatrix);
      }
      return _xmi;
    }
    
    Teuchos::RCP<XpetraCrsGraphInput > getXpetraCrsGraphInputAdapter()
    {
      if (_xgi.is_null()){
        Teuchos::RCP<TpetraCrsGraphInput> tgraph = 
          getTpetraCrsGraphInputAdapter();

        _xgi = Teuchos::rcp_implicit_cast<XpetraCrsGraphInput>(tgraph);
      }
      return _xgi;
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
    typedef Zoltan2::EpetraCrsGraphInput EpetraCrsGraphInput;   
    typedef Zoltan2::EpetraCrsMatrixInput EpetraCrsMatrixInput;  
    typedef Zoltan2::TpetraCrsGraphInput<int,int> TpetraCrsGraphInput;
    typedef Zoltan2::TpetraCrsMatrixInput<double,int,int> TpetraCrsMatrixInput;
    typedef Zoltan2::XpetraCrsGraphInput<int,int> XpetraCrsGraphInput;
    typedef Zoltan2::XpetraCrsMatrixInput<double,int,int> XpetraCrsMatrixInput;
    typedef Tpetra::CrsMatrix<double,int,int,nodeType> crsMatrix_t;
    typedef Tpetra::CrsGraph<int,int,nodeType> crsGraph_t;
    typedef Tpetra::Map<int,int,nodeType> map_t;

    std::string _fname;
    Teuchos::RCP<Teuchos::Comm<int> > _tcomm;
#ifdef HAVE_MPI
    Epetra_MpiComm *_ecomm;
#else
    Epetra_SerialComm *_ecomm;
#endif

    // TODO support the different types of nodes we want to test
    //Teuchos::RCP<> _node;

    Teuchos::RCP<Kokkos::DefaultNode::DefaultNodeType> _node;

    Teuchos::RCP<crsMatrix_t> _M; 

    Teuchos::RCP<EpetraCrsGraphInput > _egi;
    Teuchos::RCP<EpetraCrsMatrixInput > _emi;
    Teuchos::RCP<TpetraCrsGraphInput > _tgi;
    Teuchos::RCP<TpetraCrsMatrixInput > _tmi;
    Teuchos::RCP<XpetraCrsGraphInput > _xgi;
    Teuchos::RCP<XpetraCrsMatrixInput > _xmi;

    void createMatrix()
    {
      try{
        _M = Tpetra::MatrixMarket::Reader<crsMatrix_t>::readSparseFile(
                 _fname, _tcomm, _node);
      }
      catch (std::exception &e) {
        TEST_FAIL_AND_THROW(*_tcomm, 1, e.what());
      }
    }

public:
    ~TestAdapters() { if (_ecomm) delete _ecomm; }

    // Constructor  TODO not sure that Kokkos::getDefaultNode() will
    //             compile for all Node types. 
    TestAdapters(std::string s): _fname(s), _tcomm(), _ecomm(NULL),
       _node(Kokkos::DefaultNode::getDefaultNode()), _M(), 
        _egi(), _emi(),_tgi(), _tmi(), _xgi(), _xmi() 
    {
#ifdef HAVE_MPI
    _ecomm = new Epetra_MpiComm(MPI_COMM_WORLD);
    _tcomm = Zoltan2::getTeuchosMpiComm<int>(MPI_COMM_WORLD);
#else
    _ecomm = new Epetra_SerialComm;
    _tcomm = Teuchos::rcp(new Teuchos::SerialComm<int>):
#endif
    }

#ifdef HAVE_MPI
    // Must have communicator before creating adapters.
    void setMpiCommunicator(MPI_Comm comm)
    {
      _tcomm = Zoltan2::getTeuchosMpiComm<int>(comm);
      delete _ecomm;
      _ecomm = new Epetra_MpiComm(comm);
    }
#endif

    // TODO support the different types of nodes we want to test
    //void setNode(Teuchos::RCP<> n) { _node = n; }

    Teuchos::RCP<crsMatrix_t> getMatrix() 
    { 
      if (_M.is_null())
       createMatrix();
      return _M;
    }
  
    Teuchos::RCP<EpetraCrsGraphInput > getEpetraCrsGraphInputAdapter()
    {
      if (_egi.is_null()){
        if (_M.is_null())
          createMatrix();
        Teuchos::RCP<const crsGraph_t> tgraph = _M->getCrsGraph();
        Teuchos::RCP<const map_t > trowMap = tgraph->getRowMap();
        Teuchos::RCP<const map_t > tcolMap = tgraph->getColMap();

        int nElts = static_cast<int>(trowMap->getGlobalNumElements());
        int nMyElts = static_cast<int>(trowMap->getNodeNumElements());
        int base = trowMap->getIndexBase();
        Teuchos::ArrayView<const int> gids = trowMap->getNodeElementList();
        
        Epetra_BlockMap erowMap(nElts, nMyElts, 
          gids.getRawPtr(), 1, base, *_ecomm);

        Teuchos::Array<int> rowSize(nMyElts);
        for (int i=0; i < nMyElts; i++){
          rowSize[i] = static_cast<int>(_M->getNumEntriesInLocalRow(i+base));
        }

        size_t maxRow = _M->getNodeMaxNumRowEntries();
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

        _egi = Teuchos::rcp(new EpetraCrsGraphInput(egraph));
      }
      return _egi;
    }

    Teuchos::RCP<EpetraCrsMatrixInput > getEpetraCrsMatrixInputAdapter()
    {
      if (_emi.is_null()){
        Teuchos::RCP<EpetraCrsGraphInput> graphAdapter = 
          getEpetraCrsGraphInputAdapter();

        Teuchos::RCP<const Epetra_CrsGraph> egraph = graphAdapter->getGraph();

        Teuchos::RCP<Epetra_CrsMatrix> matrix = Teuchos::rcp(
          new Epetra_CrsMatrix(Copy, *egraph));

        size_t maxRow = _M->getNodeMaxNumRowEntries();
        int nrows = egraph->NumMyRows();
        int base = egraph->IndexBase();
        const Epetra_BlockMap &rowMap = egraph->RowMap();
        const Epetra_BlockMap &colMap = egraph->ColMap();
        Teuchos::Array<int> colGid(maxRow);

        for (int i=0; i < nrows; i++){
          Teuchos::ArrayView<const int> colLid;
          Teuchos::ArrayView<const double> nz;
          _M->getLocalRowView(i+base, colLid, nz);
          size_t rowSize = colLid.size();
          int rowGid = rowMap.GID(i+base);
          for (size_t j=0; j < rowSize; j++){
            colGid[j] = colMap.GID(colLid[j]);
          }
          matrix->InsertGlobalValues(
            rowGid, rowSize, nz.getRawPtr(), colGid.getRawPtr());
        }
        matrix->FillComplete();

        _emi = Teuchos::rcp(new EpetraCrsMatrixInput(matrix));
      }
      return _emi;
    }
    
    Teuchos::RCP<TpetraCrsGraphInput > getTpetraCrsGraphInputAdapter()
    {
      if (_tgi.is_null()){
        if (_M.is_null())
          createMatrix();
        Teuchos::RCP<const crsGraph_t> graph = _M->getCrsGraph();
        _tgi = Teuchos::rcp(new TpetraCrsGraphInput(graph));
      }
      return _tgi;
    }
    
    Teuchos::RCP<TpetraCrsMatrixInput > getTpetraCrsMatrixInputAdapter()
    {
      if (_tmi.is_null()){
        if (_M.is_null())
          createMatrix();
        _tmi = Teuchos::rcp(new TpetraCrsMatrixInput(_M));
      }
      return _tmi;
    }
    
    Teuchos::RCP<XpetraCrsMatrixInput > getXpetraCrsMatrixInputAdapter()
    {
      if (_xmi.is_null()){
        Teuchos::RCP<TpetraCrsMatrixInput> tmatrix = 
          getTpetraCrsMatrixInputAdapter();

        _xmi = Teuchos::rcp_implicit_cast<XpetraCrsMatrixInput>(tmatrix);
      }
      return _xmi;
    }
    
    Teuchos::RCP<XpetraCrsGraphInput > getXpetraCrsGraphInputAdapter()
    {
      if (_xgi.is_null()){
        Teuchos::RCP<TpetraCrsGraphInput> tgraph = 
          getTpetraCrsGraphInputAdapter();

        _xgi = Teuchos::rcp_implicit_cast<XpetraCrsGraphInput>(tgraph);
      }
      return _xgi;
    }
};
