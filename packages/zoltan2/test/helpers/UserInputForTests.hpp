// @HEADER
// ***********************************************************************
//
//         Zoltan2: Sandia Partitioning Ordering & Coloring Library
//
//                Copyright message goes here.   TODO
//
// ***********************************************************************
//
// Create xpetra, tpetra, or epetra graph, matrix, vector or
//   multivector objects for testing.
//
//   Two choices:
//
//   1. Read the generating matrix from a MatrixMarket file.
//   2. Build the generating matrix in-core using MueLu::Gallery.
//

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

using Teuchos::RCP;
using Teuchos::ArrayRCP;
using Teuchos::ArrayView;
using Teuchos::Array;
using Teuchos::Comm;
using Teuchos::rcp;
using Teuchos::rcp_const_cast;

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

    void readMatrixMarketFile()
    {
      try{
        M_ = Tpetra::MatrixMarket::Reader<tcrsMatrix_t>::readSparseFile(
                 fname_, tcomm_, node_);
      }
      catch (std::exception &e) {
        TEST_FAIL_AND_THROW(*tcomm_, 1, e.what());
      }
      RCP<const xcrsMatrix_t> xm = 
        Zoltan2::XpetraTraits<tcrsMatrix_t>::convertToXpetra(M_);
      xM_ = rcp_const_cast<xcrsMatrix_t>(xm);
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
  
    UserInputForTests(std::string s, const RCP<const Comm<int> > &c): 
      xdim_(0), ydim_(0), zdim_(0), mtype_(), fname_(s),
       node_(Kokkos::DefaultNode::getDefaultNode()), 
       tcomm_(c), ecomm_(),
       M_(), xM_(), eM_(), eG_()
    {
      ecomm_ = Xpetra::toEpetra(c);
    }

    // Constructor for a user object created in memory using
    // a MueLue::Gallery factory.

    UserInputForTests(int x, int y, int z, const RCP<const Comm<int> > &c,
      std::string matrixType=std::string("Laplace3D")): 
       xdim_(x), ydim_(y), zdim_(z), mtype_(matrixType), fname_(),
       node_(Kokkos::DefaultNode::getDefaultNode()),
       tcomm_(c), ecomm_(),
       M_(), xM_(), eM_(), eG_()
    {
      ecomm_ = Xpetra::toEpetra(c);
    }
#else
    // Constructor for a user object created from a Matrix
    // Market file.

    UserInputForTests(std::string s, const RCP<const Comm<int> > &c):
      xdim_(0), ydim_(0), zdim_(0), mtype_(), fname_(s), 
       node_(Kokkos::DefaultNode::getDefaultNode()),
       tcomm_(c), M_(), xM_()
    {
    }

    // Constructor for a user object created in memory using
    // a MueLue::Gallery factory.

    UserInputForTests(gno_t x, gno_t y, gno_t z, const RCP<const Comm<int> > &c,
      std::string matrixType=std::string("Laplace3D")): 
       xdim_(x), ydim_(y), zdim_(z), mtype_(matrixType), fname_(), 
       node_(Kokkos::DefaultNode::getDefaultNode()),
       tcomm_(c), M_(), xM_()
    {
    }
#endif
    

    RCP<tcrsMatrix_t> getTpetraCrsMatrix() 
    { 
      if (M_.is_null())
       createMatrix();
      return M_;
    }

    RCP<tcrsGraph_t> getTpetraCrsGraph() 
    { 
      if (M_.is_null())
       createMatrix();
      return rcp_const_cast<tcrsGraph_t>(M_->getCrsGraph());
    }

    RCP<tVector_t> getTpetraVector() 
    { 
      if (M_.is_null())
       createMatrix();
      RCP<tVector_t> V = rcp(new tVector_t(M_->getRowMap(),  1));
      V->randomize();
      
      return V;
    }

    RCP<tMVector_t> getTpetraMultiVector(int nvec) 
    { 
      if (M_.is_null())
       createMatrix();
      RCP<tMVector_t> mV = rcp(new tMVector_t(M_->getRowMap(), nvec));
      mV->randomize();
      
      return mV;
    }

    RCP<xcrsMatrix_t> getXpetraCrsMatrix() 
    { 
      if (xM_.is_null())
       createMatrix();
      return xM_;
    }

    RCP<xcrsGraph_t> getXpetraCrsGraph() 
    { 
      if (xM_.is_null())
       createMatrix();
      return rcp_const_cast<xcrsGraph_t>(xM_->getCrsGraph());
    }

    RCP<xVector_t> getXpetraVector() 
    { 
      RCP<const tVector_t> tV = getTpetraVector();
      RCP<const xVector_t> xV =
        Zoltan2::XpetraTraits<tVector_t>::convertToXpetra(tV);
      return rcp_const_cast<xVector_t>(xV);
    }

    RCP<xMVector_t> getXpetraMultiVector(int nvec) 
    { 
      RCP<const tMVector_t> tMV = getTpetraMultiVector(nvec);
      RCP<const xMVector_t> xMV =
        Zoltan2::XpetraTraits<tMVector_t>::convertToXpetra(tMV);
      return rcp_const_cast<xMVector_t>(xMV);
    }

#ifdef HAVE_EPETRA_DATA_TYPES
    RCP<Epetra_CrsGraph> getEpetraCrsGraph()
    {
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
      RCP<Epetra_CrsGraph> egraph = getEpetraCrsGraph();
      RCP<Epetra_Vector> V = 
        rcp(new Epetra_Vector(egraph->RowMap()));
      V->Random();
      return V;
    }

    RCP<Epetra_MultiVector> getEpetraMultiVector(int nvec) 
    { 
      RCP<Epetra_CrsGraph> egraph = getEpetraCrsGraph();
      RCP<Epetra_MultiVector> mV = 
        rcp(new Epetra_MultiVector(egraph->RowMap(), nvec));
      mV->Random();
      return mV;
    }
#endif
};
