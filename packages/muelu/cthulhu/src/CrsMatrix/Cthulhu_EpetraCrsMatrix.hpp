#ifndef CTHULHU_EPETRACRSMATRIX_HPP
#define CTHULHU_EPETRACRSMATRIX_HPP

#include "Cthulhu_EpetraConfigDefs.hpp"

#include <Teuchos_ArrayViewDecl.hpp>

#include <Epetra_CrsMatrix.h>
#include <Epetra_Map.h>

#include "Cthulhu_CrsMatrix.hpp"
#include "Cthulhu_EpetraCrsGraph.hpp"

#include "Cthulhu_EpetraMap.hpp"
#include "Cthulhu_EpetraMultiVector.hpp"
#include "Cthulhu_EpetraVector.hpp"
#include "Cthulhu_Trans.hpp"

namespace Cthulhu {

  class EpetraCrsMatrix
    : public CrsMatrix<double,int,int>
  {
  public:

    //! @name Constructor/Destructor Methods
    //@{ 

    //! Constructor specifying the number of non-zeros for all rows.
    EpetraCrsMatrix(const RCP<const Map<int,int> > &rowMap, size_t maxNumEntriesPerRow, Cthulhu::ProfileType pftype = Cthulhu::DynamicProfile) 
    { 
      CTHULHU_RCP_DYNAMIC_CAST(const EpetraMap, rowMap, eRowMap, "Cthulhu::EpetraCrsMatrix constructors only accept Cthulhu::EpetraMap as input arguments.");

      //TODO: test constblk and blksize=1 for the Map ?
      mtx_ = rcp(new Epetra_CrsMatrix(Copy, eRowMap->getEpetra_Map(), maxNumEntriesPerRow, false)); // TODO Copy or View by default ? // TODO: bool StaticProfile
    }

    EpetraCrsMatrix(const Teuchos::RCP<Epetra_CrsMatrix> &mtx) : mtx_(mtx) {  }

    // !Destructor.
     virtual ~EpetraCrsMatrix() { }

    //@}

    //! @name Insertion/Removal Methods
    //@{ 

    //! Insert matrix entries, using global IDs.
    void insertGlobalValues(int globalRow, const ArrayView<const int> &cols, const ArrayView<const double> &vals);

    //! Scale the current values of a matrix, this = alpha*this. 
    void scale(const double &alpha) {  mtx_->Scale(alpha); }

    //! Signal that data entry is complete, specifying domain and range maps.
    void fillComplete(const RCP<const Map<int,int> > &domainMap, const RCP<const Map<int,int> > &rangeMap, OptimizeOption os = DoOptimizeStorage) { 
       
      CTHULHU_RCP_DYNAMIC_CAST(const EpetraMap, domainMap, tDomainMap, "Cthulhu::EpetraCrsMatrix::fillComplete() only accept Cthulhu::TpetraMap as input arguments.");
      CTHULHU_RCP_DYNAMIC_CAST(const EpetraMap, rangeMap,  tRangeMap,  "Cthulhu::EpetraCrsMatrix::fillComplete() only accept Cthulhu::TpetraMap as input arguments.");
      mtx_->FillComplete(tDomainMap->getEpetra_Map(), tRangeMap->getEpetra_Map()); // TODO: os 
    }

    //! Signal that data entry is complete. 
    void fillComplete(Cthulhu::OptimizeOption os = Cthulhu::DoOptimizeStorage) {  if (os == Cthulhu::DoOptimizeStorage) mtx_->FillComplete(true); else mtx_->FillComplete(false); }

    //! leftScale
    void leftScale(const Vector<double,int,int>& x){
      CTHULHU_DYNAMIC_CAST(const EpetraVector, x, tX, "Cthulhu::EpetraCrsMatrix->leftScale() only accept Cthulhu::EpetraVector as input arguments.");
      mtx_->LeftScale(*tX.getEpetra_Vector());
    }

    //@}

    //! @name Methods implementing RowMatrix
    //@{ 

    //! Returns the communicator.
    const RCP<const Comm<int> > getComm() const {
      RCP<const Epetra_Comm> rcpComm = rcpFromRef(mtx_->Comm());
      return Epetra2Teuchos_Comm(rcpComm);
    }

    //! Returns the Map that describes the row distribution in this matrix.
    const RCP<const Cthulhu::Map<int,int> > getRowMap() const { 
      // The Epetra map is owned by the underlying Epetra matrix and freed when the matrix is deleted.
      // I have to make a copy of the map to be sure that the RCP<Map> returned by this method will remain valid even if the matrix is deleted.
      // Note that it is not a problem because the copy constructor of Epetra_Map doesn't really copy the data (there is an reference count mecanism in Epetra).

      RCP<const Epetra_Map> map = rcp(new Epetra_Map(mtx_->RowMap()));
      return rcp ( new Cthulhu::EpetraMap(map) );
    }
     
    //! \brief Returns the Map that describes the column distribution in this matrix.
    const RCP<const Cthulhu::Map<int,int> > getColMap() const { 
      // The Epetra map is owned by the underlying Epetra matrix and freed when the matrix is deleted.
      // I have to make a copy of the map to be sure that the RCP<Map> returned by this method will remain valid even if the matrix is deleted.
      // Note that it is not a problem because the copy constructor of Epetra_Map doesn't really copy the data (there is an reference count mecanism in Epetra).

      RCP<const Epetra_Map> map = rcp(new Epetra_Map(mtx_->ColMap()));
      return rcp ( new Cthulhu::EpetraMap(map) );
    }

    //! Returns the CrsGraph associated with this matrix. 
    RCP<const CrsGraph<int,int> > getCrsGraph() const {  
      RCP<const Epetra_CrsGraph> const_graph = rcp(new Epetra_CrsGraph(mtx_->Graph()));

      RCP<Epetra_CrsGraph> graph = Teuchos::rcp_const_cast<Epetra_CrsGraph>(const_graph); //TODO: can I avoid the const_cast ?
      return rcp ( new Cthulhu::EpetraCrsGraph(graph) );
    }

    //! Returns the number of global rows in this matrix.
    global_size_t getGlobalNumRows() const {  return mtx_->NumGlobalRows(); }

    //! \brief Returns the number of global columns in the matrix.
    global_size_t getGlobalNumCols() const {  return mtx_->NumGlobalCols(); }

    //! Returns the number of matrix rows owned on the calling node.
    size_t getNodeNumRows() const {  return mtx_->NumMyRows(); }

    //! Returns the global number of entries in this matrix.
    global_size_t getGlobalNumEntries() const {  return mtx_->NumGlobalNonzeros(); }

    //! Returns the local number of entries in this matrix.
    size_t getNodeNumEntries() const {  return mtx_->NumMyNonzeros(); }

    //! Returns the current number of entries on this node in the specified local row.
    size_t getNumEntriesInLocalRow(int localRow) const {  return mtx_->NumMyEntries(localRow); }

    //! \brief Returns the number of global diagonal entries, based on global row/column index comparisons. 
    global_size_t getGlobalNumDiags() const {  return mtx_->NumGlobalDiagonals(); }

    //! \brief Returns the number of local diagonal entries, based on global row/column index comparisons. 
    size_t getNodeNumDiags() const {  return mtx_->NumMyDiagonals(); }

    //! \brief Returns the maximum number of entries across all rows/columns on all nodes.
    size_t getGlobalMaxNumRowEntries() const {  return mtx_->GlobalMaxNumEntries(); }

    //! \brief Returns the maximum number of entries across all rows/columns on this node.
    size_t getNodeMaxNumRowEntries() const {  return mtx_->MaxNumEntries(); }

    //! \brief If matrix indices are in the local range, this function returns true. Otherwise, this function returns false. */
    bool isLocallyIndexed() const {  return mtx_->IndicesAreLocal(); }

    //! \brief If matrix indices are in the global range, this function returns true. Otherwise, this function returns false. */
    bool isGloballyIndexed() const {  return mtx_->IndicesAreGlobal(); }

    //! Returns \c true if fillComplete() has been called and the matrix is in compute mode.
    bool isFillComplete() const {  return mtx_->Filled(); }

    //! Extract a list of entries in a specified local row of the matrix. Put into storage allocated by calling routine.
    //TODO: throw same exception as Tpetra
    void getLocalRowCopy(int LocalRow, const ArrayView<int> &Indices, const ArrayView<double> &Values, size_t &NumEntries) const;

    //! Extract a const, non-persisting view of global indices in a specified row of the matrix.
    void getGlobalRowView(int GlobalRow, ArrayView<const int> &indices, ArrayView<const double> &values) const;

    //! Extract a const, non-persisting view of local indices in a specified row of the matrix.
    void getLocalRowView(int LocalRow, ArrayView<const int> &indices, ArrayView<const double> &values) const;

    //! \brief Get a copy of the diagonal entries owned by this node, with local row indices.
    void getLocalDiagCopy(Vector<double,int,int> &diag) const { 
      CTHULHU_DYNAMIC_CAST(EpetraVector, diag, eDiag, "Cthulhu::EpetraCrsMatrix.getLocalDiagCopy() only accept Cthulhu::EpetraVector as input arguments.");
      mtx_->ExtractDiagonalCopy(*eDiag.getEpetra_Vector()); 
    }

    //@}

    //! @name Advanced Matrix-vector multiplication and solve methods
    //@{

    //! Multiplies this matrix by a MultiVector.
    // TODO Note: Do we need to use a Tpetra::CrsMatrixMultiplyOp ?? 
    //            (Epetra Doc of multiply: it is recommended that multiply() not be called directly; instead, use the CrsMatrixMultiplyOp, as it will handle the import/exprt operations required to apply a matrix with non-trivial communication needs.)

    // TODO : templated type
    void multiply(const MultiVector<double,int,int> & X, MultiVector<double,int,int> &Y, Teuchos::ETransp trans, double alpha, double beta) const {
       

      TEST_FOR_EXCEPTION((alpha != 1) || (beta != 0), Cthulhu::Exceptions::NotImplemented, "Cthulhu::EpetraCrsMatrix.multiply() only accept alpha==1 and beta==0");
      
      CTHULHU_DYNAMIC_CAST(const EpetraMultiVector, X, eX, "Cthulhu::EpetraCrsMatrix->multiply() only accept Cthulhu::EpetraMultiVector as input arguments.");
      CTHULHU_DYNAMIC_CAST(      EpetraMultiVector, Y, eY, "Cthulhu::EpetraCrsMatrix->multiply() only accept Cthulhu::EpetraMultiVector as input arguments.");

      TEST_FOR_EXCEPTION((trans != Teuchos::NO_TRANS) && (trans == Teuchos::TRANS), Cthulhu::Exceptions::NotImplemented, "Cthulhu::EpetraCrsMatrix->multiply() only accept trans == NO_TRANS or trans == TRANS");
      bool eTrans = Teuchos2Epetra_Trans(trans);

      CTHULHU_ERR_CHECK(mtx_->Multiply(eTrans, *eX.getEpetra_MultiVector(), *eY.getEpetra_MultiVector()));
    }
    
    //@}

    //! @name Methods implementing Operator
    //@{ 

    //! \brief Computes the sparse matrix-multivector multiplication.
    void apply(const MultiVector<double,int,int> & X, MultiVector<double,int,int> &Y, 
                      Teuchos::ETransp mode = Teuchos::NO_TRANS,
                      double alpha = ScalarTraits<double>::one(),
                      double beta = ScalarTraits<double>::zero()) const { 
       

      TEST_FOR_EXCEPTION((alpha != 1) || (beta != 0), Cthulhu::Exceptions::NotImplemented, "Cthulhu::EpetraCrsMatrix.multiply() only accept alpha==1 and beta==0");
      
      CTHULHU_DYNAMIC_CAST(const EpetraMultiVector, X, eX, "Cthulhu::EpetraCrsMatrix->apply() only accept Cthulhu::EpetraMultiVector as input arguments.");
      CTHULHU_DYNAMIC_CAST(      EpetraMultiVector, Y, eY, "Cthulhu::EpetraCrsMatrix->apply() only accept Cthulhu::EpetraMultiVector as input arguments.");

      TEST_FOR_EXCEPTION((mode != Teuchos::NO_TRANS) && (mode != Teuchos::TRANS), Cthulhu::Exceptions::NotImplemented, "Cthulhu::EpetraCrsMatrix->apply() only accept mode == NO_TRANS or mode == TRANS");
      bool eTrans = Teuchos2Epetra_Trans(mode);

      // /!\ UseTranspose value
      TEST_FOR_EXCEPTION(mtx_->UseTranspose(), Cthulhu::Exceptions::NotImplemented, "An exception is throw to let you know that Cthulhu::EpetraCrsMatrix->apply() do not take into account the UseTranspose() parameter of Epetra_CrsMatrix.");
      
      CTHULHU_ERR_CHECK(mtx_->Multiply(eTrans, *eX.getEpetra_MultiVector(), *eY.getEpetra_MultiVector()));
    }

    //! \brief Returns the Map associated with the domain of this operator.
    //! This will be <tt>null</tt> until fillComplete() is called.
    const RCP<const Map<int,int> > getDomainMap() const { 
      // The Epetra map is owned by the underlying Epetra matrix and freed when the matrix is deleted.
      // I have to make a copy of the map to be sure that the RCP<Map> returned by this method will remain valid even if the matrix is deleted.
      // Note that it is not a problem because the copy constructor of Epetra_Map doesn't really copy the data (there is an reference count mecanism in Epetra).

      RCP<const Epetra_Map> map = rcp(new Epetra_Map(mtx_->DomainMap()));
      return rcp ( new Cthulhu::EpetraMap(map) );
    }
    
    //! Returns the Map associated with the domain of this operator.
    //! This will be <tt>null</tt> until fillComplete() is called.
    const RCP<const Map<int,int> > getRangeMap() const { 
      // The Epetra map is owned by the underlying Epetra matrix and freed when the matrix is deleted.
      // I have to make a copy of the map to be sure that the RCP<Map> returned by this method will remain valid even if the matrix is deleted.
      // Note that it is not a problem because the copy constructor of Epetra_Map doesn't really copy the data (there is an reference count mecanism in Epetra).

      RCP<const Epetra_Map> map = rcp(new Epetra_Map(mtx_->RangeMap()));
      return rcp ( new Cthulhu::EpetraMap(map) );
    }

    //@}

    //! @name Overridden from Teuchos::Describable 
    //@{

    //! Return a simple one-line description of this object. */
    std::string description() const;
    
    //! Print the object with some verbosity level to an FancyOStream object. */
    void describe(Teuchos::FancyOStream &out, const Teuchos::EVerbosityLevel verbLevel=Teuchos::Describable::verbLevel_default) const;    
    //@}

    RCP<const Epetra_CrsMatrix> getEpetra_CrsMatrix() const {  return mtx_; }

    RCP<Epetra_CrsMatrix> getEpetra_CrsMatrixNonConst() const {  return mtx_; }

    /** TODO : interface of Teuchos_LabeledObject.hpp **/
    void setObjectLabel (const std::string &objectLabel) {  //mtx_->setObjectLabel(objectLabel); TODO
    }

    //@{
    // Implements DistObject interface

    const Teuchos::RCP<const Map<int,int> > getMap() const { 
      RCP<const Epetra_BlockMap> map = rcp(new Epetra_BlockMap(mtx_->Map()));
      return rcp ( new Cthulhu::EpetraMap(map) );
    }
    
    void doImport(const DistObject<char, int, int> &source, const Import<int, int> &importer, CombineMode CM);

    void doExport(const DistObject<char, int, int> &dest, const Import<int, int>& importer, CombineMode CM);

    void doImport(const DistObject<char, int, int> &source, const Export<int, int>& exporter, CombineMode CM);

    void doExport(const DistObject<char, int, int> &dest, const Export<int, int>& exporter, CombineMode CM);

    //@}

  private:
    
    RCP<Epetra_CrsMatrix> mtx_;

  }; // class EpetraCrsMatrix

} // namespace Cthulhu

#endif
