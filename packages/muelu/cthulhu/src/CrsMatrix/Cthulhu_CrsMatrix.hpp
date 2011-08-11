#ifndef CTHULHU_CRSMATRIX_HPP
#define CTHULHU_CRSMATRIX_HPP

#include <Kokkos_DefaultNode.hpp>
#include <Kokkos_DefaultKernels.hpp>
#include <Kokkos_CrsMatrix.hpp>
#include "Cthulhu_ConfigDefs.hpp"
#include "Cthulhu_DistObject.hpp"
#include "Cthulhu_CrsGraph.hpp"
#include "Cthulhu_Vector.hpp"

namespace Cthulhu {

  template <class Scalar, class LocalOrdinal = int, class GlobalOrdinal = LocalOrdinal, class Node = Kokkos::DefaultNode::DefaultNodeType, class LocalMatOps = typename Kokkos::DefaultKernels<Scalar,LocalOrdinal,Node>::SparseOps>
  class CrsMatrix
    : public DistObject< char, LocalOrdinal, GlobalOrdinal, Node > {
//TODO public RowMatrix< Scalar, LocalOrdinal, GlobalOrdinal, Node >

  public:

    //! @name Constructor/Destructor Methods
    //@{ 

    // !Destructor.
    virtual ~CrsMatrix() { }

    //@}

    //! @name Insertion/Removal Methods
    //@{

    //! Insert matrix entries, using global IDs.
    virtual void insertGlobalValues(GlobalOrdinal globalRow, const ArrayView< const GlobalOrdinal > &cols, const ArrayView< const Scalar > &vals)= 0;

    //! Scale the current values of a matrix, this = alpha*this.
    virtual void scale(const Scalar &alpha)= 0;

    //@}

    //! @name Transformational Methods
    //@{

    //! Signal that data entry is complete, specifying domain and range maps.
    virtual void fillComplete(const RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &domainMap, const RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &rangeMap, OptimizeOption os=DoOptimizeStorage)= 0;

    //! Signal that data entry is complete.
    virtual void fillComplete(OptimizeOption os=DoOptimizeStorage)= 0;

    //@}

    //! @name Methods implementing RowMatrix
    //@{

    //! Returns the Map that describes the row distribution in this matrix.
    virtual const RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > >  getRowMap() const = 0;

    //! Returns the Map that describes the column distribution in this matrix.
    virtual const RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > >  getColMap() const = 0;

    //! Returns the CrsGraph associated with this matrix.
    virtual RCP< const CrsGraph< LocalOrdinal, GlobalOrdinal, Node, LocalMatOps > > getCrsGraph() const = 0;

    //! Number of global elements in the row map of this matrix.
    virtual global_size_t getGlobalNumRows() const = 0;

    //! Number of global columns in the matrix.
    virtual global_size_t getGlobalNumCols() const = 0;

    //! Returns the number of matrix rows owned on the calling node.
    virtual size_t getNodeNumRows() const = 0;

    //! Returns the global number of entries in this matrix.
    virtual global_size_t getGlobalNumEntries() const = 0;

    //! Returns the local number of entries in this matrix.
    virtual size_t getNodeNumEntries() const = 0;

    //! Returns the current number of entries on this node in the specified local row.
    virtual size_t getNumEntriesInLocalRow(LocalOrdinal localRow) const = 0;

    //! Returns the number of global diagonal entries, based on global row/column index comparisons.
    virtual global_size_t getGlobalNumDiags() const = 0;

    //! Returns the number of local diagonal entries, based on global row/column index comparisons.
    virtual size_t getNodeNumDiags() const = 0;

    //! Returns the maximum number of entries across all rows/columns on all nodes.
    virtual size_t getGlobalMaxNumRowEntries() const = 0;

    //! Returns the maximum number of entries across all rows/columns on this node.
    virtual size_t getNodeMaxNumRowEntries() const = 0;

    //! If matrix indices are in the local range, this function returns true. Otherwise, this function returns false.
    virtual bool isLocallyIndexed() const = 0;

    //! If matrix indices are in the global range, this function returns true. Otherwise, this function returns false.
    virtual bool isGloballyIndexed() const = 0;

    //! Returns true if fillComplete() has been called and the matrix is in compute mode.
    virtual bool isFillComplete() const = 0;

    //! Extract a list of entries in a specified local row of the matrix. Put into storage allocated by calling routine.
    virtual void getLocalRowCopy(LocalOrdinal LocalRow, const ArrayView< LocalOrdinal > &Indices, const ArrayView< Scalar > &Values, size_t &NumEntries) const = 0;

    //! Extract a const, non-persisting view of global indices in a specified row of the matrix.
    virtual void getGlobalRowView(GlobalOrdinal GlobalRow, ArrayView< const GlobalOrdinal > &indices, ArrayView< const Scalar > &values) const = 0;

    //! Extract a const, non-persisting view of local indices in a specified row of the matrix.
    virtual void getLocalRowView(LocalOrdinal LocalRow, ArrayView< const LocalOrdinal > &indices, ArrayView< const Scalar > &values) const = 0;

    //! Get a copy of the diagonal entries owned by this node, with local row idices.
    virtual void getLocalDiagCopy(Vector< Scalar, LocalOrdinal, GlobalOrdinal, Node > &diag) const = 0;

    //@}

    //! @name Methods implementing Operator
    //@{

    //! Computes the sparse matrix-multivector multiplication.
    virtual void apply(const MultiVector< Scalar, LocalOrdinal, GlobalOrdinal, Node > &X, MultiVector< Scalar, LocalOrdinal, GlobalOrdinal, Node > &Y, Teuchos::ETransp mode=Teuchos::NO_TRANS, Scalar alpha=ScalarTraits< Scalar >::one(), Scalar beta=ScalarTraits< Scalar >::zero()) const = 0;

    //! Returns the Map associated with the domain of this operator. This will be null until fillComplete() is called.
    virtual const RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > >  getDomainMap() const = 0;

    //! 
    virtual const RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > >  getRangeMap() const = 0;

    //@}

    //! @name Overridden from Teuchos::Describable
    //@{

    //! Return a simple one-line description of this object.
    virtual std::string description() const = 0;

    //! Print the object with some verbosity level to an FancyOStream object.
    virtual void describe(Teuchos::FancyOStream &out, const Teuchos::EVerbosityLevel verbLevel=Teuchos::Describable::verbLevel_default) const = 0;

    //@}

  }; // CrsMatrix class

} // Cthulhu namespace

#define CTHULHU_CRSMATRIX_SHORT
#endif // CTHULHU_CRSMATRIX_HPP
