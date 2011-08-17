#ifndef XPETRA_OPERATORVIEW_HPP
#define XPETRA_OPERATORVIEW_HPP

#include "Xpetra_ConfigDefs.hpp"

#include <Teuchos_Describable.hpp>

/** \file Xpetra_Operator.hpp

Declarations for the class Xpetra::Operator.
*/
namespace Xpetra {

  template <class LocalOrdinal  = int, 
            class GlobalOrdinal = LocalOrdinal, 
            class Node          = Kokkos::DefaultNode::DefaultNodeType>
  class OperatorView { // TODO : virtual public Teuchos::Describable {
    typedef Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> Map;

  public:
  
    //! @name Constructor/Destructor Methods
    //@{

    //! Constructor
    OperatorView(const RCP<const Map> &rowMap, const RCP<const Map> &colMap) 
      : rowMap_ (rowMap), colMap_(colMap)
    { }
  
    //! Destructor
    virtual ~OperatorView() {}
    
    //@}

    //! @name Map access methods
    //@{
    //! Returns the Map that describes the row distribution in this matrix.
    const RCP<const Map> & GetRowMap() const { return rowMap_; }
  
    //! \brief Returns the Map that describes the column distribution in this matrix.
    const RCP<const Map> & GetColMap() const { return colMap_; }
  
    //! Returns the Map that describes the row distribution in this matrix.
    void SetRowMap(const RCP<const Map> & rowMap) { rowMap_ = rowMap; }
  
    //! \brief Set the Map that describes the column distribution in this matrix.
    void SetColMap(const RCP<const Map> & colMap) { colMap_ = colMap; }
    //@}
  
  private:
    RCP<const Map> rowMap_;
    RCP<const Map> colMap_;

  }; // class OperatorView

} // namespace Xpetra

#define XPETRA_OPERATORVIEW_SHORT
#endif //XPETRA_OPERATOR_VIEW_DECL_HPP
