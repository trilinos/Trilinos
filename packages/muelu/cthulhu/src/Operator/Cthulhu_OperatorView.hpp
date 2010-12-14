#ifndef CTHULHU_OPERATORVIEW_HPP
#define CTHULHU_OPERATORVIEW_HPP

#include "Cthulhu_ConfigDefs.hpp"
#include "Cthulhu_BlockMap.hpp"

#include <Teuchos_Describable.hpp>

#include "Cthulhu_Debug.hpp"

/** \file Cthulhu_Operator.hpp

Declarations for the class Cthulhu::Operator.
*/
namespace Cthulhu {

  template <class LocalOrdinal  = int, 
            class GlobalOrdinal = LocalOrdinal, 
            class Node          = Kokkos::DefaultNode::DefaultNodeType>
  class OperatorView { // TODO : virtual public Teuchos::Describable {
    typedef Cthulhu::Map<LocalOrdinal, GlobalOrdinal, Node> Map;

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
    inline const RCP<const Map> & GetRowMap() const { return rowMap_; }
  
    //! \brief Returns the Map that describes the column distribution in this matrix.
    inline const RCP<const Map> & GetColMap() const { return colMap_; }
  
    //@}
  
  private:
    const RCP<const Map> rowMap_;
    const RCP<const Map> colMap_;

  }; // class OperatorView

} // namespace Cthulhu

#define CTHULHU_OPERATORVIEW_SHORT
#endif //CTHULHU_OPERATOR_VIEW_DECL_HPP
