#ifndef CTHULHU_OPERATOR_HPP
#define CTHULHU_OPERATOR_HPP

#include <Kokkos_DefaultNode.hpp>
#include <Kokkos_DefaultKernels.hpp>

#include "Cthulhu_ConfigDefs.hpp"
#include "Cthulhu_Exceptions.hpp"

#include "Cthulhu_BlockMap.hpp"
#include "Cthulhu_MultiVector.hpp"
#include "Cthulhu_CrsMatrix.hpp"
#include "Cthulhu_TpetraCrsMatrix.hpp" //TMP
#include "Cthulhu_CrsMatrixFactory.hpp"
#include "Cthulhu_OperatorView.hpp"

#include <Teuchos_SerialDenseMatrix.hpp>
#include <Teuchos_Hashtable.hpp>

#include "Cthulhu_Debug.hpp"

/** \file Cthulhu_Operator.hpp

  Declarations for the class Cthulhu::Operator.
*/
namespace Cthulhu {

  typedef std::string viewLabel_t;

template <class ScalarType, 
          class LocalOrdinal  = int, 
          class GlobalOrdinal = LocalOrdinal, 
          class Node          = Kokkos::DefaultNode::DefaultNodeType, 
          class LocalMatOps   = typename Kokkos::DefaultKernels<ScalarType,LocalOrdinal,Node>::SparseOps > //TODO: or BlockSparseOp ?
class Operator : virtual public Teuchos::Describable {

  typedef Cthulhu::Map<LocalOrdinal, GlobalOrdinal, Node> Map;
  typedef Cthulhu::CrsMatrix<ScalarType, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> CrsMatrix;
  typedef Cthulhu::TpetraCrsMatrix<ScalarType, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> TpetraCrsMatrix;
  typedef Cthulhu::CrsMatrixFactory<ScalarType, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> CrsMatrixFactory;
  typedef Cthulhu::OperatorView<LocalOrdinal, GlobalOrdinal, Node> OperatorView;

public:
  
  //! @name Constructor/Destructor Methods
  //@{

  //! Constructor
  Operator(const RCP<const Map> &rowMap, size_t maxNumEntriesPerRow, Cthulhu::ProfileType pftype = Cthulhu::DynamicProfile) 
  {
    // Set matrix data
    matrixData_ = CrsMatrixFactory::Build(rowMap, maxNumEntriesPerRow, pftype);
    
    // Create default view
    defaultViewLabel_ = "point";
    CreateView(GetDefaultViewLabel(), matrixData_->getRowMap(), matrixData_->getColMap());    

    // Set current view
    currentViewLabel_ = GetDefaultViewLabel();
    SwitchToDefaultView(); // JG: this line is useless (for the moment).
  }
  
  //! Destructor
  virtual ~Operator() {}

  //@}

  //! @name Methods
  //@{
  void CreateView(viewLabel_t viewLabel, const RCP<const Map> & rowMap, const RCP<const Map> & colMap) {
    TEST_FOR_EXCEPTION(operatorViewTable_.containsKey(viewLabel) == true, Cthulhu::Exceptions::RuntimeError, "Cthulhu::Operator.CreateView(): a view labeled '" + viewLabel + "' already exist.");
    RCP<OperatorView> view = rcp(new OperatorView(rowMap, colMap));
    operatorViewTable_.put(viewLabel, view);
  }
    
  void RemoveView(const viewLabel_t viewLabel) {
    TEST_FOR_EXCEPTION(operatorViewTable_.containsKey(viewLabel) == false, Cthulhu::Exceptions::RuntimeError, "Cthulhu::Operator.RemoveView(): view '" + viewLabel + "' does not exist.");
    TEST_FOR_EXCEPTION(viewLabel == GetDefaultViewLabel(), Cthulhu::Exceptions::RuntimeError, "Cthulhu::Operator.RemoveView(): view '" + viewLabel + "' is the default view and cannot be removed.");
    operatorViewTable_.remove(viewLabel);
  }
    
  const viewLabel_t SwitchToView(const viewLabel_t viewLabel) {
    TEST_FOR_EXCEPTION(operatorViewTable_.containsKey(viewLabel) == false, Cthulhu::Exceptions::RuntimeError, "Cthulhu::Operator.SwitchToView(): view '" + viewLabel + "' does not exist.");
    viewLabel_t oldViewLabel = GetCurrentViewLabel();
    currentViewLabel_ = viewLabel;
    return oldViewLabel;
  }

  inline const viewLabel_t SwitchToDefaultView() { return SwitchToView(GetDefaultViewLabel()); }

  inline const viewLabel_t & GetDefaultViewLabel() const { return defaultViewLabel_; }
  inline const viewLabel_t & GetCurrentViewLabel() const { return currentViewLabel_; }

  //! @name Methods working on current view
  //@{
  inline const RCP<const Map> & GetRowMap() const { return GetRowMap(GetCurrentViewLabel()); }
  inline const RCP<const Map> & GetColMap() const { return GetColMap(GetCurrentViewLabel()); }
  
  //@}

  //! @name Methods working on specific view
  //@{ 
  const RCP<const Map> & GetRowMap(viewLabel_t viewLabel) const { 
    TEST_FOR_EXCEPTION(operatorViewTable_.containsKey(viewLabel) == false, Cthulhu::Exceptions::RuntimeError, "Cthulhu::Operator.GetRowMap(): view '" + viewLabel + "' does not exist.");
    return operatorViewTable_.get(viewLabel)->GetRowMap(); 
  }

  const RCP<const Map> & GetColMap(viewLabel_t viewLabel) const { 
    TEST_FOR_EXCEPTION(operatorViewTable_.containsKey(viewLabel) == false, Cthulhu::Exceptions::RuntimeError, "Cthulhu::Operator.GetColMap(): view '" + viewLabel + "' does not exist.");
    return operatorViewTable_.get(viewLabel)->GetColMap(); 
  }
  
  //@}

  //! @name Blablabla
  //@{ 

  inline void fillComplete(Cthulhu::OptimizeOption os = Cthulhu::DoOptimizeStorage) { matrixData_->fillComplete(); } //TODO:os arg

  //@}

  //! @name Overridden from Teuchos::Describable 
  //@{
  
  /** \brief Return a simple one-line description of this object. */
  inline std::string description() const { 
    std::ostringstream oss;
    return oss.str();
  }
  
  /** \brief Print the object with some verbosity level to an FancyOStream object. */
  inline void describe(Teuchos::FancyOStream &out, const Teuchos::EVerbosityLevel verbLevel=Teuchos::Describable::verbLevel_default) const { 
    //     Teuchos::EVerbosityLevel vl = verbLevel;
    //     if (vl == VERB_DEFAULT) vl = VERB_LOW;
    //     RCP<const Comm<int> > comm = this->getComm();
    //     const int myImageID = comm->getRank(),
    //       numImages = comm->getSize();
    
    //     if (myImageID == 0) out << this->description() << std::endl; 
    
    // Teuchos::OSTab tab(out);
  }
  
  //@}
  
private:
  RCP<CrsMatrix> matrixData_;

  Teuchos::Hashtable<viewLabel_t, RCP<OperatorView> > operatorViewTable_; // hashtable storing the operator views (keys = view names, values = views).
  
  viewLabel_t defaultViewLabel_;  // label of the view associated with inital Operator construction
  viewLabel_t currentViewLabel_;  // label of the current view
  
}; //class Operator

} //namespace Cthulhu

#endif //CTHULHU_OPERATOR_DECL_HPP
