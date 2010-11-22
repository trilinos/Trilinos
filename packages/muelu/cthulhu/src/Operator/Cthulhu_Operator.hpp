#ifndef CTHULHU_OPERATOR_HPP
#define CTHULHU_OPERATOR_HPP

#include <Kokkos_DefaultNode.hpp>
#include <Kokkos_DefaultKernels.hpp>
#include <Kokkos_Operator.hpp>

#include "Cthulhu_ConfigDefs.hpp"
// #include "Cthulhu_Operator.hpp"
#include "Cthulhu_BlockMap.hpp"
#include "Cthulhu_MultiVector.hpp"
// #include "Cthulhu_BlockCrsGraph.hpp"

#include <Teuchos_SerialDenseMatrix.hpp>

#include "Cthulhu_Debug.hpp"

/** \file Cthulhu_Operator.hpp

  Declarations for the class Cthulhu::Operator.
*/
namespace Cthulhu {

template <class Scalar, 
          class LocalOrdinal  = int, 
          class GlobalOrdinal = LocalOrdinal, 
          class Node          = Kokkos::DefaultNode::DefaultNodeType, 
          class LocalMatOps   = typename Kokkos::DefaultKernels<Scalar,LocalOrdinal,Node>::BlockSparseOps >
class Operator {
 public:

  //! @name Constructor/Destructor Methods
  //@{

  //! Destructor
  virtual ~Operator();


  

};//class Operator

} //namespace Cthulhu

#endif //CTHULHU_OPERATOR_DECL_HPP
