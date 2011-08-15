#include "Cthulhu_EpetraCrsGraph.hpp"

namespace Cthulhu {

  void EpetraCrsGraph::insertGlobalIndices(int globalRow, const ArrayView<const int> &indices) { 
    int* indices_rawPtr = const_cast<int*>(indices.getRawPtr()); // there is no const in the Epetra interface :(
    CTHULHU_ERR_CHECK(graph_->InsertGlobalIndices(globalRow, indices.size(), indices_rawPtr)); 
  }

  void EpetraCrsGraph::insertLocalIndices(int localRow, const ArrayView<const int> &indices) { 
    int* indices_rawPtr = const_cast<int*>(indices.getRawPtr()); // there is no const in the Epetra interface :(
    CTHULHU_ERR_CHECK(graph_->InsertMyIndices(localRow, indices.size(), indices_rawPtr)); 
  }

  void EpetraCrsGraph::getGlobalRowView(int GlobalRow, ArrayView<const int> &Indices) const { 
    int      numEntries;
    int    * eIndices;
      
    CTHULHU_ERR_CHECK(graph_->ExtractGlobalRowView(GlobalRow, numEntries, eIndices));
    if (numEntries == 0) { eIndices = NULL; } // Cf. TEST_FOR_EXCEPT( p == 0 && size_in != 0 ) in Teuchos ArrayView constructor.

    Indices = ArrayView<const int>(eIndices, numEntries);
  }

  void EpetraCrsGraph::getLocalRowView(int LocalRow, ArrayView<const int> &indices) const {
    int      numEntries;
    int    * eIndices;
      
    CTHULHU_ERR_CHECK(graph_->ExtractMyRowView(LocalRow, numEntries, eIndices));
    if (numEntries == 0) { eIndices = NULL; } // Cf. TEST_FOR_EXCEPT( p == 0 && size_in != 0 ) in Teuchos ArrayView constructor.

    indices = ArrayView<const int>(eIndices, numEntries);
  }

} // namespace Cthulhu
