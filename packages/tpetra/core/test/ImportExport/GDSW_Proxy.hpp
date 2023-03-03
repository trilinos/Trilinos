#ifndef GDSW_PROXY_HPP
#define GDSW_PROXY_HPP

#include <vector>
#include "Tpetra_CrsGraph.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_Distributor.hpp"




template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
class TpetraFunctions
{

  using SC = Scalar;
  using LO = LocalOrdinal;
  using GO = GlobalOrdinal;
  using NT = Node;
  using Map = Tpetra::Map<LO,GO,NT>;
  using CrsMatrix = Tpetra::CrsMatrix<SC,LO,GO,NT>;
  using IndicesViewT = typename CrsMatrix::local_inds_host_view_type;
  using ValuesViewT  = typename CrsMatrix::values_host_view_type;
  
public:

    TpetraFunctions() = default;

    ~TpetraFunctions() = default;

    void importSquareMatrix(Teuchos::RCP<const CrsMatrix> inputMatrix, 
                            Teuchos::RCP<const Map> outputRowMap, 
                            Teuchos::RCP<CrsMatrix> & outputMatrix);

private:

    void communicateMatrixData(Teuchos::RCP<const CrsMatrix> inputMatrix, 
                               Teuchos::RCP<const Map> rowMap,
                               Teuchos::RCP<Tpetra::Distributor> distributor,
                               const std::vector<GO> & targetMapGIDs, 
                               const std::vector<LO> & targetMapGIDsBegin,
                               const std::vector<GO> & ownedRowGIDs,
                               const std::vector<LO> & localRowsSend,
                               const std::vector<LO> & localRowsSendBegin,
                               const std::vector<LO> & localRowsRecv,
                               const std::vector<LO> & localRowsRecvBegin,
                               Teuchos::RCP<CrsMatrix> & outputMatrix);

    void communicateRowMap(Teuchos::RCP<const Map> rowMap, 
                           Teuchos::RCP<Tpetra::Distributor> distributor, 
                           std::vector<GO> & rowMapGIDs, 
                           std::vector<LO> & rowMapGIDsBegin);

    void constructDistributor(Teuchos::RCP<const CrsMatrix> inputMatrix, 
                              Teuchos::RCP<const Map> rowMap, 
                              Teuchos::RCP<Tpetra::Distributor> & distributor,
                              std::vector<GO> & ownedRowGIDs,
                              std::vector<LO> & localRowsSend,
                              std::vector<LO> & localRowsSendBegin,
                              std::vector<LO> & localRowsRecv,
                              std::vector<LO> & localRowsRecvBegin);

    void getUniqueEntries(const std::vector<int> & vector, 
                          std::vector<int> & vectorUnique);

};


#endif // GDSW_PROXY_HPP
