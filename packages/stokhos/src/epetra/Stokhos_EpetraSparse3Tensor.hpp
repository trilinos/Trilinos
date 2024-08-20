// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef STOKHOS_EPETRA_SPARSE_3_TENSOR_HPP
#define STOKHOS_EPETRA_SPARSE_3_TENSOR_HPP

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Stokhos_OrthogPolyBasis.hpp"
#include "Stokhos_Sparse3Tensor.hpp"
#include "EpetraExt_MultiComm.h"
#include "Epetra_BlockMap.h"
#include "Epetra_CrsGraph.h"

namespace Stokhos {

  class EpetraSparse3Tensor {
  public:

    //! Constructor from a full Cijk
    EpetraSparse3Tensor(
      const Teuchos::RCP<const Stokhos::OrthogPolyBasis<int,double> >& basis,
      const Teuchos::RCP<const Stokhos::Sparse3Tensor<int,double> >& Cijk,
      const Teuchos::RCP<const EpetraExt::MultiComm>& globalMultiComm,
      int k_begin = 0, int k_end = -1);

    //! Constructor from an already parallelized Cijk
    EpetraSparse3Tensor(
      const Teuchos::RCP<const Stokhos::OrthogPolyBasis<int,double> >& basis,
      const Teuchos::RCP<const Stokhos::Sparse3Tensor<int,double> >& Cijk,
      const Teuchos::RCP<const EpetraExt::MultiComm>& globalMultiComm,
      const Teuchos::RCP<const Epetra_BlockMap>& stoch_row_map,
      const Teuchos::RCP<const Stokhos::Sparse3Tensor<int,double> >& Cijk_parallel = Teuchos::null,
      int k_begin = 0, int k_end = -1);

    //! Copy constructor with possible change in scaling and k limits
    EpetraSparse3Tensor(const EpetraSparse3Tensor& epetraCijk,
			int k_begin_ = 0, int k_end_ = -1);

    //! Destructor
    ~EpetraSparse3Tensor() {}

    //! Rebalance maps and graph using Isorropia
    void rebalance(Teuchos::ParameterList& isorropia_params);

    //! Transform Cijk to local i and j indices
    void transformToLocal();

    //! Return k_begin index
    int getKBegin() const { return k_begin; }

    //! Return k_end index
    int getKEnd() const { return k_end; }

    //! Return whether stochastic blocks are parallel distributed
    bool isStochasticParallel() const { return is_parallel; }

    //! Return global row id for given local row id
    int GRID(int lrid) const { return stoch_row_map->GID(lrid); }

    //! Return global column id for given local column id
    int GCID(int lcid) const { return stoch_col_map->GID(lcid); }

    //! Return true if global row id is on processor
    bool myGRID(int grid) const { return stoch_row_map->MyGID(grid); }

    //! Return true if global column id is on processor
    bool myGCID(int gcid) const { return stoch_col_map->MyGID(gcid); }

    //! Return number of rows on this processor
    int numMyRows() const { return stoch_row_map->NumMyElements(); }

    //! Return number of columns on this processor
    int numMyCols() const { return stoch_col_map->NumMyElements(); }

    //! Get global comm
    Teuchos::RCP<const EpetraExt::MultiComm> 
    getMultiComm() const { return globalMultiComm; }

    //! Get stochastic comm
    Teuchos::RCP<const Epetra_Comm> 
    getStochasticComm() const { return stoch_comm; }

    //! Get stochastic row map
    Teuchos::RCP<const Epetra_BlockMap> 
    getStochasticRowMap() const { return stoch_row_map; }

    //! Get stochastic column map
    Teuchos::RCP<const Epetra_BlockMap> 
    getStochasticColMap() const { return stoch_col_map; }

    //! Get Cijk
    Teuchos::RCP<const Stokhos::Sparse3Tensor<int,double> > 
    getCijk() const { return Cijk; }

    //! Get parallel Cijk
    Teuchos::RCP<const Stokhos::Sparse3Tensor<int,double> > 
    getParallelCijk() const { return Cijk_parallel; }

    //! Get stochastic graph
    Teuchos::RCP<const Epetra_CrsGraph>
    getStochasticGraph() const { return stoch_graph; }

  protected:

    typedef Stokhos::Sparse3Tensor<int,double> Cijk_type;

    //! Build parallel Cijk tensor from a parallel row map
    Teuchos::RCP<Cijk_type> buildParallelCijk() const;

  protected:

    //! Basis
    Teuchos::RCP<const Stokhos::OrthogPolyBasis<int,double> > basis;

    //! Triple product
    Teuchos::RCP<const Cijk_type> Cijk;

    //! Multi-comm
    Teuchos::RCP<const EpetraExt::MultiComm> globalMultiComm;

    //! Number of global stochastic blocks
    int num_global_stoch_blocks;

    //! Beginning of k index
    int k_begin;

    //! End of k index
    int k_end;

    //! Stochastic comm
    Teuchos::RCP<const Epetra_Comm> stoch_comm;

    //! Whether stochastic blocks are parallel
    bool is_parallel;

    //! Stochastic row-map
    Teuchos::RCP<const Epetra_BlockMap> stoch_row_map;

    //! Stochastic col-map
    Teuchos::RCP<const Epetra_BlockMap> stoch_col_map;

    //! Cijk tensor parallel over i
    Teuchos::RCP<const Cijk_type> Cijk_parallel;

    //! Stochastic operator graph
    Teuchos::RCP<const Epetra_CrsGraph> stoch_graph;

  }; // class EpetraSparse3Tensor

} // namespace Stokhos

#endif // STOKHOS_EPETRA_SPARSE_3_TENSOR_HPP
