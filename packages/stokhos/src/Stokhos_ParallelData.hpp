// @HEADER
// ***********************************************************************
// 
//                           Stokhos Package
//                 Copyright (2009) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Eric T. Phipps (etphipp@sandia.gov).
// 
// ***********************************************************************
// @HEADER

#ifndef STOKHOS_PARALLEL_DATA_HPP
#define STOKHOS_PARALLEL_DATA_HPP

#include "Teuchos_RCP.hpp"
#include "Stokhos_OrthogPolyBasis.hpp"
#include "Stokhos_Sparse3Tensor.hpp"
#include "EpetraExt_MultiComm.h"
#include "Teuchos_ParameterList.hpp"
#include "Stokhos_EpetraSparse3Tensor.hpp"

namespace Stokhos {

  class ParallelData {
  public:

    //! Constructor
    ParallelData(
      const Teuchos::RCP<const Stokhos::OrthogPolyBasis<int,double> >& basis,
      const Teuchos::RCP<const Stokhos::Sparse3Tensor<int,double> >& Cijk,
      const Teuchos::RCP<const Epetra_Comm>& globalComm,
      Teuchos::ParameterList& params);

    //! Constructor with globalMultiComm
    ParallelData(
      const Teuchos::RCP<const Stokhos::OrthogPolyBasis<int,double> >& basis,
      const Teuchos::RCP<const Stokhos::Sparse3Tensor<int,double> >& Cijk,
      const Teuchos::RCP<const EpetraExt::MultiComm>& globalMultiComm,
      Teuchos::ParameterList& params);

    //! Destructor
    ~ParallelData() {}

    //! Get global comm
    Teuchos::RCP<const EpetraExt::MultiComm> 
    getMultiComm() const { return globalMultiComm; }

    //! Get stochastic comm
    Teuchos::RCP<const Epetra_Comm> 
    getStochasticComm() const { return stoch_comm; }

    //! Get spatial comm
    Teuchos::RCP<const Epetra_Comm> 
    getSpatialComm() const { return spatial_comm; }

    //! Get Epetra Cijk
    Teuchos::RCP<const Stokhos::EpetraSparse3Tensor>
    getEpetraCijk() const { return epetraCijk; }

  protected:

    //! Multi-comm
    Teuchos::RCP<const EpetraExt::MultiComm> globalMultiComm;

    //! Stochastic comm
    Teuchos::RCP<const Epetra_Comm> stoch_comm;

    //! Spatial comm
    Teuchos::RCP<const Epetra_Comm> spatial_comm;

    //! Epetra Cijk
    Teuchos::RCP<Stokhos::EpetraSparse3Tensor> epetraCijk;

  }; // class ParallelData

  Teuchos::RCP<const EpetraExt::MultiComm> 
  buildMultiComm(const Epetra_Comm& globalComm,
		 int num_global_stochastic_blocks,
		 int num_spatial_procs = -1);

  Teuchos::RCP<const Epetra_Comm> 
  getSpatialComm(
    const Teuchos::RCP<const EpetraExt::MultiComm>& globalMultiComm);

  Teuchos::RCP<const Epetra_Comm> 
  getStochasticComm(
    const Teuchos::RCP<const EpetraExt::MultiComm>& globalMultiComm);
    

} // namespace Stokhos

#endif // STOKHOS_PARALLEL_DATA_HPP
