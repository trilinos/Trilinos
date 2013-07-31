// @HEADER
// ***********************************************************************
// 
//                           Stokhos Package
//                 Copyright (2009) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
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
