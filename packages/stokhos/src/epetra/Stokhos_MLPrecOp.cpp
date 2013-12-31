// $Id$ 
// $Source$ 
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

#include "Stokhos_MLPrecOp.hpp"

#ifdef HAVE_STOKHOS_ML

//==============================================================================
// constructor -- it's presumed that the user has constructed the ML
// object somewhere else.  Uses AMG to invert the mean stiffness matrix.
Stokhos::MLPrecOp::MLPrecOp(const Epetra_CrsMatrix& mean_op,const Teuchos::Array<double>& norms, const Epetra_Comm& Comm, const Epetra_Map& DMap,const Epetra_Map& RMap)
  : //Epetra_Operator(),
    Comm_(Comm),
    Label_(0),
    DomainMap_(DMap),
    RangeMap_(RMap),
    norms_(norms)
  {
  Label_ = "Stochos::MLPrecOp";
  

  // create a parameter list for ML options
  Teuchos::ParameterList MLList;

  // Sets default parameters for classic smoothed aggregation. After this
  // call, MLList contains the default values for the ML parameters,
  // as required by typical smoothed aggregation for symmetric systems.
  // Other sets of parameters are available for non-symmetric systems
  // ("DD" and "DD-ML"), and for the Maxwell equations ("maxwell").
  ML_Epetra::SetDefaults("SA",MLList);
  
  // overwrite some parameters. Please refer to the user's guide
  // for more information
  // some of the parameters do not differ from their default value,
  // and they are here reported for the sake of clarity
  
  // output level, 0 being silent and 10 verbose
  MLList.set("ML output", 10);
  // maximum number of levels
  MLList.set("max levels",5);
  // set finest level to 0
  MLList.set("increasing or decreasing","increasing");

  // use Uncoupled scheme to create the aggregate
  MLList.set("aggregation: type", "Uncoupled");

  // smoother is Chebyshev. Example file 
  // `ml/examples/TwoLevelDD/ml_2level_DD.cpp' shows how to use
  // AZTEC's preconditioners as smoothers

  MLList.set("smoother: type","Chebyshev");
  MLList.set("smoother: sweeps",3);

  // use both pre and post smoothing
  MLList.set("smoother: pre or post", "both");
  //MLList.set("coarse: max size", );

#ifdef HAVE_ML_AMESOS
  // solve with serial direct solver KLU
  MLList.set("coarse: type","Amesos-KLU");
#else
  // this is for testing purposes only, you should have 
  // a direct solver for the coarse problem (either Amesos, or the SuperLU/
  // SuperLU_DIST interface of ML)
  MLList.set("coarse: type","Jacobi");
#endif

  // Creates the preconditioning object. We suggest to use `new' and
  // `delete' because the destructor contains some calls to MPI (as
  // required by ML and possibly Amesos). This is an issue only if the
  // destructor is called **after** MPI_Finalize().
  
  MLPrec = new ML_Epetra::MultiLevelPreconditioner(mean_op, MLList);

  // verify unused parameters on process 0 (put -1 to print on all
  // processes)
  MLPrec->PrintUnused(0);

  
}

//==============================================================================
//Applies the preconditioner implicitly to the global stiffness matrix.
//==============================================================================

int Stokhos::MLPrecOp::ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const{
  
   if (!X.Map().SameAs(OperatorDomainMap())) 
    std::cout << "!X.Map().SameAs(OperatorDomainMap())\n";
  if (!Y.Map().SameAs(OperatorRangeMap())) 
    std::cout<< "!Y.Map().SameAs(OperatorRangeMap())\n";
  if (Y.NumVectors()!=X.NumVectors()) 
    std::cout<< "Y.NumVectors()!=X.NumVectors()\n";

  

  for(int mm = 0; mm< X.NumVectors(); mm++){
  
    
    int N_xi = norms_.size();
    int N_x = X.MyLength()/N_xi;

    // Construct a Map with NumElements and index base of 0
    //Epetra_Map::Epetra_Map Map(N_x, 0, Comm_);
    Epetra_Map Map(MLPrec->OperatorDomainMap());
    
    // Form x and y into block vectors.
    Epetra_MultiVector xBlock(Map,N_xi);
    Epetra_MultiVector yBlock(MLPrec->OperatorRangeMap(),N_xi);
 
    // get the local size of the vector
    int MyLength = xBlock.MyLength();
    for( int c=0; c<N_xi ; c++){
      for( int i=0; i<MyLength; i++){
        xBlock[c][i] = (X)[mm][c*N_x + i];
        yBlock[c][i] = 0;
      }
    }
    
    Epetra_MultiVector blockProducts(MLPrec->OperatorRangeMap(),N_xi);
    MLPrec->ApplyInverse(xBlock, blockProducts);
       
  
    for(int j = 0; j<N_xi; j++){
      (*yBlock(j)).Update(1/norms_[j],*blockProducts(j), 1.0);
    }
      
    for( int c=0; c<N_xi ; c++){
      for( int i=0; i<MyLength; i++){
        (Y)[mm][c*N_x + i] = yBlock[c][i];
      }
    }
  }
  return 1; 
}

#endif




