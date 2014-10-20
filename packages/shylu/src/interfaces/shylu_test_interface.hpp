
#ifndef INTERFACE_HPP
#define INTERFACE_HPP

/*--------------------------Trilinos Depend includes----------------*/

//Shylu
#include "shylu.h"
#include "shylu_util.h"
/*Shylu_util includes
assert, mpi, iostream, sstream,
Isorropia_config, etc
*/

//Epetra
#include "Epetra_CrsMatrix.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Vector.h"
#include "Epetra_Map.h"
//
#include "Isorropia_config.h"
#include "Isorropia_Epetra.hpp"
#include "Isorropia_EpetraProber.hpp"
#include "Isorropia_EpetraPartitioner.hpp"
#include "Isorropia_EpetraRedistributor.hpp"

#ifdef HAVE_TPETRA
#include "Tpetra_CrsMatrix_decl.hpp"
#incldue "Tpetra_CrsMatrix_def.hpp"
#endif




namespace shylu_interface{

  using std::cout;
  using std::endl;
  
  using Teuchos::ParameterList;


  int partitioning_interface(Epetra_CrsMatrix *A, Epetra_CrsMatrix *AHat, Epetra_MultiVector *b,
			     Epetra_MultiVector *bHat, ParameterList shyLUList);

 #ifdef HAVE_TPETRA
  template <class ST, class LT, class GT>
  int partitioning_interface(Tpetra::CrsMatrix<ST,LT, GT>, Tpetra::CrsMatrix<ST,LT,GT>,
			     Tpetra::MultiVector<ST,LT,GT>, Tpetra::Multiector<ST,LT,GT>,
			     ParameterList shyLUList);
#endif
  

}

#endif
