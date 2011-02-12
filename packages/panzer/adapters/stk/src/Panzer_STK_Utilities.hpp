#ifndef __Panzer_STK_Utilities_hpp__
#define __Panzer_STK_Utilities_hpp__

#include "Panzer_STK_Interface.hpp"

#include "Epetra_Vector.h"
#include "Epetra_MultiVector.h"

namespace panzer {
   template <typename LO,typename GO> class DOFManager;
}

namespace panzer_stk { 

void write_solution_data(const panzer::DOFManager<int,int> & dofMngr,panzer_stk::STK_Interface & mesh,const Epetra_MultiVector & x);
void write_solution_data(const panzer::DOFManager<int,int> & dofMngr,panzer_stk::STK_Interface & mesh,const Epetra_Vector & x);

void read_solution_data(const panzer::DOFManager<int,int> & dofMngr,const panzer_stk::STK_Interface & mesh,Epetra_MultiVector & x);
void read_solution_data(const panzer::DOFManager<int,int> & dofMngr,const panzer_stk::STK_Interface & mesh,Epetra_Vector & x);

}

#endif
