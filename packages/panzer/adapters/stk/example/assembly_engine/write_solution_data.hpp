#ifndef __write_solution_data_hpp__
#define __write_solution_data_hpp__

#include "Panzer_STK_Interface.hpp"

#include "Epetra_Vector.h"
#include "Intrepid_FieldContainer.hpp"

#include <map>
#include <vector>

namespace panzer { 
   template <typename LO,typename GO> class DOFManager;
}


void write_solution_data(const panzer::DOFManager<int,int> & dofMngr,panzer_stk::STK_Interface & mesh,const Epetra_Vector & x);

void gather_in_block(const std::string & blockId, const panzer::DOFManager<int,int> & dofMngr,
                     const Epetra_Vector & x,const std::vector<std::size_t> & localCellIds,
                     std::map<std::string,Intrepid::FieldContainer<double> > & fc);

void build_local_ids(const panzer_stk::STK_Interface & mesh,
                   std::map<std::string,Teuchos::RCP<std::vector<std::size_t> > > & localIds);

#endif
