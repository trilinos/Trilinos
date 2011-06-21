#ifndef AMESOS2_EPETRACRSMATRIX_MATRIXADAPTER_DEF_HPP
#define AMESOS2_EPETRACRSMATRIX_MATRIXADAPTER_DEF_HPP

#include <Epetra_LocalMap.h>
#include <Epetra_Import.h>

namespace Amesos {

  ConcreteMatrixAdapter<Epetra_CrsMatrix>::ConcreteMatrixAdapter(RCP<Epetra_CrsMatrix> m) 
    : AbstractConcreteMatrixAdapter<Epetra_RowMatrix,Epetra_CrsMatrix>(m) // CrsMatrix inherits from RowMatrix virtually, so a dynamic cast is necessary
    {}

  Teuchos::RCP<const MatrixAdapter<Epetra_CrsMatrix> >
  ConcreteMatrixAdapter<Epetra_CrsMatrix>::get_impl(EDistribution d) const
    {
      using Teuchos::as;
      using Teuchos::rcp;
      using Teuchos::rcpFromRef;
      
      RCP<const Epetra_Map> o_map, l_map;
      o_map = rcpFromRef(this->mat_->RowMap());
      int numrows = as<int>(this->getGlobalNumRows());
      
      if( d == Util::Rooted ){
	if( Teuchos::size(*comm_) == 1 ){ // then mat_ (and this adapter) is already "rooted"
	  return( rcpFromRef(*this) );
	}

	int num_my_elements;
	if( this->getComm()->getRank() == 0 ){
	  num_my_elements = numrows;
	} else {
	  num_my_elements = Teuchos::OrdinalTraits<int>::zero();
	}
	l_map = rcp(new Epetra_Map( -1,	// do not create a locally replicated map but compute numglobal 
				    num_my_elements,
				    o_map->IndexBase(),
				    this->mat_->Comm() ));
      } else if( d == Util::Globally_Replicated ){
	l_map = rcp(new Epetra_LocalMap( numrows,
					 o_map->IndexBase(),
					 this->mat_->Comm() ));
      }
      
      RCP<Epetra_CrsMatrix> l_mat = rcp(new Epetra_CrsMatrix(Copy, *l_map, this->getMaxRowNNZ()));
      
      Epetra_Import importer(*l_map, *o_map);
      
      l_mat->Import(*(this->mat_), importer, Insert);

      l_mat->FillComplete();	// Must be in local form for later extraction of rows
      
      return( rcp(new ConcreteMatrixAdapter<Epetra_CrsMatrix>(l_mat)) );
    }

} // end namespace Amesos

#endif	// AMESOS2_EPETRACRSMATRIX_MATRIXADAPTER_DEF_HPP
