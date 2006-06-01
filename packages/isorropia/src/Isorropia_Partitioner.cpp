//@HEADER
/*
************************************************************************

              Isorropia: Partitioning and Load Balancing Package
                Copyright (2006) Sandia Corporation

Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
license for use of this work by or on behalf of the U.S. Government.

This library is free software; you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as
published by the Free Software Foundation; either version 2.1 of the
License, or (at your option) any later version.

This library is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
USA
Questions? Contact Alan Williams (william@sandia.gov)
                or Erik Boman    (egboman@sandia.gov)

************************************************************************
*/
//@HEADER

#include <Isorropia_Partitioner.hpp>
#include <Isorropia_Zoltan_Rebalance.hpp>
#include <Isorropia_Exception.hpp>
#include <Isorropia_Epetra_utils.hpp>

#include <Teuchos_RefCountPtr.hpp>

#ifdef HAVE_EPETRA
#include <Epetra_Map.h>
#include <Epetra_Import.h>
#include <Epetra_Vector.h>
#include <Epetra_MultiVector.h>
#include <Epetra_CrsGraph.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_LinearProblem.h>
#endif

/** Isorropia is the namespace that contains isorropia's declarations
  for classes and functions.
*/
namespace Isorropia {

#ifdef HAVE_EPETRA

  Partitioner::Partitioner(Teuchos::RefCountPtr<const Epetra_CrsGraph> input_graph,
			   Teuchos::RefCountPtr<Teuchos::ParameterList> paramlist)
  : new_map_(),
    input_graph_(input_graph),
    paramlist_(paramlist)
{
}

Partitioner::~Partitioner()
{
}

void Partitioner::compute_partition()
{
  std::string bal_package_str("Balancing package");
  std::string bal_package = paramlist_->get(bal_package_str, "none_specified");
  if (bal_package == "Zoltan" || bal_package == "zoltan") {
#ifdef HAVE_EPETRAEXT_ZOLTAN
    new_map_ = Isorropia_Zoltan::create_balanced_map(*input_graph_, *paramlist_);
#else
    throw Isorropia::Exception("Zoltan requested, but epetraext-zoltan not enabled.");
#endif
  }
  else {
    Epetra_Vector* weights_nnz =
      Isorropia::Epetra_Utils::create_row_weights_nnz(*input_graph_);
    new_map_ = Isorropia::Epetra_Utils::create_balanced_map(input_graph_->RowMap(),
							    *weights_nnz);
    delete weights_nnz;
  }
}

int Partitioner::newNumMyElements() const
{
  return( new_map_->NumMyElements() );
}

void Partitioner::myNewElements(int* elementList, int len) const
{
  int err = new_map_->MyGlobalElements(elementList);
  if (err != 0) {
    throw Isorropia::Exception("error in Epetra_Map::MyGlobalElements");
  }
}

#endif //HAVE_EPETRA

}//namespace Isorropia
