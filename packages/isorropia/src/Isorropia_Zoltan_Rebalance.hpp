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

#ifndef _Isorropia_Zoltan_Rebalance_hpp_
#define _Isorropia_Zoltan_Rebalance_hpp_

#include <Isorropia_configdefs.hpp>

#ifdef HAVE_EPETRAEXT_ZOLTAN

#include <Teuchos_RefCountPtr.hpp>
#include <Teuchos_ParameterList.hpp>

#ifdef HAVE_EPETRA
class Epetra_Map;
class Epetra_BlockMap;
class Epetra_Import;
class Epetra_Vector;
class Epetra_MultiVector;
class Epetra_CrsGraph;
class Epetra_CrsMatrix;
class Epetra_RowMatrix;
class Epetra_LinearProblem;
#endif

/** Isorropia_Zoltan is the namespace that contains isorropia's
  Zoltan-specific classes and functions.
*/
namespace Isorropia_Zoltan {

#ifdef HAVE_EPETRA

/** Create a balanced copy of an input Epetra_CrsGraph.

  This function represents a basic case, not accepting weights.
  By default, it uses a self-contained internal implementation
  rather than interfacing to a third-party library such as Zoltan.

  If the ParameterList object contains a string parameter with
  name == "Balancing package" and value == "Zoltan", then the Zoltan
  library is used to perform the balancing.

  Furthmore, if the ParameterList object contains a string parameter with
  name == "Partitioning algorithm" then possible values are "graph"
  or "hypergraph" (not yet supported). If "Zoltan" has not been specified
  as the "Balancing package", then "Partitioning algorithm" is ignored.

  The rebalancing is 1-D, row-wise, and attempts to make the number
  of nonzeros equal in each partition. I.e., it is equivalent to specifying
  a weighted rebalance where the weights are the number of nonzeros in
  each row.
*/
Teuchos::RefCountPtr<Epetra_CrsGraph>
create_balanced_copy(const Epetra_CrsGraph& input_graph,
		     Teuchos::ParameterList& paramlist);

/** Create a balanced copy of an input Epetra_CrsMatrix.

  This function represents a basic case, not accepting weights.
  By default, it uses a self-contained internal implementation
  rather than interfacing to a third-party library such as Zoltan.

  If the ParameterList object contains a string parameter with
  name == "Balancing package" and value == "Zoltan", then the Zoltan
  library is used to perform the balancing.

  Furthmore, if the ParameterList object contains a string parameter with
  name == "Partitioning algorithm" then possible values are "graph"
  or "hypergraph" (not yet supported). If "Zoltan" has not been specified
  as the "Balancing package", then "Partitioning algorithm" is ignored.

  The rebalancing is 1-D, row-wise, and attempts to make the number
  of nonzeros equal in each partition. I.e., it is equivalent to specifying
  a weighted rebalance where the weights are the number of nonzeros in
  each row.
*/
Teuchos::RefCountPtr<Epetra_CrsMatrix>
create_balanced_copy(const Epetra_CrsMatrix& input_matrix,
		     Teuchos::ParameterList& paramlist);

/** Create a map which describes the layout of rows resulting from the
   repartitioning of the input graph object.
*/
Teuchos::RefCountPtr<Epetra_Map>
create_balanced_map(const Epetra_CrsGraph& input_graph,
		    Teuchos::ParameterList& paramlist);

#endif //HAVE_EPETRA

}//namespace Isorropia_Zoltan

//the following endif closes the '#ifdef HAVE_EPETRAEXT_ZOLTAN' block.
#endif

#endif

