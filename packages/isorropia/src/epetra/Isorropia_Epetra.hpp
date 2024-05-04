//@HEADER
//************************************************************************
//
//              Isorropia: Partitioning and Load Balancing Package
//                Copyright (2006) Sandia Corporation
//
//Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
//license for use of this work by or on behalf of the U.S. Government.
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
//************************************************************************
//@HEADER

#ifndef _Isorropia_Epetra_hpp_
#define _Isorropia_Epetra_hpp_

#include <Isorropia_ConfigDefs.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>

#ifdef HAVE_EPETRA
class Epetra_Comm;
class Epetra_Map;
class Epetra_BlockMap;
class Epetra_Import;
class Epetra_Vector;
class Epetra_MultiVector;
class Epetra_RowMatrix;
class Epetra_CrsGraph;
class Epetra_CrsMatrix;
class Epetra_LinearProblem;
#endif

namespace Isorropia {

namespace Epetra {

  class Partitioner;
  class CostDescriber;

#ifdef HAVE_EPETRA

/** createBalancedCopy() creates a copy with a more balanced map.
    The caller should free the copy after use.
*/
Epetra_MultiVector *
createBalancedCopy(const Epetra_MultiVector& input_vector);

/** createBalancedCopy() creates a copy with a more balanced map.
    The caller should free the copy after use.
*/
Epetra_MultiVector *
createBalancedCopy(const Epetra_MultiVector& input_vector,
                     const Teuchos::ParameterList& paramlist);

/** createBalancedCopy() creates a copy with a more balanced map.
    The caller should free the copy after use.
*/
Epetra_CrsGraph *
createBalancedCopy(const Epetra_CrsGraph& input_graph);

/** createBalancedCopy() creates a copy with a more balanced map.
    The caller should free the copy after use.
*/
Epetra_CrsGraph *
createBalancedCopy(const Epetra_CrsGraph& input_graph,
                     const Teuchos::ParameterList& paramlist);

/** createBalancedCopy() creates a copy with a more balanced map.
    The caller should free the copy after use.
*/
Epetra_CrsMatrix *
createBalancedCopy(const Epetra_CrsMatrix& input_matrix);

/** createBalancedCopy() creates a copy with a more balanced map.
    The caller should free the copy after use.
*/
Epetra_CrsMatrix *
createBalancedCopy(const Epetra_CrsMatrix& input_matrix,
                     const Teuchos::ParameterList& paramlist);

/** createBalancedCopy() creates a copy with a more balanced map.
    The caller should free the copy after use.
*/

Epetra_LinearProblem *
createBalancedCopy(const Epetra_LinearProblem & input_problem);

/** createBalancedCopy() creates a copy with a more balanced map.
    The caller should free the copy after use.
*/

Epetra_LinearProblem *
createBalancedCopy(const Epetra_LinearProblem & input_problem,
                     const Teuchos::ParameterList& paramlist);

/** create_target_map() is an internal function used by Isorropia.

   Given a Partitioner object, create a target map representing the
   new partitioning.
   This function calls partitioner.compute_partitioning() if it has not
   already been called.
*/
// Teuchos::RCP<Epetra_Map>
// create_target_map(const Epetra_Comm& comm, Partitioner& partitioner);

/** create_balanced_copy(), which is part of the Isorropia API, is used to
    create and return a balanced copy of an input Epetra_CrsMatrix.

  Isorropia will use Zoltan to perform partitioning if it has been
  configured with Zoltan support.  If Zoltan is not available, it will
  perform its own simple row partitioning.

  Zoltan will perform hypergraph partitioning using its own PHG method.  
  It will assign unit weight to each row (vertex), and unit weight to each 
  column (hyperedge).  It will attempt to balance row weights while 
  minimizing cuts in the hyperedges.

  The non-Zoltan rebalancing is 1-D, row-wise, and attempts to make the number
  of nonzeros equal in each partition. I.e., it is equivalent to specifying
  a weighted rebalance where the weights are the number of nonzeros in
  each row.

*/
__deprecated Teuchos::RCP<Epetra_CrsMatrix>
  create_balanced_copy(const Epetra_CrsMatrix& input_matrix);

/** create_balanced_copy(), which is part of the Isorropia API, is used to
    create and return a balanced copy of an input Epetra_CrsMatrix.

  Isorropia will use Zoltan to perform partitioning if it has been
  configured with Zoltan support.  If Zoltan is not available, it will
  perform its own simple row partitioning.

  The balancing algorithm (Zoltan or non-Zoltan) will use the row
  weights provided for the matrix rows.

  Zoltan will perform hypergraph partitioning using its own PHG method.  
  It will assign a unit weight to each column (hyperedge).  It will attempt 
  to balance row weights while minimizing cuts in the hyperedges.

  The non-Zoltan rebalancing is a quick 1-D, row-wise balancing.
*/
__deprecated Teuchos::RCP<Epetra_CrsMatrix>
  create_balanced_copy(const Epetra_CrsMatrix& input_matrix,
                       const Epetra_Vector &row_weights);

/** create_balanced_copy(), which is part of the Isorropia API, is used to
    create and return a balanced copy of an input Epetra_CrsMatrix.

  Isorropia can use Zoltan to perform partitioning if it has been
  configured with Zoltan support.  If Zoltan is not available, it will
  perform its own simple row partitioning.

  In addition, if the parameter "PARTITIONING_METHOD" is set to "SIMPLE_LINEAR",
  Isorropia will perform its own simple row partitioning, even if Zoltan is
  available.  This method balances the rows across processes after assigning 
  each row a weight equal to the number of non zeros it has.

  If Zoltan is called to do the partitioning, any parameters set in a
  "Zoltan" sublist of the paramlist are provided to Zoltan.

  Refer to the Zoltan users guide for specific parameters that Zoltan
  recognizes. A couple of important ones are "LB_METHOD" (valid values
  include "GRAPH", "HYPERGRAPH"), "DEBUG_LEVEL" (valid values are
  0 to 10, default is 1), etc.

  If no Zoltan parameters are provided, Zoltan will perform hypergraph
  partitioning using its own PHG method.  It will assign unit weight to
  each row (vertex), and unit weight to each column (hyperedge).  It will
  attempt to minimize cuts in the hyperedges.
*/
__deprecated Teuchos::RCP<Epetra_CrsMatrix>
  create_balanced_copy(const Epetra_CrsMatrix& input_matrix,
                     const Teuchos::ParameterList& paramlist);



/** create_balanced_copy(), which is part of the Isorropia API, is used to
    create and return a balanced copy of an input Epetra_CrsMatrix.

  Isorropia can use Zoltan to perform partitioning if it has been
  configured with Zoltan support.  If Zoltan is not available, it will
  perform its own simple row partitioning.

  In addition, if the parameter "PARTITIONING_METHOD" is set to "SIMPLE_LINEAR",
  Isorropia will perform its own simple row partitioning, even if Zoltan is
  available.  This method balances the rows across processes after assigning 
  each row a weight equal to the number of non zeros it has.  It ignores the
  costs object.

  If Zoltan is called to do the partitioning, any parameters set in a
  "Zoltan" sublist of the paramlist are provided to Zoltan.  The weights
  provided in the CostDescriber object will be supplied to Zoltan.

  Refer to the Zoltan users guide for specific parameters that Zoltan
  recognizes. A couple of important ones are "LB_METHOD" (valid values
  include "GRAPH", "HYPERGRAPH"), "DEBUG_LEVEL" (valid values are
  0 to 10, default is 1), etc.

  If no Zoltan parameters are provided, Zoltan will perform hypergraph
  partitioning using its own PHG method.  

  If an empty CostDescriber is supplied, we will assign unit weight to each 
  row (vertex), and unit weight to each column (hyperedge), in the case of 
  hypergraph partitioning.  We will assign unit weight to each row (vertex)
  and each non zero (edge) in the case of graph partitioning.
*/
__deprecated Teuchos::RCP<Epetra_CrsMatrix>
  create_balanced_copy(const Epetra_CrsMatrix& input_matrix,
                     CostDescriber &costs,
                     const Teuchos::ParameterList& paramlist);

/** create_balanced_copy(), which is part of the Isorropia API, is used to
    create and return a balanced copy of an input Epetra_CrsMatrix.


  Isorropia will use Zoltan to perform partitioning if it has been
  configured with Zoltan support.  If Zoltan is not available, it will
  perform its own simple row partitioning.

  Zoltan will perform hypergraph partitioning using its own PHG method.  
  It will assign unit weight to each row (vertex), and unit weight to each 
  column (hyperedge).  It will attempt to minimize cuts in the hyperedges.

  The non-Zoltan rebalancing is 1-D, row-wise, and attempts to make the number
  of nonzeros equal in each partition. I.e., it is equivalent to specifying
  a weighted rebalance where the weights are the number of nonzeros in
  each row.

*/
__deprecated Teuchos::RCP<Epetra_RowMatrix>
  create_balanced_copy(const Epetra_RowMatrix& input_matrix);

/** create_balanced_copy(), which is part of the Isorropia API, is used to
    create and return a balanced copy of an input Epetra_CrsMatrix.

  Isorropia will use Zoltan to perform partitioning if it has been
  configured with Zoltan support.  If Zoltan is not available, it will
  perform its own simple row partitioning.

  The balancing algorithm (Zoltan or non-Zoltan) will use the row
  weights provided for the matrix rows.

  Zoltan will perform hypergraph partitioning using its own PHG method.  
  It will assign a unit weight to each column (hyperedge).  It will attempt 
  to balance row weights while minimizing cuts in the hyperedges.

  The non-Zoltan rebalancing is a quick 1-D, row-wise balancing.
*/
__deprecated Teuchos::RCP<Epetra_RowMatrix>
  create_balanced_copy(const Epetra_RowMatrix& input_matrix,
                       const Epetra_Vector &row_weights);

/** create_balanced_copy(), which is part of the Isorropia API, is used to
    create and return a balanced copy of an input Epetra_CrsMatrix.


  Isorropia can use Zoltan to perform partitioning if it has been
  configured with Zoltan support.  If Zoltan is not available, it will
  perform its own simple row partitioning.

  In addition, if the parameter "PARTITIONING_METHOD" is set to "SIMPLE_LINEAR",
  Isorropia will perform its own simple row partitioning, even if Zoltan is
  available.  This method balances the rows across processes after assigning 
  each row a weight equal to the number of non zeros it has.

  If Zoltan is called to do the partitioning, any parameters set in a
  "Zoltan" sublist of the paramlist are provided to Zoltan.

  Refer to the Zoltan users guide for specific parameters that Zoltan
  recognizes. A couple of important ones are "LB_METHOD" (valid values
  include "GRAPH", "HYPERGRAPH"), "DEBUG_LEVEL" (valid values are
  0 to 10, default is 1), etc.

  If no Zoltan parameters are provided, Zoltan will perform hypergraph
  partitioning using its own PHG method.  It will assign unit weight to
  each row (vertex), and unit weight to each column (hyperedge).  It will
  attempt to minimize cuts in the hyperedges.
*/
__deprecated Teuchos::RCP<Epetra_RowMatrix>
  create_balanced_copy(const Epetra_RowMatrix& input_matrix,
                     const Teuchos::ParameterList& paramlist);

/** create_balanced_copy(), which is part of the Isorropia API, is used to
    create and return a balanced copy of an input Epetra_CrsMatrix.


  Isorropia can use Zoltan to perform partitioning if it has been
  configured with Zoltan support.  If Zoltan is not available, it will
  perform its own simple row partitioning.

  In addition, if the parameter "PARTITIONING_METHOD" is set to "SIMPLE_LINEAR",
  Isorropia will perform its own simple row partitioning, even if Zoltan is
  available.  This method balances the rows across processes after assigning 
  each row a weight equal to the number of non zeros it has.  The costs
  object is ignored in this case.

  If Zoltan is called to do the partitioning, any parameters set in a
  "Zoltan" sublist of the paramlist are provided to Zoltan.  The weights
  provided in the CostDescriber object will be supplied to Zoltan.

  Refer to the Zoltan users guide for specific parameters that Zoltan
  recognizes. A couple of important ones are "LB_METHOD" (valid values
  include "GRAPH", "HYPERGRAPH"), "DEBUG_LEVEL" (valid values are
  0 to 10, default is 1), etc.

  If no Zoltan parameters are provided, Zoltan will perform hypergraph
  partitioning using its own PHG method.  

  If an empty CostDescriber is supplied, we will assign unit weight to each 
  row (vertex), and unit weight to each column (hyperedge), in the case of 
  hypergraph partitioning.  We will assign unit weight to each row (vertex)
  and each non zero (edge) in the case of graph partitioning.
*/
__deprecated Teuchos::RCP<Epetra_RowMatrix>
  create_balanced_copy(const Epetra_RowMatrix& input_matrix,
                     CostDescriber &costs,
                     const Teuchos::ParameterList& paramlist);

/** create_balanced_copy(), which is part of the Isorropia API, is used to
    create and return a balanced copy of an input Epetra_CrsMatrix.


  Isorropia will use Zoltan to perform partitioning if it has been
  configured with Zoltan support.  If Zoltan is not available, it will
  perform its own simple row partitioning.

  Zoltan will perform hypergraph partitioning using its own PHG method.  
  It will assign unit weight to each row (vertex), and unit weight to each 
  column (hyperedge).  It will attempt to minimize cuts in the hyperedges.

  The non-Zoltan rebalancing is 1-D, row-wise, and attempts to make the number
  of nonzeros equal in each partition. I.e., it is equivalent to specifying
  a weighted rebalance where the weights are the number of nonzeros in
  each row.

*/
__deprecated Teuchos::RCP<Epetra_CrsGraph>
  create_balanced_copy(const Epetra_CrsGraph& input_graph);

/** create_balanced_copy(), which is part of the Isorropia API, is used to
    create and return a balanced copy of an input Epetra_CrsMatrix.


  Isorropia will use Zoltan to perform partitioning if it has been
  configured with Zoltan support.  If Zoltan is not available, it will
  perform its own simple row partitioning.

  The balancing algorithm (Zoltan or non-Zoltan) will use the row
  weights provided for the matrix rows.

  Zoltan will perform hypergraph partitioning using its own PHG method.  
  It will assign a unit weight to each column (hyperedge).  It will attempt 
  to balance row weights while minimizing cuts in the hyperedges.

  The non-Zoltan rebalancing is a quick 1-D, row-wise balancing.
*/
__deprecated Teuchos::RCP<Epetra_CrsGraph>
  create_balanced_copy(const Epetra_CrsGraph& input_matrix,
                       const Epetra_Vector &row_weights);

/** create_balanced_copy(), which is part of the Isorropia API, is used to
    create and return a balanced copy of an input Epetra_CrsMatrix.

  Isorropia can use Zoltan to perform partitioning if it has been
  configured with Zoltan support.  If Zoltan is not available, it will
  perform its own simple row partitioning.

  In addition, if the parameter "PARTITIONING_METHOD" is set to "SIMPLE_LINEAR",
  Isorropia will perform its own simple row partitioning, even if Zoltan is
  available.  This method balances the rows across processes after assigning 
  each row a weight equal to the number of non zeros it has.

  If Zoltan is called to do the partitioning, any parameters set in a
  "Zoltan" sublist of the paramlist are provided to Zoltan.

  Refer to the Zoltan users guide for specific parameters that Zoltan
  recognizes. A couple of important ones are "LB_METHOD" (valid values
  include "GRAPH", "HYPERGRAPH"), "DEBUG_LEVEL" (valid values are
  0 to 10, default is 1), etc.

  If no Zoltan parameters are provided, Zoltan will perform hypergraph
  partitioning using its own PHG method.  It will assign unit weight to
  each row (vertex), and unit weight to each column (hyperedge).  It will
  attempt to minimize cuts in the hyperedges.
*/
__deprecated Teuchos::RCP<Epetra_CrsGraph>
  create_balanced_copy(const Epetra_CrsGraph& input_graph,
                     const Teuchos::ParameterList& paramlist);

/** create_balanced_copy(), which is part of the Isorropia API, is used to
    create and return a balanced copy of an input Epetra_CrsMatrix.


  Isorropia can use Zoltan to perform partitioning if it has been
  configured with Zoltan support.  If Zoltan is not available, it will
  perform its own simple row partitioning.

  In addition, if the parameter "PARTITIONING_METHOD" is set to "SIMPLE_LINEAR",
  Isorropia will perform its own simple row partitioning, even if Zoltan is
  available.  This method balances the rows across processes after assigning 
  each row a weight equal to the number of non zeros it has.

  If Zoltan is called to do the partitioning, any parameters set in a
  "Zoltan" sublist of the paramlist are provided to Zoltan.  The weights
  provided in the CostDescriber object will be supplied to Zoltan.

  Refer to the Zoltan users guide for specific parameters that Zoltan
  recognizes. A couple of important ones are "LB_METHOD" (valid values
  include "GRAPH", "HYPERGRAPH"), "DEBUG_LEVEL" (valid values are
  0 to 10, default is 1), etc.

  If no Zoltan parameters are provided, Zoltan will perform hypergraph
  partitioning using its own PHG method.  

  If an empty CostDescriber is supplied, we will assign unit weight to each 
  row (vertex), and unit weight to each column (hyperedge), in the case of 
  hypergraph partitioning.  We will assign unit weight to each row (vertex)
  and each non zero (edge) in the case of graph partitioning.
*/
__deprecated Teuchos::RCP<Epetra_CrsGraph>
  create_balanced_copy(const Epetra_CrsGraph& input_graph,
                     CostDescriber &costs,
                     const Teuchos::ParameterList& paramlist);

/** create_balanced_copy(), which is part of the Isorropia API, is used to
    create and return a balanced copy of an input Epetra_CrsMatrix.


  Isorropia will use Zoltan to perform partitioning if it has been
  configured with Zoltan support.  If Zoltan is not available, it will
  perform its own simple row partitioning.

  Zoltan will perform hypergraph partitioning using its own PHG method.  
  It will assign unit weight to each row (vertex), and unit weight to each 
  column (hyperedge).  It will attempt to minimize cuts in the hyperedges.

  The non-Zoltan rebalancing is 1-D, row-wise, and attempts to make the number
  of nonzeros equal in each partition. I.e., it is equivalent to specifying
  a weighted rebalance where the weights are the number of nonzeros in
  each row.

*/

__deprecated Teuchos::RCP<Epetra_LinearProblem>
  create_balanced_copy(const Epetra_LinearProblem & input_problem);

/** create_balanced_copy(), which is part of the Isorropia API, is used to
    create and return a balanced copy of an input Epetra_CrsMatrix.

  Isorropia will use Zoltan to perform partitioning if it has been
  configured with Zoltan support.  If Zoltan is not available, it will
  perform its own simple row partitioning.

  The balancing algorithm (Zoltan or non-Zoltan) will use the row
  weights provided for the matrix rows.

  Zoltan will perform hypergraph partitioning using its own PHG method.  
  It will assign a unit weight to each column (hyperedge).  It will attempt 
  to balance row weights while minimizing cuts in the hyperedges.

  The non-Zoltan rebalancing is a quick 1-D, row-wise balancing.
*/
__deprecated Teuchos::RCP<Epetra_LinearProblem>
  create_balanced_copy(const Epetra_LinearProblem& input_matrix,
                       const Epetra_Vector &row_weights);

/** create_balanced_copy(), which is part of the Isorropia API, is used to
    create and return a balanced copy of an input Epetra_CrsMatrix.

  Isorropia can use Zoltan to perform partitioning if it has been
  configured with Zoltan support.  If Zoltan is not available, it will
  perform its own simple row partitioning.

  In addition, if the parameter "PARTITIONING_METHOD" is set to "SIMPLE_LINEAR",
  Isorropia will perform its own simple row partitioning, even if Zoltan is
  available.  This method balances the rows across processes after assigning 
  each row a weight equal to the number of non zeros it has.

  If Zoltan is called to do the partitioning, any parameters set in a
  "Zoltan" sublist of the paramlist are provided to Zoltan.

  Refer to the Zoltan users guide for specific parameters that Zoltan
  recognizes. A couple of important ones are "LB_METHOD" (valid values
  include "GRAPH", "HYPERGRAPH"), "DEBUG_LEVEL" (valid values are
  0 to 10, default is 1), etc.

  If no Zoltan parameters are provided, Zoltan will perform hypergraph
  partitioning using its own PHG method.  It will assign unit weight to
  each row (vertex), and unit weight to each column (hyperedge).  It will
  attempt to minimize cuts in the hyperedges.
*/

__deprecated Teuchos::RCP<Epetra_LinearProblem>
  create_balanced_copy(const Epetra_LinearProblem& input_problem,
                     const Teuchos::ParameterList& paramlist);

/** create_balanced_copy(), which is part of the Isorropia API, is used to
    create and return a balanced copy of an input Epetra_CrsMatrix.

  Isorropia can use Zoltan to perform partitioning if it has been
  configured with Zoltan support.  If Zoltan is not available, it will
  perform its own simple row partitioning.

  In addition, if the parameter "PARTITIONING_METHOD" is set to "SIMPLE_LINEAR",
  Isorropia will perform its own simple row partitioning, even if Zoltan is
  available.  This method balances the rows across processes after assigning 
  each row a weight equal to the number of non zeros it has.  This costs
  object is ignored in this case.

  If Zoltan is called to do the partitioning, any parameters set in a
  "Zoltan" sublist of the paramlist are provided to Zoltan.  The weights
  provided in the CostDescriber object will be supplied to Zoltan.

  Refer to the Zoltan users guide for specific parameters that Zoltan
  recognizes. A couple of important ones are "LB_METHOD" (valid values
  include "GRAPH", "HYPERGRAPH"), "DEBUG_LEVEL" (valid values are
  0 to 10, default is 1), etc.

  If no Zoltan parameters are provided, Zoltan will perform hypergraph
  partitioning using its own PHG method.  

  If an empty CostDescriber is supplied, we will assign unit weight to each 
  row (vertex), and unit weight to each column (hyperedge), in the case of 
  hypergraph partitioning.  We will assign unit weight to each row (vertex)
  and each non zero (edge) in the case of graph partitioning.
*/
__deprecated Teuchos::RCP<Epetra_LinearProblem>
  create_balanced_copy(const Epetra_LinearProblem& input_problem,
                     CostDescriber &costs,
                     const Teuchos::ParameterList& paramlist);

__deprecated Teuchos::RCP<Epetra_MultiVector>
create_balanced_copy(const Epetra_MultiVector &coords,
                   const Teuchos::ParameterList& paramlist);

__deprecated Teuchos::RCP<Epetra_MultiVector>
create_balanced_copy(const Epetra_MultiVector &coords,
                     const Epetra_MultiVector &weights,
                   const Teuchos::ParameterList& paramlist);


__deprecated Teuchos::RCP<Epetra_MultiVector>
create_balanced_copy(const Epetra_MultiVector &coords);

__deprecated Teuchos::RCP<Epetra_MultiVector>
create_balanced_copy(const Epetra_MultiVector &coords,
                     const Epetra_MultiVector &weights);

/** redistribute_rows() is an internal Isorropia function, not part
    of the API.

  Return a new Epetra_CrsMatrix object constructed with target_rowmap,
  and with the contents of input_matrix imported into it.

  The caller is responsible for deleting the returned object.

  param input_matrix Source/input object.

  param target_rowmap Target rowmap, required to be compatible with
     input_matrix.RowMap() in terms of number-of-elements, etc.

  param importer Optional argument. If importer is supplied, it will be
     used to perform the import operation. Otherwise, a temporary importer
     will be created and used.
*/
Teuchos::RCP<Epetra_CrsMatrix>
  redistribute_rows(const Epetra_CrsMatrix& input_matrix,
                    const Epetra_Map& target_rowmap,
                    Epetra_Import* importer=0);

/** redistribute_rows() is an internal Isorropia function, not part
    of the API.

    Return a new Epetra_CrsMatrix object constructed with target_rowmap,
  and with the contents of input_matrix imported into it.

  The caller is responsible for deleting the returned object.

  param input_matrix Source/input object.

  param target_rowmap Target rowmap, required to be compatible with
     input_matrix.RowMatrixRowMap() in terms of number-of-elements, etc.

  param importer Optional argument. If importer is supplied, it will be
     used to perform the import operation. Otherwise, a temporary importer
     will be created and used.
*/
Teuchos::RCP<Epetra_CrsMatrix>
  redistribute_rows(const Epetra_RowMatrix& input_matrix,
                    const Epetra_Map& target_rowmap,
                    Epetra_Import* importer=0);

/** Return a new Epetra_CrsGraph object constructed with target_rowmap,
  and with the contents of input_graph imported into it.

  param input_graph Source/input object.

  param target_rowmap Target rowmap, required to be compatible with
     input_graph.RowMap() in terms of number-of-elements, etc.

  param importer Optional argument. If importer is supplied, it will be
     used to perform the import operation. Otherwise, a temporary importer
     will be created and used.
*/
Teuchos::RCP<Epetra_CrsGraph>
  redistribute_rows(const Epetra_CrsGraph& input_graph,
                    const Epetra_Map& target_rowmap,
                    Epetra_Import* importer=0);

/** Return a new Epetra_MultiVector object constructed with target_map,
  and with the contents of 'input' imported into it.

  param input Source/input object.

  param target_map Target map, required to be compatible with
     input.Map() in terms of number-of-elements, etc.

  param importer Optional argument. If importer is supplied, it will be
     used to perform the import operation. Otherwise, a temporary importer
     will be created and used.
*/
Teuchos::RCP<Epetra_MultiVector>
  redistribute(const Epetra_MultiVector& input,
               const Epetra_BlockMap& target_map,
               Epetra_Import* importer=0);

/** Return a new Epetra_Vector object constructed with target_map,
  and with the contents of 'input' imported into it.

  param input Source/input object.

  param target_map Target map, required to be compatible with
     input.RowMap() in terms of number-of-elements, etc.

  param importer Optional argument. If importer is supplied, it will be
     used to perform the import operation. Otherwise, a temporary importer
     will be created and used.
*/
Teuchos::RCP<Epetra_Vector>
  redistribute(const Epetra_Vector& input,
               const Epetra_Map& target_map,
               Epetra_Import* importer=0);

#ifndef DOXYGEN_SHOULD_SKIP_THIS

/** Return a vector containing weights that are equal to the number of
  nonzeros per row in the input_matrix. The returned vector will have
  the same size and distribution as input_matrix's row-map.
*/
Epetra_MultiVector* create_row_weights_nnz(const Epetra_RowMatrix& input_matrix);

/** Return a vector containing weights that are equal to the number of
  nonzeros per row in the input_graph. The returned vector will have
  the same size and distribution as input_graph's row-map.
*/
Epetra_MultiVector* create_row_weights_nnz(const Epetra_CrsGraph& input_graph);

Epetra_MultiVector* create_unit_weights(const Epetra_MultiVector& input_coords);


/** Calculate a new partitioning, and fill output containers with new
    elements for the local partition, as well as export and import lists.
    This is a simple linear partitioning that does not use Zoltan.

    \param[in] input_map Input map describing the existing or 'old' partitioning.

    \param[in] weights Input vector giving a weight for each element in input_map.
    weights.Map() is required to be the same size and layout as input_map.

    \param[out] newPartitions contains the new partition for each element,
                    in input_map local ID order.  Partition numbers go from
                    0 to numProcs - 1
         
    \param[out] exports the number of exports, that is, the number of
                  elements in newPartitions that are not equal to my
                  process rank

    \param[out] imports the list of global IDs of the elements I will
                   import under the new partitioning

    \return Error-code, 0 if successful. This probably should be a void
    function, since a serious error will result in an exception-throw
    rather than an integer-code error-return.
*/
int
repartition(const Epetra_BlockMap& input_map,
	    const Epetra_MultiVector& weights,
	    std::vector<int>& myNewElements,
	    int& exportsSize,
	    std::vector<int>& imports);

/**  gather_all_proc_global_offsets() is an internal Isorropia function, not
     part of the API.

    Given an Epetra_BlockMap object, fill a vector of length numprocs+1
  with each processor's starting offset into the Epetra_BlockMap's global
  set of elements (the last position will contain num-global-elements).
  Gather the vector of offsets onto all processors.
*/
void gather_all_proc_global_offsets(const Epetra_BlockMap& blkmap,
                                    std::vector<int>& all_proc_offsets);


/**  compute_imbalance() is an internal Isorropia function, not
     part of the API.

     This function is used by Isorropia to compute the global imbalance
     of an initial partitioning and a new partitioning, to ensure the
     new computed partitioning is better.
*/

double compute_imbalance(int nprocs, std::vector<int> &offsets, 
                         double *wgts, double target);

#endif //DOXYGEN_SHOULD_SKIP_THIS
#endif //HAVE_EPETRA

}//namespace Epetra
}//namespace Isorropia

#endif


#if defined(Isorropia_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The Isorropia package is deprecated"
#endif
#endif

