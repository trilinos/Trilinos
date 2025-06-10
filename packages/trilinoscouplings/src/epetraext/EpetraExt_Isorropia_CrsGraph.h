// @HEADER
// *****************************************************************************
//           Trilinos: An Object-Oriented Solver Framework
//
// Copyright 2001-2024 NTESS and the Trilinos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef EPETRAEXT_ISORROPIA_CRSGRAPH_H
#define EPETRAEXT_ISORROPIA_CRSGRAPH_H

#include <EpetraExt_Transform.h>

#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>

class Epetra_CrsGraph;
class Epetra_Map;

namespace EpetraExt {

///
/** Generates an Epetra_CrsGraph based on the repartitioning algorithms of Isorropia
 */
class Isorropia_CrsGraph : public StructuralSameTypeTransform<Epetra_CrsGraph> {

  // Parameter list to specify partitioning algorithm for Isorropia.
  Teuchos::ParameterList partitionList_;

  // New graph for the repartitioned matrix.
  Teuchos::RCP<Epetra_CrsGraph> NewGraph_;

 public:

  ///
  /* Destructor
   */
  ~Isorropia_CrsGraph();

  ///
  /* Constructor
   * input param part_method - type of Isorropia partitioning to use
   */
  Isorropia_CrsGraph( Teuchos::ParameterList& partitionList )
  : partitionList_(partitionList)
  {}

  ///
  /* Generates the Isorropia partitioned Epetra_CrsGraph from the input object.
   */
  NewTypeRef operator()( OriginalTypeRef orig );

};

} //namespace EpetraExt

#endif //EPETRAEXT_ISORROPIA_CRSGRAPH_H
