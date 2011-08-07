//@HEADER
// ***********************************************************************
// 
//     EpetraExt: Epetra Extended - Linear Algebra Services Package
//                 Copyright (2001) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
//@HEADER

#ifndef EpetraExt_AMESOS_AMD_GLOBAL_CRSGRAPH_H
#define EpetraExt_AMESOS_AMD_GLOBAL_CRSGRAPH_H

#include <EpetraExt_Transform.h>
#include <Teuchos_RCP.hpp>
#include <vector>

class Epetra_CrsMatrix;
class Epetra_CrsGraph;
class Epetra_Map;
class Epetra_Import;

namespace EpetraExt {

///
/** Approximate Minimum Degree (AMD) ordering of Epetra_CrsGraph
 *
 * Uses Tim Davis' AMD algorithm to find a fill-reducing ordering from an
 * Epetra_CrsGraph.
 */

class AmesosAMDGlobal_CrsGraph : public SameTypeTransform<Epetra_CrsGraph> {

 public:

  ~AmesosAMDGlobal_CrsGraph();

  AmesosAMDGlobal_CrsGraph( bool verbose = false, bool debug = false )
  : verbose_(verbose),
    debug_(debug)
  {}

  NewTypeRef operator()( OriginalTypeRef orig );

  bool fwd();
  bool rvs();

  Teuchos::RCP<Epetra_Import> Importer() { return Importer_; }
  std::vector<int> Perm() { return perm_; }

 private:

  Teuchos::RCP<Epetra_Map> NewRowMap_;
  Teuchos::RCP<Epetra_Map> NewColMap_;
  
  Teuchos::RCP<Epetra_CrsGraph> NewGraph_;

  Teuchos::RCP<Epetra_Import> Importer_;

  std::vector<int> perm_;
  
  const bool verbose_, debug_;

};

} //namespace EpetraExt

#endif //EpetraExt_AMESOS_AMD_GLOBAL_CRSGRAPH_H
