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

#ifndef _Isorropia_Zoltan_Repartition_hpp_
#define _Isorropia_Zoltan_Repartition_hpp_

#include <Isorropia_ConfigDefs.hpp>

#ifdef HAVE_ISORROPIA_ZOLTAN

#include <Teuchos_RefCountPtr.hpp>
#include <Teuchos_ParameterList.hpp>

#include <vector>
#include <map>

#ifdef HAVE_MPI
#include <mpi.h>
#endif

#ifdef HAVE_EPETRA
class Epetra_CrsGraph;
class Epetra_RowMatrix;
#endif

#include <QueryObject.hpp>

// class Isorropia::Epetra::CostDescriber;

namespace Isorropia{

namespace Epetra {

/** The ZoltanLib namespace within the Epetra namespace contains the
    classes and functions that use the Zoltan library to partition an
    Epetra object.
*/

namespace ZoltanLib{

#ifdef HAVE_EPETRA

/** Partition an Epetra_CrsGraph using Zoltan.

    myNewElements lists the global IDs of the rows I will have
    under the new partitioning.

    exports is a map from my row global IDs to process receiving that
    row under the new partititioning.

    imports is a map from each row globalID that will be sent to me
    under the new partitioning to the process that currently owns it.

    Isorropia::Epetra::ZoltanLib::repartition() is called by the
    Isorropia::Epetra::Partitioner constructor when Zoltan is available.
*/
int
repartition(Teuchos::RefCountPtr<const Epetra_CrsGraph> input_graph,
	    Teuchos::RefCountPtr<const Isorropia::Epetra::CostDescriber> costs,
	    Teuchos::ParameterList& paramlist,
            std::vector<int>& myNewElements,
            std::map<int,int>& exports,
            std::map<int,int>& imports);


/** Partition an Epetra_RowMatrix using Zoltan.

    myNewElements lists the global IDs of the rows I will have
    under the new partitioning.

    exports is a map from my row global IDs to process receiving that
    row under the new partititioning.

    imports is a map from each row globalID that will be sent to me
    under the new partitioning to the process that currently owns it.

    Isorropia::Epetra::ZoltanLib::repartition() is called by the
    Isorropia::Epetra::Partitioner constructor when Zoltan is available.
*/
int
repartition(Teuchos::RefCountPtr<const Epetra_RowMatrix> input_matrix,
	    Teuchos::RefCountPtr<const Isorropia::Epetra::CostDescriber> costs,
	    Teuchos::ParameterList& paramlist,
            std::vector<int>& myNewElements,
            std::map<int,int>& exports,
            std::map<int,int>& imports);

#ifdef HAVE_MPI
/** load_balance() is called by Isorropia::Epetra::ZoltanLib::repartition().

    It sets up the Zoltan query functions and parameters and calls Zoltan
    to perform the partitioning.

    myNewElements lists the global IDs of the rows I will have
    under the new partitioning.

    exports is a map from my row global IDs to process receiving that
    row under the new partititioning.

    imports is a map from each row globalID that will be sent to me
    under the new partitioning to the process that currently owns it.
*/
    
int
load_balance(MPI_Comm &comm,
	     Teuchos::ParameterList& paramlist,
	     QueryObject & queryObject,
	     std::vector<int>& myNewElements,
	     std::map<int,int>& exports,
	     std::map<int,int>& imports);
#endif

#endif //HAVE_EPETRA

}//namespace ZoltanLib

}//namespace Epetra

}//namespace Isorropia

//the following endif closes the '#ifdef HAVE_ISORROPIA_ZOLTAN' block.
#endif

#endif


#if defined(Isorropia_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The Isorropia package is deprecated"
#endif
#endif

