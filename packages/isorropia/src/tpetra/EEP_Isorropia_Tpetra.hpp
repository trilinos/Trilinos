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

#ifndef _Isorropia_Tpetra_hpp_
#define _Isorropia_Tpetra_hpp_

//#include <Isorropia_ConfigDefs.hpp> // AquiToDo

#if 1 // def HAVE_TPETRA // AquiToDo
//#include <Isorropia_Exception.hpp>
//#include <Isorropia_Utils.hpp>
//#include <Isorropia_EpetraCostDescriber.hpp>
#include <EEP_Isorropia_TpetraPartitioner.hpp>
#include <EEP_Isorropia_TpetraRedistributor.hpp>
#endif

//#include <Tpetra_Map_decl.hpp>
//#include <Tpetra_Import_decl.hpp>
//#include <Tpetra_Vector_decl.hpp>
//#include <Tpetra_Comm_decl.hpp>
//#include <Tpetra_IntVector_decl.hpp>
#include <Tpetra_CrsGraph_decl.hpp>
//#include <Tpetra_CrsMatrix_decl.hpp>
//#include <Tpetra_RowMatrix_decl.hpp>
//#include <Tpetra_LinearProblem_decl.hpp>

#ifdef HAVE_MPI
#include <Teuchos_DefaultComm.hpp>
#include <mpi.h>
#endif

#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>

namespace Isorropia {

namespace Tpetra {

  //class Partitioner;
  //class CostDescriber;

#if 1 // def HAVE_TPETRA // AquiToDo

/** createBalancedCopy() creates a copy with a more balanced map.
    The caller should free the copy after use.
*/
template <class LocalOrdinal,
          class GlobalOrdinal,
          class Node>
::Tpetra::CrsGraph<LocalOrdinal, GlobalOrdinal, Node> *
createBalancedCopy(const ::Tpetra::CrsGraph<LocalOrdinal, GlobalOrdinal, Node>& input_graph, // EEP__
                   const Teuchos::ParameterList& paramlist)
{
  std::cout << "EEP Entering isorropia/src/tpetra/EEP_Isorropia_Tpetra.hpp createBalancedCopy(4)..." << std::endl;
  Teuchos::RCP< const ::Tpetra::CrsGraph<LocalOrdinal, GlobalOrdinal, Node> > rcp_input_graph = Teuchos::rcp(&(input_graph), false);

  std::cout << "EEP In isorropia/src/tpetra/EEP_Isorropia_Tpetra.hpp createBalancedCopy(4), pos 001" << std::endl; // Aqui

  Teuchos::RCP< Partitioner<LocalOrdinal, GlobalOrdinal, Node> > partitioner = Teuchos::rcp(new Partitioner<LocalOrdinal, GlobalOrdinal, Node>(rcp_input_graph, paramlist));

  std::cout << "EEP In isorropia/src/tpetra/EEP_Isorropia_Tpetra.hpp createBalancedCopy(4), pos 002" << std::endl; 

  Redistributor rd(partitioner);

  std::cout << "EEP In isorropia/src/tpetra/EEP_Isorropia_Tpetra.hpp createBalancedCopy(4), pos 003" << std::endl; 

  Teuchos::RCP<::Tpetra::CrsGraph<LocalOrdinal, GlobalOrdinal, Node>> balanced_graph = rd.redistribute(input_graph);

  std::cout << "EEP In isorropia/src/tpetra/EEP_Isorropia_Tpetra.hpp createBalancedCopy(4), pos 004" << std::endl; 
  
#if 1 // AquiToDo // EEP___
  balanced_graph.release();

  std::cout << "EEP Leaving isorropia/src/tpetra/EEP_Isorropia_Tpetra.hpp createBalancedCopy(4)" << std::endl; 
  return balanced_graph.get();
#else
  std::cout << "EEP Leaving isorropia/src/tpetra/EEP_Isorropia_Tpetra.hpp createBalancedCopy(4)" << std::endl; 
  return nullptr;
#endif // AquiToDo
}

#endif // HAVE_TPETRA

} // namespace Tpetra
} // namespace Isorropia

#endif

