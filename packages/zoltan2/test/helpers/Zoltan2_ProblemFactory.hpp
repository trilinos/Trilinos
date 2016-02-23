// @HEADER
//
// ***********************************************************************
//
//   Zoltan2: A package of combinatorial algorithms for scientific computing
//                  Copyright 2012 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Questions? Contact Karen Devine      (kddevin@sandia.gov)
//                    Erik Boman        (egboman@sandia.gov)
//                    Siva Rajamanickam (srajama@sandia.gov)
//
// ***********************************************************************
//
// @HEADER

/*! \file Zoltan2_ProblemFactory.hpp
    \brief Returns a pointer to a new problem class.
    Is not responsible for memory management!
*/
#ifndef ZOLTAN2_PROBLEM_FACTORY_HPP
#define ZOLTAN2_PROBLEM_FACTORY_HPP
#include <Zoltan2_TestHelpers.hpp>
#include <Zoltan2_Problem.hpp>
#include <Zoltan2_PartitioningProblem.hpp>
#include <Zoltan2_ColoringProblem.hpp>
#include <Zoltan2_OrderingProblem.hpp>
#include <Zoltan2_Typedefs.hpp>

using namespace Zoltan2_TestingFramework;
using namespace Zoltan2;

namespace Zoltan2_TestingFramework {
/// \brief ProblemFactory class contains 1 static factory method
  class ProblemFactory {
  public:

    /// \brif Zoltan2::Problem factory method
    ///
    /// @param kind A string equal to the type of the problem (paritioning, ordering, coloring)
    /// @param kind A string equal to the type of the adapter
    /// @param input Input adapter used to construct the problem
    /// @param params Zolta2 parameter list
    /// @param (MPI) MPI world communicator
    ///
    /// @return returns a pointer to new Zoltan2::Problem or a nullptr if kind was not known.
#ifdef HAVE_ZOLTAN2_MPI
    static Problem<basic_id_t> * newProblem(  const std::string &kind,
                                                const std::string &adapter_name, 
                                                base_adapter_t *input,
                                                ParameterList *params,
                                                MPI_Comm comm)
#else
    static Problem<basic_id_t> * newProblem(  const std::string & kind,
                                                const std::string &adapter_name, 
                                                base_adapter_t *input,
                                                ParameterList *params)
#endif
    {
      if(kind == "partitioning") {
#ifdef HAVE_ZOLTAN2_MPI
        if (adapter_name == "BasicIdentifier")
          return reinterpret_cast< Problem<basic_id_t> *>( new Zoltan2::PartitioningProblem<basic_vector_adapter>(reinterpret_cast<basic_vector_adapter *>(input), params, comm));
        else if (adapter_name == "XpetraMultiVector")
          return reinterpret_cast< Problem<basic_id_t> *>( new Zoltan2::PartitioningProblem<xpetra_mv_adapter>(reinterpret_cast<xpetra_mv_adapter *>(input), params, comm));
        else if (adapter_name == "XpetraCrsGraph")
          return reinterpret_cast< Problem<basic_id_t> *>( new Zoltan2::PartitioningProblem<xcrsGraph_adapter>(reinterpret_cast<xcrsGraph_adapter *>(input), params, comm));
        else if (adapter_name == "XpetraCrsMatrix")
          return reinterpret_cast< Problem<basic_id_t> *>( new Zoltan2::PartitioningProblem<xcrsMatrix_adapter>(reinterpret_cast<xcrsMatrix_adapter *>(input), params, comm));
        else if (adapter_name == "BasicVector")
          return reinterpret_cast< Problem<basic_id_t> *>( new Zoltan2::PartitioningProblem<basic_vector_adapter>(reinterpret_cast<basic_vector_adapter *>(input), params, comm));
        else if (adapter_name == "PamgenMesh")
          return reinterpret_cast< Problem<basic_id_t> *>( new Zoltan2::PartitioningProblem<pamgen_adapter_t>(reinterpret_cast<pamgen_adapter_t*>(input), params, comm));
#else
        if (adapter_name == "BasicIdentifier")
          return reinterpret_cast< Problem<basic_id_t> *>( new Zoltan2::PartitioningProblem<basic_vector_adapter>(reinterpret_cast<basic_vector_adapter *>(input), params));
        else if (adapter_name == "XpetraMultiVector")
          return reinterpret_cast< Problem<basic_id_t> *>( new Zoltan2::PartitioningProblem<xpetra_mv_adapter>(reinterpret_cast<xpetra_mv_adapter *>(input), params));
        else if (adapter_name == "XpetraCrsGraph")
          return reinterpret_cast< Problem<basic_id_t> *>( new Zoltan2::PartitioningProblem<xcrsGraph_adapter>(reinterpret_cast<xcrsGraph_adapter *>(input), params));
        else if (adapter_name == "XpetraCrsMatrix")
          return reinterpret_cast< Problem<basic_id_t> *>( new Zoltan2::PartitioningProblem<xcrsMatrix_adapter>(reinterpret_cast<xcrsMatrix_adapter *>(input), params));
        else if (adapter_name == "BasicVector")
          return reinterpret_cast< Problem<basic_id_t> *>( new Zoltan2::PartitioningProblem<basic_vector_adapter>(reinterpret_cast<basic_vector_adapter *>(input), params));
        else if (adapter_name == "PamgenMesh")
          return reinterpret_cast< Problem<basic_id_t> *>( new Zoltan2::PartitioningProblem<pamgen_adapter_t>(reinterpret_cast<pamgen_adapter_t*>(input), params));
#endif
      } else if (kind == "ordering") {
#ifdef HAVE_ZOLTAN2_MPI
        if (adapter_name == "BasicIdentifier")
          return reinterpret_cast< Problem<basic_id_t> *>( new Zoltan2::OrderingProblem<basic_vector_adapter>(reinterpret_cast<basic_vector_adapter *>(input), params, comm));
        else if (adapter_name == "XpetraMultiVector")
          return reinterpret_cast< Problem<basic_id_t> *>( new Zoltan2::OrderingProblem<xpetra_mv_adapter>(reinterpret_cast<xpetra_mv_adapter *>(input), params, comm));
        else if (adapter_name == "XpetraCrsGraph")
          return reinterpret_cast< Problem<basic_id_t> *>( new Zoltan2::OrderingProblem<xcrsGraph_adapter>(reinterpret_cast<xcrsGraph_adapter *>(input), params, comm));
        else if (adapter_name == "XpetraCrsMatrix")
          return reinterpret_cast< Problem<basic_id_t> *>( new Zoltan2::OrderingProblem<xcrsMatrix_adapter>(reinterpret_cast<xcrsMatrix_adapter *>(input), params, comm));
        else if (adapter_name == "BasicVector")
          return reinterpret_cast< Problem<basic_id_t> *>( new Zoltan2::OrderingProblem<basic_vector_adapter>(reinterpret_cast<basic_vector_adapter *>(input), params, comm));
        else if (adapter_name == "PamgenMesh")
          return reinterpret_cast< Problem<basic_id_t> *>( new Zoltan2::OrderingProblem<pamgen_adapter_t>(reinterpret_cast<pamgen_adapter_t*>(input), params, comm));
#else
        if (adapter_name == "BasicIdentifier")
          return reinterpret_cast< Problem<basic_id_t> *>( new Zoltan2::OrderingProblem<basic_vector_adapter>(reinterpret_cast<basic_vector_adapter *>(input), params));
        else if (adapter_name == "XpetraMultiVector")
          return reinterpret_cast< Problem<basic_id_t> *>( new Zoltan2::OrderingProblem<xpetra_mv_adapter>(reinterpret_cast<xpetra_mv_adapter *>(input), params));
        else if (adapter_name == "XpetraCrsGraph")
          return reinterpret_cast< Problem<basic_id_t> *>( new Zoltan2::OrderingProblem<xcrsGraph_adapter>(reinterpret_cast<xcrsGraph_adapter *>(input), params));
        else if (adapter_name == "XpetraCrsMatrix")
          return reinterpret_cast< Problem<basic_id_t> *>( new Zoltan2::OrderingProblem<xcrsMatrix_adapter>(reinterpret_cast<xcrsMatrix_adapter *>(input), params));
        else if (adapter_name == "BasicVector")
          return reinterpret_cast< Problem<basic_id_t> *>( new Zoltan2::OrderingProblem<basic_vector_adapter>(reinterpret_cast<basic_vector_adapter *>(input), params));
        else if (adapter_name == "PamgenMesh")
          return reinterpret_cast< Problem<basic_id_t> *>( new Zoltan2::OrderingProblem<pamgen_adapter_t>(reinterpret_cast<pamgen_adapter_t*>(input), params));
#endif
      } else if (kind == "coloring") {
#ifdef HAVE_ZOLTAN2_MPI
        if (adapter_name == "BasicIdentifier")
          return reinterpret_cast< Problem<basic_id_t> *>( new Zoltan2::ColoringProblem<basic_vector_adapter>(reinterpret_cast<basic_vector_adapter *>(input), params, comm));
        else if (adapter_name == "XpetraMultiVector")
          return reinterpret_cast< Problem<basic_id_t> *>( new Zoltan2::ColoringProblem<xpetra_mv_adapter>(reinterpret_cast<xpetra_mv_adapter *>(input), params, comm));
        else if (adapter_name == "XpetraCrsGraph")
          return reinterpret_cast< Problem<basic_id_t> *>( new Zoltan2::ColoringProblem<xcrsGraph_adapter>(reinterpret_cast<xcrsGraph_adapter *>(input), params, comm));
        else if (adapter_name == "XpetraCrsMatrix")
          return reinterpret_cast< Problem<basic_id_t> *>( new Zoltan2::ColoringProblem<xcrsMatrix_adapter>(reinterpret_cast<xcrsMatrix_adapter *>(input), params, comm));
        else if (adapter_name == "BasicVector")
          return reinterpret_cast< Problem<basic_id_t> *>( new Zoltan2::ColoringProblem<basic_vector_adapter>(reinterpret_cast<basic_vector_adapter *>(input), params, comm));
        else if (adapter_name == "PamgenMesh")
          return reinterpret_cast< Problem<basic_id_t> *>( new Zoltan2::ColoringProblem<pamgen_adapter_t>(reinterpret_cast<pamgen_adapter_t*>(input), params, comm));
#else
        if (adapter_name == "BasicIdentifier")
          return reinterpret_cast< Problem<basic_id_t> *>( new Zoltan2::ColoringProblem<basic_vector_adapter>(reinterpret_cast<basic_vector_adapter *>(input), params));
        else if (adapter_name == "XpetraMultiVector")
          return reinterpret_cast< Problem<basic_id_t> *>( new Zoltan2::ColoringProblem<xpetra_mv_adapter>(reinterpret_cast<xpetra_mv_adapter *>(input), params));
        else if (adapter_name == "XpetraCrsGraph")
          return reinterpret_cast< Problem<basic_id_t> *>( new Zoltan2::ColoringProblem<xcrsGraph_adapter>(reinterpret_cast<xcrsGraph_adapter *>(input), params));
        else if (adapter_name == "XpetraCrsMatrix")
          return reinterpret_cast< Problem<basic_id_t> *>( new Zoltan2::ColoringProblem<xcrsMatrix_adapter>(reinterpret_cast<xcrsMatrix_adapter *>(input), params));
        else if (adapter_name == "BasicVector")
          return reinterpret_cast< Problem<basic_id_t> *>( new Zoltan2::ColoringProblem<basic_vector_adapter>(reinterpret_cast<basic_vector_adapter *>(input), params));
        else if (adapter_name == "PamgenMesh")
          return reinterpret_cast< Problem<basic_id_t> *>( new Zoltan2::ColoringProblem<pamgen_adapter_t>(reinterpret_cast<pamgen_adapter_t*>(input), params));
#endif
      }
     return nullptr; // problem type not known 
    }
  };
}
#endif // ZOLTAN2_PROBLEM_FACTORY_HPP

