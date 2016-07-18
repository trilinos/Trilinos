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

/*! \file Zoltan2_EvaluatePartitionFactory.hpp
    \brief Returns a pointer to a new evaluate partition class.
    Is not responsible for memory management!
*/
#ifndef ZOLTAN2_EVALUATE_PARTITION_FACTORY_HPP
#define ZOLTAN2_EVALUATE_PARTITION_FACTORY_HPP
namespace Zoltan2_TestingFramework {
///brief EvaluatePartitionFActory class contains 1 static factory method
  class EvaluatePartitionFactory {
  public:

    /// \brief Zoltan2::EvaluatePartition factory method
    ///
    /// @param problem partitioning problem
    /// @param adapter_name a string equal to the type of the adapter
    /// @param input Input adapter used to construct the partitioning problem
    ///
    /// @return returns a pointer to new Zoltan2::EvaluatePartition or a nullptr if adapter type was not known.
    static EvaluatePartition<basic_id_t>*newEvaluatePartition
    (partitioning_problem_t *problem, const std::string &adapter_name,
     base_adapter_t *input, ParameterList *params)
    {
      RCP<const Teuchos::Comm<int> > CommT = problem->getComm();

      if (adapter_name == "BasicIdentifier") {
	return reinterpret_cast<Zoltan2::EvaluatePartition<basic_id_t> *>
	  (new Zoltan2::EvaluatePartition<basic_id_t>
	   (reinterpret_cast<const basic_id_t *>(input),
	    params, CommT, reinterpret_cast
	    <const Zoltan2::PartitioningSolution<basic_id_t> *>
	    (&problem->getSolution())));
      } else if (adapter_name == "XpetraMultiVector") {
	return reinterpret_cast<Zoltan2::EvaluatePartition<basic_id_t> *>
	  (new Zoltan2::EvaluatePartition<xpetra_mv_adapter>
	   (reinterpret_cast<const xpetra_mv_adapter *>(input),
	    params, CommT, reinterpret_cast
	    <const Zoltan2::PartitioningSolution<xpetra_mv_adapter> *>
	    (&problem->getSolution())));
      } else if (adapter_name == "XpetraCrsGraph") {
	return reinterpret_cast<Zoltan2::EvaluatePartition<basic_id_t> *>
	  (new Zoltan2::EvaluatePartition<xcrsGraph_adapter>
	   (reinterpret_cast<const xcrsGraph_adapter *>(input),
	    params, CommT, reinterpret_cast
	    <const Zoltan2::PartitioningSolution<xcrsGraph_adapter> *>
	    (&problem->getSolution())));
      } else if (adapter_name == "XpetraCrsMatrix") {
	return reinterpret_cast<Zoltan2::EvaluatePartition<basic_id_t> *>
	  (new Zoltan2::EvaluatePartition<xcrsMatrix_adapter>
	   (reinterpret_cast<const xcrsMatrix_adapter *>(input),
	    params, CommT, reinterpret_cast
	    <const Zoltan2::PartitioningSolution<xcrsMatrix_adapter> *>
	    (&problem->getSolution())));
      } else if (adapter_name == "BasicVector") {
	return reinterpret_cast<Zoltan2::EvaluatePartition<basic_id_t> *>
	  (new Zoltan2::EvaluatePartition<basic_vector_adapter>
	   (reinterpret_cast<const basic_vector_adapter *>(input),
	    params, CommT, reinterpret_cast
	    <const Zoltan2::PartitioningSolution<basic_vector_adapter> *>
	    (&problem->getSolution())));
      } else if (adapter_name == "PamgenMesh") {
	return reinterpret_cast<Zoltan2::EvaluatePartition<basic_id_t> *>
	  (new Zoltan2::EvaluatePartition<pamgen_adapter_t>
	   (reinterpret_cast<const pamgen_adapter_t *>(input),
	    params, CommT, reinterpret_cast
	    <const Zoltan2::PartitioningSolution<pamgen_adapter_t> *>
	    (&problem->getSolution())));
      }
      return nullptr; // adapter type not known
    }
  };
}
#endif // ZOLTAN2_EVALUATE_PARTITION_FACTORY_HPP
