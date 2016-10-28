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

/*! \file Zoltan2_EvaluateOrderingFactory.hpp
    \brief Returns a pointer to a new evaluate ordering class.
    Is not responsible for memory management!
*/
#ifndef ZOLTAN2_EVALUATE_ORDERING_FACTORY_HPP
#define ZOLTAN2_EVALUATE_ORDERING_FACTORY_HPP
namespace Zoltan2_TestingFramework {

///brief EvaluateOrderingFactory class contains 1 static factory method
class EvaluateOrderingFactory {
  public:
    /// \brief Zoltan2::EvaluateOrdering factory method
    ///
    /// @param problem ordering problem
    /// @param adapter_name a string equal to the type of the adapter
    /// @param input Input adapter used to construct the ordering problem
    /// @param input params Zoltan2 parameters
    ///
    /// @return returns a pointer to new Zoltan2::EvaluateOrdering or a nullptr
    // if adapter type was not known.

    static EvaluateOrdering<basic_id_t>* newEvaluateOrdering
    (ordering_problem_t *problem, const std::string &adapter_name,
      base_adapter_t *input, ParameterList *params) {

      #define ZOLTAN2_EVAL_LOCAL_ORDERING(stringName, adapterClass)              \
      if (adapter_name == stringName) {                                          \
        return reinterpret_cast<EvaluateOrdering<basic_id_t> *>                  \
          (new EvaluateLocalOrdering<adapterClass>                               \
            (reinterpret_cast<const adapterClass *>(input),                      \
            params, problem->getComm(),                                          \
            reinterpret_cast<const LocalOrderingSolution<adapterClass::lno_t> *> \
              (problem->getLocalOrderingSolution())));                           \
      }

      // EvaluateGlobalOrdering not tested/implemented yet
      #define ZOLTAN2_EVAL_GLOBAL_ORDERING(stringName, adapterClass)             \
      if (adapter_name == stringName) {                                          \
        return reinterpret_cast<EvaluateOrdering<basic_id_t> *>                  \
          (new EvaluateGlobalOrdering<adapterClass>                              \
            (reinterpret_cast<const adapterClass *>(input),                      \
            params, problem->getComm(),                                          \
            reinterpret_cast<const GlobalOrderingSolution<adapterClass::gno_t> *>\
              (problem->getGlobalOrderingSolution())));                          \
      }

      ZOLTAN2_EVAL_LOCAL_ORDERING("BasicIdentifier",   basic_id_t)
      ZOLTAN2_EVAL_LOCAL_ORDERING("XpetraMultiVector", xpetra_mv_adapter)
      ZOLTAN2_EVAL_LOCAL_ORDERING("XpetraCrsGraph",    xcrsGraph_adapter)
      ZOLTAN2_EVAL_LOCAL_ORDERING("XpetraCrsMatrix",   xcrsMatrix_adapter)
      ZOLTAN2_EVAL_LOCAL_ORDERING("BasicVector",       basic_vector_adapter)
      ZOLTAN2_EVAL_LOCAL_ORDERING("PamgenMesh",        pamgen_adapter_t)
      return nullptr;
    }
};

} // Zoltan2_TestingFramework

#endif // ZOLTAN2_EVALUATE_ORDERING_FACTORY_HPP
