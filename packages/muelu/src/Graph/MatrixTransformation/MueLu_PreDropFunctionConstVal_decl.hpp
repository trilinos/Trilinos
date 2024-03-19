// @HEADER
//
// ***********************************************************************
//
//        MueLu: A package for multigrid based preconditioning
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
// Questions? Contact
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Andrey Prokopenko (aprokop@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
#ifndef MUELU_PREDROPFUNCTIONCONSTVAL_DECL_HPP
#define MUELU_PREDROPFUNCTIONCONSTVAL_DECL_HPP

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_PreDropFunctionBaseClass.hpp"
#include "MueLu_PreDropFunctionConstVal_fwd.hpp"

namespace MueLu {

/*!
 * Example implementation for dropping values smaller then a constant threshold
 *
 */
template <class Scalar        = DefaultScalar,
          class LocalOrdinal  = DefaultLocalOrdinal,
          class GlobalOrdinal = DefaultGlobalOrdinal,
          class Node          = DefaultNode>
class PreDropFunctionConstVal : public MueLu::PreDropFunctionBaseClass<Scalar, LocalOrdinal, GlobalOrdinal, Node> {
#undef MUELU_PREDROPFUNCTIONCONSTVAL_SHORT
#include "MueLu_UseShortNames.hpp"

 public:
  //! Constructor
  explicit PreDropFunctionConstVal(const Scalar threshold = 0.0);

  //! Destructor
  virtual ~PreDropFunctionConstVal() {}

  /*! Drop
   * @param lrow (size_t): local row index (=lrowid)
   * @param grow (GlobalOrdinal: global row id
   * @param k    (size_t): local column iterator
   * @param lcid (LocalOrdinal): local column id (=indices[k])
   * @param gcid (GlobalOrdinal): global column id
   * @param indices (ArrrayView): array of local column ids in current row (lrow)
   * @param vals (ArrayView): array of corresponding values in current row (lrow)
   * @return bool: false, if value in (lrow, lcid) shall be kept, true if it should be dropped
   */
  bool Drop(size_t lrow, GlobalOrdinal grow, size_t k, LocalOrdinal lcid, GlobalOrdinal gcid, const Teuchos::ArrayView<const LocalOrdinal>& indices, const Teuchos::ArrayView<const Scalar>& vals);

  //! Return threshold value.
  Scalar GetThreshold() const;

  //! @name Overridden from Teuchos::Describable
  //@{

  //! Return a simple one-line description of this object.
  std::string description() const;

  //! Print the object with some verbosity level to an FancyOStream object.
  // using MueLu::Describable::describe; // overloading, not hiding
  // void describe(Teuchos::FancyOStream &out, const VerbLevel verbLevel = Default) const;

  //@}

 private:
  Scalar threshold_;
};

}  // namespace MueLu

#define MUELU_PREDROPFUNCTIONCONSTVAL_SHORT
#endif  // MUELU_PREDROPFUNCTIONCONSTVAL_DECL_HPP
