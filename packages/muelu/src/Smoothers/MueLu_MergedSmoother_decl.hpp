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
#ifndef MUELU_MERGEDSMOOTHER_DECL_HPP
#define MUELU_MERGEDSMOOTHER_DECL_HPP

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_MergedSmoother_fwd.hpp"
#include "MueLu_SmootherPrototype.hpp"

/*! @class MergedSmoother
    @ingroup MueLuSmootherClasses 
*/

namespace MueLu {

  class Level;

  template <class Scalar = SmootherPrototype<>::scalar_type,
            class LocalOrdinal = typename SmootherPrototype<Scalar>::local_ordinal_type,
            class GlobalOrdinal = typename SmootherPrototype<Scalar, LocalOrdinal>::global_ordinal_type,
            class Node = typename SmootherPrototype<Scalar, LocalOrdinal, GlobalOrdinal>::node_type>
  class MergedSmoother : public SmootherPrototype<Scalar,LocalOrdinal,GlobalOrdinal,Node> {
#undef MUELU_MERGEDSMOOTHER_SHORT
#include "MueLu_UseShortNames.hpp"

  public:
    //! @name Constructors / destructors
    //@{

    //! Constructor
    MergedSmoother(ArrayRCP<RCP<SmootherPrototype> >& smootherList, bool verbose = false);

    //! Copy constructor (performs a deep copy of input object)
    MergedSmoother(const MergedSmoother& src);

    //! Copy method (performs a deep copy of input object)
    RCP<SmootherPrototype> Copy() const;

    //! Destructor
    virtual ~MergedSmoother() { }
    //@}

    //! @name Set/Get methods
    //@{

    void StandardOrder() { reverseOrder_ = false; }
    void ReverseOrder()  { reverseOrder_ = true;  }

    // TODO: GetOrder() is a better name (+ enum to define order)
    bool GetReverseOrder() const { return reverseOrder_; }

    //  TODO  const ArrayRCP<const RCP<const SmootherPrototype> > & GetSmootherList() const;
    const ArrayRCP<const RCP<SmootherPrototype> > GetSmootherList() const { return smootherList_; }

    //@}

    void DeclareInput(Level &currentLevel) const;

    //! @name Setup and Apply methods.
    //@{

    /*! @brief Set up. */
    void Setup(Level& level);

    /*! @brief Apply

    Solves the linear system <tt>AX=B</tt> using the smoothers of the list.

    @param X initial guess
    @param B right-hand side
    @param InitialGuessIsZero
    */
    void Apply(MultiVector& X, const MultiVector& B, bool InitialGuessIsZero = false) const;

    //@}

    //! Custom SetFactory
    void SetFactory(const std::string& varName, const RCP<const FactoryBase>& factory);

    void print(Teuchos::FancyOStream& out, const VerbLevel verbLevel = Default) const;

    void CopyParameters(RCP<SmootherPrototype> src); // TODO: wrong prototype. We do not need an RCP here.

    ArrayRCP<RCP<SmootherPrototype> > SmootherListDeepCopy(const ArrayRCP<const RCP<SmootherPrototype> >& srcSmootherList);

    //! Get a rough estimate of cost per iteration
    size_t getNodeSmootherComplexity() const;

    //@}

  private:
    // List of smoothers. It is an ArrayRCP of RCP because:
    //  1) I need a vector of pointers (to avoid slicing problems)
    //  2) I can use an std::vector insead of an ArrayRCP but then the constructor will do a copy of user input
    ArrayRCP<RCP<SmootherPrototype> > smootherList_;

    //
    bool reverseOrder_;

    // tmp, for debug
    bool verbose_;

  }; //class MergedSmoother

} //namespace MueLu

#define MUELU_MERGEDSMOOTHER_SHORT
#endif // MUELU_MERGEDSMOOTHER_DECL_HPP
