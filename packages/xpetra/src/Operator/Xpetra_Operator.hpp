// @HEADER
//
// ***********************************************************************
//
//             Xpetra: A linear algebra interface package
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
#ifndef XPETRA_OPERATOR_HPP
#define XPETRA_OPERATOR_HPP

#include "Xpetra_ConfigDefs.hpp"

#include <Teuchos_Describable.hpp>
#include <Teuchos_BLAS_types.hpp>
#include <Teuchos_ScalarTraits.hpp>

#include "Xpetra_Map.hpp"
#include "Xpetra_MultiVector.hpp"

namespace Xpetra {

  template<class Scalar         = MultiVector<>::scalar_type,
           class LocalOrdinal   = typename MultiVector<Scalar>::local_ordinal_type,
           class GlobalOrdinal  = typename MultiVector<Scalar, LocalOrdinal>::global_ordinal_type,
           class Node           = typename MultiVector<Scalar, LocalOrdinal, GlobalOrdinal>::node_type>
  class Operator : virtual public Teuchos::Describable {
    typedef Xpetra::Map<LocalOrdinal,GlobalOrdinal,Node> Map;
    typedef Xpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> MultiVector;
  public:
    virtual ~Operator() { }
    /** \name Typedefs that give access to the template parameters. */
    //@{

    //! The type of the entries of the input and output multivectors.
    typedef Scalar scalar_type;

    //! The local index type.
    typedef LocalOrdinal local_ordinal_type;

    //! The global index type.
    typedef GlobalOrdinal global_ordinal_type;

    //! The Kokkos Node type.
    typedef Node node_type;

    //@}
    /** \name Pure virtual functions to be overridden by subclasses. */
    //@{

    //! The Map associated with the domain of this operator, which must be compatible with X.getMap().
    virtual Teuchos::RCP<const Map> getDomainMap() const = 0;

    //! The Map associated with the range of this operator, which must be compatible with Y.getMap().
    virtual Teuchos::RCP<const Map> getRangeMap() const = 0;

    //! \brief Computes the operator-multivector application.
    /*! Loosely, performs \f$Y = \alpha \cdot A^{\textrm{mode}} \cdot X + \beta \cdot Y\f$. However, the details of operation
        vary according to the values of \c alpha and \c beta. Specifically
        - if <tt>beta == 0</tt>, apply() <b>must</b> overwrite \c Y, so that any values in \c Y (including NaNs) are ignored.
        - if <tt>alpha == 0</tt>, apply() <b>may</b> short-circuit the operator, so that any values in \c X (including NaNs) are ignored.
     */
    virtual void
    apply (const MultiVector& X, MultiVector& Y,
           Teuchos::ETransp mode = Teuchos::NO_TRANS,
           Scalar alpha = Teuchos::ScalarTraits<Scalar>::one(),
           Scalar beta = Teuchos::ScalarTraits<Scalar>::zero()) const = 0;

    /// \brief Whether this operator supports applying the transpose or conjugate transpose.
    ///
    /// By default, this returns false.  Subclasses must override this
    /// method if they can support apply() with
    /// <tt>mode=Teuchos::TRANS</tt> or
    /// <tt>mode=Teuchos::CONJ_TRANS</tt>.
    virtual bool hasTransposeApply() const { return false; }

    //@}

    virtual void removeEmptyProcessesInPlace(const RCP<const Map>& newMap) { }
  };

} // Xpetra namespace

#define XPETRA_OPERATOR_SHORT
#endif // XPETRA_OPERATOR_HPP
