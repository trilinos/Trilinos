// @HEADER
// ***********************************************************************
//
//                           Stokhos Package
//                 Copyright (2009) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// Questions? Contact Eric T. Phipps (etphipp@sandia.gov).
//
// ***********************************************************************
// @HEADER

#ifndef STOKHOS_SGMODELEVALUATORBASE_HPP
#define STOKHOS_SGMODELEVALUATORBASE_HPP

#include "EpetraExt_ModelEvaluator.h"
#include "EpetraExt_BlockVector.h"

#include "Teuchos_RCP.hpp"
#include "Teuchos_Array.hpp"
#include "Stokhos_EpetraVectorOrthogPoly.hpp"

namespace Stokhos {

  //! Base class for stochastic Galerkin model evaluators
  class SGModelEvaluatorBase : public virtual EpetraExt::ModelEvaluator {
  public:

    // Constructor
    SGModelEvaluatorBase() {}

    //! Destructor
    virtual ~SGModelEvaluatorBase() {}

    //! Set initial solution polynomial
    virtual void set_x_sg_init(const Stokhos::EpetraVectorOrthogPoly& x_sg_in) = 0;

    //! Return initial SG x
    virtual Teuchos::RCP<const Stokhos::EpetraVectorOrthogPoly> get_x_sg_init() const = 0;

    //! Set initial parameter polynomial
    virtual void set_p_sg_init(int i, const Stokhos::EpetraVectorOrthogPoly& p_sg_in) = 0;

    //! Return initial SG parameters
    virtual Teuchos::RCP<const Stokhos::EpetraVectorOrthogPoly> get_p_sg_init(int l) const = 0;

    //! Get indices of SG parameters
    /*!
     * These indices determine which parameter vectors support SG
     */
    virtual Teuchos::Array<int> get_p_sg_map_indices() const = 0;

    //! Get indices of SG responses
    /*!
     * These indices determine which response vectors support SG
     */
    virtual Teuchos::Array<int> get_g_sg_map_indices() const = 0;

    //! Get base maps of SG responses
    virtual Teuchos::Array< Teuchos::RCP<const Epetra_Map> > get_g_sg_base_maps() const = 0;

    //! Return overlap stochastic map
    virtual Teuchos::RCP<const Epetra_BlockMap> get_overlap_stochastic_map() const = 0;

    //! Return x sg overlap map
    virtual Teuchos::RCP<const Epetra_BlockMap> get_x_sg_overlap_map() const = 0;

    //! Return x sg importer
    virtual Teuchos::RCP<const Epetra_Import> get_x_sg_importer() const = 0;

    //! Create vector orthog poly using x map and owned sg map
    virtual Teuchos::RCP<Stokhos::EpetraVectorOrthogPoly>
    create_x_sg(Epetra_DataAccess CV = Copy,
                const Epetra_Vector* v = NULL) const = 0;

    //! Create vector orthog poly using x map and overlap sg map
    virtual Teuchos::RCP<Stokhos::EpetraVectorOrthogPoly>
    create_x_sg_overlap(Epetra_DataAccess CV = Copy,
                        const Epetra_Vector* v = NULL) const = 0;

    //! Create vector orthog poly using x map and owned sg map
    virtual Teuchos::RCP<Stokhos::EpetraMultiVectorOrthogPoly>
    create_x_mv_sg(int num_vecs,
                Epetra_DataAccess CV = Copy,
                const Epetra_MultiVector* v = NULL) const = 0;

    //! Create vector orthog poly using x map and overlap sg map
    virtual Teuchos::RCP<Stokhos::EpetraMultiVectorOrthogPoly>
    create_x_mv_sg_overlap(int num_vecs,
                           Epetra_DataAccess CV = Copy,
                           const Epetra_MultiVector* v = NULL) const = 0;

    //! Create vector orthog poly using p map
    virtual Teuchos::RCP<Stokhos::EpetraVectorOrthogPoly>
    create_p_sg(int l, Epetra_DataAccess CV = Copy,
                const Epetra_Vector* v = NULL) const = 0;

    //! Create vector orthog poly using f map and owned sg map
    virtual Teuchos::RCP<Stokhos::EpetraVectorOrthogPoly>
    create_f_sg(Epetra_DataAccess CV = Copy,
                const Epetra_Vector* v = NULL) const = 0;

    //! Create vector orthog poly using f map and overlap sg map
    virtual Teuchos::RCP<Stokhos::EpetraVectorOrthogPoly>
    create_f_sg_overlap(Epetra_DataAccess CV = Copy,
                        const Epetra_Vector* v = NULL) const = 0;

    //! Create multi-vector orthog poly using f map and owned sg map
    virtual Teuchos::RCP<Stokhos::EpetraMultiVectorOrthogPoly>
    create_f_mv_sg(int num_vecs, Epetra_DataAccess CV = Copy,
                   const Epetra_MultiVector* v = NULL) const = 0;

    //! Create multi-vector orthog poly using f map and overlap sg map
    virtual Teuchos::RCP<Stokhos::EpetraMultiVectorOrthogPoly>
    create_f_mv_sg_overlap(int num_vecs, Epetra_DataAccess CV = Copy,
                           const Epetra_MultiVector* v = NULL) const = 0;

    //! Create vector orthog poly using g map
    virtual Teuchos::RCP<Stokhos::EpetraVectorOrthogPoly>
    create_g_sg(int l, Epetra_DataAccess CV = Copy,
                const Epetra_Vector* v = NULL) const = 0;

    //! Create multi-vector orthog poly using g map
    virtual Teuchos::RCP<Stokhos::EpetraMultiVectorOrthogPoly>
    create_g_mv_sg(int l, int num_vecs, Epetra_DataAccess CV = Copy,
                const Epetra_MultiVector* v = NULL) const = 0;

  };

}

#endif // STOKHOS_SGMODELEVALUATORBASE_HPP
