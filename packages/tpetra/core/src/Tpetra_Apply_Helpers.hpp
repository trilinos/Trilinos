// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TPETRA_APPLY_HELPERS_HPP
#define TPETRA_APPLY_HELPERS_HPP
#include <type_traits>
#include "Teuchos_ScalarTraits.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "Tpetra_Details_Profiling.hpp"

namespace Tpetra {
  //! @name Methods performing multiple matrix-vector products at once
  //@{
  /// \brief Does multiply matrix apply() calls with a single X vector
  ///
  /// Computes Y[i] = beta*Y[i] + alpha * Matrices[i] * X, for i=1,...,N
  ///
  /// This routine only communicates the interprocessor portions of X once,
  /// if such a thing is possible (aka the ColMap's of the matrices match). 
  ///
  /// If the Y's are replicated this will not use the reduced communication path.
  ///
  /// If X is aliased to any Y but the last one, this routine will throw
  ///
  /// \param Matrices [in] - [std::vector|Teuchos::Array|Teuchos::ArrayRCP] of Tpetra::CrsMatrix objects.
  ///                        These matrices can different numbers of rows. 
  /// \param X [in]        - Tpetra::MultiVector or Tpetra::Vector object.
  /// \param Y [out]       - [std::vector|Teuchos::Array|Teuchos::ArrayRCP] of Tpetra::MultiVector or Tpetra::Vector objects.
  ///                        These must have the same number of vectors as X.
  /// \param alpha [in]    - alpha parameter.  Defaults to one.
  /// \param beta [in]    -  beta parameter.  Defaults to zero.
  /// \param params [in/out] - The "can batch" parameter can either be unset, be true or be false on input.  If it is unset,
  ///                        maps will be checked with isSameAs() for compatibility.  If true, the map check will be skipped and batching 
  ///                        will be used if none of the cheap checks fail.  If false, batching will not be used.  This parameter will
  ///                        be set on output to either true or false depending on if batching was used during this call. Defaults to NULL.

  template <class MatrixArray, class MultiVectorArray> 
  void batchedApply(const MatrixArray &Matrices, 
                    const typename std::remove_pointer<typename MultiVectorArray::value_type>::type &X,
                    MultiVectorArray &Y,
                    typename std::remove_pointer<typename MatrixArray::value_type>::type::scalar_type alpha = Teuchos::ScalarTraits< typename std::remove_pointer<typename MatrixArray::value_type>::type::scalar_type>::one(),
                    typename std::remove_pointer<typename MatrixArray::value_type>::type::scalar_type beta  = Teuchos::ScalarTraits< typename std::remove_pointer<typename MatrixArray::value_type>::type::scalar_type>::zero(),
                    Teuchos::RCP<Teuchos::ParameterList> params = Teuchos::null) {
    Tpetra::Details::ProfilingRegion regionTotal ("Tpetra::batchedApply");

    using size_type   = typename MatrixArray::size_type;
    using matrix_type = typename std::remove_pointer<typename MatrixArray::value_type>::type;
    using map_type    = typename matrix_type::map_type;
    using import_type = typename matrix_type::import_type;
    using export_type = typename matrix_type::export_type;
    using MV          = typename matrix_type::MV;
    using scalar_type = typename matrix_type::scalar_type;

    using Teuchos::RCP;
    using Teuchos::rcp_const_cast;

    const scalar_type ONE  = Teuchos::ScalarTraits<scalar_type>::one();
    const scalar_type ZERO = Teuchos::ScalarTraits<scalar_type>::zero();

    size_type N = Matrices.size();
    if(N==0) return;
    int numRanks = X.getMap()->getComm()->getSize();

    // If X is aliased to any Y but the last one, throw
    for(size_type i=0; i<N-1; i++) {
      TEUCHOS_TEST_FOR_EXCEPTION( &X == Y[i], std::runtime_error, "Tpetra::batchedApply(): X cannot be aliased to any Y except the final one.");
    }

    /* Checks for whether the reduced communication path can be used */
    RCP<const map_type> compare_colMap    = Matrices[0]->getColMap();
    RCP<const import_type> importer       = Matrices[0]->getGraph()->getImporter();

    
    bool can_batch, check_maps;
    if(params.is_null() || !params->isParameter("can batch")) {
      can_batch  = (importer.is_null() || N==1) ? false :true;
      check_maps = true;
    }
    else {
      can_batch  = (!params->get<bool>("can batch") || importer.is_null() || N==1) ? false : true;
      check_maps = false;
    }
    // We can't batch with replicated Y's
    if(numRanks > 1) {
      for(size_type i=0; i<N && can_batch; i++) {
        if(!Y[i]->isDistributed())
          can_batch = false;
      }
    }

    // Do the domain/column maps all match?
    for(size_type i=1; i<N && check_maps && can_batch; i++) {
      if(!Matrices[i]->getColMap()->isSameAs(*compare_colMap)) {
        can_batch=false;
      }
    }

    if(can_batch) {
      /* Batching path: Guarantees an existing importer and N>1 */

      // Special case for alpha = 0
      if (alpha == ZERO) {
        if (beta == ZERO) {
          for(size_type i=0; i<N; i++) Y[i]->putScalar(ZERO);
        } else if (beta != ONE) {
          for(size_type i=0; i<N; i++) Y[i]->scale(beta);
        }
        if(!params.is_null()) params->set("can batch",true);
        return;
      }
      
      const bool Y_is_overwritten = (beta == ZERO);

      // Start by importing X to Matrices[0]'s temporary
      RCP<const MV> X_colMap;
      {
        Tpetra::Details::ProfilingRegion regionImport ("Tpetra::batchedApply: Import");
        RCP<MV> X_colMapNonConst = Matrices[0]->getColumnMapMultiVector(X);

        // Import from the domain Map MV to the column Map MV.
        X_colMapNonConst->doImport(X, *importer, INSERT);
        X_colMap = rcp_const_cast<const MV>(X_colMapNonConst);
      }

      for(size_type i=0; i<N; i++) {
        RCP<const export_type> exporter = Matrices[i]->getGraph()->getExporter();

        // Temporary MV for doExport (if needed),
        RCP<MV> Y_rowMap = Matrices[i]->getRowMapMultiVector(*Y[i]);
        if (!exporter.is_null()) {
          Matrices[i]->localApply(*X_colMap, *Y_rowMap, Teuchos::NO_TRANS, alpha, ZERO);
          {
            Tpetra::Details::ProfilingRegion regionExport ("Tpetra::batchedApply: Export");             
            if (Y_is_overwritten) {
              Y[i]->putScalar(ZERO);
            }
            else {
              Y[i]->scale (beta);
            }
            Y[i]->doExport(*Y_rowMap, *exporter, ADD_ASSIGN);
          }
        }
        else { // Don't do an Export: row Map and range Map are the same.
          // Check for aliasing
          if (! Y[i]->isConstantStride() || X_colMap.getRawPtr() == Y[i]) {
            Y_rowMap = Matrices[i]->getRowMapMultiVector(*Y[i], true);
            if (beta != ZERO) {
              Tpetra::deep_copy (*Y_rowMap, *Y[i]);
            }

            Matrices[i]->localApply(*X_colMap, *Y_rowMap, Teuchos::NO_TRANS, alpha, beta);
            Tpetra::deep_copy(*Y[i], *Y_rowMap);
          }
          else {
            Matrices[i]->localApply(*X_colMap, *Y[i], Teuchos::NO_TRANS, alpha, beta);
          }
        }        
      }
      if(!params.is_null()) params->set("can batch",true);
    } 
    else {
      /* Non-batching path */
      for(size_type i=0; i<N; i++) {
        Matrices[i]->apply(X,*Y[i],Teuchos::NO_TRANS, alpha, beta);
      }
      if(!params.is_null()) params->set("can batch",false);
    }
  }
  //@}

}// namespace Tpetra

#endif // TPETRA_APPLY_HELPERS_HPP

