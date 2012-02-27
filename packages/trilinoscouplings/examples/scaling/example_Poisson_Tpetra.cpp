// @HEADER
// ************************************************************************
//
//                           Intrepid Package
//                 Copyright (2007) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Pavel Bochev  (pbboche@sandia.gov),
//                    Denis Ridzal  (dridzal@sandia.gov),
//                    Kara Peterson (kjpeter@sandia.gov).
//
// ************************************************************************
// @HEADER

/** \file   example_Poisson.cpp
    \brief  Example solution of a Poisson equation on a hexahedral mesh using
            nodal (Hgrad) elements.

           This example uses the following Trilinos packages:
    \li     Pamgen to generate a Hexahedral mesh.
    \li     Sacado to form the source term from user-specified manufactured solution.
    \li     Intrepid to build the discretization matrix and right-hand side.
    \li     Epetra to handle the global matrix and vector.
    \li     Isorropia to partition the matrix. (Optional)
    \li     ML to solve the linear system.


    \verbatim

     Poisson system:
 
            div A grad u = f in Omega
                       u = g on Gamma

      where
             A is a symmetric, positive definite material tensor
             f is a given source term

     Corresponding discrete linear system for nodal coefficients(x):

                 Kx = b

            K - HGrad stiffness matrix
            b - right hand side vector

    \endverbatim

    \author Created by P. Bochev, D. Ridzal, K. Peterson, D. Hensinger, C. Siefert.

    \remark Usage:
    \code   ./TrilinosCouplings_examples_scaling_example_Poisson.exe \endcode

    \remark Example driver requires input file named Poisson.xml with Pamgen 
            formatted mesh description and settings for Isorropia (a version 
            is included in the Trilinos repository with this driver).

    \remark The exact solution (u) and material tensor (A) are set in the
            functions "exactSolution" and "materialTensor" and may be
            modified by the user.
            
*/

/*** Uncomment if you would like output data for plotting ***/
//#define DUMP_DATA

/**************************************************************/
/*                          Includes                          */
/**************************************************************/

// TrilinosCouplings includes
#include "TrilinosCouplings_config.h"

// Intrepid includes
#include "Intrepid_FunctionSpaceTools.hpp"
#include "Intrepid_CellTools.hpp"
#include "Intrepid_ArrayTools.hpp"
#include "Intrepid_HGRAD_HEX_C1_FEM.hpp"
#include "Intrepid_RealSpaceTools.hpp"
#include "Intrepid_DefaultCubatureFactory.hpp"
#include "Intrepid_Utils.hpp"

// Epetra includes
#include "Epetra_Time.h"
#include "Epetra_Map.h"
#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_FECrsMatrix.h"
#include "Epetra_FEVector.h"
#include "Epetra_Import.h"

// Teuchos includes
#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_BLAS.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"

// Tpetra includes
#include "Tpetra_MultiVector.hpp"
#include "Teuchos_CommHelpers.hpp"
#include <iterator>

// Shards includes
#include "Shards_CellTopology.hpp"

// EpetraExt includes
#include "EpetraExt_RowMatrixOut.h"
#include "EpetraExt_MultiVectorOut.h"

// Pamgen includes
#include "create_inline_mesh.h"
#include "im_exodusII_l.h"
#include "im_ne_nemesisI_l.h"
#include "pamgen_extras.h"

// AztecOO includes
#include "AztecOO.h"

// ML Includes
#include "ml_MultiLevelPreconditioner.h"
#include "ml_epetra_utils.h"

#ifdef TrilinosCouplings_ENABLE_Isorropia
#define TC_HAVE_ISORROPIA
#endif

#ifdef TC_HAVE_ISORROPIA
// Isorropia includes
#include "Isorropia_Epetra.hpp"
#include "Isorropia_EpetraRedistributor.hpp"
#include "Isorropia_EpetraPartitioner.hpp"
#endif

// Sacado includes
#include "Sacado.hpp"

using namespace Intrepid;

namespace {
  /// \class MultiVectorFillerData
  /// \brief Implementation of fill and local assembly for \c MultiVectorFiller.
  /// \author Mark Hoemmen
  ///
  /// \tparam MV Specialization of \c Tpetra::MultiVector.
  template<class MV>
  class MultiVectorFillerData {
  public:
    typedef typename MV::scalar_type scalar_type;
    typedef typename MV::local_ordinal_type local_ordinal_type;
    typedef typename MV::global_ordinal_type global_ordinal_type;
    typedef typename MV::node_type node_type;

    /// \brief Default constructor (sets number of columns to zero).
    ///
    /// Before using this object, you should call \c setNumColumns()
    /// to set the number of columns in the output multivector.
    /// Otherwise, the two-argument version of \c
    /// sumIntoGlobalValues() won't actually do anything.
    MultiVectorFillerData () : 
      numCols_ (0)
    {}

    /// \brief Constructor.
    ///
    /// \param numColumns [in] The (expected) number of columns in the
    ///   output multivector.  You can always change this later by
    ///   calling \c setNumColumns().
    ///
    /// \note If the number of columns given here is not the same as
    ///   the number of columns in the output multivector, you should
    ///   call \c setNumColumns() first before inserting any data.
    ///   Otherwise, the two-argument version of \c
    ///   sumIntoGlobalValues() won't do the right thing.
    MultiVectorFillerData (const size_t numColumns) : 
      numCols_ (numColumns),
      sourceIndices_ (numColumns),
      sourceValues_ (numColumns)
    {}

    //! Set the number of columns in the output multivector.
    void
    setNumColumns (const size_t newNumColumns) 
    {
      const size_t oldNumColumns = getNumColumns();
      if (newNumColumns >= oldNumColumns) {
	for (size_t j = oldNumColumns; j < newNumColumns; ++j) {
	  sourceIndices_.push_back (Teuchos::Array<const global_ordinal_type> ());
	  sourceValues_.push_back (Teuchos::Array<const scalar_type> ());
	}
      } 
      else {
	// This may not necessarily deallocate any data, but that's OK.
	sourceIndices_.resize (newNumColumns);
	sourceValues_.resize (newNumColumns);
      }
      numCols_ = oldNumColumns;
    }

    void
    sumIntoGlobalValues (Teuchos::ArrayView<const global_ordinal_type> rows, 
			 size_t column,
			 Teuchos::ArrayView<const scalar_type> values)
    {
      if (column >= getNumColumns()) {
	for (size_t j = column; j < getNumColumns(); ++j) {
	  sourceIndices_.push_back (Teuchos::Array<const global_ordinal_type> ());
	  sourceValues_.push_back (Teuchos::Array<const scalar_type> ());
	}
      }
      std::copy (rows.begin(), rows.end(), std::back_inserter (sourceIndices_[column]));
      std::copy (values.begin(), values.end(), std::back_inserter (sourceValues_[column]));
    }

    /// Data for each column are stored contiguously in rows and in
    /// values.  Thus, rows and values are in rowwise order, even
    /// though they may be stored in columnwise order in the
    /// multivector.
    ///
    /// Be sure that the number of columns is set correctly before
    /// calling this.
    void
    sumIntoGlobalValues (Teuchos::ArrayView<const global_ordinal_type> rows, 
			 Teuchos::ArrayView<const scalar_type> values)
    {
      typedef global_ordinal_type GO;
      typedef scalar_type ST;
      typedef typename Teuchos::ArrayView<const GO>::const_iterator GoIter;
      typedef typename Teuchos::ArrayView<const ST>::const_iterator StIter;

      const size_t numColumns = getNumColumns();
      GoIter rowIter = rows.begin();
      StIter valIter = values.begin();
      for (size_t j = 0; j < numColumns; ++j) {
	GoIter rowIterNext = rowIter + numColumns;
	StIter valIterNext = valIter + numColumns;
	std::copy (rowIter, rowIterNext, std::back_inserter (sourceIndices_[j]));
	std::copy (valIter, valIterNext, std::back_inserter (sourceValues_[j]));
	rowIter = rowIterNext;
	valIter = valIterNext;
      }
    }

    void
    locallyAssemble (MV& X)
    {
      using Teuchos::Array;
      using Teuchos::ArrayRCP;
      using Teuchos::ArrayView;
      using Teuchos::RCP;
      typedef local_ordinal_type LO;
      typedef global_ordinal_type GO;
      typedef scalar_type ST;
      typedef node_type NT;
      typedef Tpetra::Map<LO, GO, NT> map_type;

      RCP<const map_type> map = X->getMap();
      Array<LO> localIndices;
      const size_t numColumns = getNumColumns();
      for (size_t j = 0; j < numColumns; ++j) {
	const typename Array<const GO>::size_type numIndices = sourceIndices_[j].size();
	// Precompute all the local indices before saving to the
	// vector.  Hopefully this gives us a little bit of locality
	// in the global->local conversion, at the expense of a little
	// more storage.
	if (localIndices.size() < numIndices) {
	  localIndices.resize (numIndices);
	}
	ArrayView<LO> localIndicesView = localIndices.view (0, numIndices);
	ArrayView<const GO> globalIndicesView = sourceIndices_[j].view (0, numIndices);
	for (typename ArrayView<const GO>::size_type i = 0; i < numIndices; ++i) {
	  localIndices[i] = map->getLocalElement (globalIndicesView[i]);
	}

	ArrayRCP<ST> X_j = X->getDataNonConst (j);
	ArrayView<const ST> localValues = sourceValues_[j].view (0, numIndices);
	for (typename ArrayView<const GO>::size_type i = 0; i < numIndices; ++i) {
	  X_j[localIndices[i]] = localValues[i];
	}
      }
    }

    //! Clear the contents of the vector, making it implicitly a vector of zeros.
    void clear() {
      Teuchos::Array<Teuchos::Array<global_ordinal_type> > newSourceIndices;
      Teuchos::Array<Teuchos::Array<scalar_type> > newSourceValues;
      // The standard STL idiom for clearing the contents of a vector completely.
      std::swap (sourceIndices_, newSourceIndices);
      std::swap (sourceValues_, newSourceValues);
    }

  private:
    size_t numCols_;
    Teuchos::Array<Teuchos::Array<global_ordinal_type> > sourceIndices_;
    Teuchos::Array<Teuchos::Array<scalar_type> > sourceValues_;

    size_t getNumColumns() const { return numCols_; }
  };



  /// \class MultiVectorFillerData2
  /// \brief Second implementation of fill and local assembly for \c MultiVectorFiller.
  /// \author Mark Hoemmen
  ///
  /// \tparam MV Specialization of \c Tpetra::MultiVector.
  template<class MV>
  class MultiVectorFillerData2 {
  public:
    typedef typename MV::scalar_type scalar_type;
    typedef typename MV::local_ordinal_type local_ordinal_type;
    typedef typename MV::global_ordinal_type global_ordinal_type;
    typedef typename MV::node_type node_type;

    typedef Tpetra::Map<local_ordinal_type, global_ordinal_type, node_type> map_type;

    /// \brief Default constructor (sets number of columns to zero).
    ///
    /// \param map [in] Map over which to distribute the initial fill.
    ///
    /// Before using this object, you should call \c setNumColumns()
    /// to set the number of columns in the output multivector.
    /// Otherwise, the two-argument version of \c
    /// sumIntoGlobalValues() won't actually do anything.
    MultiVectorFillerData2 (const Teuchos::RCP<const map_type>& map) : 
      map_ (map),
      numCols_ (0)
    {}

    /// \brief Constructor.
    ///
    /// \param map [in] Map over which to distribute the initial fill.
    ///
    /// \param numColumns [in] The (expected) number of columns in the
    ///   output multivector.  You can always change this later by
    ///   calling \c setNumColumns().
    ///
    /// \note If the number of columns given here is not the same as
    ///   the number of columns in the output multivector, you should
    ///   call \c setNumColumns() first before inserting any data.
    ///   Otherwise, the two-argument version of \c
    ///   sumIntoGlobalValues() won't do the right thing.
    MultiVectorFillerData2 (const Teuchos::RCP<const map_type>& map,
			    const size_t numColumns) : 
      map_ (map),
      numCols_ (numColumns),
      localVec_ (new MV (map, numColumns)),
      nonlocalIndices_ (numColumns),
      nonlocalValues_ (numColumns)
    {}

    //! Set the number of columns in the output multivector.
    void
    setNumColumns (const size_t newNumColumns) 
    {
      using Teuchos::Array;
      using Teuchos::Range1D;
      using Teuchos::RCP;
      typedef global_ordinal_type GO;
      typedef scalar_type ST;

      const size_t oldNumColumns = numCols_;
      if (newNumColumns == oldNumColumns) {
	return; // No side effects if no change.
      }

      RCP<MV> newLocalVec;
      if (newNumColumns > oldNumColumns) {
	newLocalVec = (new MV (map_, newNumColumns));
	// Assign the contents of the old local multivector to the
	// first oldNumColumns columns of the new local multivector,
	// then get rid of the old local multivector.
	RCP<MV> newLocalVecView = 
	  newLocalVec.subViewNonConst (Range1D (0, oldNumColumns-1));
	*newLocalVecView = *localVec_;
      } 
      else {
	if (newNumColumns == 0) {
	  // Tpetra::MultiVector doesn't let you construct a
	  // multivector with zero columns.
	  newLocalVec = Teuchos::null;
	}
	else {
	  newLocalVec = 
	    localVec_.subViewNonConst (Range1D (0, newNumColumns-1));
	}
      }

      // Leave most side effects until the end, for exception safety.
      nonlocalIndices_.resize (newNumColumns);
      nonlocalValues_.resize (newNumColumns);
      localVec_ = newLocalVec;
      numCols_ = newNumColumns;
    }

    void
    sumIntoGlobalValues (Teuchos::ArrayView<const global_ordinal_type> rows, 
			 size_t columnIndex,
			 Teuchos::ArrayView<const scalar_type> values)
    {
      using Teuchos::ArrayRCP;
      using Teuchos::ArrayView;
      typedef local_ordinal_type LO;
      typedef global_ordinal_type GO;
      typedef scalar_type ST;

      if (columnIndex >= getNumColumns()) {
	// Automatically expand the number of columns.  This
	// implicitly ensures that localVec_ is not null.
	setNumColumns (columnIndex + 1);
      }
      
      typename ArrayView<const GO>::const_iterator rowIter = rows.begin();
      typename ArrayView<const ST>::const_iterator valIter = values.begin();
      for ( ; rowIter != rows.end() && valIter != values.end(); ++rowIter, ++valIter) {
	const GO globalRowIndex = *rowIter;
	// Converting from global to local index could be logarithmic
	// in the number of global indices that this process owns,
	// depending on the Map implementation.  However, the lookup
	// allows us to store data in the local multivector, rather
	// than in a separate data structure.
	const LO localRowIndex = map_->getLocalElement (globalRowIndex);
	if (localRowIndex == Teuchos::OrdinalTraits<LO>::invalid()) {
	  nonlocalIndices_[columnIndex].push_back (globalRowIndex);
	  nonlocalValues_[columnIndex].push_back (*valIter);
	}
	else {
	  // FIXME (mfh 27 Feb 2012) This will be very slow for GPU
	  // Node types.  In that case, we should hold on to the view
	  // of localVec_ as long as the number of columns doesn't
	  // change, and make modifications to the view until
	  // localAssemble() is called.
	  ArrayRCP<ST> X_j = localVec_->getDataNonConst (columnIndex);
	  // FIXME (mfh 27 Feb 2012) Allow different combine modes.
	  // The current combine mode just adds to the current value
	  // at that location.
	  X_j[localRowIndex] += *valIter;
	}
      }
    }

    /// Data for each column are stored contiguously in rows and in
    /// values.  Thus, rows and values are in rowwise order, even
    /// though they may be stored in columnwise order in the
    /// multivector.
    ///
    /// Be sure that the number of columns is set correctly before
    /// calling this.
    void
    sumIntoGlobalValues (Teuchos::ArrayView<const global_ordinal_type> rows, 
			 Teuchos::ArrayView<const scalar_type> values)
    {
      using Teuchos::ArrayView;
      typedef typename ArrayView<const global_ordinal_type>::size_type size_type;

      const size_t numCols = getNumColumns();
      for (size_t j = 0; j < numCols; ++j) {
	const size_type offset = numCols*j;
	const size_type len = numCols;
	sumIntoGlobalValues (rows.view (offset, len), j, values.view (offset, len));
      }
    }

    /// \brief Locally assemble into X.
    ///
    /// \param X [in/out] Multivector (overlapping source distribution).
    ///
    /// \param f [in/out] Binary function that defines the combine
    ///   mode.  It must define scalar_type operator (const
    ///   scalar_type&, const scalar_type&).  It need not necessarily
    ///   be commutative or even associative, but it should be
    ///   thread-safe in case we decide to parallelize local assembly.
    ///   We call it via X(i,j) = f(X(i,j), Y(i,j)), so write your
    ///   possibly nonassociative or noncommutative operation
    ///   accordingly.
    ///
    /// X is distributed by the source Map (with possible overlap) of
    /// the Export operation.  The source Map of the Export includes
    /// both the elements owned by this object's constructor's input
    /// Map, and the indices inserted by \c sumIntoGlobalValues().
    ///
    /// Precondition: The set of global indices in X's Map equals the
    /// union of the global indices in map_ and (the union of the
    /// entries of nonlocalIndices_[j] for all valid columns j).
    ///
    /// \note You can get the usual ADD combine mode by supplying f =
    ///   std::plus<scalar_type>.
    template<class BinaryFunction>
    void 
    locallyAssemble (MV& X, BinaryFunction& f)
    {
      using Teuchos::ArrayRCP;
      using Teuchos::ArrayView;
      using Teuchos::RCP;
      typedef local_ordinal_type LO;
      typedef global_ordinal_type GO;
      typedef scalar_type ST;

      RCP<const map_type> srcMap = X.getMap();
      ArrayView<const GO> localIndices = map_->getNodeElementList ();

      for (size_t j = 0; j < X.getNumVectors(); ++j) {
	ArrayRCP<ST> X_j = X.getDataNonConst (j);

	// First add all the local data into X_j.
	ArrayRCP<const ST> local_j = localVec_.getDataNonConst (j);
	for (typename ArrayView<const GO>::const_iterator it = localIndices.begin(); 
	     it != localIndices.end(); ++it) {
	  const LO rowIndLocal = map_->getLocalElement (*it);
	  const LO rowIndX = srcMap->getLocalElement (*it);

	  TEUCHOS_TEST_FOR_EXCEPTION(rowIndX == Teuchos::OrdinalTraits<LO>::invalid(), 
            std::invalid_argument, "locallyAssemble(): Input multivector X does "
            "not own the global index " << *it << ".  This probably means that "
            "X was not constructed with the right Map.");
	  // FIXME (mfh 27 Feb 2012) We hard-code the ADD combine mode
	  // for now.  Later, accept other combine modes.
	  X_j[rowIndX] = f (X_j[rowIndX], local_j[rowIndLocal]);
	}

	// Now add the nonlocal data into X_j.
	ArrayView<const GO> nonlocalIndices = nonlocalIndices_[j].view();
	typename ArrayView<const GO>::const_iterator indexIter = nonlocalIndices.begin(); 
	ArrayView<const ST> nonlocalValues = nonlocalValues_[j].view();
	typename ArrayView<const ST>::const_iterator valueIter = nonlocalValues.begin();
	for ( ; indexIter != nonlocalIndices.end() && valueIter != nonlocalValues.end();
	      ++indexIter, ++valueIter) {
	  const LO rowIndX = srcMap->getLocalElement (*indexIter);
	  X_j[rowIndX] = f (X_j[rowIndX], *valueIter);
	}
      }
    }

    //! \c locallyAssemble() for the usual ADD combine mode.
    void 
    locallyAssemble (MV& X)
    {
      std::plus<double> f;
      locallyAssemble<std::plus<scalar_type> > (X, f);
    }

    /// \brief Clear the contents of the vector.
    /// 
    /// This fills the vector with zeros, and also removes nonlocal data.
    void clear() {
      Teuchos::Array<Teuchos::Array<global_ordinal_type> > newNonlocalIndices;
      Teuchos::Array<Teuchos::Array<scalar_type> > newNonlocalValues;
      // The standard STL idiom for clearing the contents of a vector
      // completely.  Setting the size to zero may not necessarily
      // deallocate data.
      std::swap (nonlocalIndices_, newNonlocalIndices);
      std::swap (nonlocalValues_, newNonlocalValues);

      // Don't actually deallocate the multivector of local entries.
      // Just fill it with zero.  This is because the caller hasn't
      // reset the number of columns.
      if (! localVec_.is_null()) {
	localVec_.putScalar (Teuchos::ScalarTraits<scalar_type>::zero());
      }
    }

  private:
    Teuchos::RCP<const map_type> map_;
    size_t numCols_;
    Teuchos::RCP<MV> localVec_;
    Teuchos::Array<Teuchos::Array<global_ordinal_type> > nonlocalIndices_;
    Teuchos::Array<Teuchos::Array<scalar_type> > nonlocalValues_;

    size_t getNumColumns() const { return numCols_; }
  };



  template<class MV>
  class MultiVectorFiller {
  public:
    typedef typename MV::scalar_type scalar_type;
    typedef typename MV::local_ordinal_type local_ordinal_type;
    typedef typename MV::global_ordinal_type global_ordinal_type;
    typedef typename MV::node_type node_type;
    typedef Tpetra::Map<local_ordinal_type, global_ordinal_type, node_type> map_type;

    /// \brief Constructor
    ///
    /// \param map [in] A Map with the same communicator and Kokkos
    ///   Node as the output multivector of \c globalAssemble().
    ///
    /// \param numCols [in] Expected number of columns in the output
    ///   multivector of \c globalAssemble().
    ///
    /// \note If the number of columns given here is not the same as
    ///   the number of columns in the output multivector, the
    ///   two-argument version of \c sumIntoGlobalValues() won't do
    ///   the right thing.  Use the three-argument version of that
    ///   method if you don't know how many columns there will be in
    ///   advance.
    ///
    /// \note Not providing the output multivector in the constructor
    ///   gives \c globalAssemble() more flexibility.  For example,
    ///   its output multivector may have any distribution with the
    ///   same global number of rows.  Furthermore, decoupling
    ///   insertion from the output multivector lets this object store
    ///   preassembled data in whatever format it likes.  It doesn't
    ///   have to insert directly into the output multivector until
    ///   assembly time.  (This may improve performance by amortizing
    ///   the cost of global->local index conversions, for example.)
    MultiVectorFiller (const Teuchos::RCP<const map_type>& map, 
		       const size_t numCols);

    /// \brief Assemble all the data (local and nonlocal) into X_out.
    ///
    /// You can call this method multiple times with different
    /// multivector output arguments.  If those arguments have the
    /// same Map, this method will attempt to reuse the Export object
    /// each time.  It will only reuse if the new target Map is the
    /// same as (in the sense of \c isSameAs()) the previous target
    /// Map, unless you force reuse with the second argument (that
    /// saves a few global reductions for the check).
    ///
    /// \param X_out [in/out] MultiVector with a nonoverlapping Map.
    ///   That Map need not be the same as the one with which this
    ///   object was created.
    ///
    /// \param forceReuseMap [in] Assume that X_out has the same Map
    ///   (in the sense of \c isSameAs()) as the target Map in the
    ///   previous call to this method.  If this method was not called
    ///   before, then assume that X_out has the same Map as the
    ///   argument to the constructor of this object.
    void globalAssemble (MV& X_out, const bool forceReuseMap = false);

    /// \brief Sum data into the multivector.
    ///
    /// \param rows [in] Array of global rows for which to insert
    ///   values.  Must have the same length as values.
    ///
    /// \param column [in] Index of the column in which to insert.
    ///
    /// \param values [in] Array of values to insert.  Must have the
    ///   same length as rows.  rows[i] is the row in which values[i]
    ///   is to be inserted.
    void
    sumIntoGlobalValues (Teuchos::ArrayView<const global_ordinal_type> rows, 
			 size_t column,
			 Teuchos::ArrayView<const scalar_type> values)
    {
      data_->sumIntoGlobalValues (rows, column, values);
    }

    /// \brief Sum data into the multivector.
    ///
    /// In rows and values, the data for each column are stored
    /// contiguously.  Thus, each array is really a matrix in rowwise
    /// order.  (The multivector may use columnwise order internally.)
    ///
    /// Be sure that the number of columns is set correctly before
    /// calling this.
    ///
    /// \param rows [in] Array of global rows for which to insert
    ///   values.  Must have the same length as values.
    ///
    /// \param values [in] Array of values to insert.  Must have the
    ///   same length as rows.  rows[i] is the row in which values[i]
    ///   is to be inserted.
    void
    sumIntoGlobalValues (Teuchos::ArrayView<const global_ordinal_type> rows, 
			 Teuchos::ArrayView<const scalar_type> values)
    {
      data_->sumIntoGlobalValues (rows, values);
    }

  private:
    //! Map with which this object was created ("ctor" == "constructor").
    Teuchos::RCP<const map_type> ctorMap_;

    /// \brief Source Map of the last call to \c globalAssemble().
    ///
    /// A possibly overlapping Map that describes the distribution of
    /// input data, and is the source of the Export.
    Teuchos::RCP<const map_type> sourceMap_;

    /// \brief The target Map of the last call to \c globalAssemble().
    ///
    /// A nonoverlapping Map which is the target of the Export, and
    /// describes the distribution of the output multivector of \c
    /// globalAssemble().
    Teuchos::RCP<const map_type> targetMap_;

    /// \brief The source MV of the Export operation in \c globalAssemble().
    /// 
    /// The \c globalAssemble() method uses this multivector as the
    /// source of the Export operation.  The Export redistributes the
    /// data from a possibly overlapping distribution (reflecting how
    /// elements were inserted) to a nonoverlapping distribution (that
    /// of the output multivector \c X_out of \c globalAssemble()).
    /// This is the simplest way to implement redistribution
    ///
    /// We avoid resizing and reallocating \c sourceVec_ by using a
    /// contiguous subview as the source of the Export, if \c X_out
    /// has fewer columns than \c sourceVec_.
    Teuchos::RCP<MV> sourceVec_;

    /// \brief Storage for inserted indices and values.
    ///
    /// We separate this to facilitate experimentation with different
    /// storage formats.
    MultiVectorFillerData<MV> data_;

    typedef Tpetra::Export<local_ordinal_type, global_ordinal_type, node_type> export_type;
    Teuchos::RCP<export_type> exporter_;

    /// \brief Assemble the local data into \c X_in.
    ///
    /// This method is called by \c globalAssemble(), in which \c X_in
    /// is the multivector with the (possibly overlapping) source
    /// distribution.
    void locallyAssemble (MV& X_in) {
      data_->locallyAssemble (X_in);
    }

    /// \brief Compute the Map from the source indices, if necessary.
    ///
    /// indexBase, comm, and node are the same as the input arguments
    /// of the Map constructors.
    ///
    /// \param initialMap [in] A previously computed Map, if
    ///   available.  If you pass in \c Teuchos::null, this method
    ///   will always compute the Map.  If forceReuseMap is true and
    ///   initialMap is not null, this method will not compute the
    ///   Map, even if that results in an incorrect Map (that does not
    ///   match the source indices).
    ///
    /// \param sourceIndices [in] Global indices, some of which may
    ///   not be owned by sourceOwnedMap.  We assume that they are
    ///   sorted and unique (this is fair because the order of indices
    ///   is irrelevant in finite element assembly, once the local
    ///   assembly is done).
    ///
    /// \param forceReuseMap [in] If true and if initialMap is not
    ///   null, this method will not compute the Map, even if that
    ///   results in an incorrect Map (that does not match the source
    ///   indices).
    ///
    /// \note This is a collective operation.
    ///
    /// \return (the nonlocal Map, whether we had to compute the Map).
    static std::pair<Teuchos::RCP<const map_type>, bool> 
    computeMap (const global_ordinal_type indexBase,
		const Teuchos::RCP<const Teuchos::Comm<int> >& comm, 
		const Teuchos::RCP<node_type>& node,
		const Teuchos::RCP<const map_type>& initialMap, 
		Teuchos::ArrayView<const global_ordinal_type> sourceIndices,
		const bool forceReuseMap);

    /// \brief The inserted row indices, possibly unsorted and with duplicates.
    ///
    /// If \c indices is the view returned by this method, and \c
    /// values is the view returned by \c getSourceValues(), then
    /// indices[i] is the row index of values[i].
    ///
    /// \note This method is "conceptually const."  It may do caching,
    ///   so we didn't declare this const, to avoid "mutable" data.
    Teuchos::ArrayView<const global_ordinal_type> 
    getSourceIndices() 
    {
      return data_->getSourceIndices();
    }

    /// \brief The inserted values.
    ///  
    /// If \c indices is the view returned by \c getSourceIndices(),
    /// and \c values is the view returned by this method, then
    /// indices[i] is the row index of values[i].
    ///
    /// \note This method is "conceptually const."  It may do caching,
    ///   so we didn't declare this const, to avoid "mutable" data.
    Teuchos::ArrayView<const scalar_type> 
    getSourceValues() 
    {
      return data_->getSourceValues();
    }

    /// \brief Get a copy of all the inserted row indices, sorted and made unique.
    ///
    /// This method is suited for assembling the right-hand sides of
    /// the element systems in the finite element method.  In that
    /// case, the order of insertions doesn't matter, and duplicate
    /// insertions get merged addditively.
    ///
    /// \note This method is "conceptually const."  It may do caching,
    ///   so we didn't declare this const, to avoid "mutable" data.
    Teuchos::Array<global_ordinal_type> getSortedUniqueSourceIndices();
  };

  template<class MV>
  MultiVectorFiller<MV>::MultiVectorFiller (const Teuchos::RCP<const typename MultiVectorFiller<MV>::map_type>& map, const size_t numCols) 
    : ctorMap_ (map), sourceMap_ (Teuchos::null), 
      targetMap_ (Teuchos::null), data_ (numCols), 
      exporter_ (Teuchos::null)
  {}

  template<class MV>
  Teuchos::Array<typename MultiVectorFiller<MV>::global_ordinal_type>
  MultiVectorFiller<MV>::getSortedUniqueSourceIndices() 
  {
    using Teuchos::Array;
    using Teuchos::ArrayView;
    typedef global_ordinal_type GO;

    // We've chosen for now to decouple sourceIndices from the
    // actual representation of data.  This is why
    // getSourceIndices() returns an Array rather than an ArrayView.
    // Later, this method may which to cache the Array internally,
    // to avoid copying out of the internal representation each
    // time.  The current getSourceIndices() interface requires
    // copying sourceIndices twice (possibly), so that it can return
    // an Array instead of an ArrayView.
    ArrayView<const GO> sourceIndices = getSourceIndices();
    Array<GO> src (sourceIndices.size());
    std::copy (sourceIndices.begin(), sourceIndices.end(), src.begin());
    //if (! std::is_sorted (src.begin(), src.end())) {
    std::sort (src.begin(), src.end());
    //}
    typename Array<GO>::const_iterator it = std::unique (src.begin(), src.end());
    const size_t newSize = Teuchos::as<size_t> (it - sourceIndices.begin());
    src.resize (newSize);
    return src;
  }


  template<class MV>
  std::pair<Teuchos::RCP<const typename MultiVectorFiller<MV>::map_type>, bool> 
  MultiVectorFiller<MV>::
  computeMap (const global_ordinal_type indexBase,
	      const Teuchos::RCP<const Teuchos::Comm<int> >& comm, 
	      const Teuchos::RCP<node_type>& node,
	      const Teuchos::RCP<const map_type>& initialMap, 
	      Teuchos::ArrayView<const global_ordinal_type> sourceIndices,
	      const bool forceReuseMap)
  {
    using Teuchos::Array;
    using Teuchos::ArrayView;
    using Teuchos::RCP;
    using Teuchos::rcp;
    using Tpetra::global_size_t;
    typedef global_ordinal_type GO;

    // Passing "invalid" for the numGlobalElements argument of the
    // Map constructor tells the Map to compute the global number of
    // elements itself.
    const global_size_t invalid = Teuchos::OrdinalTraits<global_size_t>::invalid();
#if 0
    //
    // Split source indices into owned and nonowned.
    //
    ArrayView<const GO> ownedIndices, nonownedIndices;
    {
      std::pair<Array<GO>, Array<GO> > result = 
	splitSourceIndices (sourceIndices, ownedIndices);
      // Swapping avoids array copies.
      std::swap (ownedIndices, result.first);
      std::swap (nonownedIndices, result.second);
    }
#endif // 0
    RCP<const map_type> map;
    bool computedMap = false;
    if (forceReuseMap && ! initialMap.is_null()) {
      map = initialMap;
      computedMap = false;
    }
    else {
      // Future optimization: Check whether the input Map has the same
      // global indices and index base on all processes.  (We assume
      // that the input Map has the same communicator and Kokkos
      // Node.)  If so, we can reuse the input Map.  This check would
      // look very much like Map::isSameAs(), except that we only have
      // one Map here.
      map = rcp (new map_type (invalid, sourceIndices, indexBase, comm, node));
      computedMap = true;
    }
    return std::make_pair (map, computedMap);
  }

  template<class MV>
  void 
  MultiVectorFiller<MV>::globalAssemble (MV& X_out, const bool forceReuseMap)
  {
    using Teuchos::ArrayView;
    using Teuchos::Array;
    using Teuchos::Range1D;
    using Teuchos::RCP;
    using Teuchos::rcp;
    typedef global_ordinal_type GO;

    const size_t numVecs = X_out.getNumVectors();

    if (numVecs == 0) {
      // Nothing to do!  Of course, this does not check for whether
      // X_out has the right number of rows.  That's OK, though.
      // Compare to the definition of the BLAS' _AXPY for an input
      // vector containing NaNs, but multiplied by alpha=0.
      return;
    }
    //
    // Get the target Map of the Export.  If X_out's Map is the same
    // as the target Map from last time, then we may be able to
    // recycle the Export from last time, if the source Map is also
    // the same.
    //
    RCP<const map_type> targetMap;
    bool assumeSameTargetMap = false;
    if (targetMap_.is_null()) {
      targetMap_ = X_out.getMap();
      targetMap = targetMap_;
      assumeSameTargetMap = false;
    }
    else {
      if (forceReuseMap) {
	targetMap = targetMap_;
	assumeSameTargetMap = true;
      }
      else  {
	// If X_out's Map is the same as targetMap_, we may be able to
	// reuse the Export.  Constructing the Export may be more
	// expensive than calling isSameAs() (which requires just a
	// few reductions and reading through the lists of owned
	// global indices), so it's worth checking.
	if (targetMap_.isSameAs (X_out.getMap())) {
	  assumeSameTargetMap = true;
	  targetMap = targetMap_;
	}
      }
    }
    //
    // Get the source Map of the Export.  If the source Map of the
    // Export is the same as last time, then we may be able to recycle
    // the Export from last time, if the target Map is also the same.
    //
    RCP<const map_type> sourceMap;
    bool computedSourceMap = false;
    {
      if (forceReuseMap && ! sourceMap_.is_null()) {
	sourceMap = sourceMap_;
      }
      else {
	// This reads all of the index data and forms a sorted unique
	// Array as a copy.  That makes the operation expensive, which
	// is why we check whether we want to and can reuse the Map
	// before asking for the list.  We want a deep copy because
	// the original order of the indices matters for local
	// assembly.
	Array<const GO> sourceIndices = getSortedUniqueSourceIndices();
	std::pair<RCP<const map_type>, bool> result = 
	  computeMap (ctorMap_.getIndexBase(), ctorMap_.getComm(), 
		      ctorMap_.getNode(), sourceMap_, sourceIndices, 
		      forceReuseMap);
	sourceMap = result.first;
	computedSourceMap = result.second;
      }
    }
    if (computedSourceMap) {
      sourceMap_ = sourceMap;
    }
    // 
    // Now that we have the source and target Maps of the Export, we
    // can check whether we can recycle the Export from last time.
    //
    const bool mustComputeExport = 
      (exporter_.is_null() || (assumeSameTargetMap && ! computedSourceMap));
    if (mustComputeExport) {
      exporter_ = rcp (new export_type (sourceMap_, targetMap_));
    }

    // Source multivector for the Export.
    RCP<MV> X_in;
    const bool mustReallocateVec = sourceVec_.is_null() || 
      sourceVec_->getNumVectors() < numVecs || ! assumeSameTargetMap;

    if (mustReallocateVec) {
      // Reallocating nonlocalVec_ ensures that it has the right Map.
      sourceVec_ = rcp (new MV (X_out, numVecs));
      X_in = sourceVec_;
    } else {
      if (sourceVec_->getNumVectors() == numVecs) {
	X_in = sourceVec_;
      } else { // sourceVec_ has more vectors than needed.
	X_in = sourceVec_->subView (Range1D (0, numVecs-1));
      }
    }

    // "Locally assemble" the data into X_in by summing together
    // entries with the same indices.
    locallyAssemble (X_in);
    
    // Do the Export.
    const Tpetra::CombineMode combineMode = Tpetra::ADD;
    X_in->doExport (X_out, exporter_, combineMode);
  }

  /// \class MultiVectorFillerTester
  /// \brief Tests for \c MultiVectorFiller
  /// \author Mark Hoemmen
  ///
  /// \tparam MV A specialization of \c Tpetra::MultiVector.
  template<class MV>
  class MultiVectorFillerTester {
  public:
    typedef typename MV::scalar_type scalar_type;
    typedef typename MV::local_ordinal_type local_ordinal_type;
    typedef typename MV::global_ordinal_type global_ordinal_type;
    typedef typename MV::node_type node_type;
    typedef Tpetra::Map<local_ordinal_type, global_ordinal_type, node_type> map_type;

    /// \brief Test global assembly when constructor Map = target Map.
    ///
    /// Constructor Map = target Map is a common case for finite
    /// element assembly.  This method current only tests the version
    /// of \c sumIntoGlobalValues() that works on one column at a
    /// time.
    ///
    /// If any test fails, this method throws an exception.
    static void 
    testSameMap (const Teuchos::RCP<const map_type>& targetMap, 
		 const global_ordinal_type eltSize, // Must be odd
		 const size_t numCols)
    {
      using Teuchos::Array;
      using Teuchos::ArrayRCP;
      using Teuchos::ArrayView;
      using Teuchos::as;
      using Teuchos::Comm;
      using Teuchos::ptr;
      using Teuchos::RCP;
      using Teuchos::rcp;
      using Teuchos::REDUCE_SUM;
      using Teuchos::reduceAll;
      using std::cerr;
      using std::endl;

      typedef local_ordinal_type LO;
      typedef global_ordinal_type GO;
      typedef scalar_type ST;
      typedef Teuchos::ScalarTraits<ST> STS;

      TEUCHOS_TEST_FOR_EXCEPTION(eltSize % 2 == 0, std::invalid_argument,
        "Element size (eltSize) argument must be odd.");
      TEUCHOS_TEST_FOR_EXCEPTION(numCols == 0, std::invalid_argument,
        "Number of columns (numCols) argument must be nonzero.");

      //RCP<MV> X = rcp (new MV (targetMap, numCols));

      Array<GO> rows (eltSize);
      Array<ST> values (eltSize);
      std::fill (values.begin(), values.end(), STS::one());

      // Make this a pointer so we can free its contents, in case
      // those contents depend on the input to globalAssemble().
      RCP<MultiVectorFiller<MV> > filler = 
	rcp (new MultiVectorFiller<MV> (targetMap, numCols));

      const GO minGlobalIndex = targetMap->getMinGlobalIndex();
      const GO maxGlobalIndex = targetMap->getMaxGlobalIndex();
      const GO minAllGlobalIndex = targetMap->getMinAllGlobalIndex();
      const GO maxAllGlobalIndex = targetMap->getMaxAllGlobalIndex();
      for (size_t j = 0; j < numCols; ++j) {
	for (GO i = minGlobalIndex; i <= maxGlobalIndex; ++i) {
	  // Overlap over processes, without running out of bounds.
	  const GO start = std::max (i - eltSize/2, minAllGlobalIndex);
	  const GO end = std::min (i + eltSize/2, maxAllGlobalIndex);
	  const GO len = end - start + 1;

	  ArrayView<GO> rowsView = rows.view (start, len);
	  for (GO k = 0; k < len; ++k) {
	    rowsView[k] = start + k;
	  }
	  ArrayView<const ST> valuesView = values.view (start, len);
	  // Convert rowsView to const view.
	  filler->sumIntoGlobalValues (rowsView.view(), j, valuesView);
	}
      }

      MV X_out (targetMap, numCols);
      filler->globalAssemble (X_out);
      filler = Teuchos::null;

      const GO indexBase = targetMap->getIndexBase();
      Array<GO> errorLocations;
      for (size_t j = 0; j < numCols; ++j) {
	ArrayRCP<const ST> X_j = X_out.getData (j);

	// Each entry of the column should have the value eltSize,
	// except for the first and last few entries in the whole
	// column (globally, not locally).
	for (GO i = minGlobalIndex; i <= maxGlobalIndex; ++i) {
	  const LO localIndex = targetMap->getLocalElement (i);
	  TEUCHOS_TEST_FOR_EXCEPTION(i == Teuchos::OrdinalTraits<LO>::invalid(),
            std::logic_error, "Global index " << i << " is not in the multi"
            "vector's Map.");

	  if (i <= minAllGlobalIndex + eltSize/2) {
	    if (X_j[localIndex] != STS::one() + as<ST>(i) - as<ST>(indexBase)) {
	      errorLocations.push_back (i);
	    }
	  } 
	  else if (i >= maxAllGlobalIndex - eltSize/2) {
	    if (X_j[localIndex] != STS::one() + as<ST>(maxAllGlobalIndex) - as<ST>(i)) {
	      errorLocations.push_back (i);
	    }
	  }
	  else {
	    if (X_j[localIndex] != as<ST>(eltSize)) {
	      errorLocations.push_back (i);
	    }
	  }
	} // for each global index which my process owns

	const typename Array<GO>::size_type localNumErrors = errorLocations.size();
	typename Array<GO>::size_type globalNumErrors = 0;
	RCP<Comm<int> > comm = targetMap->getComm();
	reduceAll (*comm, REDUCE_SUM, localNumErrors, ptr (&globalNumErrors));

	if (globalNumErrors != 0) {
	  std::ostringstream os;
	  os << "Proc " << comm->getRank() << ": " << localNumErrors 
	     << " incorrect value" << (localNumErrors != 1 ? "s" : "") 
	     << ".  Error locations: [ ";
	  std::copy (errorLocations.begin(), errorLocations.end(),
		     std::ostream_iterator<GO> (os, " "));
	  os << "].";
	  // Iterate through all processes in the communicator,
	  // printing out each process' local errors.
	  for (int p = 0; p < comm->getSize(); ++p) {
	    if (p = comm->getRank()) {
	      cerr << os.str() << endl;
	    }
	    // Barriers to let output finish.
	    comm->barrier();
	    comm->barrier();
	    comm->barrier();
	  }
	  TEUCHOS_TEST_FOR_EXCEPTION(globalNumErrors != 0, std::logic_error, 
            "Over all procs: " << globalNumErrors << " total error" 
            << (globalNumErrors != 1 ? "s" : "") << ".");
	} // if there were any errors in column j
      } // for each column j
    }
  };

} // namespace (anonymous)



/*********************************************************/
/*                     Typedefs                          */
/*********************************************************/
typedef Sacado::Fad::SFad<double,3>      Fad3; //# ind. vars fixed at 3
typedef Intrepid::FunctionSpaceTools     IntrepidFSTools;
typedef Intrepid::RealSpaceTools<double> IntrepidRSTools;
typedef Intrepid::CellTools<double>      IntrepidCTools;

/**********************************************************************************/
/***************** FUNCTION DECLARATION FOR ML PRECONDITIONER *********************/
/**********************************************************************************/

/** \brief Solve the linear system using AztecOO (CG) and ML.

  \return The number of iterations that the solve took.

  \param  ProblemType        [in]    problem type
  \param  MLList             [in]    ML parameter list 
  \param  A                  [in]    discrete operator matrix
  \param  xexact             [in]    exact solution
  \param  b                  [in]    right-hand-side vector
  \param  uh                 [out]   solution vector
  \param  TotalErrorResidual [out]   error residual
  \param  TotalErrorExactSol [out]   error in uh
 */
int TestMultiLevelPreconditionerLaplace(char ProblemType[],
				 Teuchos::ParameterList   & MLList,
                                 Epetra_CrsMatrix   & A,
                                 const Epetra_MultiVector & xexact,
                                 Epetra_MultiVector & b,
                                 Epetra_MultiVector & uh,
                                 double & TotalErrorResidual,
                                 double & TotalErrorExactSol);



/**********************************************************************************/
/******** FUNCTION DECLARATIONS FOR EXACT SOLUTION AND SOURCE TERMS ***************/
/**********************************************************************************/

/** \brief  User-defined exact solution.

    \param  x           [in]    x-coordinate of the evaluation point
    \param  y           [in]    y-coordinate of the evaluation point
    \param  z           [in]    z-coordinate of the evaluation point

    \return Value of the exact solution at (x,y,z)
 */
template<typename Scalar>
Scalar exactSolution(const Scalar& x, const Scalar& y, const Scalar& z);

/** \brief  User-defined material tensor.

    \param  material    [out]   3 x 3 material tensor evaluated at (x,y,z)
    \param  x           [in]    x-coordinate of the evaluation point
    \param  y           [in]    y-coordinate of the evaluation point
    \param  z           [in]    z-coordinate of the evaluation point

    \warning Symmetric and positive definite tensor is required for every (x,y,z).
*/
template<typename Scalar>
void materialTensor(Scalar material[][3], const Scalar&  x, const Scalar&  y, const Scalar&  z);

/** \brief  Computes gradient of the exact solution. Requires user-defined exact solution.

    \param  gradExact  [out]   gradient of the exact solution evaluated at (x,y,z)
    \param  x          [in]    x-coordinate of the evaluation point
    \param  y          [in]    y-coordinate of the evaluation point
    \param  z          [in]    z-coordinate of the evaluation point
 */
template<typename Scalar>
void exactSolutionGrad(Scalar gradExact[3], const Scalar& x, const Scalar& y, const Scalar& z);


/** \brief Computes source term: f = -div(A.grad u).  Requires user-defined exact solution
           and material tensor.

    \param  x          [in]    x-coordinate of the evaluation point
    \param  y          [in]    y-coordinate of the evaluation point
    \param  z          [in]    z-coordinate of the evaluation point

    \return Source term corresponding to the user-defined exact solution evaluated at (x,y,z)
 */
template<typename Scalar>
Scalar sourceTerm(Scalar& x, Scalar& y, Scalar& z);


/** \brief Computation of the material tensor at array of points in physical space.

    \param worksetMaterialValues      [out]     Rank-2, 3 or 4 array with dimensions (C,P), (C,P,D) or (C,P,D,D)
                                                with the values of the material tensor
    \param evaluationPoints           [in]      Rank-3 (C,P,D) array with the evaluation points in physical frame
*/
template<class ArrayOut, class ArrayIn>
void evaluateMaterialTensor(ArrayOut &        worksetMaterialValues,
                            const ArrayIn &   evaluationPoints);


/** \brief Computation of the source term at array of points in physical space.

    \param sourceTermValues           [out]     Rank-2 (C,P) array with the values of the source term
    \param evaluationPoints           [in]      Rank-3 (C,P,D) array with the evaluation points in physical frame
*/
template<class ArrayOut, class ArrayIn>
void evaluateSourceTerm(ArrayOut &       sourceTermValues,
                        const ArrayIn &  evaluationPoints);

/** \brief Evaluate the exact solution at the given array of points in
  physical space.

  \param exactSolutionValues [out] Rank-2 (C,P) array with the values
    of the exact solution.

  \param evaluationPoints [in] Rank-3 (C,P,D) array with the
    evaluation points in physical frame.
*/
template<class ArrayOut, class ArrayIn>
void evaluateExactSolution(ArrayOut &       exactSolutionValues,
                           const ArrayIn &  evaluationPoints);


/** Evaluate the gradient of the exact solution at the given array of
    points in physical space.

  \param exactSolutionGradValues [out] Rank-3 (C,P,D) array with the
    values of the gradient of the exact solution.

  \param evaluationPoints [in] Rank-3 (C,P,D) array with the
    evaluation points in physical frame.
*/
template<class ArrayOut, class ArrayIn>
void evaluateExactSolutionGrad(ArrayOut &       exactSolutionGradValues,
                               const ArrayIn &  evaluationPoints);


/**********************************************************************************/
/******************************** MAIN ********************************************/
/**********************************************************************************/

int 
main(int argc, char *argv[]) 
{
  using Teuchos::ptr;
  using std::cerr;
  using std::cout;
  using std::endl;

  int error = 0;

#ifdef HAVE_MPI
  Teuchos::GlobalMPISession mpiSession(&argc, &argv,0);
  const int rank = mpiSession.getRank();
  const int numProcs = mpiSession.getNProc();
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
  const int rank = 0;
  const int numProcs = 1;
#endif

  const int MyPID = Comm.MyPID();
  Epetra_Time Time(Comm);

   //Check number of arguments
  if (argc > 3) {
      cerr << "\n>>> ERROR: Invalid number of arguments.\n\n";
      cerr << "Usage:\n\n";
      cerr << "  ./TrilinosCouplings_examples_scaling_Example_Poisson_Tpetra.exe [meshfile.xml] [solver.xml]\n\n";
      cerr << "   meshfile.xml(optional) - xml file with description of Pamgen mesh\n\n";
      cerr << "   solver.xml(optional) - xml file with ML solver options\n\n";
      return EXIT_FAILURE;
   }

 if (MyPID == 0){
  cout \
    << "===============================================================================\n" \
    << "|                                                                             |\n" \
    << "|          Example: Solve Poisson Equation on Hexahedral Mesh, with Tpetra    |\n" \
    << "|                                                                             |\n" \
    << "|  Questions? Contact  Pavel Bochev  (pbboche@sandia.gov),                    |\n" \
    << "|                      Denis Ridzal  (dridzal@sandia.gov),                    |\n" \
    << "|                      Kara Peterson (kjpeter@sandia.gov).                    |\n" \
    << "|                                                                             |\n" \
    << "|  Intrepid's website:   http://trilinos.sandia.gov/packages/intrepid         |\n" \
    << "|  Pamgen's website:     http://trilinos.sandia.gov/packages/pamgen           |\n" \
    << "|  ML's website:         http://trilinos.sandia.gov/packages/ml               |\n" \
    << "|  Isorropia's website:  http://trilinos.sandia.gov/packages/isorropia        |\n" \
    << "|  Trilinos website:     http://trilinos.sandia.gov                           |\n" \
    << "|                                                                             |\n" \
    << "===============================================================================\n";
  }

  if (MyPID == 0) {
#ifdef HAVE_MPI
    cout << "PARALLEL executable \n"; 
#else
    cout << "SERIAL executable \n";
#endif
  }

  /**********************************************************************************/
  /********************************** GET XML INPUTS ********************************/
  /**********************************************************************************/

  // Command line for xml file, otherwise use default
  std::string xmlMeshInFileName, xmlSolverInFileName;
  if (argc >= 2) {
    xmlMeshInFileName = argv[1];
  }
  else {
    xmlMeshInFileName = "Poisson.xml";
  }
  if (argc >= 3) {
    xmlSolverInFileName = argv[2];
  }

  // Read xml file into parameter list
  Teuchos::ParameterList inputMeshList;
  Teuchos::ParameterList inputSolverList;

  if (xmlMeshInFileName.length()) {
    if (MyPID == 0) {
      cout << "\nReading parameter list from the XML file \""<<xmlMeshInFileName<<"\" ...\n\n";
    }
    Teuchos::updateParametersFromXmlFile (xmlMeshInFileName, ptr (&inputMeshList));
    if (MyPID == 0) {
      inputMeshList.print(cout,2,true,true);
      cout << "\n";
    }
  }
  else {
    cerr << "Cannot read input file: " << xmlMeshInFileName << "\n";
    return EXIT_FAILURE;
  }

  if (xmlSolverInFileName.length()) {
    if (MyPID == 0) {
      cout << "\nReading parameter list from the XML file \""<<xmlSolverInFileName<<"\" ...\n\n";
    }
    Teuchos::updateParametersFromXmlFile (xmlSolverInFileName, ptr (&inputSolverList));
  } 
  else if (MyPID == 0) {
    cout << "Using default solver values ..." << endl;
  }

  // Get pamgen mesh definition
  std::string meshInput = Teuchos::getParameter<std::string> (inputMeshList, "meshInput");

  // Get Isorropia and Zoltan parameters.
  Teuchos::ParameterList iso_paramlist = inputMeshList.sublist ("Isorropia Input");
  if (MyPID == 0) {
    cout << "Isorropia/Zoltan parameters" << endl;
    iso_paramlist.print (cout, 2, true, true);
  }


  /**********************************************************************************/
  /***************************** GET CELL TOPOLOGY **********************************/
  /**********************************************************************************/

  // Get cell topology for base hexahedron
  shards::CellTopology cellType(shards::getCellTopologyData<shards::Hexahedron<8> >() );

  // Get dimensions 
  int numNodesPerElem = cellType.getNodeCount();
  int spaceDim = cellType.getDimension();
  int dim = 3;

  /**********************************************************************************/
  /******************************* GENERATE MESH ************************************/
  /**********************************************************************************/

  if (MyPID == 0) {
    cout << "Generating mesh ... \n\n";
  }

  long long *  node_comm_proc_ids   = NULL;
  long long *  node_cmap_node_cnts  = NULL;
  long long *  node_cmap_ids        = NULL;
  long long ** comm_node_ids        = NULL;
  long long ** comm_node_proc_ids   = NULL;

  // Generate mesh with Pamgen
  long long maxInt = 9223372036854775807LL;
  Create_Pamgen_Mesh (meshInput.c_str(), dim, rank, numProcs, maxInt);

  std::string msg ("Poisson: ");
  if (MyPID == 0) {
    cout << msg << "Pamgen Setup     = " << Time.ElapsedTime() << endl; 
    Time.ResetStartTime();
  }
    
   // Get mesh size info
  char title[100];
  long long numDim;
  long long numNodes;
  long long numElems;
  long long numElemBlk;
  long long numNodeSets;
  long long numSideSets;
  int id = 0;

  im_ex_get_init_l (id, title, &numDim, &numNodes, &numElems, 
		    &numElemBlk, &numNodeSets, &numSideSets);

  long long numNodesGlobal;
  long long numElemsGlobal;
  long long numElemBlkGlobal;
  long long numNodeSetsGlobal;
  long long numSideSetsGlobal;

  im_ne_get_init_global_l (id, &numNodesGlobal, &numElemsGlobal, 
			   &numElemBlkGlobal, &numNodeSetsGlobal,
			   &numSideSetsGlobal);

  // Print mesh information
  if (MyPID == 0){
    cout << " Number of Global Elements: " << numElemsGlobal << " \n";
    cout << "    Number of Global Nodes: " << numNodesGlobal << " \n\n";
  }

  long long * block_ids = new long long [numElemBlk];
  error += im_ex_get_elem_blk_ids_l (id, block_ids);

  long long  *nodes_per_element   = new long long [numElemBlk];
  long long  *element_attributes  = new long long [numElemBlk];
  long long  *elements            = new long long [numElemBlk];
  char      **element_types       = new char * [numElemBlk];
  long long **elmt_node_linkage   = new long long * [numElemBlk];

  for (long long i = 0; i < numElemBlk; ++i) {
    element_types[i] = new char [MAX_STR_LENGTH + 1];
    error += im_ex_get_elem_block_l(id, 
				    block_ids[i], 
				    element_types[i],
				    (long long*)&(elements[i]),
				    (long long*)&(nodes_per_element[i]), 
				    (long long*)&(element_attributes[i]));
  }

  // connectivity
  for (long long b = 0; b < numElemBlk; b++) {
    elmt_node_linkage[b] =  new long long [nodes_per_element[b]*elements[b]];
    error += im_ex_get_elem_conn_l(id,block_ids[b],elmt_node_linkage[b]);
  }

  // Get node-element connectivity
  int telct = 0;
  FieldContainer<int> elemToNode(numElems,numNodesPerElem);
  for(long long b = 0; b < numElemBlk; b++){
    for(long long el = 0; el < elements[b]; el++){
      for (int j=0; j<numNodesPerElem; j++) {
	elemToNode(telct,j) = elmt_node_linkage[b][el*numNodesPerElem + j]-1;
      }
      ++telct;
    }
  }
 
   // Read node coordinates and place in field container
  FieldContainer<double> nodeCoord(numNodes,dim);
  double * nodeCoordx = new double [numNodes];
  double * nodeCoordy = new double [numNodes];
  double * nodeCoordz = new double [numNodes];
  im_ex_get_coord_l (id, nodeCoordx, nodeCoordy, nodeCoordz);
  for (int i=0; i<numNodes; i++) {          
    nodeCoord(i,0)=nodeCoordx[i];
    nodeCoord(i,1)=nodeCoordy[i];
    nodeCoord(i,2)=nodeCoordz[i];
  }
  delete [] nodeCoordx;
  delete [] nodeCoordy;
  delete [] nodeCoordz;

  // parallel info
  long long num_internal_nodes;
  long long num_border_nodes;
  long long num_external_nodes;
  long long num_internal_elems;
  long long num_border_elems;
  long long num_node_comm_maps;
  long long num_elem_comm_maps;
  im_ne_get_loadbal_param_l( id, 
			     &num_internal_nodes,
			     &num_border_nodes, 
			     &num_external_nodes,
			     &num_internal_elems, 
			     &num_border_elems,
			     &num_node_comm_maps,
			     &num_elem_comm_maps,
			     0/*unused*/ );

  if (num_node_comm_maps > 0) {
    node_comm_proc_ids   = new long long  [num_node_comm_maps];
    node_cmap_node_cnts  = new long long  [num_node_comm_maps];
    node_cmap_ids        = new long long  [num_node_comm_maps];
    comm_node_ids        = new long long* [num_node_comm_maps];
    comm_node_proc_ids   = new long long* [num_node_comm_maps];
  
    long long *  elem_cmap_ids        = new long long [num_elem_comm_maps];
    long long *  elem_cmap_elem_cnts  = new long long [num_elem_comm_maps];


    if (im_ne_get_cmap_params_l (id, 
				 node_cmap_ids,
				 (long long*)node_cmap_node_cnts, 
				 elem_cmap_ids,
				 (long long*)elem_cmap_elem_cnts, 
				 0/*not used proc_id*/ ) < 0 ) {
      ++error;
    }
      
    for (long long j = 0; j < num_node_comm_maps; ++j) {
      comm_node_ids[j]       = new long long [node_cmap_node_cnts[j]];
      comm_node_proc_ids[j]  = new long long [node_cmap_node_cnts[j]];
      if (im_ne_get_node_cmap_l (id, 
				 node_cmap_ids[j], 
				 comm_node_ids[j], 
				 comm_node_proc_ids[j],
				 0/*not used proc_id*/ ) < 0 ) {
	++error;
      }
      node_comm_proc_ids[j] = comm_node_proc_ids[j][0];
    }

    delete [] elem_cmap_ids;
    delete [] elem_cmap_elem_cnts;      
  }

  if (!Comm.MyPID()) {
    cout << msg << "Mesh Queries     = " << Time.ElapsedTime() << endl; 
    Time.ResetStartTime();
  }

  // Calculate global node ids
  long long * globalNodeIds = new long long[numNodes];
  bool * nodeIsOwned = new bool[numNodes];

  calc_global_node_ids (globalNodeIds,
			nodeIsOwned,
			numNodes,
			num_node_comm_maps,
			node_cmap_node_cnts,
			node_comm_proc_ids,
			comm_node_ids,
			rank);    

  if (MyPID == 0) {
    cout << msg << "Global Node Nums = " << Time.ElapsedTime() << endl; 
    Time.ResetStartTime();
  }

  // Container indicating whether a node is on the boundary (1-yes 0-no)
  FieldContainer<int> nodeOnBoundary(numNodes);

  // Get boundary (side set) information
  long long * sideSetIds = new long long [numSideSets];
  long long numSidesInSet;
  long long numDFinSet;
  im_ex_get_side_set_ids_l(id,sideSetIds);
  for (int i=0; i<numSideSets; i++) {
    im_ex_get_side_set_param_l(id,sideSetIds[i],&numSidesInSet,&numDFinSet);
    if (numSidesInSet > 0){
      long long * sideSetElemList = new long long [numSidesInSet];
      long long * sideSetSideList = new long long [numSidesInSet];
      im_ex_get_side_set_l(id,sideSetIds[i],sideSetElemList,sideSetSideList);
      for (int j=0; j<numSidesInSet; j++) {
             
	int sideNode0 = cellType.getNodeMap(2,sideSetSideList[j]-1,0);
	int sideNode1 = cellType.getNodeMap(2,sideSetSideList[j]-1,1);
	int sideNode2 = cellType.getNodeMap(2,sideSetSideList[j]-1,2);
	int sideNode3 = cellType.getNodeMap(2,sideSetSideList[j]-1,3);
             
	nodeOnBoundary(elemToNode(sideSetElemList[j]-1,sideNode0))=1;
	nodeOnBoundary(elemToNode(sideSetElemList[j]-1,sideNode1))=1;
	nodeOnBoundary(elemToNode(sideSetElemList[j]-1,sideNode2))=1;
	nodeOnBoundary(elemToNode(sideSetElemList[j]-1,sideNode3))=1;
      }
      delete [] sideSetElemList;
      delete [] sideSetSideList;
    }
  }
  delete [] sideSetIds;

  if (MyPID == 0) {
    cout << msg << "Boundary Conds   = " << Time.ElapsedTime() << endl; 
    Time.ResetStartTime();
  }
 

  /**********************************************************************************/
  /********************************* GET CUBATURE ***********************************/
  /**********************************************************************************/

  // Get numerical integration points and weights
  DefaultCubatureFactory<double>  cubFactory;                                   
  int cubDegree = 2;
  Teuchos::RCP<Cubature<double> > hexCub = cubFactory.create(cellType, cubDegree); 
  
  int cubDim       = hexCub->getDimension();
  int numCubPoints = hexCub->getNumPoints();
  
  FieldContainer<double> cubPoints(numCubPoints, cubDim);
  FieldContainer<double> cubWeights(numCubPoints);
  
  hexCub->getCubature(cubPoints, cubWeights);
  
  if (MyPID==0) {
    cout << "Getting cubature                            "
	 << Time.ElapsedTime() << " sec \n"  ; 
    Time.ResetStartTime();
  }

  /**********************************************************************************/
  /*********************************** GET BASIS ************************************/
  /**********************************************************************************/

  // Define basis 
  Basis_HGRAD_HEX_C1_FEM<double, FieldContainer<double> > hexHGradBasis;
  int numFieldsG = hexHGradBasis.getCardinality();
  FieldContainer<double> HGBValues(numFieldsG, numCubPoints); 
  FieldContainer<double> HGBGrads(numFieldsG, numCubPoints, spaceDim); 

  // Evaluate basis values and gradients at cubature points
  hexHGradBasis.getValues(HGBValues, cubPoints, OPERATOR_VALUE);
  hexHGradBasis.getValues(HGBGrads, cubPoints, OPERATOR_GRAD);

  if (MyPID == 0) {
    cout << "Getting basis                               "
	 << Time.ElapsedTime() << " sec \n"  ; 
    Time.ResetStartTime();
  }

  /**********************************************************************************/
  /********************* BUILD MAPS FOR GLOBAL SOLUTION *****************************/
  /**********************************************************************************/

  // Count owned nodes
  int ownedNodes=0;
  for (int i = 0; i < numNodes; ++i) {
    if (nodeIsOwned[i]) {
      ownedNodes++;
    }
  }

  // Build a list of the OWNED global ids...
  // NTS: will need to switch back to long long
  int *ownedGIDs = new int[ownedNodes];    
  int oidx = 0;
  for(int i = 0; i < numNodes; ++i) {
    if (nodeIsOwned[i]) {
      ownedGIDs[oidx] = static_cast<int> (globalNodeIds[i]);
      ++oidx;
    }
  }
  // Generate Epetra Map for nodes
  Epetra_Map globalMapG (-1, ownedNodes, ownedGIDs, 0, Comm);

  // Global arrays in Epetra format
  Epetra_FECrsMatrix StiffMatrix (Copy, globalMapG, 20*numFieldsG);
  Epetra_FEVector rhsVector (globalMapG);

  if (MyPID == 0) {
    cout << msg << "Build global maps                           "
	 << Time.ElapsedTime() << " sec \n";  
    Time.ResetStartTime();
  }

#ifdef DUMP_DATA
  /**********************************************************************************/
  /**** PUT COORDINATES AND NODAL VALUES IN ARRAYS FOR OUTPUT (FOR PLOTTING ONLY) ***/
  /**********************************************************************************/

  // Put coordinates in multivector for output
  Epetra_MultiVector nCoord(globalMapG,3);
  Epetra_MultiVector nBound(globalMapG,1);

  int indOwned = 0;
  for (int inode = 0; inode < numNodes; ++inode) {
    if (nodeIsOwned[inode]) {
      nCoord[0][indOwned]=nodeCoord(inode,0);
      nCoord[1][indOwned]=nodeCoord(inode,1);
      nCoord[2][indOwned]=nodeCoord(inode,2);
      nBound[0][indOwned]=nodeOnBoundary(inode);
      ++indOwned;
    }
  }
  EpetraExt::MultiVectorToMatrixMarketFile("coords.dat",nCoord,0,0,false);
  EpetraExt::MultiVectorToMatrixMarketFile("nodeOnBound.dat",nBound,0,0,false);

  // Put element to node mapping in multivector for output
  Epetra_Map globalMapElem (numElemsGlobal, numElems, 0, Comm);
  Epetra_MultiVector elem2node (globalMapElem, numNodesPerElem);
  for (int ielem = 0; ielem < numElems; ++ielem) {
    for (int inode = 0; inode < numNodesPerElem; ++inode) {
      elem2node[inode][ielem] = globalNodeIds[elemToNode(ielem,inode)];
    }
  }
  EpetraExt::MultiVectorToMatrixMarketFile("elem2node.dat",elem2node,0,0,false);

  if (MyPID == 0) {
    Time.ResetStartTime();
  }
#endif // DUMP_DATA

  /**********************************************************************************/
  /************************** DIRICHLET BC SETUP ************************************/
  /**********************************************************************************/

  int numBCNodes = 0;
  for (int inode = 0; inode < numNodes; ++inode) {
    if (nodeOnBoundary(inode) && nodeIsOwned[inode]) {
      numBCNodes++;
    }
  }

  // Vector for use in applying BCs
  Epetra_MultiVector v (globalMapG, true);
  v.PutScalar(0.0);

  // Set v to boundary values on Dirichlet nodes
  int * BCNodes = new int [numBCNodes];
  int indbc=0;
  int iOwned=0;
  for (int inode = 0; inode < numNodes; ++inode) {
    if (nodeIsOwned[inode]) {
      if (nodeOnBoundary(inode)) {
	BCNodes[indbc]=iOwned;
	indbc++;
	double x = nodeCoord(inode, 0);
	double y = nodeCoord(inode, 1);
	double z = nodeCoord(inode, 2);
	v[0][iOwned] = exactSolution(x, y, z);
      }
      ++iOwned;
    }
  }

  if (MyPID == 0) {
    cout << msg << "Get Dirichlet boundary values               "
	 << Time.ElapsedTime() << " sec \n\n"; 
    Time.ResetStartTime();
  }

  /**********************************************************************************/
  /******************** DEFINE WORKSETS AND LOOP OVER THEM **************************/
  /**********************************************************************************/

  // Define desired workset size and count how many worksets there are
  // on this processor's mesh block.
  int desiredWorksetSize = numElems; // change to desired workset size!
  //int desiredWorksetSize = 100;    // change to desired workset size!
  int numWorksets        = numElems / desiredWorksetSize;

  // When numElems is not divisible by desiredWorksetSize, increase
  // workset count by 1.
  if (numWorksets*desiredWorksetSize < numElems) {
    numWorksets += 1;
  }

 if (MyPID == 0) {
    cout << "Building discretization matrix and right hand side... \n\n";
    cout << "\tDesired workset size:                 " << desiredWorksetSize <<"\n";
    cout << "\tNumber of worksets (per processor):   " << numWorksets <<"\n\n";
    Time.ResetStartTime();
  }

  for (int workset = 0; workset < numWorksets; ++workset) {
    // Compute cell numbers where the workset starts and ends
    int worksetSize  = 0;
    int worksetBegin = (workset + 0)*desiredWorksetSize;
    int worksetEnd   = (workset + 1)*desiredWorksetSize;

    // When numElems is not divisible by desiredWorksetSize, the last
    // workset ends at numElems
    worksetEnd = (worksetEnd <= numElems) ? worksetEnd : numElems;

    // Now we know the actual workset size and can allocate the array
    // for the cell nodes.
    worksetSize = worksetEnd - worksetBegin;
    FieldContainer<double> cellWorkset (worksetSize, numNodesPerElem, spaceDim);

    // Copy coordinates into cell workset.
    int cellCounter = 0;
    for (int cell = worksetBegin; cell < worksetEnd; ++cell) {
      for (int node = 0; node < numNodesPerElem; ++node) {
        cellWorkset(cellCounter, node, 0) = nodeCoord(elemToNode(cell, node), 0);
        cellWorkset(cellCounter, node, 1) = nodeCoord(elemToNode(cell, node), 1);
        cellWorkset(cellCounter, node, 2) = nodeCoord(elemToNode(cell, node), 2);
      }
      ++cellCounter;
    }

    /**********************************************************************************/
    /*                                Allocate arrays                                 */
    /**********************************************************************************/

    // Containers for Jacobians, integration measure & cubature points in workset cells
    FieldContainer<double> worksetJacobian  (worksetSize, numCubPoints, spaceDim, spaceDim);
    FieldContainer<double> worksetJacobInv  (worksetSize, numCubPoints, spaceDim, spaceDim);
    FieldContainer<double> worksetJacobDet  (worksetSize, numCubPoints);
    FieldContainer<double> worksetCubWeights(worksetSize, numCubPoints);
    FieldContainer<double> worksetCubPoints (worksetSize, numCubPoints, cubDim);

    // Containers for basis values transformed to workset cells and them multiplied by cubature weights
    FieldContainer<double> worksetHGBValues        (worksetSize, numFieldsG, numCubPoints);
    FieldContainer<double> worksetHGBValuesWeighted(worksetSize, numFieldsG, numCubPoints);
    FieldContainer<double> worksetHGBGrads         (worksetSize, numFieldsG, numCubPoints, spaceDim);
    FieldContainer<double> worksetHGBGradsWeighted (worksetSize, numFieldsG, numCubPoints, spaceDim);

    // Containers for diffusive & advective fluxes & non-conservative adv. term and reactive terms
    FieldContainer<double> worksetDiffusiveFlux(worksetSize, numFieldsG, numCubPoints, spaceDim);

    // Containers for material values and source term. Require user-defined functions
    FieldContainer<double> worksetMaterialVals (worksetSize, numCubPoints, spaceDim, spaceDim);
    FieldContainer<double> worksetSourceTerm   (worksetSize, numCubPoints);

    // Containers for workset contributions to the discretization matrix and the right hand side
    FieldContainer<double> worksetStiffMatrix (worksetSize, numFieldsG, numFieldsG);
    FieldContainer<double> worksetRHS         (worksetSize, numFieldsG);

    if (MyPID == 0) {
      cout << msg << "Allocate arrays                             "
	   << Time.ElapsedTime() << " sec \n"; 
      Time.ResetStartTime();
    }

    /**********************************************************************************/
    /*                                Calculate Jacobians                             */
    /**********************************************************************************/

    IntrepidCTools::setJacobian(worksetJacobian, cubPoints, cellWorkset, cellType);
    IntrepidCTools::setJacobianInv(worksetJacobInv, worksetJacobian );
    IntrepidCTools::setJacobianDet(worksetJacobDet, worksetJacobian );

    if (MyPID == 0) {
      cout << msg << "Calculate Jacobians                         "
	   << Time.ElapsedTime() << " sec \n"; 
      Time.ResetStartTime();
    }

    /**********************************************************************************/
    /*          Cubature Points to Physical Frame and Compute Data                    */
    /**********************************************************************************/

    // map cubature points to physical frame
    IntrepidCTools::mapToPhysicalFrame (worksetCubPoints, cubPoints, cellWorkset, cellType);

    // get A at cubature points
    evaluateMaterialTensor (worksetMaterialVals, worksetCubPoints);

    // get source term at cubature points
    evaluateSourceTerm (worksetSourceTerm, worksetCubPoints);

    if (MyPID == 0) {
      cout << msg << "Map to physical frame and get source term   "
	   << Time.ElapsedTime() << " sec \n"; 
      Time.ResetStartTime();
    }

    /**********************************************************************************/
    /*                         Compute Stiffness Matrix                               */
    /**********************************************************************************/

    // Transform basis gradients to physical frame:
    IntrepidFSTools::HGRADtransformGRAD<double>(worksetHGBGrads,                // DF^{-T}(grad u)
                                                worksetJacobInv,   HGBGrads);

    // Compute integration measure for workset cells:
    IntrepidFSTools::computeCellMeasure<double>(worksetCubWeights,              // Det(DF)*w = J*w
                                                worksetJacobDet, cubWeights);


    // Multiply transformed (workset) gradients with weighted measure
    IntrepidFSTools::multiplyMeasure<double>(worksetHGBGradsWeighted,           // DF^{-T}(grad u)*J*w
                                             worksetCubWeights, worksetHGBGrads);


    // Compute the diffusive flux:
    IntrepidFSTools::tensorMultiplyDataField<double>(worksetDiffusiveFlux,      //  A*(DF^{-T}(grad u)
                                                     worksetMaterialVals,
                                                     worksetHGBGrads);

    // Integrate to compute workset diffusion contribution to global matrix:
    IntrepidFSTools::integrate<double>(worksetStiffMatrix,                      // (DF^{-T}(grad u)*J*w)*(A*DF^{-T}(grad u))
                                       worksetHGBGradsWeighted,
                                       worksetDiffusiveFlux, COMP_BLAS);

    if(MyPID==0) {
      cout << msg << "Compute stiffness matrix                    "
	   << Time.ElapsedTime() << " sec \n"; 
      Time.ResetStartTime();}

    /**********************************************************************************/
    /*                                   Compute RHS                                  */
    /**********************************************************************************/

    // Transform basis values to physical frame:
    IntrepidFSTools::HGRADtransformVALUE<double>(worksetHGBValues,              // clones basis values (u)
                                                 HGBValues);

    // Multiply transformed (workset) values with weighted measure
    IntrepidFSTools::multiplyMeasure<double>(worksetHGBValuesWeighted,          // (u)*J*w
                                             worksetCubWeights, worksetHGBValues);

    // Integrate worksetSourceTerm against weighted basis function set
    IntrepidFSTools::integrate<double>(worksetRHS,                             // f.(u)*J*w
                                       worksetSourceTerm,
                                       worksetHGBValuesWeighted,  COMP_BLAS);

    if (MyPID == 0) {
      cout << msg << "Compute right-hand side                     "
	   << Time.ElapsedTime() << " sec \n"; 
      Time.ResetStartTime();
    }

    /**********************************************************************************/
    /*                         Assemble into Global Matrix                            */
    /**********************************************************************************/

    //"WORKSET CELL" loop: local cell ordinal is relative to numElems
    for (int cell = worksetBegin; cell < worksetEnd; ++cell) {

      // Compute cell ordinal relative to the current workset
      int worksetCellOrdinal = cell - worksetBegin;

      // "CELL EQUATION" loop for the workset cell: cellRow is
      // relative to the cell DoF numbering.
      for (int cellRow = 0; cellRow < numFieldsG; ++cellRow) {

        int localRow  = elemToNode(cell, cellRow);
        int globalRow = globalNodeIds[localRow];
        double sourceTermContribution =  worksetRHS(worksetCellOrdinal, cellRow);

        rhsVector.SumIntoGlobalValues(1, &globalRow, &sourceTermContribution);

        // "CELL VARIABLE" loop for the workset cell: cellCol is
        // relative to the cell DoF numbering.
        for (int cellCol = 0; cellCol < numFieldsG; cellCol++) {

          int localCol  = elemToNode(cell, cellCol);
          int globalCol = globalNodeIds[localCol];
          double operatorMatrixContribution = worksetStiffMatrix(worksetCellOrdinal, cellRow, cellCol);

          StiffMatrix.InsertGlobalValues(1, &globalRow, 1, &globalCol, &operatorMatrixContribution);

        } // *** cell col loop ***
      } // *** cell row loop ***
    } // *** workset cell loop **
  } // *** workset loop ***

  /**********************************************************************************/
  /********************* ASSEMBLE OVER MULTIPLE PROCESSORS **************************/
  /**********************************************************************************/

  StiffMatrix.GlobalAssemble();
  StiffMatrix.FillComplete();
  rhsVector.GlobalAssemble();

  if (MyPID==0) {
    cout << msg << "Global assembly                             "
	 << Time.ElapsedTime() << " sec \n"; 
    Time.ResetStartTime();
  }

  /**********************************************************************************/
  /******************************* ADJUST MATRIX DUE TO BC **************************/
  /**********************************************************************************/

  // Apply stiffness matrix to v
  Epetra_MultiVector rhsDir(globalMapG,true);
  StiffMatrix.Apply(v,rhsDir);

  // Update right-hand side
  rhsVector.Update(-1.0,rhsDir,1.0);

  // Adjust rhs due to Dirichlet boundary conditions
  iOwned=0;
  for (int inode=0; inode<numNodes; inode++){
    if (nodeIsOwned[inode]){
      if (nodeOnBoundary(inode)){
	rhsVector[0][iOwned]=v[0][iOwned];
      }
      iOwned++;
    }
  }

  // Zero out rows and columns of stiffness matrix corresponding to Dirichlet edges
  //  and add one to diagonal.
  ML_Epetra::Apply_OAZToMatrix(BCNodes, numBCNodes, StiffMatrix);

  delete [] BCNodes;

  if (MyPID == 0) {
    cout << msg << "Adjust global matrix and rhs due to BCs     " << Time.ElapsedTime()
	 << " sec \n"; Time.ResetStartTime();
  }

  /**********************************************************************************/
  /************************** PARTITION MATRIX WITH ISORROPIA ***********************/
  /**********************************************************************************/

#ifdef TC_HAVE_ISORROPIA
  if (MyPID == 0) {
    cout << msg << "Adjust Matrix = " << Time.ElapsedTime() << endl;
    cout << "rows = " << StiffMatrix.NumGlobalRows() << endl ;
    cout << "cols = " << StiffMatrix.NumGlobalCols() << endl ;
    cout << "nnzs = " << StiffMatrix.NumGlobalNonzeros() << endl ;
    Time.ResetStartTime();
  }

  /* ****************** Call Isorropia to partition the matrix *********** */
  Teuchos::RCP <const Epetra_RowMatrix> rowmat = Teuchos::rcpFromRef
    (StiffMatrix) ;

  // Create the partitioner.
  Isorropia::Epetra::Partitioner iso_part (rowmat, iso_paramlist) ;
  Teuchos::RCP<Isorropia::Epetra::Partitioner> partitioner =
    Teuchos::rcpFromRef (iso_part) ;
  if (MyPID == 0) {
    cout << msg << "Partition Time = " << Time.ElapsedTime() << endl;
    Time.ResetStartTime();
  }

  // Create the redistributor
  Isorropia::Epetra::Redistributor redist (partitioner) ;

  Teuchos::RCP<Epetra_CrsMatrix> bal_matrix ;
  Teuchos::RCP<Epetra_MultiVector> bal_rhs ;
  Teuchos::RCP<Epetra_Map> iso_bal_map ;
  try {
    // Redistribute the matrix and the rhs.
    bal_matrix = redist.redistribute(StiffMatrix) ;
    bal_rhs = redist.redistribute(rhsVector) ;
  }
  catch (std::exception& exc) {
    cout << "example_Poisson: Isorropia exception ' "
	 << exc.what() << " ' on proc " << MyPID  << endl ;
    return  1;
  }

  if (MyPID == 0) {
    cout << msg << "Redistribute Time = " << Time.ElapsedTime() << endl;
    Time.ResetStartTime();
  }

  // Get the map from Isorropia for lhs.
  iso_bal_map = iso_part.createNewMap() ;
  globalMapG = *iso_bal_map ;

  if (MyPID == 0) {
    cout << msg << "Isorropia create new map = " << Time.ElapsedTime()
	 << endl;
    Time.ResetStartTime();
  }
#endif // TC_HAVE_ISORROPIA

#ifdef DUMP_DATA
  // Dump matrices to disk
  EpetraExt::RowMatrixToMatlabFile("stiff_matrix.dat",StiffMatrix);
  EpetraExt::MultiVectorToMatrixMarketFile("rhs_vector.dat",rhsVector,0,0,false);
#endif // DUMP_DATA

  /**********************************************************************************/
  /*********************************** SOLVE ****************************************/
  /**********************************************************************************/

  // Run the solver
  Teuchos::ParameterList MLList = inputSolverList;
  ML_Epetra::SetDefaults("SA", MLList, 0, 0, false);
  Epetra_FEVector exactNodalVals(globalMapG);
  Epetra_FEVector femCoefficients(globalMapG);
  double TotalErrorResidual = 0.0;
  double TotalErrorExactSol = 0.0;

  // Get exact solution at nodes
  for (int i = 0; i < numNodes; ++i) {
    if (nodeIsOwned[i]) {
      double x = nodeCoord(i,0);
      double y = nodeCoord(i,1);
      double z = nodeCoord(i,2);
      double exactu = exactSolution(x, y, z);

      int rowindex=globalNodeIds[i];
      exactNodalVals.SumIntoGlobalValues(1, &rowindex, &exactu);
    }
  }
  exactNodalVals.GlobalAssemble();

  char probType[10] = "laplace";
   
#ifdef TC_HAVE_ISORROPIA
  TestMultiLevelPreconditionerLaplace(probType,           MLList, 
				      *bal_matrix,        exactNodalVals,
				      *bal_rhs,           femCoefficients, 
				      TotalErrorResidual, TotalErrorExactSol);
                                        
#else
  TestMultiLevelPreconditionerLaplace(probType,             MLList,
				      StiffMatrix,          exactNodalVals,
				      rhsVector,            femCoefficients,
				      TotalErrorResidual,   TotalErrorExactSol);
#endif // TC_HAVE_ISORROPIA

  /**********************************************************************************/
  /**************************** CALCULATE ERROR *************************************/
  /**********************************************************************************/

  if (MyPID == 0) {
    Time.ResetStartTime();
  }

  double L2err = 0.0;
  double L2errTot = 0.0;
  double H1err = 0.0;
  double H1errTot = 0.0;
  double Linferr = 0.0;
  double LinferrTot = 0.0;

#ifdef HAVE_MPI
  // Import solution onto current processor
  Epetra_Map    solnMap (numNodesGlobal, numNodesGlobal, 0, Comm);
  Epetra_Import solnImporter (solnMap, globalMapG);
  Epetra_Vector uCoeff (solnMap);
  uCoeff.Import (femCoefficients, solnImporter, Insert);
#endif

  // Define desired workset size
  desiredWorksetSize = numElems;
  int numWorksetsErr    = numElems/desiredWorksetSize;

  // When numElems is not divisible by desiredWorksetSize, increase workset count by 1
  if(numWorksetsErr*desiredWorksetSize < numElems) numWorksetsErr += 1;

  // Get cubature points and weights for error calc (may be different from previous)
  Intrepid::DefaultCubatureFactory<double>  cubFactoryErr;
  int cubDegErr = 3;
  Teuchos::RCP<Intrepid::Cubature<double> > cellCubatureErr = cubFactoryErr.create(cellType, cubDegErr);
  int cubDimErr       = cellCubatureErr->getDimension();
  int numCubPointsErr = cellCubatureErr->getNumPoints();
  Intrepid::FieldContainer<double> cubPointsErr(numCubPointsErr, cubDimErr);
  Intrepid::FieldContainer<double> cubWeightsErr(numCubPointsErr);
  cellCubatureErr->getCubature(cubPointsErr, cubWeightsErr);

  // Evaluate basis values and gradients at cubature points
  Intrepid::FieldContainer<double> uhGVals(numFieldsG, numCubPointsErr);
  Intrepid::FieldContainer<double> uhGrads(numFieldsG, numCubPointsErr, spaceDim);
  hexHGradBasis.getValues(uhGVals, cubPointsErr, Intrepid::OPERATOR_VALUE);
  hexHGradBasis.getValues(uhGrads, cubPointsErr, Intrepid::OPERATOR_GRAD);

   // Loop over worksets
  for (int workset = 0; workset < numWorksetsErr; ++workset) {

    // Compute cell numbers where the workset starts and ends
    int worksetSize  = 0;
    int worksetBegin = (workset + 0)*desiredWorksetSize;
    int worksetEnd   = (workset + 1)*desiredWorksetSize;

    // when numElems is not divisible by desiredWorksetSize, the last
    // workset ends at numElems.
    worksetEnd = (worksetEnd <= numElems) ? worksetEnd : numElems;

    // now we know the actual workset size and can allocate the array
    // for the cell nodes.
    worksetSize = worksetEnd - worksetBegin;
    Intrepid::FieldContainer<double> cellWorksetEr(worksetSize, numNodesPerElem, spaceDim);
    Intrepid::FieldContainer<double> worksetApproxSolnCoef(worksetSize, numNodesPerElem);

    // Loop over cells to fill arrays with coordinates and discrete
    // solution coefficients.
    int cellCounter = 0;
    for (int cell = worksetBegin; cell < worksetEnd; ++cell) {

      for (int node = 0; node < numNodesPerElem; node++) {
	cellWorksetEr(cellCounter, node, 0) = nodeCoord( elemToNode(cell, node), 0);
	cellWorksetEr(cellCounter, node, 1) = nodeCoord( elemToNode(cell, node), 1);
	cellWorksetEr(cellCounter, node, 2) = nodeCoord( elemToNode(cell, node), 2);

	int rowIndex  = globalNodeIds[elemToNode(cell, node)];
#ifdef HAVE_MPI
	worksetApproxSolnCoef(cellCounter, node) = uCoeff.Values()[rowIndex];
#else
	worksetApproxSolnCoef(cellCounter, node) = femCoefficients.Values()[rowIndex];
#endif // HAVE_MPI
      }

      ++cellCounter;
    } 

    // Containers for Jacobian
    Intrepid::FieldContainer<double> worksetJacobianE(worksetSize, numCubPointsErr, spaceDim, spaceDim);
    Intrepid::FieldContainer<double> worksetJacobInvE(worksetSize, numCubPointsErr, spaceDim, spaceDim);
    Intrepid::FieldContainer<double> worksetJacobDetE(worksetSize, numCubPointsErr);
    Intrepid::FieldContainer<double> worksetCubWeightsE(worksetSize, numCubPointsErr);

    // Containers for basis values and gradients in physical space
    Intrepid::FieldContainer<double> uhGValsTrans(worksetSize,numFieldsG, numCubPointsErr);
    Intrepid::FieldContainer<double> uhGradsTrans(worksetSize, numFieldsG, numCubPointsErr, spaceDim);

    // compute cell Jacobians, their inverses and their determinants
    IntrepidCTools::setJacobian(worksetJacobianE, cubPointsErr, cellWorksetEr, cellType);
    IntrepidCTools::setJacobianInv(worksetJacobInvE, worksetJacobianE );
    IntrepidCTools::setJacobianDet(worksetJacobDetE, worksetJacobianE );

    // map cubature points to physical frame
    Intrepid::FieldContainer<double> worksetCubPoints(worksetSize, numCubPointsErr, cubDimErr);
    IntrepidCTools::mapToPhysicalFrame(worksetCubPoints, cubPointsErr, cellWorksetEr, cellType);

    // evaluate exact solution and gradient at cubature points
    Intrepid::FieldContainer<double> worksetExactSoln(worksetSize, numCubPointsErr);
    Intrepid::FieldContainer<double> worksetExactSolnGrad(worksetSize, numCubPointsErr, spaceDim);
    evaluateExactSolution(worksetExactSoln, worksetCubPoints);
    evaluateExactSolutionGrad(worksetExactSolnGrad, worksetCubPoints);

    // transform basis values to physical coordinates
    IntrepidFSTools::HGRADtransformVALUE<double>(uhGValsTrans, uhGVals);
    IntrepidFSTools::HGRADtransformGRAD<double>(uhGradsTrans, worksetJacobInvE, uhGrads);

    // compute weighted measure
    IntrepidFSTools::computeCellMeasure<double>(worksetCubWeightsE, worksetJacobDetE, cubWeightsErr);

    // evaluate the approximate solution and gradient at cubature points
    Intrepid::FieldContainer<double> worksetApproxSoln(worksetSize, numCubPointsErr);
    Intrepid::FieldContainer<double> worksetApproxSolnGrad(worksetSize, numCubPointsErr, spaceDim);
    IntrepidFSTools::evaluate<double>(worksetApproxSoln, worksetApproxSolnCoef, uhGValsTrans);
    IntrepidFSTools::evaluate<double>(worksetApproxSolnGrad, worksetApproxSolnCoef, uhGradsTrans);

    // get difference between approximate and exact solutions
    Intrepid::FieldContainer<double> worksetDeltaSoln(worksetSize, numCubPointsErr);
    Intrepid::FieldContainer<double> worksetDeltaSolnGrad(worksetSize, numCubPointsErr, spaceDim);
    IntrepidRSTools::subtract(worksetDeltaSoln, worksetApproxSoln, worksetExactSoln);
    IntrepidRSTools::subtract(worksetDeltaSolnGrad, worksetApproxSolnGrad, worksetExactSolnGrad);

    // take absolute values
    IntrepidRSTools::absval(worksetDeltaSoln);
    IntrepidRSTools::absval(worksetDeltaSolnGrad);
    // apply cubature weights to differences in values and grads for use in integration
    Intrepid::FieldContainer<double> worksetDeltaSolnWeighted(worksetSize, numCubPointsErr);
    Intrepid::FieldContainer<double> worksetDeltaSolnGradWeighted(worksetSize, numCubPointsErr, spaceDim);
    IntrepidFSTools::scalarMultiplyDataData<double>(worksetDeltaSolnWeighted,
						    worksetCubWeightsE, worksetDeltaSoln);
    IntrepidFSTools::scalarMultiplyDataData<double>(worksetDeltaSolnGradWeighted,
						    worksetCubWeightsE, worksetDeltaSolnGrad);

    // integrate to get errors on each element
    Intrepid::FieldContainer<double> worksetL2err(worksetSize);
    Intrepid::FieldContainer<double> worksetH1err(worksetSize);
    IntrepidFSTools::integrate<double>(worksetL2err, worksetDeltaSoln,
				       worksetDeltaSolnWeighted, Intrepid::COMP_BLAS);
    IntrepidFSTools::integrate<double>(worksetH1err, worksetDeltaSolnGrad,
				       worksetDeltaSolnGradWeighted, Intrepid::COMP_BLAS);

    // loop over cells to get errors for total workset
    cellCounter = 0;
    for (int cell = worksetBegin; cell < worksetEnd; ++cell) {
      // loop over cubature points
      for(int nPt = 0; nPt < numCubPointsErr; nPt++){
	Linferr = std::max(Linferr, worksetDeltaSoln(cellCounter,nPt));
      }

      L2err += worksetL2err(cellCounter);
      H1err += worksetH1err(cellCounter);

      cellCounter++;
    } // end cell loop
  } // end loop over worksets

#ifdef HAVE_MPI
  // sum over all processors
  Comm.SumAll(&L2err,&L2errTot,1);
  Comm.SumAll(&H1err,&H1errTot,1);
  Comm.MaxAll(&Linferr,&LinferrTot,1);
#else
  L2errTot = L2err;
  H1errTot = H1err;
  LinferrTot = Linferr;
#endif // HAVE_MPI

  if (MyPID == 0) {
    cout << "\n" << "L2 Error:  " << sqrt(L2errTot) <<"\n";
    cout << "H1 Error:  " << sqrt(H1errTot) <<"\n";
    cout << "LInf Error:  " << LinferrTot <<"\n\n";
  }

  if (MyPID == 0) {
    cout << msg << "Calculate error                             "
	 << Time.ElapsedTime() << " s \n"; 
    Time.ResetStartTime();
  }

  // Cleanup
   for (long long b = 0; b < numElemBlk; ++b) {     
     delete [] elmt_node_linkage[b];
     delete [] element_types[b];
   }
   delete [] block_ids;
   delete [] nodes_per_element;
   delete [] element_attributes;
   delete [] element_types;
   delete [] elmt_node_linkage;
   // delete [] ownedGIDs;
   delete [] elements;
   delete [] globalNodeIds;
   delete [] nodeIsOwned;
   if (num_node_comm_maps > 0) {
     delete [] node_comm_proc_ids;
     delete [] node_cmap_node_cnts;
     delete [] node_cmap_ids;
     for (long long i = 0; i < num_node_comm_maps; ++i) {
       delete [] comm_node_ids[i];
       delete [] comm_node_proc_ids[i];
     }
      
     delete [] comm_node_ids;
     delete [] comm_node_proc_ids;
   }

   // delete mesh
   Delete_Pamgen_Mesh();
   
   return EXIT_SUCCESS;
}

/**********************************************************************************/
/********************************* END MAIN ***************************************/
/**********************************************************************************/

/**********************************************************************************/
/************ USER DEFINED FUNCTIONS FOR EXACT SOLUTION ***************************/
/**********************************************************************************/

template<typename Scalar>
Scalar 
exactSolution (const Scalar& x, const Scalar& y, const Scalar& z) 
{
  // Patch test: tri-linear function is in the FE space and should be recovered
  return 1. + x + y + z + x*y + x*z + y*z + x*y*z;

  // Analytic solution with homogeneous Dirichlet boundary data
  // return sin(M_PI*x)*sin(M_PI*y)*sin(M_PI*z)*exp(x+y+z);

  // Analytic solution with inhomogeneous Dirichlet boundary data
  // return exp(x + y + z)/(1. + x*y + y*z + x*y*z);
}

template<typename Scalar>
void 
materialTensor (Scalar material[][3], const Scalar& x, const Scalar& y, const Scalar& z) 
{
  const Scalar zero = Teuchos::ScalarTraits<Scalar>::zero();
  const Scalar one = Teuchos::ScalarTraits<Scalar>::one();

  material[0][0] = one;
  material[0][1] = zero;
  material[0][2] = zero;

  material[1][0] = zero;
  material[1][1] = one;
  material[1][2] = zero;

  material[2][0] = zero;
  material[2][1] = zero;
  material[2][2] = one;
}

/**********************************************************************************/
/************** AUXILIARY FUNCTIONS FROM EXACT SOLUTION ***************************/
/**********************************************************************************/

/************ Grad of Exact Solution ****************/
template<typename Scalar>
void 
exactSolutionGrad (Scalar gradExact[3], const Scalar& x, const Scalar& y, const Scalar& z) 
{
  // To enable derivatives of the gradient (i.e., 2nd derivatives of
  // the exact solution), we need 2 levels of fad types.
  Sacado::Fad::SFad<Scalar,3> fad_x = x;
  Sacado::Fad::SFad<Scalar,3> fad_y = y;
  Sacado::Fad::SFad<Scalar,3> fad_z = z;
  Sacado::Fad::SFad<Scalar,3> u;

  // Indicate the independent variables
  fad_x.diff(0,3);
  fad_y.diff(1,3);
  fad_z.diff(2,3);

  u = exactSolution(fad_x, fad_y, fad_z);

  gradExact[0] = u.dx(0);
  gradExact[1] = u.dx(1);
  gradExact[2] = u.dx(2);
}

/************ Source Term (RHS) ****************/
template<typename Scalar>
Scalar 
sourceTerm (Scalar& x, Scalar& y, Scalar& z)
{
  Scalar u;
  Scalar grad_u[3];
  Scalar flux[3];
  Scalar material[3][3];
  Scalar f = 0.;

  // Indicate the independent variables
  x.diff(0,3);
  y.diff(1,3);
  z.diff(2,3);

  // Get exact solution and its gradient
  u = exactSolution(x, y, z);
  exactSolutionGrad(grad_u, x, y, z);

  // Get material tensor
  materialTensor<Scalar>(material, x, y, z);

  // Compute total flux = (A.grad u)
  for (int i = 0; i < 3; ++i) {
    // Add diffusive flux
    for (int j = 0; j < 3; ++j) {
      flux[i] += material[i][j]*grad_u[j];
    }
  }

  // Compute source term (right hand side): f = -div(A.grad u)
  f = -(flux[0].dx(0) + flux[1].dx(1) + flux[2].dx(2));

  return f;
}

/**********************************************************************************/
/*************************** EVALUATION METHODS ***********************************/
/**********************************************************************************/

/************ Material Tensor ****************/
template<class ArrayOut, class ArrayIn>
void 
evaluateMaterialTensor (ArrayOut&        matTensorValues,
			const ArrayIn&   evaluationPoints)
{
  int numWorksetCells  = evaluationPoints.dimension(0);
  int numPoints        = evaluationPoints.dimension(1);
  int spaceDim         = evaluationPoints.dimension(2);

  double material[3][3];

  for (int cell = 0; cell < numWorksetCells; ++cell) {
    for (int pt = 0; pt < numPoints; ++pt) {
      double x = evaluationPoints(cell, pt, 0);
      double y = evaluationPoints(cell, pt, 1);
      double z = evaluationPoints(cell, pt, 2);

      materialTensor<double>(material, x, y, z);

      for (int row = 0; row < spaceDim; ++row){
        for (int col = 0; col < spaceDim; ++col){
          matTensorValues(cell, pt, row, col) = material[row][col];
        }
      }
    }
  }
}

/************ Source Term (RHS) ****************/
template<class ArrayOut, class ArrayIn>
void 
evaluateSourceTerm (ArrayOut&       sourceTermValues,
		    const ArrayIn&  evaluationPoints)
{
  int numWorksetCells  = evaluationPoints.dimension(0);
  int numPoints = evaluationPoints.dimension(1);

  for (int cell = 0; cell < numWorksetCells; cell++){
    for (int pt = 0; pt < numPoints; pt++){
      Sacado::Fad::SFad<double,3> x = evaluationPoints(cell, pt, 0);
      Sacado::Fad::SFad<double,3> y = evaluationPoints(cell, pt, 1);
      Sacado::Fad::SFad<double,3> z = evaluationPoints(cell, pt, 2);

      sourceTermValues(cell, pt) = sourceTerm<Sacado::Fad::SFad<double,3> >(x, y, z).val();
    }
  }
}


template<class ArrayOut, class ArrayIn>
void 
evaluateExactSolution (ArrayOut&       exactSolutionValues,
		       const ArrayIn&  evaluationPoints)
{
  int numWorksetCells  = evaluationPoints.dimension(0);
  int numPoints = evaluationPoints.dimension(1);

  for (int cell = 0; cell < numWorksetCells; cell++) {
    for (int pt = 0; pt < numPoints; pt++) {
      double x = evaluationPoints(cell, pt, 0);
      double y = evaluationPoints(cell, pt, 1);
      double z = evaluationPoints(cell, pt, 2);

      exactSolutionValues(cell, pt) = exactSolution<double>(x, y, z);
    }
  }
}


template<class ArrayOut, class ArrayIn>
void 
evaluateExactSolutionGrad (ArrayOut&       exactSolutionGradValues,
			   const ArrayIn&  evaluationPoints)
{
  int numWorksetCells  = evaluationPoints.dimension(0);
  int numPoints = evaluationPoints.dimension(1);
  int spaceDim  = evaluationPoints.dimension(2);

  double gradient[3];

  for (int cell = 0; cell < numWorksetCells; cell++){
    for (int pt = 0; pt < numPoints; pt++){
      double x = evaluationPoints(cell, pt, 0);
      double y = evaluationPoints(cell, pt, 1);
      double z = evaluationPoints(cell, pt, 2);

      exactSolutionGrad<double>(gradient, x, y, z);

      for(int row = 0; row < spaceDim; row++){
        exactSolutionGradValues(cell, pt, row) = gradient[row];
      }
    }
  }
}


int 
TestMultiLevelPreconditionerLaplace (char ProblemType[],
				     Teuchos::ParameterList   & MLList,
				     Epetra_CrsMatrix   & A,
				     const Epetra_MultiVector & xexact,
				     Epetra_MultiVector & b,
				     Epetra_MultiVector & uh,
				     double & TotalErrorResidual,
				     double & TotalErrorExactSol)
{
  Epetra_MultiVector x (xexact); 
  x.PutScalar (0.0); // Initial guess for the iterative solve starts as zero.
  
  Epetra_LinearProblem Problem(&A,&x,&b); 
  Epetra_MultiVector* lhs = Problem.GetLHS();
  Epetra_MultiVector* rhs = Problem.GetRHS();
  
  Epetra_Time Time(A.Comm());
  
  // =================== //
  // call ML and AztecOO //
  // =================== //
  
  AztecOO solver(Problem);  
  ML_Epetra::MultiLevelPreconditioner *MLPrec = new ML_Epetra::MultiLevelPreconditioner(A, MLList, true);
  
  // Tell AztecOO to use this preconditioner, then solve.
  solver.SetPrecOperator(MLPrec);
  solver.SetAztecOption(AZ_solver, AZ_cg);
  solver.SetAztecOption(AZ_output, 1);

  solver.Iterate(200, 1e-10);
  
  delete MLPrec;

  uh = *lhs;
  
  // ==================================================== //
  // compute difference between exact solution and ML one //
  // ==================================================== //  
  double d = 0.0, d_tot = 0.0;  
  for (int i=0 ; i<lhs->Map().NumMyElements(); ++i) {
    d += ((*lhs)[0][i] - xexact[0][i]) * ((*lhs)[0][i] - xexact[0][i]);
  }
  
  A.Comm().SumAll(&d,&d_tot,1);
  
  // ================== //
  // compute ||Ax - b|| //
  // ================== //
  double Norm;
  Epetra_Vector Ax(rhs->Map());
  A.Multiply(false, *lhs, Ax);
  Ax.Update(1.0, *rhs, -1.0);
  Ax.Norm2(&Norm);
  
  string msg = ProblemType;
  
  if (A.Comm().MyPID() == 0) {
    cout << msg << endl << "......Using " << A.Comm().NumProc() << " processes" << endl;
    cout << msg << "......||A x - b||_2 = " << Norm << endl;
    cout << msg << "......||x_exact - x||_2 = " << sqrt(d_tot) << endl;
    cout << msg << "......Total Time = " << Time.ElapsedTime() << endl;
  }
  
  TotalErrorExactSol += sqrt(d_tot);
  TotalErrorResidual += Norm;
  
  return solver.NumIters();
}


