// @HEADER
// ***********************************************************************
// 
//          Tpetra: Templated Linear Algebra Services Package
//                 Copyright (2008) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ************************************************************************
// @HEADER

#ifndef TPETRA_RTIOP_HPP
#define TPETRA_RTIOP_HPP

#include "Tpetra_RTI.hpp"
#include "Tpetra_RTI_detail.hpp"

namespace Tpetra {

  namespace RTI {

    /// \class KernelOp
    /// \brief Operator wrapping a Kokkos (Classic) kernel using RTI.
    ///
    /// This Tpetra::Operator subclass wraps a Kokkos (Classic) kernel
    /// using the Tpetra Reduction/Transformation Interface (RTI).
    /// The first four template parameters are the same (and in the
    /// same order) as those of Tpetra::Operator.  The fifth template
    /// parameter is the type of the Kokkos (Classic) kernel.
    template <class S, class LO, class GO, class Node, class Kernel>
    class KernelOp : public Tpetra::Operator<S,LO,GO,Node> {
    protected:
      RCP<const Import<LO,GO,Node> > _importer;
      RCP<const Export<LO,GO,Node> > _exporter;
      mutable RCP<MultiVector<S,LO,GO,Node> > _importMV, _exportMV;
      RCP<const Map<LO,GO,Node> > _rangeMap, _domainMap;
      mutable Kernel _kernel;
      
    public:
      //! Constructor.
      KernelOp (Kernel kernel,
		const RCP<const Map<LO,GO,Node> > & domainMap,
		const RCP<const Map<LO,GO,Node> > & rangeMap,
		const RCP<const Import<LO,GO,Node> > & importer,
		const RCP<const Export<LO,GO,Node> > & exporter) 
        : _importer (importer), _exporter (exporter)
        , _rangeMap (rangeMap), _domainMap (domainMap)
        , _kernel (kernel)
      {
	const char tfecfFuncName[] = "KernelOp(kernel,domainMap,rangeMap,importer,exporter)";
	if (_rangeMap == null) {
	  _rangeMap = _domainMap;
	}
	TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
          _domainMap == null || _rangeMap == null, std::runtime_error,
	  ": neither domainMap nor rangeMap may be specified null:\ndomainMap: " 
	  << _domainMap << "\nrangeMap: " << _rangeMap << "\n");
#ifdef HAVE_TPETRA_DEBUG
        TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
          _rangeMap->getNode() != _domainMap->getNode(), std::runtime_error, 
	  ": all specified maps must have the same Node instance.");
	if (_importer != null) {
	  TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC( 
            ! _importer->getSourceMap ()->isSameAs (*_domainMap), 
	    std::runtime_error, 
	    ": domain Map is not consistent with importer.");
	  TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
            _importer->getSourceMap ()->getNode () != _domainMap->getNode (),
	    std::runtime_error, 
	    ": all specified Maps must have the same Node instance.");
	}
	if (_exporter != null) {
	  TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
            ! _exporter->getTargetMap ()->isSameAs (*_rangeMap), 
	    std::runtime_error, ": range Map is not consistent with importer.");
          TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
            _exporter->getTargetMap ()->getNode () != _domainMap->getNode (), 
	    std::runtime_error, 
	    ": all specified Maps must have the same Node instance.");
	}
#endif // HAVE_TPETRA_DEBUG
      }

      RCP<const Map<LO,GO,Node> > getDomainMap () const { return _domainMap; }
      RCP<const Map<LO,GO,Node> > getRangeMap ()  const { return _rangeMap; }

      void
      apply (const MultiVector<S,LO,GO,Node> &X, 
	     MultiVector<S,LO,GO,Node> &Y,
	     Teuchos::ETransp mode = Teuchos::NO_TRANS, 
	     S alpha = Teuchos::ScalarTraits<S>::one (), 
	     S beta = Teuchos::ScalarTraits<S>::zero ()) const
      {
	const size_t numVectors = X.getNumVectors ();
	RCP<MultiVector<S,LO,GO,Node> > mvec_inout;
	RCP<const MultiVector<S,LO,GO,Node> > mvec_in2;

	if (_importer != null) {
	  if (_importMV != null && _importMV->getNumVectors () != numVectors) {
	    _importMV = null;
	  }
	  if (_importMV == null) {
	    _importMV = createMultiVector<S> (_importer->getTargetMap (), numVectors);
	  }
	  _importMV->doImport (X, *_importer, INSERT);
	  mvec_in2 = _importMV;
	}
	else {
	  mvec_in2 = rcpFromRef(X);
	}

	if (_exporter != null) {
	  if (_exportMV != null && _exportMV->getNumVectors () != numVectors) {
	    _exportMV = null;
	  }
	  if (_exportMV == null) {
	    _exportMV = createMultiVector<S> (_exporter->getSourceMap (), numVectors);
	  }
	  mvec_inout = _exportMV;
	}
	else {
	  mvec_inout = rcpFromRef (Y);
	}
	_kernel.setAlphaBeta (alpha, beta);
	//
	for (size_t j=0; j < numVectors; ++j) {
	  RCP<       Vector<S,LO,GO,Node> > vec_inout = mvec_inout->getVectorNonConst(j);
	  RCP< const Vector<S,LO,GO,Node> > vec_in2   = mvec_in2->getVector(j);
	  Tpetra::RTI::detail::binary_transform( *vec_inout, *vec_in2, _kernel );
	}
	// export
	if (_exporter != null) {
	  Y.doExport (*_exportMV, *_exporter, ADD);
	}
      }
    };

    //! Non-member constructor for a Tpetra::RTI::KernelOp object.
    template <class S, class LO, class GO, class Node, class Kernel> 
    RCP<const KernelOp<S,LO,GO,Node,Kernel> >
    kernelOp (Kernel kernel, 
	      const RCP<const Map<LO,GO,Node> > & domainMap,
	      const RCP<const Map<LO,GO,Node> > & rangeMap = null,
	      const RCP<const Import<LO,GO,Node> > & importer = null,
	      const RCP<const Export<LO,GO,Node> > & exporter = null) 
    {
      return Teuchos::rcp (new KernelOp<S,LO,GO,Node,Kernel> (kernel, domainMap, rangeMap, 
							      importer, exporter) );
    }

    //! Tpetra::Operator wrapping a binary functor using the Tpetra Reduction/Transformation Interface
    template <class S, class LO, class GO, class Node, class Op> 
    class BinaryOp : 
      public KernelOp<S,LO,GO,Node,Tpetra::RTI::detail::BinaryFunctorAdapterWithAlphaBeta<Op,S> >
    {
    public:
      BinaryOp (Op op, 
		const RCP<const Map<LO,GO,Node> > & domainMap,
		const RCP<const Map<LO,GO,Node> > & rangeMap,
		const RCP<const Import<LO,GO,Node> > & importer,
		const RCP<const Export<LO,GO,Node> > & exporter) 
        : KernelOp<S,LO,GO,Node,Tpetra::RTI::detail::BinaryFunctorAdapterWithAlphaBeta<Op,S> >( Tpetra::RTI::detail::BinaryFunctorAdapterWithAlphaBeta<Op,S>(op), domainMap, rangeMap, importer, exporter ) {}
    };

    //! Non-member constructor for a Tpetra::RTI::BinaryOp object.
    template <class S, class LO, class GO, class Node, class Op> 
    RCP<const BinaryOp<S,LO,GO,Node,Op> >
    binaryOp (Op op, 
	      const RCP<const Map<LO,GO,Node> > & domainMap,
	      const RCP<const Map<LO,GO,Node> > & rangeMap = null,
	      const RCP<const Import<LO,GO,Node> > & importer = null,
	      const RCP<const Export<LO,GO,Node> > & exporter = null) 
    {
      return Teuchos::rcp (new BinaryOp<S,LO,GO,Node,Op> (op, domainMap, rangeMap, 
							  importer, exporter) );
    }
  } // namespace RTI
} // namespace Tpetra

#endif // TPETRA_RTIOP_HPP
