//@HEADER
// ************************************************************************
// 
//               Tpetra: Templated Linear Algebra Services Package 
//                 Copyright (2008) Sandia Corporation
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
// Lesser General Public License for more detail.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ************************************************************************
//@HEADER

#ifndef TPETRA_RTIOP_HPP
#define TPETRA_RTIOP_HPP

#include "Tpetra_RTI.hpp"
#include "Tpetra_RTI_detail.hpp"

namespace Tpetra {

  namespace RTI {

    //! Tpetra::Operator wrapping a Kokkos kernel using the Tpetra Reduction/Transformation Interface
    template <class S, class LO, class GO, class Node, class Kernel>
    class KernelOp : public Tpetra::Operator<S,LO,GO,Node> 
    {
      protected:
        RCP< const Import<LO,GO,Node> > _importer;
        RCP< const Export<LO,GO,Node> > _exporter;
        mutable RCP< MultiVector<S,LO,GO,Node> >  _importMV,  _exportMV;
        RCP< const Map<LO,GO,Node> > _rangeMap, _domainMap;
        mutable Kernel _kernel;
      public:
        KernelOp(Kernel kernel,
                 const RCP<const Map<LO,GO,Node> > & domainMap,
                 const RCP<const Map<LO,GO,Node> > & rangeMap,
                 const RCP<const Import<LO,GO,Node> > & importer,
                 const RCP<const Export<LO,GO,Node> > & exporter) 
        : _importer(importer), _exporter(exporter)
        , _rangeMap(rangeMap), _domainMap(domainMap)
        , _kernel(kernel)
        {
          std::string tfecfFuncName("KernelOp(kernel,domainMap,rangeMap,importer,exporter)");
          if (_rangeMap == null) _rangeMap = _domainMap;
          TEST_FOR_EXCEPTION_CLASS_FUNC( _domainMap == null || _rangeMap == null, std::runtime_error,
              ":KernelOp(): neither domainMap nor rangeMap may be specified null:\ndomainMap: " << _domainMap << "\nrangeMap: " << _rangeMap << "\n");
#ifdef HAVE_TPETRA_DEBUG
          TEST_FOR_EXCEPTION_CLASS_FUNC( _rangeMap->getNode() != _domainMap->getNode(), std::runtime_error, ": all specified maps must have the same node.");
          if (_importer != null) {
            TEST_FOR_EXCEPTION_CLASS_FUNC( !_importer->getSourceMap()->isSameAs(*_domainMap), std::runtime_error, ": domain map is not consistent with importer.");
            TEST_FOR_EXCEPTION_CLASS_FUNC( _importer->getSourceMap()->getNode() != _domainMap->getNode(), std::runtime_error, ": all specified maps must have the same node.");
          }
          if (_exporter != null) {
            TEST_FOR_EXCEPTION_CLASS_FUNC( !_exporter->getTargetMap()->isSameAs(*_rangeMap), std::runtime_error, ": range map is not consistent with importer.");
            TEST_FOR_EXCEPTION_CLASS_FUNC( _exporter->getTargetMap()->getNode() != _domainMap->getNode(), std::runtime_error, ": all specified maps must have the same node.");
          }
#endif
        }
        //
        const RCP<const Map<LO,GO,Node> > & getDomainMap() const { return _domainMap; }
        //
        const RCP<const Map<LO,GO,Node> > & getRangeMap()  const { return _rangeMap; }
        //
        void apply(const MultiVector<S,LO,GO,Node> &X, 
                         MultiVector<S,LO,GO,Node> &Y,
                   Teuchos::ETransp mode = Teuchos::NO_TRANS, 
                   S alpha = Teuchos::ScalarTraits<S>::one(), 
                   S beta = Teuchos::ScalarTraits<S>::zero()) const
        {
          const size_t numVectors = X.getNumVectors();
          RCP< MultiVector<S,LO,GO,Node> > mvec_inout;
          RCP< const MultiVector<S,LO,GO,Node> > mvec_in2;
          //
          if (_importer != null) {
            if (_importMV != null && _importMV->getNumVectors() != numVectors) _importMV = null;
            if (_importMV == null) _importMV = createMultiVector<S>(_importer->getTargetMap(), numVectors);
            _importMV->doImport( X, *_importer, INSERT );
            mvec_in2 = _importMV;
          }
          else {
            mvec_in2 = rcpFromRef(X);
          }
          //
          if (_exporter != null) {
            if (_exportMV != null && _exportMV->getNumVectors() != numVectors) _exportMV = null;
            if (_exportMV == null) _exportMV = createMultiVector<S>(_exporter->getSourceMap(), numVectors);
            mvec_inout = _exportMV;
          }
          else {
            mvec_inout = rcpFromRef(Y);
          }
          _kernel.setAlphaBeta(alpha,beta);
          //
          for (size_t j=0; j < numVectors; ++j) 
          {
            RCP<       Vector<S,LO,GO,Node> > vec_inout = mvec_inout->getVectorNonConst(j);
            RCP< const Vector<S,LO,GO,Node> > vec_in2   = mvec_in2->getVector(j);
            Tpetra::RTI::detail::binary_transform( *vec_inout, *vec_in2, _kernel );
          }
          // export
          if (_exporter != null) {
            Y.doExport(*_exportMV, *_exporter, ADD);
          }
        }
    };

    //! Non-member constructor for a Tpetra::RTI::KernelOp object.
    template <class S, class LO, class GO, class Node, class Kernel> 
    RCP< const KernelOp<S,LO,GO,Node,Kernel> >
    kernelOp(Kernel kernel, 
          const RCP<const Map<LO,GO,Node> > & domainMap,
          const RCP<const Map<LO,GO,Node> > & rangeMap = null,
          const RCP<const Import<LO,GO,Node> > & importer = null,
          const RCP<const Export<LO,GO,Node> > & exporter = null) 
    {
      return Teuchos::rcp(new KernelOp<S,LO,GO,Node,Kernel>(kernel,domainMap,rangeMap,importer,exporter) );
    }

    //! Tpetra::Operator wrapping a binary functor using the Tpetra Reduction/Transformation Interface
    template <class S, class LO, class GO, class Node, class Op> 
    class BinaryOp : public KernelOp<S,LO,GO,Node,Tpetra::RTI::detail::BinaryFunctorAdapterWithAlphaBeta<Op,S> >
    {
      public:
        BinaryOp(Op op, 
              const RCP<const Map<LO,GO,Node> > & domainMap,
              const RCP<const Map<LO,GO,Node> > & rangeMap,
              const RCP<const Import<LO,GO,Node> > & importer,
              const RCP<const Export<LO,GO,Node> > & exporter) 
        : KernelOp<S,LO,GO,Node,Tpetra::RTI::detail::BinaryFunctorAdapterWithAlphaBeta<Op,S> >( Tpetra::RTI::detail::BinaryFunctorAdapterWithAlphaBeta<Op,S>(op), domainMap, rangeMap, importer, exporter ) {}
    };

    //! Non-member constructor for a Tpetra::RTI::BinaryOp object.
    template <class S, class LO, class GO, class Node, class Op> 
    RCP< const BinaryOp<S,LO,GO,Node,Op> >
    binaryOp(Op op, 
             const RCP<const Map<LO,GO,Node> > & domainMap,
             const RCP<const Map<LO,GO,Node> > & rangeMap = null,
             const RCP<const Import<LO,GO,Node> > & importer = null,
             const RCP<const Export<LO,GO,Node> > & exporter = null) 
    {
      return Teuchos::rcp(new BinaryOp<S,LO,GO,Node,Op>(op,domainMap,rangeMap,importer,exporter) );
    }

  } // end of namespace RTI

} // end of namespace Tpetra

#endif // TPETRA_RTIOP_HPP
