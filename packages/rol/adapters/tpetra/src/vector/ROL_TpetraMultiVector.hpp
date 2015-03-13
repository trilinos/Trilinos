// @HEADER
// ************************************************************************
//
//               Rapid Optimization Library (ROL) Package
//                 Copyright (2014) Sandia Corporation
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
// Questions? Contact lead developers:
//              Drew Kouri   (dpkouri@sandia.gov) and
//              Denis Ridzal (dridzal@sandia.gov)
//
// ************************************************************************
// @HEADER

#ifndef ROL_TPETRAMULTIVECTOR_HPP
#define ROL_TPETRAMULTIVECTOR_HPP

/** \class ROL::EpetraMultiVector
    \brief Implements the ROL::Vector interface for a Tpetra_MultiVector.
    \author Created by Greg von Winckel
*/

#include "ROL_Vector.hpp"
#include "Tpetra_MultiVector.hpp"

namespace ROL {

// Template on the Real/Scalar type, Local Ordinal, Global Ordinal, and Node
template <class Real, class LO, class GO, class Node>
class TpetraMultiVector : public Vector<Real> {

    typedef Tpetra::MultiVector<Real,LO,GO,Node> MV;
    typedef Teuchos::RCP<MV> MVP;
    typedef TpetraMultiVector TMV;
    typedef Teuchos::RCP<TMV> TMVP; 
 
    private:
        MVP tpetra_vec_;

    public:
        virtual ~TpetraMultiVector() {}

        TpetraMultiVector(const MVP &tpetra_vec) : tpetra_vec_(tpetra_vec) {}
       
        /** \brief Compute \f$y \leftarrow x + y\f$ where \f$y = \mbox{*this}\f$.
        */
        void plus(const Vector<Real> &x) {
            const TMV &ex = Teuchos::dyn_cast<const TMV>(x);
            tpetra_vec_->update(1.0,*ex.getVector(),1.0); 
        }

        /** \brief Compute \f$y \leftarrow \alpha y\f$ where \f$y = \mbox{*this}\f$.
        */
        void scale( const Real alpha ) {
            tpetra_vec_->scale( alpha ); 
        }

        /** \brief Returns \f$ \langle y,x \rangle \f$ where \f$y = \mbox{*this}\f$.
        */
        Real dot( const Vector<Real> &x ) const {
            Real v; // Need this to make a 1-element ArrayView
            Teuchos::ArrayView<Real> val(&v,1);
            const TMV &ex = Teuchos::dyn_cast<const TMV>(x);
            tpetra_vec_->dot( *ex.getVector(), val );
            return val[0];
        }

        /** \brief Returns \f$ \| y \| \f$ where \f$y = \mbox{*this}\f$.
        */
        Real norm() const {
            Real v; // Need this to make a 1-element ArrayView
            Teuchos::ArrayView<Real> val(&v,1);
            tpetra_vec_->norm2(val);
            return val[0];
        } 

        /** \brief Clone to make a new (uninitialized) vector.
        */
        Teuchos::RCP<Vector<Real> > clone() const {
            using Teuchos::rcp; 
            size_t n = tpetra_vec_->getNumVectors();  
            return rcp( new TMV(rcp( new MV(tpetra_vec_->getMap(),n,false)) ));
        }

        /** \brief Compute \f$y \leftarrow \alpha x + y\f$ where \f$y = \mbox{*this}\f$.
        */
        virtual void axpy( const Real alpha, const Vector<Real> &x ) {
            const TMV &ex = Teuchos::dyn_cast<const TMV>(x);
            tpetra_vec_->update( alpha, *ex.getVector(), 1.0 );
        }

        /**  \brief Set to zero vector.
        */
        virtual void zero() {
            tpetra_vec_->putScalar(0.0);
        }

        /**  \brief Set \f$y \leftarrow x\f$ where \f$y = \mbox{*this}\f$.
        */
        virtual void set( const Vector<Real> &x ) {
            const TMV &ex = Teuchos::dyn_cast<const TMV>(x);
            tpetra_vec_->scale(1.0,*ex.getVector());
        }

        MVP getVector() const {
            return this->tpetra_vec_;
        }

        MVP getVector() {
            return this->tpetra_vec_;
        }

        Teuchos::RCP<Vector<Real> > basis( const int i ) const {
            using Teuchos::rcp; 
            size_t n = tpetra_vec_->getNumVectors();  
             
            MVP e = rcp( new MV(tpetra_vec_->getMap(),n,true) );  

            e->replaceLocalValue(i,0,1.0);

            return rcp(new TMV(e) );  
        }

        int dimension() const { return tpetra_vec_->getGlobalLength(); }

}; // class TpetraMultiVector

} // namespace ROL

#endif
