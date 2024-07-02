// @HEADER
// *****************************************************************************
//                 Anasazi: Block Eigensolvers Package
//
// Copyright 2004 NTESS and the Anasazi contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file AnasaziThyraDebugAdapter.hpp
  \brief Declarations of Anasazi multi-vector and operator classes using Thyra_MultiVectorBase and Thyra_LinearOpBase classes
*/

#ifndef ANASAZI_THYRA_DEBUG_ADAPTER_HPP
#define ANASAZI_THYRA_DEBUG_ADAPTER_HPP

#include "AnasaziConfigDefs.hpp"
#include "AnasaziTypes.hpp"
#include "AnasaziMultiVec.hpp"
#include "AnasaziOperator.hpp"

#include "AnasaziThyraAdapter.hpp"
#include <Thyra_DetachedMultiVectorView.hpp>
#include <Thyra_MultiVectorBase.hpp>
#include <Thyra_MultiVectorStdOps.hpp>

#include "Teuchos_Assert.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_TimeMonitor.hpp"

namespace Anasazi {

  //@}

  ///////////////////////////////////////////////////////////////
  //
  //--------template class AnasaziThyraMultiVec-----------------
  //
  ///////////////////////////////////////////////////////////////

  /*!
    \brief Basic adapter class for Anasazi::MultiVec that uses Thyra::MultiVectorBase<ScalarType>.

    \note This adapter is only to be used for debugging purposes.  For production use Anasazi::MultiVecTraits templated
    on Thyra::MultiVectorBase<> and Thyra::LinearOpBase<>.

  */
  template<class ScalarType>
  class ThyraMultiVec : public MultiVec<ScalarType> {
  public:

    typedef MultiVecTraits<ScalarType,Thyra::MultiVectorBase<ScalarType> > MVT;
    typedef Teuchos::ScalarTraits<ScalarType> SCT;
    typedef typename Teuchos::ScalarTraits<ScalarType>::magnitudeType MagnitudeType;
    typedef Teuchos::ScalarTraits<MagnitudeType> MT;

    //! @name Constructors/Destructors
    //@{

    //! Basic ThyraMultiVec constructor (wraps Thyra::MultiVectorBase<> object).
    /*! @param mv [in] a reference-counted pointer to a Thyra::MultiVectorBase<> object.

      \returns Pointer to a ThyraMultiVec object.
    */
    ThyraMultiVec( const Teuchos::RCP<Thyra::MultiVectorBase< ScalarType > > & mv ) :
      _timerCreate(Teuchos::TimeMonitor::getNewTimer("ThyraMultiVec::create")),
      _timerClone(Teuchos::TimeMonitor::getNewTimer("ThyraMultiVec::clone")),
      _timerDestroy(Teuchos::TimeMonitor::getNewTimer("ThyraMultiVec::destroy")),
      _timerMvTimesMatAddMv(Teuchos::TimeMonitor::getNewTimer("ThyraMultiVec::mvtimesmataddmv")),
      _timerMvTransMv(Teuchos::TimeMonitor::getNewTimer("ThyraMultiVec::mvtransmv")),
      _timerMvAddMv(Teuchos::TimeMonitor::getNewTimer("ThyraMultiVec::mvaddmv")),
      _timerMvDot(Teuchos::TimeMonitor::getNewTimer("ThyraMultiVec::mvdot")),
      _timerMvNorm(Teuchos::TimeMonitor::getNewTimer("ThyraMultiVec::mvnorm")),
      _timerMvScale(Teuchos::TimeMonitor::getNewTimer("ThyraMultiVec::mvscale")),
      _timerSetBlock(Teuchos::TimeMonitor::getNewTimer("ThyraMultiVec::setblock")),
      _timerMvInit(Teuchos::TimeMonitor::getNewTimer("ThyraMultiVec::mvinit")),
      _timerMvRandom(Teuchos::TimeMonitor::getNewTimer("ThyraMultiVec::mvrandom"))
    {
      Teuchos::TimeMonitor timer(*_timerCreate);
      Thyra_MV = mv;
    }

    //! Basic ThyraMultiVec constructor (wraps Thyra::MultiVectorBase<> object).
    /*! @param mv [in] a reference-counted pointer to a Thyra::MultiVectorBase<> object.
        @param timers [in] a vector containing timers for this wrapper to use.

      \returns Pointer to a ThyraMultiVec object.
    */
    ThyraMultiVec( const Teuchos::RCP<Thyra::MultiVectorBase< ScalarType > > & mv, std::vector<Teuchos::RCP<Teuchos::Time> >& timers )
    {
      copyTimers( timers );
      Teuchos::TimeMonitor timer(*_timerCreate);
      Thyra_MV = mv;
    }

    //! Copy constructor.
    /*! @param mv [in] a ThyraMultiVec object.

      \returns Pointer to a ThyraMultiVec object, where the underlying Thyra::MultiVectorBase<> object has been deep copied.
    */
    ThyraMultiVec( const ThyraMultiVec<ScalarType> & mv )
    {
      copyTimers( mv.getTimers() );
      Teuchos::TimeMonitor timer(*_timerCreate);
      Thyra_MV = MVT::CloneCopy( *(mv.getRCP()) );
    }

    //! Destructor
    virtual ~ThyraMultiVec() { Teuchos::TimeMonitor timer(*_timerDestroy); }

    //@}

    //! @name Creation methods
    //@{

    /*! \brief Creates a new empty ThyraMultiVec containing \c numvecs columns.

    \returns Pointer to an ThyraMultiVec
    */
    MultiVec<ScalarType> * Clone ( const int numvecs ) const
    {
      Teuchos::TimeMonitor timer(*_timerClone);
      std::vector<Teuchos::RCP<Teuchos::Time> >  myTimers = getTimers();
      return new ThyraMultiVec<ScalarType>( MVT::Clone( *Thyra_MV, numvecs ), myTimers ); }

    /*! \brief Creates a new ThyraMultiVec and copies contents of \c *this into
      the new vector (deep copy).

      \returns Pointer to an ThyraMultiVec
    */
    MultiVec<ScalarType> * CloneCopy () const
    {
      Teuchos::TimeMonitor timer(*_timerClone);
      std::vector<Teuchos::RCP<Teuchos::Time> >  myTimers = getTimers();
      return new ThyraMultiVec<ScalarType>( MVT::CloneCopy( *Thyra_MV ), myTimers );
    }

    /*! \brief Creates a new ThyraMultiVec and copies the selected contents of \c *this
      into the new vector (deep copy).

      The copied vectors from \c *this are indicated by the \c index.size() indices in \c index.

      \returns Pointer to an ThyraMultiVec
    */
    MultiVec<ScalarType> * CloneCopy ( const std::vector<int>& index ) const
    {
      Teuchos::TimeMonitor timer(*_timerClone);
      std::vector<Teuchos::RCP<Teuchos::Time> >  myTimers = getTimers();
      return new ThyraMultiVec<ScalarType>( MVT::CloneCopy( *Thyra_MV, index ), myTimers );
    }

    /*! \brief Creates a new ThyraMultiVec that shares the selected contents of \c *this.

    The index of the \c numvecs vectors shallow copied from \c *this are indicated by the
    indices given in \c index.

    \returns Pointer to an ThyraMultiVec
    */
    MultiVec<ScalarType> * CloneViewNonConst ( const std::vector<int>& index )
    {
      Teuchos::TimeMonitor timer(*_timerClone);
      std::vector<Teuchos::RCP<Teuchos::Time> >  myTimers = getTimers();
      return new ThyraMultiVec<ScalarType>( MVT::CloneViewNonConst( *Thyra_MV, index ), myTimers );
    }

    /*! \brief Creates a new ThyraMultiVec that shares the selected contents of \c *this.

    The index of the \c numvecs vectors shallow copied from \c *this are indicated by the
    indices given in \c index.

    \returns Pointer to an ThyraMultiVec
    */
    const MultiVec<ScalarType> * CloneView ( const std::vector<int>& index ) const
    {
      Teuchos::TimeMonitor timer(*_timerClone);
      std::vector<Teuchos::RCP<Teuchos::Time> >  myTimers = getTimers();
      Teuchos::RCP<Thyra::MultiVectorBase<ScalarType> > nonconst_ptr_to_const_view = Teuchos::rcp_const_cast<Thyra::MultiVectorBase<ScalarType> >( MVT::CloneView(*Thyra_MV,index) );
      const MultiVec<ScalarType> * const_ret = new ThyraMultiVec<ScalarType>( nonconst_ptr_to_const_view, myTimers );
      return const_ret;
    }

    //@}

    //! @name Attribute methods
    //@{

    //! Obtain the vector length of *this.
    int GetNumberVecs () const { return MVT::GetNumberVecs( *Thyra_MV ); }

    //! Obtain the number of vectors in *this.
    ptrdiff_t GetGlobalLength () const { return MVT::GetGlobalLength( *Thyra_MV ); }

    //@}

    //! @name Update methods
    //@{
    /*! \brief Update \c *this with \f$\alpha AB + \beta (*this)\f$.
     */
    void MvTimesMatAddMv ( ScalarType alpha, const MultiVec<ScalarType>& A,
                           const Teuchos::SerialDenseMatrix<int,ScalarType>& B,
                           ScalarType beta )
    {
      Teuchos::TimeMonitor timer(*_timerMvTimesMatAddMv);
      const Anasazi::ThyraMultiVec<ScalarType>* vec_A = dynamic_cast<const Anasazi::ThyraMultiVec<ScalarType>* >(&A);
      MVT::MvTimesMatAddMv( alpha, *(vec_A->getRCP()), B, beta, *Thyra_MV );
    }

    /*! \brief Replace \c *this with \f$\alpha A + \beta B\f$.
     */
    void MvAddMv ( ScalarType alpha, const MultiVec<ScalarType>& A,
                   ScalarType beta, const MultiVec<ScalarType>& B)
    {
      Teuchos::TimeMonitor timer(*_timerMvAddMv);
      const Anasazi::ThyraMultiVec<ScalarType>* vec_A = dynamic_cast<const Anasazi::ThyraMultiVec<ScalarType>* >(&A);
      const Anasazi::ThyraMultiVec<ScalarType>* vec_B = dynamic_cast<const Anasazi::ThyraMultiVec<ScalarType>* >(&B);
      MVT::MvAddMv( alpha, *(vec_A->getRCP()), beta, *(vec_B->getRCP()), *Thyra_MV );
    }

    /*! \brief Compute a dense matrix \c B through the matrix-matrix multiply \f$\alpha A^T(*this)\f$.
    */
    void MvTransMv ( ScalarType alpha, const MultiVec<ScalarType>& A, Teuchos::SerialDenseMatrix<int,ScalarType>& B
#ifdef HAVE_ANASAZI_EXPERIMENTAL
        , ConjType conj = Anasazi::CONJ
#endif
        ) const
    {
      Teuchos::TimeMonitor timer(*_timerMvTransMv);
      const Anasazi::ThyraMultiVec<ScalarType>* vec_A = dynamic_cast<const Anasazi::ThyraMultiVec<ScalarType>* >(&A);
      MVT::MvTransMv( alpha, *(vec_A->getRCP()), *Thyra_MV, B );
    }

    /*! \brief Compute a vector \c b where the components are the individual dot-products, i.e. \f$ b[i] = A[i]^H(this[i])\f$ where \c A[i] is the i-th column of \c A.
    */
    void MvDot ( const MultiVec<ScalarType>& A, std::vector<ScalarType> &b
#ifdef HAVE_ANASAZI_EXPERIMENTAL
        , ConjType conj = Anasazi::CONJ
#endif
        ) const
    {
      Teuchos::TimeMonitor timer(*_timerMvDot);
      const Anasazi::ThyraMultiVec<ScalarType>* vec_A = dynamic_cast<const Anasazi::ThyraMultiVec<ScalarType>* >(&A);
      MVT::MvDot( *Thyra_MV, *(vec_A->getRCP()), b );
    }

    /*! \brief Scale each element of the vectors in \c *this with \c alpha.
     */
    void MvScale ( ScalarType alpha ) { Teuchos::TimeMonitor timer(*_timerMvScale); MVT::MvScale( *Thyra_MV, alpha ); }

    /*! \brief Scale each element of the \c i-th vector in \c *this with \c alpha[i].
     */
    void MvScale ( const std::vector<ScalarType>& alpha ) { Teuchos::TimeMonitor timer(*_timerMvScale); MVT::MvScale( *Thyra_MV, alpha ); }

    //@}
    //! @name Norm method
    //@{

    /*! \brief Compute the 2-norm of each individual vector of \c *this.
      Upon return, \c normvec[i] holds the 2-norm of the \c i-th vector of \c *this
    */
    void MvNorm ( std::vector<typename Teuchos::ScalarTraits<ScalarType>::magnitudeType> &normvec ) const { Teuchos::TimeMonitor timer(*_timerMvNorm); MVT::MvNorm( *Thyra_MV, normvec ); }
    //@}

    //! @name Initialization methods
    //@{
    /*! \brief Copy the vectors in \c A to a set of vectors in \c *this.

    The \c numvecs vectors in \c A are copied to a subset of vectors in \c *this
    indicated by the indices given in \c index.
    */
    void SetBlock ( const MultiVec<ScalarType>& A, const std::vector<int>& index )
    {
      Teuchos::TimeMonitor timer(*_timerSetBlock);
      const Anasazi::ThyraMultiVec<ScalarType>* vec_A = dynamic_cast<const Anasazi::ThyraMultiVec<ScalarType>* >(&A);
      MVT::SetBlock( *(vec_A->getRCP()), index, *Thyra_MV );
    }

    /*! \brief Fill the vectors in \c *this with random numbers.
     */
    void MvRandom() { Teuchos::TimeMonitor timer(*_timerMvRandom); MVT::MvRandom( *Thyra_MV ); }

    /*! \brief Replace each element of the vectors in \c *this with \c alpha.
     */
    void MvInit ( ScalarType alpha ) { Teuchos::TimeMonitor timer(*_timerMvInit);  MVT::MvInit( *Thyra_MV, alpha ); }

    //!{ @name Accessor methods
    //@{
    /*! \brief Return the reference-counted pointer held by this object
     */
    Teuchos::RCP< Thyra::MultiVectorBase<ScalarType> > getRCP() { return Thyra_MV; }

    /*! \brief Return the const reference-counted pointer held by this object
     */
    Teuchos::RCP< const Thyra::MultiVectorBase<ScalarType> > getRCP() const { return Thyra_MV; }

    /*! \brief Return a std::vector<> of timers held by this object
     */
    std::vector<Teuchos::RCP<Teuchos::Time> > getTimers() const {
     std::vector<Teuchos::RCP<Teuchos::Time> > timers;
     timers.push_back( _timerCreate );
     timers.push_back( _timerClone );
     timers.push_back( _timerDestroy );
     timers.push_back( _timerMvTimesMatAddMv );
     timers.push_back( _timerMvTransMv );
     timers.push_back( _timerMvAddMv );
     timers.push_back( _timerMvDot );
     timers.push_back( _timerMvNorm );
     timers.push_back( _timerMvScale );
     timers.push_back( _timerSetBlock );
     timers.push_back( _timerMvInit );
     timers.push_back( _timerMvRandom );

     return timers;
   }

   /*! \brief Copy a std::vector<> of timers into this object
    */
   void copyTimers( const std::vector<Teuchos::RCP<Teuchos::Time> >& timers ) {
     _timerCreate = timers[0];
     _timerClone = timers[1];
     _timerDestroy = timers[2];
     _timerMvTimesMatAddMv = timers[3];
     _timerMvTransMv = timers[4];
     _timerMvAddMv = timers[5];
     _timerMvDot = timers[6];
     _timerMvNorm = timers[7];
     _timerMvScale = timers[8];
     _timerSetBlock = timers[9];
     _timerMvInit = timers[10];
     _timerMvRandom = timers[11];
   }
    //@}

    //@}
    //! @name Print method
    //@{
    /*! \brief Print \c *this ThyraMultiVec.
     */
    void MvPrint( std::ostream& os ) const { MVT::MvPrint( *Thyra_MV, os ); }
    //@}

  private:

    Teuchos::RCP<Thyra::MultiVectorBase<ScalarType> > Thyra_MV;
    Teuchos::RCP<Teuchos::Time> _timerCreate, _timerClone, _timerDestroy;
    Teuchos::RCP<Teuchos::Time> _timerMvTimesMatAddMv, _timerMvTransMv, _timerMvAddMv, _timerMvDot;
    Teuchos::RCP<Teuchos::Time> _timerMvNorm, _timerMvScale, _timerSetBlock, _timerMvInit, _timerMvRandom;
  };
  //-------------------------------------------------------------

  ///////////////////////////////////////////////////////////////
  //
  //--------template class AnasaziThyraOp---------------------
  //
  ///////////////////////////////////////////////////////////////

  /*!
    \brief Basic adapter class for Anasazi::Operator that uses Thyra_Operator.

    \note The Thyra package performs double-precision arithmetic, so the use of Thyra with Anasazi will
    only provide a double-precision eigensolver.
  */
  template<class ScalarType>
  class ThyraOp : public virtual Operator<ScalarType> {
  public:

    typedef OperatorTraits<ScalarType,Thyra::MultiVectorBase<ScalarType>,Thyra::LinearOpBase<ScalarType> > OPT;

    //! @name Constructor/Destructor
    //@{

    //! Basic constructor.  Accepts reference-counted pointer to an Thyra_Operator.
    ThyraOp(const Teuchos::RCP<const Thyra::LinearOpBase<ScalarType> > &Op ) { Thyra_Op = Op; }

    //! Destructor
    ~ThyraOp() {}
    //@}

    //! @name Operator application method
    //@{

    /*! \brief This method takes the Anasazi::MultiVec \c X and
      applies the operator to it resulting in the Anasazi::MultiVec \c Y.
    */
    void Apply ( const MultiVec<ScalarType>& X, MultiVec<ScalarType>& Y ) const
    {
      const Anasazi::ThyraMultiVec<ScalarType>* vec_X = dynamic_cast<const Anasazi::ThyraMultiVec<ScalarType>* >(&X);
      Anasazi::ThyraMultiVec<ScalarType>* vec_Y = dynamic_cast<Anasazi::ThyraMultiVec<ScalarType>* >(&Y);
      OPT::Apply( *Thyra_Op, *(vec_X->getRCP()), *(vec_Y->getRCP()) );
    }

    //@}

  private:
    Teuchos::RCP<const Thyra::LinearOpBase<ScalarType> > Thyra_Op;
  };

} // end of Anasazi namespace

#endif
// end of file ANASAZI_THYRA_DEBUG_ADAPTER_HPP
