/**
  \file   Amesos2_TpetraMultiVecAdapter_decl.hpp
  \author Eric T Bavier <etbavier@sandia.gov>
  \date   Wed May 26 19:49:10 CDT 2010
  
  \brief  Amesos2::MultiVecAdapter specialization for the
          Tpetra::MultiVector class.
*/

#ifndef AMESOS2_TPETRA_MULTIVEC_ADAPTER_DECL_HPP
#define AMESOS2_TPETRA_MULTIVEC_ADAPTER_DECL_HPP

#include <Teuchos_RCP.hpp>
#include <Teuchos_as.hpp>
#include <Tpetra_MultiVector.hpp>

#include "Amesos2_MultiVecAdapter_decl.hpp"

namespace Amesos {


template< typename Scalar,
          typename LocalOrdinal,
          typename GlobalOrdinal,
          class    Node >
class MultiVecAdapter<Tpetra::MultiVector<Scalar,
                                          LocalOrdinal,
                                          GlobalOrdinal,
                                          Node> >
{
public:
  // /* \brief localize the data in this MultiVec
  //  * 
  //  * Initialize \c data_ which is a serial local contiguous copy of the data
  //  * contained in \c mv_
  //  */
  // void localize()
  //   {
  //     mv_->get1dCopy(data_, lda_);
  //     localized_ = true;
  //   }
  
public:
  // public type definitions
  typedef Tpetra::MultiVector<Scalar,
                              LocalOrdinal,
                              GlobalOrdinal,
                              Node>       multivec_type;
  typedef Scalar                          scalar_type;
  typedef LocalOrdinal                    local_ordinal_type;
  typedef GlobalOrdinal                   global_ordinal_type;
  typedef Node                            node_type;
  typedef typename Tpetra::global_size_t  global_size_type;

  static const char* name;

  MultiVecAdapter()
    : mv_(Teuchos::null)
    , l_mv_(Teuchos::null)
//    , l_l_mv_(Teuchos::null)
    { }


  /// Copy constructor
  MultiVecAdapter(const MultiVecAdapter<multivec_type>& adapter)
    : mv_(adapter.mv_)
    , l_mv_(adapter.l_mv_)
    , l_l_mv_(adapter.l_l_mv_)
    { }

  
  // MultiVecAdapter(multivec_type* const m)
  //   : mv_(Teuchos::rcp(m))
  //   , data_(mv_->getNumVectors() * mv_->getGlobalLength())
  //   , lda_(mv_->getGlobalLength())
  //   { }


  MultiVecAdapter(const Teuchos::RCP<multivec_type>& m)
    : mv_(m)
    , l_mv_(Teuchos::null)
//    , l_l_mv_(Teuchos::null)
    { }


  ~MultiVecAdapter()
    {
      // TODO: Should the destructor export changes to the serial
      // version before it destroys itself, or should we leave the
      // user responsible for calling the updateValues method?
    }


  /** \brief Scales values of \c this by a factor of \c alpha
   *
   * \f$ this = alpha * this\f$
   * 
   * \param [in] alpha scalar factor
   */
  MultiVecAdapter& scale(const Scalar alpha)
    {
      mv_->scale(alpha);

      return *this;
    }


  
  /** \brief Updates the values of \c this to \f$this = alpha*this + beta*B\f$
   *
   * \param [in] beta scalar coefficient of \c B
   * \param [in] B additive MultiVector
   * \param [in] alpha scalar coefficient of \c this
   */
  MultiVecAdapter& update(
    const Scalar beta,
    const MultiVecAdapter<multivec_type>& B,
    const Scalar alpha)
    {
      mv_->update(beta, B.mv_, alpha);
    
      return *this;
    }

  
  /// Checks whether this multivector is local to the calling node.
  bool isLocal() const 
    {
      if(getComm()->getSize() == 1){
        return true;
      } // There may be other conditions to check
    }

  
  Teuchos::RCP<Teuchos::Comm<int> >& getComm() const
    {
      return mv_->getMap()->getComm();
    }

  /// Get the length of vectors local to the calling node
  size_t getLocalLength() const
    {
      return Teuchos::as<LocalOrdinal>(mv_->getLocalLength());
    }


  /// Get the number of vectors on this node
  size_t getLocalNumVectors() const
    {
      return mv_->getNumVectors();
    }


  /// Get the length of vectors in the global space
  global_size_type getGlobalLength() const
    {
      return mv_->getGlobalLength();
    }


  /// Get the number of global vectors
  size_t getGlobalNumVectors() const
    {
      return mv_->getNumVectors();
    }


  /// Return the stride between vectors on this node
  size_t getStride() const
    {
      return mv_->getStride();
    }


  /// Return \c true if this MV has constant stride between vectors on this node
  bool isConstantStride() const
    {
      return mv_->isConstantStride();
    }


  /// Const vector access
  Teuchos::RCP<const Tpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> >
  getVector( size_t j ) const
    {
      return mv_->getVector(j);
    }
  
  
  /// Nonconst vector access
  Teuchos::RCP<Tpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> >
  getVectorNonConst(size_t j)
    {
      return mv_->getVectorNonConst(j);
    }
    

  /** \brief Copies the multivector's data into the user-provided vector.
   *  
   *  Each multivector is \lda apart in memory.
   */
  void get1dCopy(const Teuchos::ArrayView<Scalar>& A, size_t lda) const;
  

  /** \brief Extracts a 1 dimensional view of this MultiVector's data
   *
   * Guarantees that the view returned will reside in contiguous storage.
   * 
   * \warning
   * It is recommended to use the \c get1dCopy function, from a
   * data-hiding perspective. Use if you know what you are doing.
   *
   * \param local if \c true , each node will get a view of the vectors it is
   * in possession of.  The default, \c false , will give each calling node a
   * view of the global multivector.
   *
   * \note This function is not declared \c const as it normally would be,
   * since it must modify local copies of the vector data before retrning the
   * result.
   */
  Teuchos::ArrayRCP<Scalar> get1dViewNonConst(bool local = false);


  void get2dCopy(Teuchos::ArrayView<Teuchos::ArrayView<Scalar> > A) const;

  
  Teuchos::ArrayRCP<Teuchos::ArrayRCP<Scalar> > get2dView() const;
  

  /** Exports any local changes to this MultiVector to the global space.
   */
  void globalize();


  /** \brief Export \c newVals into the global MultiVector space.
   * 
   * \note we assume the vectors in newVals have the same leading dimension as
   * those in \c this
   * 
   * \param newVals The values to be exported into the global space.
   */
  template<typename Value_t>
  void globalize(const Teuchos::ArrayView<Value_t>& newVals);


  /// Get a short description of this adapter class
  std::string description() const;
  

  /// Print a description of this adapter to the Fancy Output Stream.
  void describe(Teuchos::FancyOStream& os,
    const Teuchos::EVerbosityLevel verbLevel) const;
  

private:
  /* This function "localizes" the wrapped \c mv_ . It defines the private
   * maps \c o_map_ and \c l_map_ and imports global data into the local node.
   * If \c mv_ is not distributed, this method does nothing.
   *
   * It is intended to set things up properly for calls to \c get1dCopy() and \c get1dView().
   *
   * \sa get1dCopy(), get1dView()
   */
  void localize();

  Teuchos::RCP<
  Tpetra::MultiVector<
    Scalar,LocalOrdinal,GlobalOrdinal,Node
    >
  > mv_;
Teuchos::RCP<
  Tpetra::MultiVector<
    Scalar,LocalOrdinal,GlobalOrdinal,Node
    >
  > l_mv_; ///< local multivector
Teuchos::RCP<Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node > > l_l_mv_;   ///< local-local multivector

/// used for transferring between local and global multivectors
Teuchos::RCP<Tpetra::Import<LocalOrdinal,GlobalOrdinal,Node> > importer_; 

Teuchos::RCP<
  Tpetra::Map<
    LocalOrdinal,GlobalOrdinal,Node
    >
  > l_map_; ///< local map
  Teuchos::RCP<
    const Tpetra::Map<
      LocalOrdinal,GlobalOrdinal,Node
      >
    > o_map_; ///< original map

};                              // end class MultiVector<Tpetra::CrsMatrix>

} // end namespace Amesos


#endif // AMESOS2_TPETRA_MULTIVEC_ADAPTER_DECL_HPP
