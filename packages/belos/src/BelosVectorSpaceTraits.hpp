//@HEADER
// ************************************************************************
//
//                 Belos: Block Linear Solvers Package
//                  Copyright 2004 Sandia Corporation
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
//@HEADER

#ifndef __Belos_VectorSpaceTraits_hpp
#define __Belos_VectorSpaceTraits_hpp

namespace Belos {

  /// \class UndefinedVectorSpaceTraits
  ///
  /// UndefinedVectorSpaceTraits exists only so that VectorSpaceTraits
  /// has a "default" unspecialized implementation which will result
  /// in compile-time errors.
  template<class VectorSpaceType>
  class UndefinedVectorSpaceTraits {
  public:
    /// \brief Signal compile-time error when attempting to instantiate.
    ///
    /// \note Any attempt to compile this function results in a
    /// compile time error.  This means that the template
    /// specialization of Belos::VectorSpaceTraits does not exist for
    /// type VectorSpaceType, or is not complete.
    static inline void notDefined() { VectorSpaceType::this_type_is_missing_a_specialization(); };

    typedef VectorSpaceType vector_space_type;

  private:
    UndefinedVectorSpaceTraits();
    UndefinedVectorSpaceTraits(const UndefinedVectorSpaceTraits&);
    UndefinedVectorSpaceTraits& operator= (const UndefinedVectorSpaceTraits&);
  };

  /// \class VectorSpaceTraits
  /// \brief Traits class for Belos-supported vector space objects.
  ///
  /// Belos needs to reason about vector spaces, because it is both
  /// possible and allowed in preconditioned iterative solvers for
  /// operators to have different domains and ranges.  For example,
  /// consider left-preconditioned GMRES, which involves solving
  /// \f$M^{-1} A x = M^{-1} b\f$ for operators \f$M^{-1}\f$ and
  /// \f$A\f$.  The solution vector x (and all approximate solution
  /// vectors computed by the iterative method) is in the domain of A,
  /// and the right-hand side b (and all residual vectors computed by
  /// the iterative method) is in the range of A.  The range of A is
  /// the domain of \f$M^{-1}\f$, and the range of \f$M^{-1}\f$ is the
  /// domain of A, but there is no requirement that those two spaces
  /// be identical.  One could even imagine performance reasons why
  /// the two spaces should be different (e.g., cotuning the data
  /// distributions of the two operators).  Split preconditioning may
  /// involve up to three different vector spaces.
  ///
  /// Of course, it's possible to write iterative solvers so that
  /// preconditioned residual vectors and nonpreconditioned residual
  /// vectors don't mix.  However, lazy initialization and caching
  /// vectors for later use are important optimizations for
  /// TsqrOrthoManager (for example), and it's important to know that
  /// the cached vectors are compatible with any new input vectors.
  template<class VectorSpaceType>
  class VectorSpaceTraits {
  public:
    typedef VectorSpaceType vector_space_type;

    /// \brief Are vectors from the two vector spaces compatible?
    ///
    /// If this function returns true, then
    /// - You may combine (e.g., add) vectors from the two spaces
    /// - If first resp. second is the domain of an operator_type,
    ///   you may give a multivector_type in the second resp. first
    ///   vector space to the operator_type object as input
    /// - If first resp. second is the range of an operator_type,
    ///   you may give a multivector_type in the second resp. first
    ///   vector space to the operator_type object as output
    /// Otherwise, you have no guarantee that any of these operations
    /// will succeed.
    static bool 
    compatible (const vector_space_type& first, 
		const vector_space_type& second) {
      // This will result in a compile-time error, on purpose, since
      // this traits class is only meaningful when specialized.
      UndefinedVectorSpaceTraits<VectorSpaceType>::not_defined();
      return false;
    }

    /// \brief Return a persistent view of the given vector space.
    ///
    /// A "persistent view" means that the vector space is guaranteed
    /// to survive the life of the distributed object which was
    /// queried in order to get the vector space.  In the case of
    /// Epetra, this involves making a deep copy, if space.strength()
    /// == Teuchos::RCP_WEAK.
    static Teuchos::RCP<const vector_space_type>
    persistentView (const Teuchos::RCP<const vector_space_type>& space) {
      // This will result in a compile-time error, on purpose, since
      // this traits class is only meaningful when specialized.
      UndefinedVectorSpaceTraits<VectorSpaceType>::not_defined();
      return Teuchos::null;
    }
  };


  /// \class DefaultVectorSpace
  /// \brief A default implementation of VectorSpaceType.
  ///
  /// For Belos::MultiVec and Belos::Operator, we have provided a
  /// trivial definition of the various vector space properties.  This
  /// ensures that for those two types, all vector spaces are
  /// identical, and all vectors belong to the same vector space.
  /// This might change in the future, so don't rely on this behavior,
  /// or on the implementation details of DefaultVectorSpace.
  class DefaultVectorSpace {
  public:
    bool compatible (const DefaultVectorSpace& other) const {
      return index_ == other.index_;
    }

    /// \brief Return the "default" vector space.
    /// 
    /// \warning This method is not reentrant!
    static Teuchos::RCP<const DefaultVectorSpace> 
    getDefaultVectorSpace () 
    {
      static Teuchos::RCP<const DefaultVectorSpace> theSpace;
      if (theSpace.is_null())
	{
	  const int defaultIndex = 0;
	  theSpace = Teuchos::rcp (new DefaultVectorSpace (defaultIndex));
	}
      return theSpace;
    }
      
  private:
    //! Default construction is forbidden.
    DefaultVectorSpace ();

    //! Construction itself is an implementation detail and therefore private.
    DefaultVectorSpace (const int index) : index_ (index) {}

    //! Implementation detail (the index of the vector space).
    int index_;
  };

  //! Full specialization of VectorSpaceTraits for DefaultVectorSpace.
  template<>
  class VectorSpaceTraits<DefaultVectorSpace> {
  public:
    typedef DefaultVectorSpace vector_space_type;

    /// \brief Are vectors from the two vector spaces compatible?
    static bool 
    compatible (const vector_space_type& first, 
		const vector_space_type& second) {
      return first.compatible (second);
    }

    /// \brief Return a persistent view of the given vector space.
    ///
    /// A "persistent view" means that the vector space is guaranteed
    /// to survive the life of the distributed object which was
    /// queried in order to get the vector space.  
    static Teuchos::RCP<const vector_space_type>
    persistentView (const Teuchos::RCP<const vector_space_type>& space) {
      return space;
    }
  };

  
} // namespace Belos

#endif // __Belos_VectorSpaceTraits_hpp
