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

#ifndef TPETRA_FEMULTIVECTOR_DECL_HPP
#define TPETRA_FEMULTIVECTOR_DECL_HPP

/// \file Tpetra_FEMultiVector_decl.hpp
/// \brief Declaration of the Tpetra::MultiVector class
///
/// If you want to use Tpetra::FEMultiVector, include "Tpetra_FEMultiVector.hpp"
/// (a file which CMake generates and installs for you).  If you only want
/// the declaration of Tpetra::FEMultiVector, include this file
/// (Tpetra_FEMultiVector_decl.hpp).
///

#include "Tpetra_MultiVector_decl.hpp"

namespace Tpetra {

#ifndef DOXYGEN_SHOULD_SKIP_THIS
  // forward declaration of MultiVector (declared later in this file)
  template<class S, class LO, class GO, class N> class FEMultiVector;
#endif // DOXYGEN_SHOULD_SKIP_THIS


  template <class Scalar = ::Tpetra::Details::DefaultTypes::scalar_type,
            class LocalOrdinal = ::Tpetra::Details::DefaultTypes::local_ordinal_type,
            class GlobalOrdinal = ::Tpetra::Details::DefaultTypes::global_ordinal_type,
            class Node = ::Tpetra::Details::DefaultTypes::node_type>
  class FEMultiVector :
    public MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>
  {
  private:
    friend class MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
  public:
    //! @name Typedefs to facilitate template metaprogramming.
    //@{

    //! This class' first template parameter; the Scalar type.
    typedef Scalar scalar_type;
    //! This class' second template parameter; the type of local indices.
    typedef LocalOrdinal local_ordinal_type;
    //! This class' third template parameter; the type of global indices.
    typedef GlobalOrdinal global_ordinal_type;
    //! This class' fourth template parameter; the Kokkos Node type.
    typedef Node node_type;

    //! The dual_view_type picked up from MultiVector
    typedef typename MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::dual_view_type dual_view_type;

    //! The type of the Map specialization used by this class.
    typedef typename MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::map_type map_type;

    //! Grab impl_scalar_type from superclass
    typedef typename MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::impl_scalar_type impl_scalar_type;

    //! Grab dot_type from superclass
    typedef typename MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::dot_type dot_type;

    //! Grab mag_type from superclass
    typedef typename MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::mag_type mag_type;

    //! Grab device_type from superclass
    typedef typename MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::device_type device_type;

    //! Grab execution_space from superclass
    typedef typename MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::execution_space execution_space;

    //! The type of the base class of this class.
    typedef MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> base_type;

    //@}
    //! @name Constructors and destructor
    //@{
    /// \brief Basic constructor.
    /// \param map [in] Map describing the distribution of rows of the
    ///   resulting MultiVector.  If the importer is not null, this must be the same as importer->getSourceMap()
    /// \param importer [in] Import describing the distribution of rows of the the multivector in two separate modes.
    ///   In the case of finite element assembly, importer->getSourceMap() should correspond to the uniquely owned entries in the domain map
    ///   of the finite element problem.  importer->getTargetMap() should correspond to the overlapped entries needed for assembly.  The 
    ///   cannonical way to get this importer is to take the Import object associated with the CrsGraph for your finite element matrix.
    /// \param numVectors [in] Number of columns of the resulting MultiVector.
    // NOTE: A map AND an importer need to be arguments because in serial, the importer will be null
    // This will default to importer->getTargetMap() being the active MultiVector
    FEMultiVector(const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > & map,
                  const Teuchos::RCP<const Import<LocalOrdinal,GlobalOrdinal,Node> >& importer,
                  const size_t numVecs,
                  const bool zeroOut = true);

    //! Destructor (virtual for memory safety of derived classes).
    virtual ~FEMultiVector () {}
    
    //@}
    //! @name Post-construction modification routines
    //@{

    // ! Calls endFill()
    void globalAssemble() {endFill();}

    //! Calls doTargetToSource() and then activates the source map mode
    void endFill()
    {
      if(*activeMultiVector_ == FE_ACTIVE_TARGET) {
        doTargetToSource(Tpetra::ADD);
        switchActiveMultiVector();
      }
      else
        throw std::runtime_error("FEMultiVector: Source MultiVector already active.  Cannot endFill()");
    }

    //! Activates the target map mode
    void beginFill()
    {
      // Note: This does not throw an error since the on construction, the FEMV is in target mode.  Ergo, calling beginFill(),
      // like one should expect to do in a rational universe, should not cause an error.
      if(*activeMultiVector_ == FE_ACTIVE_SOURCE) {
        switchActiveMultiVector();
      }
    }


    //! Migrate data from the target to the source map
    // Since this is non-unique -> unique, we need a combine mode.
    // Precondition: Target MultiVector must be active
    void doTargetToSource(const CombineMode CM=Tpetra::ADD);


    //! Migrate data from the source to the target map
    // Precondition: Source MultiVector must be active
    void doSourceToTarget(const CombineMode CM=Tpetra::ADD);

    //! Switches which Multivector is active (without migrating data)
    void switchActiveMultiVector();

  private:
    /// \brief Default c'tor (private so it does not get used)
    ///
    FEMultiVector() { }

  protected:
    //@{
    //! @name Internal routines and data structures

    /// \brief Replace the underlying Map in place.  NOTE: FEMultiVector does not allow this and will throw
    void replaceMap (const Teuchos::RCP<const map_type>& map);

    // Enum for activity
    enum FEWhichActive
    {
      FE_ACTIVE_TARGET,
      FE_ACTIVE_SOURCE
    };


    // This is whichever multivector isn't currently active
    Teuchos::RCP< MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> > inactiveMultiVector_;
    // This is in RCP to make shallow copies of the FEMultiVector work correctly
    Teuchos::RCP<FEWhichActive> activeMultiVector_;

    // Importer
    Teuchos::RCP<const Import<LocalOrdinal,GlobalOrdinal,Node> > importer_;

    //@}
  }; // class FEMultiVector

} // namespace Tpetra

#endif // TPETRA_FEMULTIVECTOR_DECL_HPP
