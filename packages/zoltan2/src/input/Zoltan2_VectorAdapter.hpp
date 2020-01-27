// @HEADER
//
// ***********************************************************************
//
//   Zoltan2: A package of combinatorial algorithms for scientific computing
//                  Copyright 2012 Sandia Corporation
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
// Questions? Contact Karen Devine      (kddevin@sandia.gov)
//                    Erik Boman        (egboman@sandia.gov)
//                    Siva Rajamanickam (srajama@sandia.gov)
//
// ***********************************************************************
//
// @HEADER

/*! \file Zoltan2_VectorAdapter.hpp
    \brief Defines the VectorAdapter interface.
*/


#ifndef _ZOLTAN2_VECTORADAPTER_HPP_
#define _ZOLTAN2_VECTORADAPTER_HPP_

#include <Zoltan2_Adapter.hpp>

namespace Zoltan2 {

  /*!  \brief VectorAdapter defines the interface for vector input.

    Adapter objects provide access for Zoltan2 to the user's data.
    Many built-in adapters are already defined for common data structures,
    such as Tpetra and Epetra objects and C-language pointers to arrays.

    Data types:
    \li \c scalar_t weights and vector element values
    \li \c lno_t    local indices and local counts
    \li \c gno_t    global indices and global counts
    \li \c node_t is a Kokkos Node type

    The Kokkos node type can be safely ignored.

    The template parameter \c User is a user-defined data type
    which, through a traits mechanism, provides the actual data types
    with which the Zoltan2 library will be compiled.
    \c User may be the actual class or structure used by application to
    represent a vector, or it may be the helper class BasicUserTypes.
    See InputTraits for more information.

    The \c scalar_t type, representing use data such as matrix values, is
    used by Zoltan2 for weights, coordinates, part sizes and
    quality metrics.
    Some User types (like Tpetra::CrsMatrix) have an inherent scalar type,
    and some
    (like Tpetra::CrsGraph) do not.  For such objects, the scalar type is
    set by Zoltan2 to \c float.  If you wish to change it to double, set
    the second template parameter to \c double.


    VectorAdapter may be a single vector or a set of corresponding vectors
    which have with the same global identifiers and the same distribution
    across processes.

  \todo We can make global Ids optional.  We don't need them.

*/

template <typename User>
  class VectorAdapter : public BaseAdapter<User> {
public:

#ifndef DOXYGEN_SHOULD_SKIP_THIS
  typedef typename InputTraits<User>::scalar_t scalar_t;
  typedef typename InputTraits<User>::lno_t    lno_t;
  typedef typename InputTraits<User>::gno_t    gno_t;
  typedef typename InputTraits<User>::part_t   part_t;
  typedef typename InputTraits<User>::node_t   node_t;
  typedef typename InputTraits<User>::offset_t offset_t;
  typedef User user_t;
  typedef User userCoord_t;
  typedef VectorAdapter<User> base_adapter_t;
#endif

  /*! \brief Destructor
   */
  virtual ~VectorAdapter() {};

  ////////////////////////////////////////////////////
  // The Adapter interface.
  ////////////////////////////////////////////////////

  enum BaseAdapterType adapterType() const {return VectorAdapterType;}

  ///////////////////////////////////////////////////////////////
  // User's adapter interface:
  // The user must implement these methods in his VectorAdapter
  ///////////////////////////////////////////////////////////////

  /*! \brief Return the number of vectors.
   */
  virtual int getNumEntriesPerID() const = 0;

  /*! \brief Provide a pointer to the elements of the specified vector.
      \param elements will on return point to the vector values
        corresponding to the global Ids.
      \param stride the k'th element is located at elements[stride*k]
      \param idx ranges from zero to one less than getNumEntriesPerID(), and
         represents the vector for which data is being requested.
   */
  virtual void getEntriesView(const scalar_t *&elements,
    int &stride, int idx = 0) const {
    // If adapter does not define getEntriesView, getEntriesKokkosView is called.
    // If adapter does not define getEntriesKokkosView, getEntriesView is called.
    // Allows forward and backwards compatibility.
    // coordinates in MJ are LayoutLeft since Tpetra Multivector gives LayoutLeft
    Kokkos::View<scalar_t **, Kokkos::LayoutLeft,
      typename node_t::device_type> kokkosEntries;
    getEntriesKokkosView(kokkosEntries);
    elements = kokkosEntries.data();
    stride = 1;
  }

  /*! \brief Provide a Kokkos view to the elements of the specified vector.
      \param elements will on return point to the vector values
        corresponding to the global Ids.
   */
  virtual void getEntriesKokkosView(
    // coordinates in MJ are LayoutLeft since Tpetra Multivector gives LayoutLeft
    Kokkos::View<scalar_t **, Kokkos::LayoutLeft,
      typename node_t::device_type> & elements) const {
    // If adapter does not define getEntriesKokkosView, getEntriesView is called.
    // If adapter does not define getEntriesView, getEntriesKokkosView is called.
    // Allows forward and backwards compatibility.
    // coordinates in MJ are LayoutLeft since Tpetra Multivector gives LayoutLeft
    typedef Kokkos::View<scalar_t **, Kokkos::LayoutLeft,
      typename node_t::device_type> kokkos_entries_view_t;
    elements = kokkos_entries_view_t("entries", this->getLocalNumIDs(),
      this->getNumEntriesPerID());
    typename kokkos_entries_view_t::HostMirror host_elements =
      Kokkos::create_mirror_view(elements);
    for(int j = 0; j < this->getNumEntriesPerID(); ++j) {
      const scalar_t * ptr_elements;
      int stride;
      getEntriesView(ptr_elements, stride, j);
      size_t i = 0;
      for(size_t n = 0; n < this->getLocalNumIDs() * stride; n += stride) {
        host_elements(i++,j) = ptr_elements[n];
      }
    }
    Kokkos::deep_copy(elements, host_elements);
  }

  /*! \brief Write files that can be used as input to Zoltan or Zoltan2 driver
   *  Creates chaco-formatted input files for coordinates and weights that
   *  can be used as input for Zoltan or Zoltan2 drivers.
   *  This routine is SERIAL and can be quite slow.
   *  It is meant as a debugging tool only, to allow Zoltan developers to 
   *  replicate performance that applications are seeing using the applicatios'
   *  input.
   */
  void generateFiles(
    const char *fileprefix, 
    const Teuchos::Comm<int> &comm
  ) const 
  {
    // Generate the graph file with weights using the base adapter method
    this->generateWeightFileOnly(fileprefix, comm);

    //  Generate the coords file with local method
    this->generateCoordsFileOnly(fileprefix, comm);
  }

  ////////////////////////////////////////////////////////////////
  // Handy pseudonyms, since vectors are often used as coordinates
  // User should not implement these methods.
  ////////////////////////////////////////////////////////////////

  inline int getDimension() const {return getNumEntriesPerID();}

  inline void getCoordinatesView(const scalar_t *&elements, int &stride,
                                 int idx = 0) const
  {
    getEntriesView(elements, stride, idx);
  }

  inline void getCoordinatesKokkosView(
    // coordinates in MJ are LayoutLeft since Tpetra Multivector gives LayoutLeft
    Kokkos::View<scalar_t **, Kokkos::LayoutLeft, typename node_t::device_type> & elements) const
  {
    getEntriesKokkosView(elements);
  }

private:

  void generateCoordsFileOnly(
    const char* fileprefix, 
    const Teuchos::Comm<int> &comm) const;

};

template <typename User>
void VectorAdapter<User>::generateCoordsFileOnly(
  const char *fileprefix, 
  const Teuchos::Comm<int> &comm
) const
{
  // Writes a chaco-formatted coordinates file
  // This function is SERIAL and can be quite slow.  Use it for debugging only.

  int np = comm.getSize();
  int me = comm.getRank();

  // append suffix to filename
  
  std::string filenamestr = fileprefix;
  filenamestr = filenamestr + ".coords";
  const char *filename = filenamestr.c_str();

  for (int p = 0; p < np; p++) {

    // Is it this processor's turn to write to files?
    if (me == p) {

      std::ofstream fp;
      if (me == 0) {
        // open file for writing
        fp.open(filename, std::ios::out);
      }
      else {
        // open file for appending
        fp.open(filename, std::ios::app);
      }
    
      // Get the vector entries
      size_t len = this->getLocalNumIDs();
      int nvec = this->getNumEntriesPerID();
      const scalar_t **values = new const scalar_t *[nvec];
      int *strides = new int[nvec];
      for (int n = 0; n < nvec; n++)
        getEntriesView(values[n], strides[n], n);

      // write vector entries to coordinates file

      for (size_t i = 0; i < len; i++) {
        for (int n = 0; n < nvec; n++)
          fp << values[n][i*strides[n]] << " ";
        fp << "\n";
      }

      // clean up and close the file
      delete [] strides;
      delete [] values;
      fp.close();
    }
    comm.barrier();
  }
}


}  //namespace Zoltan2

#endif
