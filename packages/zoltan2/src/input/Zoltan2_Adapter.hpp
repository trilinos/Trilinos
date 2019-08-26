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

/*! \file Zoltan2_BaseAdapter.hpp
    \brief Defines the Adapter interface for accessing user data.
*/

#ifndef _ZOLTAN2_ADAPTER_HPP_
#define _ZOLTAN2_ADAPTER_HPP_

#include <Kokkos_Core.hpp>
#include <Zoltan2_Standards.hpp>
#include <Zoltan2_InputTraits.hpp>
#include <Zoltan2_PartitioningSolution.hpp>

namespace Zoltan2 {

/*! \brief An enum to identify general types of adapters.
 *
 */
enum BaseAdapterType {
  InvalidAdapterType = 0, /*!< \brief unused value */
  IdentifierAdapterType,  /*!< \brief identifier data, just a list of IDs*/
  VectorAdapterType,      /*!< \brief vector data */
  MatrixAdapterType,      /*!< \brief matrix data */
  GraphAdapterType,       /*!< \brief graph data */
  MeshAdapterType         /*!< \brief mesh data */
};

/*! \brief BaseAdapter defines methods required by all Adapters

    Adapters provide access from Zoltan2 to the user's data.  The
    methods in the interface must be defined by users.  Many built-in
    adapters are already defined for common data structures, such as
    Tpetra and Epetra objects and C-language pointers to arrays.

 */

class BaseAdapterRoot {
public:
  virtual ~BaseAdapterRoot() {}; // required virtual declaration

  /*! \brief Returns the number of objects on this process
   *
   *  Objects may be coordinates, graph vertices, matrix rows, etc.
   *  They are the objects to be partitioned, ordered, or colored.
   */
  virtual size_t getLocalNumIDs() const = 0;

  /*! \brief Returns the number of weights per object.
   *   Number of weights per object should be zero or greater.  If
   *   zero, then it is assumed that all objects are equally weighted.
   *   Default is zero weights per ID.
   */
  virtual int getNumWeightsPerID() const { return 0; };
};

template <typename User>
  class BaseAdapter : public BaseAdapterRoot {

public:
  typedef typename InputTraits<User>::lno_t lno_t;
  typedef typename InputTraits<User>::gno_t gno_t;
  typedef typename InputTraits<User>::scalar_t scalar_t;
  typedef typename InputTraits<User>::part_t part_t;  
  typedef typename InputTraits<User>::offset_t offset_t;

  /*! \brief Returns the type of adapter.
   */
  virtual enum BaseAdapterType adapterType() const = 0;

  /*! \brief Destructor
   */
  virtual ~BaseAdapter() {};

  /*! \brief Provide a pointer to this process' identifiers.

      \param ids will on return point to the list of the global Ids for 
        this process.
   */
  virtual void getIDsView(const gno_t *&ids) const {
    Kokkos::View<gno_t *> kokkosIds;
    getIDsKokkosView(kokkosIds);
    ids = kokkosIds.data();
  }

  /*! \brief Provide a pointer to this process' identifiers.

      \param ids will on return point to the list of the global Ids for 
        this process.
   */
  virtual void getIDsKokkosView(Kokkos::View<gno_t *> &/* ids */) const {
    Z2_THROW_NOT_IMPLEMENTED
  }

  ///*! \brief Provide pointer to a weight array with stride.
  // *    \param wgt on return a pointer to the weights for this idx
  // *    \param stride on return, the value such that
  // *       the \t nth weight should be found at <tt> wgt[n*stride] </tt>.
  // *    \param idx  the weight index, zero or greater
  // *   This function must be implemented in derived adapter if
  // *   getNumWeightsPerID > 0.
  // *   This function should not be called if getNumWeightsPerID is zero.
  // */ 
  virtual void getWeightsView(const scalar_t *&wgt, int &stride,
                              int idx = 0) const {
    Kokkos::View<scalar_t *> tempWeightsView;
    getWeightsKokkosView(tempWeightsView, idx);
    wgt = tempWeightsView.data();
    stride = 1;
  }

  ///*! \brief Provide pointer to a weight View.
  // *    \param wgt on return a pointer to the weights for this idx
  // *    \param idx  the weight index, zero or greater
  // *   This function must be implemented in derived adapter if
  // *   getNumWeightsPerID > 0.
  // *   This function should not be called if getNumWeightsPerID is zero.
  // */ 
  virtual void getWeightsKokkosView(Kokkos::View<scalar_t *> &/* wgt */, 
                              int /* idx */ = 0) const {
    Z2_THROW_NOT_IMPLEMENTED
  }

  /*! \brief Provide pointer to an array containing the input part 
   *         assignment for each ID.
   *         The input part information may be used for re-partitioning
   *         to reduce data movement, or for mapping parts to processes.
   *         Adapters may return NULL for this pointer (the default
   *         behavior); if NULL is returned, algorithms will assume
   *         the rank 
   *    \param inputPart on return a pointer to input part numbers
   */ 
  void getPartsView(const part_t *&inputPart) const {
    // Default behavior:  return NULL for inputPart array;
    // assume input part == rank
    inputPart = NULL;
  }

  /*! \brief Apply a PartitioningSolution to an input.
   *
   *  This is not a required part of the InputAdapter interface. However
   *  if the Caller calls a Problem method to redistribute data, it needs
   *  this method to perform the redistribution.
   *
   *  \param in  An input object with a structure and assignment of
   *             of global Ids to processes that matches that of the input
   *             data that instantiated this Adapter.
   *  \param out On return this should point to a newly created object
   *             with the specified partitioning.
   *  \param solution  The Solution object created by a Problem should
   *      be supplied as the third argument.  It must have been templated
   *      on user data that has the same global ID distribution as this
   *      user data.
   *  \return   Returns the number of local Ids in the new partition.
   */
  template <typename Adapter>
    void applyPartitioningSolution(const User &in, User *&out,
      const PartitioningSolution<Adapter> &solution) const {
    Z2_THROW_NOT_IMPLEMENTED
  }

protected:

  // Write Chaco-formatted graph and assign files echoing adapter input
  // This routine is serial and may be slow; use it only for debugging
  // This function does not write edge info to the graph file, as the
  // BaseAdapter does not know about edge info; it writes
  // only the Chaco header and vertex weights (if applicable).
  void generateWeightFileOnly(const char* fileprefix, 
                              const Teuchos::Comm<int> &comm) const;

};

template <typename User>
void BaseAdapter<User>::generateWeightFileOnly(
  const char *fileprefix,
  const Teuchos::Comm<int> &comm
) const
{
  int np = comm.getSize();
  int me = comm.getRank();

  size_t nLocalIDs = this->getLocalNumIDs();

  // Write .graph file:  header and weights only (no edges)
  // Adapters with edges have to implement their own generateFiles function
  // to provide edge info.
  {
    // append suffix to filename
    std::string filenamestr = fileprefix;
    filenamestr = filenamestr + ".graph";
    const char *filename = filenamestr.c_str();

    size_t nGlobalIDs;
    Teuchos::reduceAll(comm, Teuchos::REDUCE_SUM, 1, &nLocalIDs, &nGlobalIDs);

    int nWgts = this->getNumWeightsPerID();

    for (int p = 0; p < np; p++) {

      // Is it this processor's turn to write to files?
      if (me == p) {

        std::ofstream fp;

        if (me == 0) {
          // open file for writing
          fp.open(filename, std::ios::out);
          // write Chaco header info
          // this function assumes no edges
          fp << nGlobalIDs << " " << 0 << " " 
             << (nWgts ? "010" : "000") << " "
             << (nWgts > 1 ? std::to_string(nWgts) : " ") << std::endl;
        }
        else {
          // open file for appending
          fp.open(filename, std::ios::app);
        }

        if (nWgts) {

          // get weight data
          const scalar_t **wgts = new const scalar_t *[nWgts];
          int *strides = new int[nWgts];
          for (int n = 0; n < nWgts; n++)
            getWeightsView(wgts[n], strides[n], n);
  
          // write weights to file
          for (size_t i = 0; i < nLocalIDs; i++) {
            for (int n = 0; n < nWgts; n++)
              fp << wgts[n][i*strides[n]] << " ";
            fp << "\n";
          }
  
          delete [] strides;
          delete [] wgts;
        }
  
        fp.close();
      }
  
      comm.barrier();
    }
  }

  // write assignments file
  {
    std::string filenamestr = fileprefix;
    filenamestr = filenamestr + ".assign";
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
  
        const part_t *parts;
        this->getPartsView(parts);
       
        for (size_t i = 0; i < nLocalIDs; i++) {
          fp << (parts != NULL ? parts[i] : me) << "\n";
        }
        fp.close();
      }
  
      comm.barrier();
    }
  }
}

}  //namespace Zoltan2
  
#endif
