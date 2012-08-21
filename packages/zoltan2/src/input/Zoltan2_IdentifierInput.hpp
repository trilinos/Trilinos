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

/*! \file Zoltan2_IdentifierInput.hpp
    \brief Defines the IdentifierInput interface.
*/


#ifndef _ZOLTAN2_IDENTIFIERINPUT_HPP_
#define _ZOLTAN2_IDENTIFIERINPUT_HPP_

#include <Zoltan2_InputAdapter.hpp>
#include <Zoltan2_PartitioningSolution.hpp>

#include <string>

namespace Zoltan2 {

/*!  \brief IdentifierInput defines the interface for identifiers.

    Zoltan2 can partition a simple list of weighted identifiers 
    with no geometry or topology provided.  IdentifierInput defines
    the interface for input adapters of this type.

    InputAdapter objects provide access for Zoltan2 to the user's data.
    Many built-in adapters are already defined for common data structures,
    such as Tpetra and Epetra objects and C-language pointers to arrays.

    Data types:
    \li \c scalar_t object weights
    \li \c lno_t    local indices and local counts
    \li \c gno_t    global indices and global counts
    \li \c gid_t    application global Ids
    \li \c node_t is a sub class of Kokkos::StandardNodeMemoryModel

    See IdentifierTraits to understand why the user's global ID type (\c gid_t)
    may differ from that used by Zoltan2 (\c gno_t).

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

*/

template <typename User>
  class IdentifierInput : public InputAdapter<User> {

public:

#ifndef DOXYGEN_SHOULD_SKIP_THIS
  typedef typename InputTraits<User>::scalar_t    scalar_t;
  typedef typename InputTraits<User>::lno_t    lno_t;
  typedef typename InputTraits<User>::gno_t    gno_t;
  typedef typename InputTraits<User>::gid_t    gid_t;
  typedef typename InputTraits<User>::node_t   node_t;
  typedef User user_t;
#endif

  /*! \brief Destructor 
   */
  virtual ~IdentifierInput() {};

  ////////////////////////////////////////////////////
  // The InputAdapter interface.
  ////////////////////////////////////////////////////

  enum InputAdapterType inputAdapterType() const {return IdentifierAdapterType;}

  ////////////////////////////////////////////////////
  // My interface.
  ////////////////////////////////////////////////////

  /*! \brief Return the number of identifiers on this process.
   */
  virtual size_t getLocalNumberOfIdentifiers() const = 0;

  /*! \brief Return the number of weights associated with each identifier.
   *   If the number of weights is zero, then we assume that the
   *   identifiers are equally weighted.
   */
  virtual int getNumberOfWeights() const = 0;

  /*! \brief Provide a pointer to this process' identifiers.

      \param Ids will on return point to the list of the global Ids for 
        this process.

       \return The number of ids in the Ids list.
   */

  virtual size_t getIdentifierList(gid_t const *&Ids) const = 0;

  /*! \brief Provide a pointer to one of the dimensions of this process' 
                optional weights.

      \param dimension is a value ranging from zero to one less than 
                   getNumberOfWeights()
      \param weights on return will contain a list of the weights for
               the dimension specified.  If weights for
           this dimension are to be uniform for all identifierse in the
           global problem, the \c weights should be a NULL pointer.

      \param stride on return will indicate the stride of the weights list.


       If stride is \c k then the weight 
       corresponding to the identifier Ids[n] (returned in getIdentifierList)
       should be found at weights[k*n].

       \return The number of values in the weights list.  This may be greater
          than the number of identifiers, because the stride may be greater
          than one.
   */

  virtual size_t getIdentifierWeights(int dimension,
     const scalar_t *&weights, int &stride) const = 0;


 /*! \brief Apply a PartitioningSolution to an input.
   *
   *  This is not a required part of the IdentifierInput interface. However
   *  if the Caller calls a Problem method to redistribute data, it needs
   *  this method to perform the redistribution.
   *
   *  \param in  An input object with a structure and assignment of
   *           of global Ids to processes that matches that of the input
   *           data that instantiated this InputAdapter.
   *  \param out On return this should point to a newly created object
   *            with the specified partitioning.
   *  \param solution  The Solution object created by a Problem should
   *      be supplied as the third argument.  It must have been templated
   *      on user data that has the same global ID distribution as this
   *      user data.
   *  \return   Returns the number of local Ids in the new partitioning.
   */

  template <typename Adapter>
    size_t applyPartitioningSolution(User &in, User *&out,
      const PartitioningSolution<Adapter> &solution) const
  {
    return 0;
  } 
};
  
  
}  //namespace Zoltan2
  
#endif
