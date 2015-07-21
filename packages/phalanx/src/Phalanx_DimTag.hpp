/*
//@HEADER
//@HEADER
// 
// The Array Dims are a copy of shards array dims.
//
// ************************************************************************
//
//                Shards : Shared Discretization Tools
//                 Copyright 2008 Sandia Corporation
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
// Questions? Contact Carter Edwards (hcedwar@sandia.gov),
//                    Pavel Bochev (pbboche@sandia.gov), or
//                    Denis Ridzal (dridzal@sandia.gov).
//
// ************************************************************************
*/

#ifndef PHALANX_DIM_TAG_HPP
#define PHALANX_DIM_TAG_HPP

//----------------------------------------------------------------------

#include "Phalanx_KokkosDeviceTypes.hpp"
#include <vector>
#include <string>

//----------------------------------------------------------------------

namespace PHX {

//----------------------------------------------------------------------

/** \class  DimTag
 *  \brief  Abstract base class for array dimension tags supplied to
 *          the Array template class.
 *  \sa Array
 *
 *  A derived array dimension tag class must provide the
 *  <b> name </b> method and <b> tag </b> singleton method
 *  as in the following example.
 *  <PRE>
 *  struct MyTag : public PHX::DimTag {
 *    const char * name() const ;
 *    static const MyTag & tag();
 *  };
 *  </PRE>
 *  An example implementation of these methods is as follows.
 *  <PRE>
 *  const char * MyTag::name() const
 *  { static const char my_name[] = "MyTag" ; return my_name ; }
 *
 *  const MyTag & MyTag::tag()
 *  { static const MyTag my_tag ; return my_tag ; }
 *  </PRE>
 */
class DimTag {
public:

  typedef PHX::index_size_type size_type;

  /** \brief Name of the tag, typically the name of the derived class. */
  virtual const char * name() const = 0 ;

  /** \brief  Given a dimension and index produce a string for output.
   *
   *          Default to converting <b> index </b> to a string.
   */
  virtual std::string to_string( size_type dimension ,
                                 size_type index ) const ;

  /** \brief Given a dimension and input strige produce an index.
   *
   *          Default to converting <b> label </b> to an integer.
   */
  virtual size_type to_index( size_type dimension ,
                              const std::string & label ) const ; 
 
protected:
  virtual ~DimTag();
  DimTag() {}
  
private:
  DimTag( const DimTag & );
  DimTag & operator = ( const DimTag & );
};

/** \brief  Macro for declaration of a simple ArrayDimTag
 *  \param ADT  name of the tag.
 */
#define PHX_DIM_TAG_DECLARATION( ADT ) \
  class ADT : public PHX::DimTag { \
  public: \
    const char * name() const ; \
    static const ADT & tag(); \
  private: \
    ~ADT(); \
    ADT(); \
    ADT( const ADT & ); \
    ADT & operator = ( const ADT & ); \
  };

/** \brief  Macro for implementing the body of a simple ArrayDimTag
 *  \param ADT  name of the tag.
 */
#define PHX_DIM_TAG_IMPLEMENTATION( ADT ) \
  ADT::ADT() {} \
  ADT::~ADT() {} \
  const char * ADT::name() const { static const char n[] = # ADT; return n; } \
  const ADT & ADT::tag() { static const ADT self ; return self ; }

} // end namespace PHX

#endif
