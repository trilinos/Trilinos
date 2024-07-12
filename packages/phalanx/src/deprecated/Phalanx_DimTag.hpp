//@HEADER
// *****************************************************************************
//        Phalanx: A Partial Differential Equation Field Evaluation 
//       Kernel for Flexible Management of Complex Dependency Chains
//
// Copyright 2008 NTESS and the Phalanx contributors.
//
//                Shards : Shared Discretization Tools
//            The Array Dims are a copy of shards array dims.
//
// Copyright 2008 NTESS and the Shards contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef PHALANX_DIM_TAG_HPP
#define PHALANX_DIM_TAG_HPP

//----------------------------------------------------------------------

#include "Phalanx_KokkosDeviceTypes.hpp"
#include "Phalanx_MDField_ExtentTraits.hpp"
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

/** \brief  Macro for implementing the body of a simple ArrayDimTag
 *  \param ADT  name of the tag.
 */
#define PHX_EXTENT_SPECIALIZATION( ADT ) \
  namespace PHX { \
    template<> struct is_extent<ADT> : std::true_type {}; \
  }

/** \brief  Macro for implementing the body of a simple ArrayDimTag
 *  \param ADT  name of the tag.
 */
#define PHX_EXTENT_PRINT( ADT ) \
  namespace PHX { \
    template<> print<ADT> {"ADT"}; \
  }

} // end namespace PHX

#endif
