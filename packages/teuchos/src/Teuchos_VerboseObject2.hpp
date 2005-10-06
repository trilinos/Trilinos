// @HEADER
// ***********************************************************************
// 
//                    Teuchos: Common Tools Package
//                 Copyright (2004) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#ifndef TEUCHOS_VERBOSE_OBJECT_HPP
#define TEUCHOS_VERBOSE_OBJECT_HPP

#include "Teuchos_RefCountPtr.hpp"
#include "Teuchos_FancyOStream2.hpp"
#include "Teuchos_VerbosityLevel.hpp"

namespace Teuchos {

/** \brief Base class for objects that can print their activities to a stream.
 */
class VerboseObject {
public:

  /** \name Public static member functions */
  //@{

  /** \brief Set the default output stream object. */
  static void setDefaultOStream( const RefCountPtr<FancyOStream> &defaultOStream );

  /** \brief Get the default output stream object. */
  static RefCountPtr<FancyOStream> getDefaultOStream();

  /** \brief Set the default verbosity level. */
  static void setDefaultVerbLevel( const EVerbosityLevel defaultVerbLevel);

  /** \brief Get the default verbosity level. */
  static EVerbosityLevel getDefaultVerbLevel();

  //@}

  /** \name Constructors/Initializers */
  //@{
  
  /** \brief . */
  virtual ~VerboseObject() {}
  
  /** \brief Calls <tt>initializeVerboseObject()</tt>.
   */
  explicit
  VerboseObject(
    const EVerbosityLevel              verbLevel = VERB_DEFAULT
    ,const RefCountPtr<FancyOStream>   &oStream  = Teuchos::null
    );
  
  /** \brief Calls <tt>initializeVerboseObject()</tt>.
   */
  virtual void initializeVerboseObject(
    const EVerbosityLevel              verbLevel = VERB_DEFAULT
    ,const RefCountPtr<FancyOStream>   &oStream  = Teuchos::null
    );

  /** \brief Set output stream */
  virtual VerboseObject& setOStream(const RefCountPtr<FancyOStream> &oStream);

  /** \brief Set the verbosity level */
  virtual VerboseObject& setVerbLevel(const EVerbosityLevel verbLevel);

  /** \brief Set the verbosity level */
  virtual VerboseObject& setLinePrefix(const std::string &linePrefix);

  //@}

  /** \name Query functions */
  //@{

  /** \brief Return the output stream to be used.
   */
  virtual RefCountPtr<FancyOStream> getOStream() const;

  /** \brief Get the verbosity level */
  virtual EVerbosityLevel getVerbLevel() const;

  //@}

  /** \name Utilities */
  //@{

  /** \brief Create a tab object with this' tab tag.
   *
   * Returns <tt>OSTab(this->getOStream(),tabs)<tt>
   */
  virtual OSTab getOSTab(const int tabs = 1, const std::string &linePrefix = "") const;

  //@}
  
private:

  static RefCountPtr<FancyOStream>   defaultOStream_;
  static EVerbosityLevel             defaultVerbLevel_;

  RefCountPtr<FancyOStream>   thisOStream_;
  EVerbosityLevel             thisVerbLevel_;
  std::string                 thisLinePrefix_;

};

} // namespace Teuchos

#endif // TEUCHOS_VERBOSE_OBJECT_HPP
