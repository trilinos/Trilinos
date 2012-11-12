// @HEADER
// ***********************************************************************
//
//                    Teuchos: Common Tools Package
//                 Copyright (2004) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// ***********************************************************************
// @HEADER

#ifndef TEUCHOS_LABELED_OBJECT_HPP
#define TEUCHOS_LABELED_OBJECT_HPP

#include "Teuchos_ConfigDefs.hpp"


namespace Teuchos {


/** \brief Base class for objects that contain a std::string label.
 *
 * The object label std::string <tt>objectLabel</tt> set in
 * <tt>setObjectLabel()</tt> should be a simple one-line label given to an
 * object to differentiate it from all other objects.  A subclass
 * implementation can define a default label in some cases but typically this
 * label is designed for end users to set to give the object a name that is
 * meaningful to the user.  The label should not contain any information about
 * the actual type of the object.  Adding type information is appropriate in
 * the <tt>Describable</tt> interface, which inherits from this interface.
 *
 * This base class provides a default implementation for the functions
 * <tt>setObjectLabel()</tt> and <tt>getObjectLabel()</tt> as well as private
 * data to hold the label.  Subclasses can override these functions but
 * general, there should be no need to do so.
 *
 * \ingroup teuchos_outputting_grp
 */
class TEUCHOS_LIB_DLL_EXPORT LabeledObject {
public:
  /** \brief Construct with an empty label. */
  LabeledObject();
  /** \brief . */
  virtual ~LabeledObject();
  /** \brief Set the object label (see LabeledObject). */
  virtual void setObjectLabel( const std::string &objectLabel );
  /** \brief Get the object label (see LabeledObject). */
  virtual std::string getObjectLabel() const;
private:
  std::string objectLabel_;
};


} // namespace Teuchos


#endif // TEUCHOS_LABELED_OBJECT_HPP
