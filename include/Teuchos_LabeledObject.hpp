// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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
class TEUCHOSCORE_LIB_DLL_EXPORT LabeledObject {
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
