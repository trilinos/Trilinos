// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TEUCHOS_EXCEPTIONS_HPP
#define TEUCHOS_EXCEPTIONS_HPP


#include "Teuchos_ConfigDefs.hpp"


namespace Teuchos {


/** \brief Base exception class for Teuchos
 *
 * \ingroup teuchos_mem_mng_grp
 */
class ExceptionBase : public std::logic_error
{public:ExceptionBase(const std::string& what_arg) : std::logic_error(what_arg) {}};
// 2007/11/07: rabartl: Above, I had to change the name from Exception to
// ExceptionBase because Marzio did a 'using namespace Teuchos' and then he
// declared his own Exception class.  The file Laplacian3D.cpp failed to
// compile.  STOP DOING USING NAMESPACE BLAH!!!!!!


/** \brief Thrown if a duplicate owning RCP is creatd the the same object.
 *
 * \ingroup teuchos_mem_mng_grp
 */
class DuplicateOwningRCPError : public ExceptionBase
{public:DuplicateOwningRCPError(const std::string& what_arg) : ExceptionBase(what_arg) {}};


/** \brief Null reference error exception class.
 *
 * \ingroup teuchos_mem_mng_grp
 */
class NullReferenceError : public ExceptionBase
{public:NullReferenceError(const std::string& what_arg) : ExceptionBase(what_arg) {}};


/** \brief Null reference error exception class.
 *
 * \ingroup teuchos_mem_mng_grp
 */
class NonconstAccessError : public ExceptionBase
{public:NonconstAccessError(const std::string& what_arg) : ExceptionBase(what_arg) {}};


/** \brief Range error exception class.
 *
 * \ingroup teuchos_mem_mng_grp
 */
class RangeError : public ExceptionBase
{public:RangeError(const std::string& what_arg) : ExceptionBase(what_arg) {}};


/** \brief Dangling reference error exception class.
 *
 * \ingroup teuchos_mem_mng_grp
 */
class DanglingReferenceError : public ExceptionBase
{public:DanglingReferenceError(const std::string& what_arg) : ExceptionBase(what_arg) {}};


/** \brief Incompatiable iterators error exception class.
 *
 * \ingroup teuchos_mem_mng_grp
 */
class IncompatibleIteratorsError : public ExceptionBase
{public:IncompatibleIteratorsError(const std::string& what_arg) : ExceptionBase(what_arg) {}};

/** \brief Optionally thrown when a sublist is set twice by either
 * updateParametersFromXmlFile(), updateParametersFromXmlFileAndUpdate() or
 * updateParametersFromXmlString()
 *
 * \relates \c ParameterList
 */
class DuplicateParameterSublist : public ExceptionBase {

public:
  DuplicateParameterSublist(const std::string& what_arg):
    ExceptionBase(what_arg){}

};

/** \brief Thrown when a Parameter Entry that is already being tracked
 * is attempted to be inserted again into the masterParameterEntryMap
 * and masterIDMap
 *
 * \relates \c ParameterEntry
 */
class DuplicateParameterEntryException : public ExceptionBase {

public:
  DuplicateParameterEntryException(const std::string& what_arg):
    ExceptionBase(what_arg){}

};

/** \brief Thrown when a Parameter Entry ID that is already being used
 * is attempted to be reused again.
 *
 * \relates \c ParameterEntry
 */
class DuplicateParameterEntryIDException : public ExceptionBase {

public:
  DuplicateParameterEntryIDException(const std::string& what_arg):
    ExceptionBase(what_arg){}

};

/** \brief Thrown when a ParameterEntryValidatorID that
 * is already being used is attempted to be reused again.
 *
 * \relates ParameterEntryValidator
 */
class DuplicateValidatorIDException : public ExceptionBase {

public:
  DuplicateValidatorIDException(const std::string& what_arg):
    ExceptionBase(what_arg){}

};

/**
 * @brief Exception class for non-printable parameter types,
 * such as enum class/std::vector and many more
 * which don't define an operator<<.
 * Thrown during runtime when trying to print a parameter list
 * with a non-printable parameter entry.
 *
 * \relates ParameterEntry
 */
class NonprintableTypeException : public ExceptionBase {

public:
    NonprintableTypeException(const std::string& what_arg) :
            ExceptionBase(what_arg) {}

};



} // end namespace Teuchos


#endif	// TEUCHOS_EXCEPTIONS_HPP
