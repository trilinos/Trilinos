// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef Teuchos_FILEINPUTSOURCE_H
#define Teuchos_FILEINPUTSOURCE_H

/*! \file Teuchos_FileInputSource.hpp
    \brief Definition of XMLInputSource derived class for reading XML from a file
*/

#include "Teuchos_ConfigDefs.hpp"
#include "Teuchos_XMLInputSource.hpp"


namespace Teuchos
{
  /** \ingroup XML
   * \brief Instantiation of XMLInputSource class for reading XML from a file.
   */
  class TEUCHOSPARAMETERLIST_LIB_DLL_EXPORT FileInputSource : public XMLInputSource
    {
    public:
      /** \brief Constructor */
      FileInputSource(const std::string& filename);

      /** \brief Destructor */
      virtual ~FileInputSource(){;}

      /** \brief Create a FileInputStream */
      virtual RCP<XMLInputStream> stream() const;

    private:
      std::string filename_;
    };

}
#endif

