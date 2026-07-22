// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef Teuchos_PARAMETERXMLFILEREADER_H
#define Teuchos_PARAMETERXMLFILEREADER_H

#include "Teuchos_ConfigDefs.hpp"
#include "Teuchos_FileInputSource.hpp"
#include "Teuchos_XMLParameterListReader.hpp"


namespace Teuchos
{
  /**
   * Reader for getting parameter lists from XML files
   */
  class ParameterXMLFileReader
    {
    public:
      /** \brief Constructor */
      ParameterXMLFileReader(const std::string& filename);

      /** */
      ParameterList getParameters() const ;
    private:
      FileInputSource fis_;
    };

}
#endif

