/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/
#ifndef stk_percept_OptionMask_hpp
#define stk_percept_OptionMask_hpp

#include <iostream>
#include <iomanip>
#include <fstream>

#include <stdio.h>
#include <stdlib.h>      

#if defined( STK_HAS_MPI )
#include <mpi.h>
#endif

#include <stk_util/diag/Writer.hpp>
#include <stk_util/diag/PrintTimer.hpp>
#include <stk_util/util/Bootstrap.hpp>
#include <stk_util/util/IndentStreambuf.hpp>

#include <stk_percept/RunEnvironment.hpp>
#include <stk_percept/Util.hpp>


/// copied and edited from stk_util/use_cases/UseCaseEnvironment 

namespace {

  // Parse command line bit masks and produce -h documentation. (Probably moved to Util at some point)
  typedef unsigned long OptionMask;

  struct OptionMaskName
  {
    OptionMaskName()
      : m_name(""),
        m_mask(0),
        m_description("")
    {}

    OptionMaskName(const std::string &name, const OptionMask &mask, const std::string &description = "No description available")
      : m_name(name),
        m_mask(mask),
        m_description(description)
    {}

    virtual ~OptionMaskName()
    {}

    std::string		m_name;
    OptionMask		m_mask;
    std::string		m_description;
  };


  class OptionMaskNameMap: public std::map<std::string, OptionMaskName>
  {
  public:
    void mask(const std::string &name, const OptionMask mask, const std::string &description) {
      iterator it = find(name);
      if (it == end())
        insert(std::make_pair(name, OptionMaskName(name, mask, description)));
      else {
        (*it).second.m_mask = mask;
        (*it).second.m_description = description;
      }
    }
  };

  class OptionMaskParser
  {
  public:
    typedef OptionMask Mask;		///< Mask for this option

  public:
    /**
     * Creates a new <b>OptionMaskParser</b> instance.
     *
     */
    OptionMaskParser(const std::string &description)
      : m_optionMaskNameMap(),
        m_description(description),
        m_optionMask(0),
        m_status(true)
    {}

    virtual ~OptionMaskParser()
    {}

    Mask parse(const char *mask) const;

    virtual void parseArg(const std::string &name) const;

    std::string describe() const {
      std::ostringstream strout;
      strout << m_description << std::endl;
      for (OptionMaskNameMap::const_iterator it = m_optionMaskNameMap.begin(); it != m_optionMaskNameMap.end(); ++it)
        strout << "  " << (*it).first << std::setw(14 - (*it).first.size()) << " " << (*it).second.m_description << std::endl;
      return strout.str();
    }

    void mask(const std::string &name, const Mask mask, const std::string &description) {
      m_optionMaskNameMap.mask(name, mask, description);
    }

  protected:
    OptionMaskNameMap		m_optionMaskNameMap;	///< Mask name vector
    std::string                   m_description;          ///< Help description
    mutable OptionMask		m_optionMask;		///< Most recently parsed mask
    mutable bool			m_status;		///< Result of most recent parse
  };


  OptionMaskParser::Mask
  OptionMaskParser::parse(
                          const char *          mask) const
  {
    if (mask) {
      const std::string mask_string(mask);

      m_status = true;

      std::string::const_iterator it0 = mask_string.begin();
      std::string::const_iterator it1;
      std::string::const_iterator it2;
      std::string::const_iterator it3;
      do {
        // Trim preceeding spaces
        while (it0 != mask_string.end() && *it0 == ' ')
          it0++;

        if (it0 == mask_string.end())
          break;

        for (it1 = it0; it1 != mask_string.end(); ++it1) {
          if (*it1 == '(' || *it1 == ':' || *it1 == ',')
            break;
        }

        // Trim trailing spaces
        it2 = it1;
        while (it2 != it0 && *(it2 - 1) == ' ')
          --it2;

        std::string name(it0, it2);

        // Get argument list
        if (*it1 == '(') {
          it2 = it1 + 1;

          // Trim preceeding spaces
          while (it2 != mask_string.end() && *it2 == ' ')
            ++it2;

          int paren_count = 0;

          for (; it1 != mask_string.end(); ++it1) {
            if (*it1 == '(')
              ++paren_count;
            else if (*it1 == ')') {
              --paren_count;
              if (paren_count == 0)
                break;
            }
          }
          it3 = it1;

          // Trim trailing spaces
          while (it3 != it2 && *(it3 - 1) == ' ')
            --it3;

          // Find next argument start
          for (; it1 != mask_string.end(); ++it1)
            if (*it1 == ':' || *it1 == ',')
              break;
        }
        else
          it2 = it3 = it1;

        const std::string arg(it2, it3);

        parseArg(name);

        it0 = it1 + 1;
      } while (it1 != mask_string.end());
    }

    return m_optionMask;
  }


  void
  OptionMaskParser::parseArg(
                             const std::string &	name) const
  {
    OptionMaskNameMap::const_iterator mask_entry = m_optionMaskNameMap.find(name);

    if (mask_entry != m_optionMaskNameMap.end()) m_optionMask |= (*mask_entry).second.m_mask;
    else {
      Mask	mask_hex = 0;
      std::istringstream mask_hex_stream(name.c_str());
      if (mask_hex_stream >> std::resetiosflags(std::ios::basefield) >> mask_hex)
        m_optionMask |= mask_hex;
      else
        m_status = false;
    }
  }
}

#endif
