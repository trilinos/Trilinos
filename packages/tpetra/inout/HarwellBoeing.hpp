#ifndef __HarwellBoeing_hpp
#define __HarwellBoeing_hpp

#include <algorithm>
#include <map>
#include <string>

namespace HarwellBoeing {

  typedef std::map<std::string, std::string> dict_type;

  /// Is the string s a valid Fortran format specifier letter?
  ///
  /// Fortran format specifier strings contain a letter which
  /// indicates the type of data to read.  For example, "A" indicates
  /// a character string, "I" an integer, "X" a space, and "D" a real
  /// double-precision floating-point value.  If \c s is in the list
  /// of letters we recognize, return true, else return false.
  bool validFormatLetter (const std::string& s) {
    const char* validValues[] = {"A", "D", "E", "F", "G", "I", "X"};
    const int numValidValues = 7;
    return (validValues + numValidValues != 
	    std::find (s, validValues, validValues + numValidValues));
  }

  /// \brief Parse Fortran format specifier.
  ///
  /// Given a Fortran format specifier stored as a string, return a
  /// dictionary of (key,value) pairs containing the information
  /// extracted from the format specifier.
  ///
  /// \param fmt [in] Fortran format specifier string
  /// \param tolerant [in] Whether to parse tolerantly
  ///
  /// \return Dictionary containing the parsed information
  dict_type
  parseFormat (const std::string& fmt, const bool tolerant) 
  {
    const std::vector<std::string> tokens = tokenizeFormat (fmt);
    std::map<std::string, std::string> theDict;
    // Fill in the dictionary with (key,value) pairs retrieved from
    // the token string.
    parseFormatHelper (tokens.begin(), tokens.end(), theDict, tolerant);
    return theDict;
  }

  std::vector<std::string>
  tokenizeFormat (const std::string& fmt)
  {
    std::vector<std::string> tokens;
    tokenizeFormatHelper (tokens, fmt, 0, tokens.size());
  }

  void
  tokenizeFormatHelper (std::vector<std::string>& tokens, 
			const std::string& fmt, 
			const size_t start, 
			const size_t end) // exclusive range
  {
    if (start == end || start == std::string::npos)
      return;
    const char curChar = fmt[start];
    if (curChar == ' ' || curChar == '\t')
      tokenizeFormatHelper (tokens, fmt, start+1, end);
    else if 
  }


  bool
  nextTokenIs (const std::string& token, const std::string& target)
  {
    // FIXME This might need revision, depending on the token.
    return token == target;
  }

  template<class ConstIterOfString>
  std::pair<ConstIterOfString, bool>
  parseInteger (ConstIterOfString iter,
		ConstIterOfString end,
		const std::string& key,
		dict_type& dict,
		const bool optional)
  {
    if (optional)
      {
	if (iter == end || ! isInteger (*iter))
	  return std::make_pair (iter, false);
	else
	  return std::make_pair (takeToken (key, dict), true);
      }
    else if (! isInteger (*iter))
      throw std::invalid_argument ("Format format specifier: expected "
				   "integer, but didn't get it.");
    else 
      return std::make_pair (takeToken (key, dict), true);
  }

  template<class ConstIterOfString>
  ConstIterOfString
  parseIntegerFormat (ConstIterOfString iter,
		      ConstIterOfString end,
		      dict_type& dict)
  {
    iter = parseInteger (iter, end, "field-width", dict, false).first;
    if (nextTokenIs (*iter, "."))
      {
	++iter;
	iter = parseInteger (iter, end, "min-num-digits", dict, false).first;
      }
    return iter;
  }

  template<class ConstIterOfString>
  ConstIterOfString
  parseFixedFormat (ConstIterOfString iter,
		    ConstIterOfString end,
		    dict_type& dict)
  {
    iter = parseInteger (iter, end, "field-width", dict, false).first;
    if (nextTokenIs (*iter, "."))
      {
	++iter;
	iter = parseInteger (iter, end, "digits-after-decimal-point", dict, false).first;
      }
    return iter;
  }


  template<class ConstIterOfString>
  ConstIterOfString
  parseFloatFormat (ConstIterOfString iter,
		    ConstIterOfString end,
		    dict_type& dict)
  {
    iter = parseInteger (iter, end, "field-width", dict, false).first;
    if (! nextTokenIs (*iter, "."))
      throw std::invalid_argument("Missing \".\" after <field-width> in floating-point format specifier.");
    else
      {
	++iter; // take the "."
	iter = parseInteger (iter, end, "significand-length", dict, false).first;
	if (nextTokenIs (*iter, "E"))
	  {
	    ++iter; // take the "E"
	    return parseInteger (iter, end, "exponent-length", dict, false).first;
	  }
      }
  }


  template<class ConstIterOfString>
  takeToken (const std::string& key, dict_type& theDict)
  {
    theDict[key] = *iter;
    ++iter;
    return iter;
  }

  template<class ConstIterOfString>
  ConstIterOfString
  parseFormatHelper (ConstIterOfString iter,
		     ConstIterOfString end,
		     dict_type& theDict,
		     const bool tolerant)
  {
    if (iter == end)
      throw std::invalid_argument("Fortran format specifier is empty.");
    // The Fortran format specifier may optionally be enclosed in
    // parentheses.
    else if (nextTokenIs (*iter, "("))
      {
	++iter;
	iter = parseUnenclosedFormat (iter, end, theDict);
	if (iter == end || ! nextTokenIs (*iter, ")"))
	  {
	    if (! tolerant)
	      throw std::invalid_argument ("Fortran format specifier begins "
					   "with an open parenthesis, but does "
					   "not end with a close parenthesis.");
	    // In tolerant mode, just leave any leftover tokens.
	    else
	      return iter;
	  }
      }
    else
      return parseUnenclosedFormat (iter, end, theDict, tolerant);
  }

} // namespace HarwellBoeing

#endif // __HarwellBoeing_hpp
