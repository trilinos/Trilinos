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

// //////////////////////////////////////////////////
// Teuchos_CommandLineProcessor.hpp

#ifndef TEUCHOS_COMMAND_LINE_PROCESSOR_HPP
#define TEUCHOS_COMMAND_LINE_PROCESSOR_HPP

/*! \file Teuchos_CommandLineProcessor.hpp
	\brief Basic command line parser for input from <tt>(argc,argv[])</tt> 
*/
#include "Teuchos_map.hpp"

namespace Teuchos {

///
/** \class CommandLineProcessor
 *
 *   \brief Class that helps parse command line input arguments from <tt>(argc,argv[])</tt> and set options.
 *
 * This class will process command-line arguments in the form of <tt>(argc,argv[])</tt>
 * and set user-defined options.  This class can also work in a number of modes.
 * This processor can require that all options be recognized or not.
 *
*/

/* ToDo: RAB: 2003/10/02: Add support for required options as per KRL's suggestion
 *
 * ToDo:  Finish documentation.
 */

class CommandLineProcessor {
public:

	/** \enum EParseCommandLineReturn
            \brief Return value for <tt>CommandLineProcessor::parse()</tt>.
	 */
	enum EParseCommandLineReturn {
		PARSE_SUCCESSFUL              =  0 /*!< Parsing the command line was successful. */
		,PARSE_HELP_PRINTED            =  1 /*!< The help statement was printed for the command line parser. */
		,PARSE_UNRECOGNIZED_OPTION     = -1 /*!< The command line parser encountered an unrecognized option. */
	};

	///
	//@{ \name Constructors
	/** \brief Default Constructor
	 *
	 * @param  throwExceptions
	 *               [in] If <tt>true</tt> then <tt>this->parse()</tt> with throw
	 *               exceptions instead of returning <tt>!=PARSE_SUCCESSFUL</tt>.
	 * @param  recogniseAllOptions
	 *               [in] If <tt>true</tt> then <tt>this->parse()</tt> with simply
	 *               ignore options that it does not recognize.
	 */
	CommandLineProcessor(
		bool    throwExceptions      = true
		,bool   recogniseAllOptions  = true
		);

	//@}

	//@{ \name Set up options

	///
	/** \brief Set a boolean option.
	 *
	 * @param  option_true    [in] (null terminated string) If this option is found then
	 *                        <tt>*option_val = true</tt> will be set.
	 * @param  option_false   [in] (null terminated string) If this option is found then
	 *                        <tt>*option_val = false</tt> will be set.
	 * @param  option_val     [in/out] On input, <tt>*option_val</tt> gives the default value
	 *                        of the option (used for printing in --help).  On output,
	 *                        will be set according to <tt>(argc,argv[])</tt>.
	 * @param  documentation  [in] If <tt>!=NULL</tt>, then this null terminated string
	 *                        gives the documentation for the option.
	 */
	void setOption(
		const char     option_true[]
		,const char    option_false[]
		,bool          *option_val
		,const char    documentation[] = NULL
		);

	///
	/** \brief Set an integer option.
	 *
	 * @param  option_name    [in] (null terminated string) The name of the option
	 *                        (without the leading '--' or trailing '=').
	 * @param  option_val     [in/out] On input, <tt>*option_val</tt> gives the default value
	 *                        of the option (used for printing in --help).  On output,
	 *                        will be set according to <tt>(argc,argv[])</tt>.
	 * @param  documentation  [in] If <tt>!=NULL</tt>, then this null terminated string
	 *                        gives the documentation for the option.
	 */
	void setOption(
		const char     option_name[]
		,int           *option_val
		,const char    documentation[] = NULL
		);

	///
	/** \brief Set a floating-point option.
	 *
	 * @param  option_name    [in] (null terminated string) The name of the option
	 *                        (without the leading '--' or trailing '=').
	 * @param  option_val     [in/out] On input, <tt>*option_val</tt> gives the default value
	 *                        of the option (used for printing in --help).  On output,
	 *                        will be set according to <tt>(argc,argv[])</tt>.
	 * @param  documentation  [in] If <tt>!=NULL</tt>, then this null terminated string
	 *                        gives the documentation for the option.
	 */
	void setOption(
		const char     option_name[]
		,double        *option_val
		,const char    documentation[] = NULL
		);

	///
	/** \brief Set a string option.
	 *
	 * @param  option_name    [in] (null terminated string) The name of the option
	 *                        (without the leading '--' or trailing '=').
	 * @param  option_val     [in/out] On input, <tt>*option_val</tt> gives the default value
	 *                        of the option (used for printing in --help).  On output,
	 *                        will be set according to <tt>(argc,argv[])</tt>.
	 * @param  documentation  [in] If <tt>!=NULL</tt>, then this null terminated string
	 *                        gives the documentation for the option.
	 */
	void setOption(
		const char     option_name[]
		,std::string   *option_val
		,const char    documentation[] = NULL
		);

	//@}

	//@{ \name Parse methods

	///
	/** \brief Parse a command line.
	 *
	 * @param  argc    [in] number of entries in argv[]
	 * @param  argv    [in/out] array (length argc) of command line arguments.
	 *                 argv[0] should be the name of the program on the shell as usual.
	 * @param  errout  [out] If <tt>!=NULL</tt> then error and help messages are sent here.
	 *                 The default is set to <tt>&std::cerr</tt>.
	 *
	 * Postconditions:<ul>
	 * <li>If an unrecognized option if found
	 *     <ul>
	 *     <li>If <tt>this->throwExceptions()==true</tt>
	 *         <ul><li>This method will throw an <tt>UnrecognizedOption</tt> exception</ul>
	 *     <li>else
	 *         <ul><li>This method will return <tt>PARSE_UNRECOGNIZED_OPTION</tt></ul>
	 *     <li>endif
	 *     </ul>
	 * <li>Else if the option <tt>--help</tt> is found
	 *     <ul>
	 *     <li>If <tt>this->throwExceptions()==true</tt>
	 *         <ul><li>This method will throw a <tt>HelpPrinted</tt> exception</ul>
	 *     <li>else
	 *         <ul><li>This method will return <tt>PARSE_HELP_PRINTED</tt></ul>
	 *     <li>endif
	 *     </ul>
	 * <li>Else
	 *     <ul><li>This method will return <tt>PARSE_SUCCESSFUL</tt></ul>
	 * </ul>
	 */
	EParseCommandLineReturn  parse(
		int             argc
		,char*          argv[]
		,std::ostream   *errout    = &std::cerr
		) const;

	//@}

	//@{ \name Behavior modes

	/// Set if an exception is thrown, there is a parse error, or help is printed.
	void throwExceptions ( const bool & throwExceptions ) { throwExceptions_ = throwExceptions; };
	
	/// Return true if an exception is thrown, there is a parse error, or help is printed.
	const bool& throwExceptions() const { return throwExceptions_; };

	/// Set if all options must be recognized or not.
	void recogniseAllOptions ( const bool & recogniseAllOptions ) { recogniseAllOptions_ = recogniseAllOptions; };

	/// Return true if all options are being recognized by the parser.
	const bool& recogniseAllOptions() const { return recogniseAllOptions_; };

	//@}

	//@{ \name Exception classes

	/// Thrown if a parse exception occurs and  throwExceptions==true
	class ParseError : public std::logic_error
	{public: ParseError(const std::string& what_arg) : std::logic_error(what_arg) {}};

	/// Thrown if --help was specified and throwExceptions==true
	class HelpPrinted : public ParseError
	{public: HelpPrinted(const std::string& what_arg) : ParseError(what_arg) {}};

	/// Thrown if an unrecognized option was found and throwExceptions==true
	class UnrecognizedOption : public ParseError
	{public: UnrecognizedOption(const std::string& what_arg) : ParseError(what_arg) {}};

	//@}

public:
	//
	enum EOptType { OPT_NONE, OPT_BOOL_TRUE, OPT_BOOL_FALSE, OPT_INT, OPT_DOUBLE, OPT_STRING };
	// RAB: 2003/10/10: Note: I had to move this out of the private section since
	// the sun compiler (version 7) complained (rightly it now appears after looking
	// up what the ISO/ANSI C++ standard says) about the declaration for opt_val_val_t
	// not being able to access a private member of CommandLineProcessor.

private:

	// /////////////////////////////////
	// Private types

	//
	struct opt_val_val_t {
		opt_val_val_t()
			:opt_type(OPT_NONE),opt_val(NULL)
			{}
		opt_val_val_t( EOptType opt_type_in, void* opt_val_in )
			:opt_type(opt_type_in),opt_val(opt_val_in)
			{}
		EOptType     opt_type;
		void         *opt_val; // Will be bool*, int* or double*
	};

	//
	typedef Teuchos::map<std::string,opt_val_val_t>   options_list_t;

	//
	struct opt_doc_t {
		opt_doc_t()
			:opt_type(OPT_NONE)
			{}
		opt_doc_t(EOptType opt_type_in, const std::string& opt_name_in, const std::string& opt_name_false_in
					  ,const std::string &documentation_in, void* default_val_in )
			:opt_type(opt_type_in),opt_name(opt_name_in),opt_name_false(opt_name_false_in)
			,documentation(documentation_in),default_val(default_val_in)
			{}
		EOptType     opt_type;
		std::string  opt_name;
		std::string  opt_name_false; // only for bool
		std::string  documentation;
		void         *default_val;
	};
	
	//
	typedef std::vector<opt_doc_t>   options_documentation_list_t;

	// /////////////////////////////////
	// Private data members

	bool 				throwExceptions_;
	bool				recogniseAllOptions_;
	options_list_t                   options_list_;
	options_documentation_list_t     options_documentation_list_;

	// /////////////////////////////////
	// Private member functions

	// Get the option and the value from an entry in argv[].
	// Will return false if entry is not formated properly.
	bool get_opt_val(
		const char     str[]
		,std::string   *opt_name
		,std::string   *opt_val_str // May be empty on return
		) const;

	// Print help message
	void print_help_msg(
		int             argc
		,char*          argv[]
		,std::ostream   *errout
		) const;

	// String for option type
	std::string opt_type_str( EOptType ) const;

	// Print bad option
	void print_bad_opt(
		int             argv_i
		,char*          argv[]
		,std::ostream   *errout
		) const;
	

}; // end class CommandLineProcessor

// /////////////////////////
// Inline members

inline
std::string CommandLineProcessor::opt_type_str( EOptType opt_type ) const
{
	std::string str;
	switch( opt_type ) {
		case OPT_BOOL_TRUE:
			str = "bool";
			break;
		case OPT_INT:
			str = "int";
			break;
		case OPT_DOUBLE:
			str = "double";
			break;
		case OPT_STRING:
			str = "string";
			break;
		default:
			assert(0); // Local programming error only
	} 
	return str;
}

} // end namespace Teuchos

#endif // TEUCHOS_COMMAND_LINE_PROCESSOR_HPP
