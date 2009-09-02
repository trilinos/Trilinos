/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#ifndef _fei_ParameterSet_hpp_
#define _fei_ParameterSet_hpp_

#include "fei_macros.hpp"
#include "fei_Param.hpp"

#include <vector>

namespace fei {

/** Container that functions as a database of named parameters, 
intended to store control parameters such as solver tolerances and
other named options.

Individual parameters are instances of fei::Param. This is a
simplistic way of simulating a std::map that can store values of varying
type.

Example creation and initialization:
<pre>
  fei::ParameterSet paramset;

  paramset.add(fei::Param("tolerance", 1.e-8));
  paramset.add(fei::Param("FEI_OUTPUT_PATH", "/home/me/tmp"));
  paramset.add(fei::Param("BLOCK_MATRIX", true));
</pre>

Note that various fei interfaces used to accept parameters in the form
of a length and a list of pointers-to-char-pointer (i.e., 
const char*const* paramStrings). A couple of utility functions exist
for converting to and from fei::ParameterSet. See the following
functions which are declared in fei_utils.hpp:<br>
fei::utils::wrap_strings()<br>
fei::utils::parse_strings()<br>
fei::utils::convert_ParameterSet_to_strings()<br>
fei::utils::strings_to_char_ptrs()
*/
class ParameterSet {
 public:
  /** Constructor */
  ParameterSet();

  /** Destructor */
  virtual ~ParameterSet();

  /** Add a new named parameter object to this parameter-set. If a parameter
      with the same name is already present, it is replaced by this new one
      if the 'maintain_unique_keys' argument is true. Otherwise, the duplicate
      parameter is simply appended to the internal list of parameters.

  @param param Named parameter to be added. A copy of this parameter object
  is stored internally. The caller is free to destroy the input argument as
  soon as this method returns.

  @param maintain_unique_keys Optional argument, defaults to true. If this
  argument is false, then multiple parameters with the same name may be stored. 
  */
  void add(const Param& param,
	   bool maintain_unique_keys=true);

  /** Query for pointer to named parameter. The internal list of named
      parameters is searched using a linear search. If multiple
      parameters have the same name (i.e., the add() method has been used with
      the optional 'maintain_unique_keys' argument specified as false), then
      the parameter returned is an arbitrary one of the duplicates. In other
      words, parameters are not stored in any particular order.
  */
  const Param* get(const char* name) const;

  /** Query for the number of named parameters currently stored. */
  int size() const;

  /** Constant iterator for visiting each named parameter in a parameter-set.*/
  class const_iterator {
  public:
    /** default constructor */
    const_iterator() : paramarray_(NULL),
                       dummyParam_((const char*)NULL, (const char*)NULL),
                       offset_(0) {}
    /** constructor */
    const_iterator(int offset, std::vector<const Param*>* params)
     :  paramarray_(params), dummyParam_((const char*)NULL, (const char*)NULL),
       offset_(offset)
      {
	if (params != NULL) {
	  if (params->empty()) {paramarray_ = NULL; offset_ = 0;}
	}
      }

    /** destructor */
    ~const_iterator(){}

    /** operator* */
    const Param& operator*() const
    {
      if (paramarray_ == NULL) return(dummyParam_);
      if (offset_ >= paramarray_->size()) return(dummyParam_);
      return( *((*paramarray_)[offset_]) );
    }

    /** operator++ */
    const_iterator& operator++()
    {
      if (paramarray_ != NULL) {
        if (offset_ < paramarray_->size()) {
          ++offset_;
        }
        if (offset_ == paramarray_->size()) {
          offset_ = 0;
          paramarray_ = NULL;
        }
      }
      return( *this );
    }

    /** operator= */
    const_iterator& operator=(const const_iterator& src)
    {
      paramarray_ = src.paramarray_;
      offset_ = src.offset_;
      return( *this );
    }

    /** operator== */
    bool operator==(const const_iterator& rhs)
    {
      return( paramarray_ == rhs.paramarray_ && offset_ == rhs.offset_ );
    }

    /** operator!= */
    bool operator!=(const const_iterator& rhs)
    {
      return( !(*this == rhs) );
    }

  private:
    std::vector<const fei::Param*>* paramarray_;
    const fei::Param dummyParam_;
    unsigned offset_;
  };

  /** Return an iterator pointing to the beginning of the list of parameters */
  const_iterator begin() const;

  /** Return an iterator pointing just past the end of the list of parameters */
  const_iterator end() const;

  /** Check whether a fei::Param with 'name' is present and has type
      fei::Param::INT. If so, set paramValue and return 0. Otherwise return -1.
  */
  int getIntParamValue(const char* name,
		       int& paramValue) const;

  /** Check whether a fei::Param with 'name' is present and has type
      fei::Param::INT or fei::Param::DOUBLE.
      If so, set paramValue and return 0. Otherwise return -1.
  */
  int getDoubleParamValue(const char* name,
			  double& paramValue) const;

  /** Check whether a fei::Param with 'name' is present and has type
      fei::Param::STRING. If so, set paramValue and return 0. Otherwise return -1.
  */
  int getStringParamValue(const char* name,
			  std::string& paramValue) const;

  /** Check whether a fei::Param with 'name' is present and has type
      fei::Param::BOOL. If so, set paramValue and return 0. Otherwise return -1.
  */
  int getBoolParamValue(const char* name,
			bool& paramValue) const;

  /** Check whether a fei::Param with 'name' is present and has type
      fei::Param::VOID. If so, set paramValue and return 0. Otherwise return -1.
  */
  int getVoidParamValue(const char* name,
			const void*& paramValue) const;

 private:
  int findOffset(const fei::Param* param) const;
  int findOffset(const char* name) const;
  std::vector<const Param*>* params_;
};//class ParameterSet

}//namespace fei

inline
int fei::ParameterSet::findOffset(const char* name) const
{
  if (params_->empty() || name == NULL) return( -1 );

  std::vector<const Param*>::const_iterator
    p_iter = params_->begin(),
    p_end = params_->end();

  int i = 0;
  for(; p_iter != p_end; ++p_iter, ++i) {
    const Param* prm = *p_iter;
    if (prm->getName() == name) {
      return(i);
    }
  }
  return(-1);
}

inline
int fei::ParameterSet::findOffset(const fei::Param* param) const
{
  if (param == NULL) return( -1 );

  return( findOffset( param->getName().c_str() ) );
}

inline
fei::ParameterSet::const_iterator fei::ParameterSet::begin() const
{
  return( const_iterator(0, params_) );
}

inline
fei::ParameterSet::const_iterator fei::ParameterSet::end() const
{
  return( const_iterator(0, NULL) );
}

inline
void fei::ParameterSet::add(const fei::Param& param, bool maintain_unique_keys)
{
  int index = findOffset(&param);
  const fei::Param* newparam = new fei::Param(param);
  if (index < 0) {
    params_->push_back(newparam);
  }
  else {
    if (maintain_unique_keys) {
      delete (*params_)[index];
      (*params_)[index] = newparam;
    }
    else {
      params_->push_back(newparam);
    }
  }
}

inline
int fei::ParameterSet::size() const
{
  return( params_->size() );
}

inline
const fei::Param* fei::ParameterSet::get(const char* name) const
{
  if (params_ == NULL) return(NULL);

  int index = findOffset(name);
  if (index < 0) return(NULL);
  return( (*params_)[index] );
}

#endif
