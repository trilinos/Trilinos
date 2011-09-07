// @HEADER
// ***********************************************************************
// 
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
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

#ifndef THYRA_PARAMETER_DRIVEN_MULTI_VECTOR_INPUT_HPP
#define THYRA_PARAMETER_DRIVEN_MULTI_VECTOR_INPUT_HPP

#include "Thyra_MultiVectorFileIOBase.hpp"
#include "Thyra_DetachedVectorView.hpp"
#include "Thyra_MultiVectorStdOps.hpp"
#include "Teuchos_ParameterListAcceptor.hpp"
#include "Teuchos_VerboseObject.hpp"
#include "Teuchos_StandardCompositionMacros.hpp"
#include "Teuchos_implicit_cast.hpp"

namespace Thyra {

/** \brief Concrete utility class that an ANA can use for reading in a
 * (multi)vector as directed by a parameter sublist.
 *
 * This class is made one-hundred percent general for all ANAs by accepting an
 * abstract <tt>MultiVectorFileIOBase</tt> object that actually reads the
 * (multi)vectors from file(s) (or whatever implementation is used for the
 * storage of the (multi)vectors).
 *
 * This class can also read in small (multi)vectors directly from the
 * parameter sublist (see <tt>getValidParameters()</tt>).
 *
 * Note that in order to use objects of this type, the client must minimally
 * set the vector space using <tt>set_vecSpc()</tt> before any vectors can be
 * extracted and if file IO is performed, then a
 * <tt>MultiVectorFileIOBase</tt> object must be set using
 * <tt>set_fileIO()</tt>.  Note that the parameter sublist can be set using
 * <tt>setParameterList()</tt> without first setting any of these two objects.
 * In other words, this object can accept a parameter list without knowing
 * anything about the actually type of the (multi)vector that will be read
 * using the <tt>readMultiVector()</tt> or <tt>readVector()</tt> functions
 * called later.
 *
 * This simple utility class is not meant to be too fancy so please just study
 * the actual implementation to see what this class does.
 *
 * ToDo: When needed, implement a function that can read a multi-vector of any
 * number of columns based on the number of columns stored in the file.
 *
 * ToDo: When needed, add the ability to read in a multi-vector with more than
 * one column directly from the array parameter.
 *
 * \ingroup Thyra_Op_Vec_ANA_Development_grp
 */
template<class Scalar>
class ParameterDrivenMultiVectorInput
  : public Teuchos::ParameterListAcceptor
  , public Teuchos::VerboseObject<ParameterDrivenMultiVectorInput<Scalar> >
{
public:

  /** \name Constructors/Initializers */
  //@{

  /** \brief . */
  ParameterDrivenMultiVectorInput();

  /** \brief Set the vector space used to create the (multi)vectors that are
   * read in.
   */
  STANDARD_CONST_COMPOSITION_MEMBERS( VectorSpaceBase<Scalar>, vecSpc );

  /** \brief Set the MultiVectorFileIOBase object that will be used to read
   * the vector from file(s).
   */
  STANDARD_COMPOSITION_MEMBERS( MultiVectorFileIOBase<Scalar>, fileIO );

  //@}

  /** \name Overridden from ParameterListAcceptor */
  //@{

  /** \brief . */
  void setParameterList(Teuchos::RCP<Teuchos::ParameterList> const& paramList);
  /** \brief . */
  Teuchos::RCP<Teuchos::ParameterList> getNonconstParameterList();
  /** \brief . */
  Teuchos::RCP<Teuchos::ParameterList> unsetParameterList();
  /** \brief . */
  Teuchos::RCP<const Teuchos::ParameterList> getParameterList() const;
  /** \brief . */
  Teuchos::RCP<const Teuchos::ParameterList> getValidParameters() const;

  //@}

  /** \name Informational */
  //@{

  /** \brief Return the value of the parameter "File Name Base" that was read in
   * from the <tt>setParameterList()</tt> function.
   */
  const std::string& readinFileNameBase() const;

  /** \brief Return the value of the parameter "Explicit Array" that was read
   * in from the <tt>setParameterList()</tt> function.
   */
  const Teuchos::Array<Scalar>& readinExplicitArray() const;

  /** \brief Return the value of the parameter "Explicit Array" that was read
   * in from the <tt>setParameterList()</tt> function.
   */
  Scalar readinScaleBy() const;

  //@}

  /** \name (Multi)Vector Readers */
  //@{
  
  /** \brief Read a MultiVector that has already been allocated, as directed
   * by the set parameter sublist.
   *
   * \param  mvName
   *           [in] The name of the multi-vector being read.  This string name
   *           is added to any exception messages that are thrown and is printed
   *           to <tt>getDefaultOStream()</tt> if verbosity is turned on.
   * \param  mv
   *           [in/out] The multi-vector to be read.
   *
   * <b>Preconditions:</b><ul>
   * <li><tt>this->getParameterList().get()!=NULL</tt>
   * <li><tt>this->get_vecSpc().get()!=NULL</tt>
   * <li>[<tt>this->readinFileNameBase().length() > 0</tt>]
   *     <tt>this->get_fileIO().get()!=NULL</tt>
   * <li><tt>mv!=NULL</tt>
   * <li><tt>this->vecSpc().isCompatible(*mv->range())==true</tt>
   * </ul>
   *
   * \returns <tt>true</tt> if a vector was read, <tt>false</tt> otherwise.
   */
  bool readMultiVector(
    const std::string                 &mvName
    ,Thyra::MultiVectorBase<Scalar>   *mv
    ) const;
  
  /** \brief Read a Vector as directed by the set parameter sublist,
   * allocating the Vector object if it has not already been allocated.
   *
   * \param  vName
   *           [in] The name of the vector being read.  This string name
   *           is added to any exception messages that are thrown and is printed
   *           to <tt>getDefaultOStream()</tt> if verbosity is turned on.
   * \param  v
   *           [in/out] RCP to the vector to be read.  If
   *           <tt>v->get()==NULL</tt> before this function is called, then a
   *           new vector object will be allocated and set on return.
   *           On return <tt>*(*v)</tt> will contain the values read as
   *           specified by the parameter sublist.
   *
   * <b>Preconditions:</b><ul>
   * <li><tt>this->getParameterList().get()!=NULL</tt>
   * <li><tt>this->get_vecSpc().get()!=NULL</tt>
   * <li>[<tt>this->readinFileNameBase().length() > 0</tt>]
   *     <tt>this->get_fileIO().get()!=NULL</tt>
   * <li><tt>v!=NULL</tt>
   * <li>[<tt>v->get()!=NULL</tt>]
   *     <tt>this->vecSpc().isCompatible(*(*v)->space())==true</tt>
   * </ul>
   *
   * <b>Postconditions:</b><ul>
   * <li>If <tt>v->get()==NULL</tt> on input no vector was read,
   *     then <tt>v->get()==NULL</tt> on output.
   * <li>If <tt>v->get()==NULL</tt> on input a vector was read,
   *     then <tt>v->get()!=NULL</tt> on output.
   * </ul>
   *
   * \returns <tt>true</tt> if a vector was read, <tt>false</tt> otherwise.
   *
   * This function helps to avoid reallocations on multiple reads so that the
   * same memory can be used over and over again.
   *
   * This function simply calls the above <tt>readMultiVector()</tt> function
   * but it allocates a <tt>Vector</tt> object instead of a multi-vector
   * object.
   */
  bool readVector(
    const std::string &vName
    ,Teuchos::RCP<Thyra::VectorBase<Scalar> > *v
    ) const;
  
  /** \brief Read a newly allocated Vector as directed by the set parameter
   * sublist.
   *
   * This function just calls the above <tt>readVector()</tt> function (see
   * its preconditions).
   *
   * \returns <tt>returnVal.get()!=NULL</tt> if a vector was read and
   * <tt>returnVal.get()==NULL</tt> if no vector was read.
   */
  Teuchos::RCP<Thyra::VectorBase<Scalar> >
  readVector( const std::string &vName ) const;
  
  //@}

private:

  mutable Teuchos::RCP<const Teuchos::ParameterList> validParamList_;
  Teuchos::RCP<Teuchos::ParameterList> paramList_;

  std::string fileNameBase_;
  Teuchos::Array<Scalar> explicitArray_;
  Scalar scaleBy_;
  Scalar addScalar_;

  static const std::string FileNameBase_name_;
  static const std::string FileNameBase_default_;

  static const std::string ExplicitArray_name_;
  static const std::string ExplicitArray_default_;

  static const std::string ScaleBy_name_;
  static const double ScaleBy_default_;

  static const std::string AddScalar_name_;
  static const double AddScalar_default_;

};


/** \brief Read a vector and override if one is read.
 *
 * \relates ParameterDrivenMultiVectorInput
 */
template<class Scalar>
RCP<const VectorBase<Scalar> >
readVectorOverride(
  const ParameterDrivenMultiVectorInput<Scalar> &pdmvi,
  const std::string &vName,
  const RCP<const VectorBase<Scalar> > &defaultVector
  )
{
  RCP<const VectorBase<Scalar> >
    vector = pdmvi.readVector(vName);
  if (!is_null(vector))
    return vector;
  return defaultVector;
}


// //////////////////////////////////////////
// Inline functions

template<class Scalar>
inline
const std::string&
ParameterDrivenMultiVectorInput<Scalar>::readinFileNameBase() const
{
  return fileNameBase_;
}

template<class Scalar>
inline
const Teuchos::Array<Scalar>&
ParameterDrivenMultiVectorInput<Scalar>::readinExplicitArray() const
{
  return explicitArray_;
}

template<class Scalar>
inline
Scalar
ParameterDrivenMultiVectorInput<Scalar>::readinScaleBy() const
{
  return scaleBy_;
}

// //////////////////////////////////////////
// Implementations

namespace PDMVIUtilityPack {

template<class Scalar>
void copy(
  const Teuchos::Array<Scalar>  &array
  ,VectorBase<Scalar>           *vec
  )
{
  using Teuchos::implicit_cast;
  TEST_FOR_EXCEPT(vec==0);
  DetachedVectorView<Scalar> dVec(*vec);
  TEST_FOR_EXCEPT(implicit_cast<int>(dVec.subDim())!=implicit_cast<int>(array.size())); // ToDo: Give a very good error message!
  for( Ordinal i = 0; i < dVec.subDim(); ++i ) {
    dVec[i] = array[i];
  }
}

} // namespace PDMVIUtilityPack

// Static data members

template<class Scalar>
const std::string
ParameterDrivenMultiVectorInput<Scalar>::FileNameBase_name_ = "File Name Base";
template<class Scalar>
const std::string
ParameterDrivenMultiVectorInput<Scalar>::FileNameBase_default_ = "";

template<class Scalar>
const std::string
ParameterDrivenMultiVectorInput<Scalar>::ExplicitArray_name_ = "Explicit Array";
template<class Scalar>
const std::string
ParameterDrivenMultiVectorInput<Scalar>::ExplicitArray_default_ = "{}";

template<class Scalar>
const std::string
ParameterDrivenMultiVectorInput<Scalar>::ScaleBy_name_ = "Scale By";
template<class Scalar>
const double
ParameterDrivenMultiVectorInput<Scalar>::ScaleBy_default_ = 1.0;

template<class Scalar>
const std::string
ParameterDrivenMultiVectorInput<Scalar>::AddScalar_name_ = "Add Scalar";
template<class Scalar>
const double
ParameterDrivenMultiVectorInput<Scalar>::AddScalar_default_ = 0.0;

// Constructors/Initializers

template<class Scalar>
ParameterDrivenMultiVectorInput<Scalar>::ParameterDrivenMultiVectorInput()
  :fileNameBase_(FileNameBase_default_),
   scaleBy_(ScaleBy_default_),
   addScalar_(AddScalar_default_)
{}

// Overridden from ParameterListAcceptor

template<class Scalar>
void ParameterDrivenMultiVectorInput<Scalar>::setParameterList(
  Teuchos::RCP<Teuchos::ParameterList> const& paramList
  )
{
  TEST_FOR_EXCEPT(0==paramList.get());
  paramList->validateParameters(*getValidParameters());
  paramList_ = paramList;
  fileNameBase_ = paramList_->get(
    FileNameBase_name_,FileNameBase_default_ );
  explicitArray_ = Teuchos::getArrayFromStringParameter<Scalar>(
    *paramList_,ExplicitArray_name_
    ,-1     // An array of any size will do here
    ,false  // The parameter does not need to exist
    );
  scaleBy_ = paramList_->get(ScaleBy_name_,ScaleBy_default_);
  addScalar_ = paramList_->get(AddScalar_name_,AddScalar_default_);
#ifdef TEUCHOS_DEBUG
  paramList_->validateParameters(*getValidParameters(),0);
#endif // TEUCHOS_DEBUG
}

template<class Scalar>
Teuchos::RCP<Teuchos::ParameterList>
ParameterDrivenMultiVectorInput<Scalar>::getNonconstParameterList()
{
  return paramList_;
}

template<class Scalar>
Teuchos::RCP<Teuchos::ParameterList>
ParameterDrivenMultiVectorInput<Scalar>::unsetParameterList()
{
  Teuchos::RCP<Teuchos::ParameterList>
    _paramList = paramList_;
  paramList_ = Teuchos::null;
  return _paramList;
}

template<class Scalar>
Teuchos::RCP<const Teuchos::ParameterList>
ParameterDrivenMultiVectorInput<Scalar>::getParameterList() const
{
  return paramList_;
}

template<class Scalar>
Teuchos::RCP<const Teuchos::ParameterList>
ParameterDrivenMultiVectorInput<Scalar>::getValidParameters() const
{
  if(!validParamList_.get()) {
    Teuchos::RCP<Teuchos::ParameterList>
      pl = Teuchos::rcp(new Teuchos::ParameterList);
    pl->set(
      FileNameBase_name_,FileNameBase_default_
      ,"Base-name of file(s) that will be used to read in the vector.\n"
      "If this parameter is empty \"\", no file(s) will be read.\n"
      "Note that a MultiVectorFileIOBase object and a VectorSpaceBase object\n"
      "must be set internally for this to work."
      );
    pl->set(
      ExplicitArray_name_,ExplicitArray_default_
      ,"The vector specified explicitly as a string interpreted as a Teuchos::Array\n"
      "object.  If this array is set, it will override the vector specified\n"
      "by the above \"" + FileNameBase_name_ + "\" parameter.\n"
      "Note that a VectorSpaceBase object\n"
      "must be set internally for this to work."
      );
    pl->set(
      ScaleBy_name_,ScaleBy_default_,
      "A factor by which the read in vector will be scaled by."
      );
    pl->set(
      AddScalar_name_, AddScalar_default_,
      "A scalar that will added to the read in vector after it\n"
      "optionally scaled."
      );
    validParamList_ = pl;
  }
  return validParamList_;
}

// (Multi)Vector Readers

template<class Scalar>
bool ParameterDrivenMultiVectorInput<Scalar>::readMultiVector(
  const std::string                 &mvName
  ,Thyra::MultiVectorBase<Scalar>   *mv
  ) const
{
  using Teuchos::implicit_cast;
  TEST_FOR_EXCEPT(0==mv);
  typedef Teuchos::ScalarTraits<Scalar> ST;
  const Teuchos::EVerbosityLevel verbLevel = this->getVerbLevel();
  Teuchos::RCP<Teuchos::FancyOStream>
    out = this->getOStream();
  const bool trace = ( verbLevel >= implicit_cast<int>(Teuchos::VERB_LOW) );
  Teuchos::OSTab tab(out);
  bool vectorWasRead = false;
  if(fileNameBase_.length()) {
    if( out.get() && trace )
      *out << "\nReading \"" << mvName << "\" from the file(s) with base name \""
           << fileNameBase_ << "\" ...\n";
    fileIO().readMultiVectorFromFile(fileNameBase_,mv);
    vectorWasRead = true;
  }
  if(explicitArray_.size()) {
    if( implicit_cast<Ordinal>(explicitArray_.size()) != vecSpc().dim() ) {
      // Call back to throw an exception with a better erro message!
      Teuchos::getArrayFromStringParameter<Scalar>(
        *paramList_,ExplicitArray_name_,vecSpc().dim(),false);
      TEST_FOR_EXCEPT(!"Should never get here!");
    }
    if( out.get() && trace )
      *out << "\nSetting \"" << mvName << "\" directly from the parameter array "
           << explicitArray_ << " ...\n";
    TEST_FOR_EXCEPTION(
      mv->domain()->dim()!=implicit_cast<Ordinal>(1), std::logic_error
      ,"Error! We can not handle reading in multi-vectors directly from"
      " the parameter list yet!"
      );
    PDMVIUtilityPack::copy(explicitArray_,&*mv->col(0));
    // ToDo: Find a way to read a matrix from a file (perhaps a nested
    // Array<Array<Scalar> > or something!)
    vectorWasRead = true;
  }
  if( scaleBy_ != ST::one() && vectorWasRead ) {
    if( out.get() && trace )
      *out << "\nScaling \"" << mvName << "\" by " << scaleBy_ << " ...\n";
    Vt_S(&*mv,scaleBy_);
  }
  if( addScalar_ != ST::zero() && vectorWasRead ) {
    if( out.get() && trace )
      *out << "\nAdding scalar " << addScalar_ << " to \"" << mvName << "\" ...\n";
    Vp_S(&*mv,addScalar_);
  }
  return vectorWasRead;
}

template<class Scalar>
bool ParameterDrivenMultiVectorInput<Scalar>::readVector(
  const std::string                                   &vName
  ,Teuchos::RCP<Thyra::VectorBase<Scalar> >   *v
  ) const
{
  TEST_FOR_EXCEPT(0==v);
  bool vectorWasRead = false;
  if( fileNameBase_.length() || explicitArray_.size() ) {
    if(!(*v).get())
      (*v) = createMember(this->vecSpc());
    vectorWasRead = this->readMultiVector(vName,&*(*v));
  }
  return vectorWasRead;
}

template<class Scalar>
Teuchos::RCP<Thyra::VectorBase<Scalar> >
ParameterDrivenMultiVectorInput<Scalar>::readVector(
  const std::string &vName
  ) const
{
  Teuchos::RCP<Thyra::VectorBase<Scalar> > v;
  const bool vectorWasRead = readVector(vName,&v);
  if(!vectorWasRead)
    v = Teuchos::null;
  return v;
}

} // namespace Thyra

#endif // THYRA_PARAMETER_DRIVEN_MULTI_VECTOR_INPUT_HPP
