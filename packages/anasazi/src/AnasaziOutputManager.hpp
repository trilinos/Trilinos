// @HEADER
// *****************************************************************************
//                 Anasazi: Block Eigensolvers Package
//
// Copyright 2004 NTESS and the Anasazi contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ANASAZI_OUTPUT_MANAGER_HPP
#define ANASAZI_OUTPUT_MANAGER_HPP

/*!     \file AnasaziOutputManager.hpp
        \brief Abstract class definition for Anasazi Output Managers.
*/

#include "AnasaziConfigDefs.hpp"
#include "AnasaziTypes.hpp"

#include "Teuchos_FancyOStream.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_oblackholestream.hpp"

/*!  \class Anasazi::OutputManager

  \brief Output managers remove the need for the eigensolver to know any information
  about the required output.  Calling isVerbosity( MsgType type ) informs the solver if
  it is supposed to output the information corresponding to the message type.

  \author Chris Baker, Ulrich Hetmaniuk, Rich Lehoucq, and Heidi Thornquist
*/

namespace Anasazi {

template <class ScalarType>
class OutputManager {

  public:

    //!@name Constructors/Destructor 
  //@{ 

  //! Default constructor
  OutputManager( int vb = Anasazi::Errors,
                 const Teuchos::RCP<Teuchos::FancyOStream> &fos = Teuchos::getFancyOStream(Teuchos::rcpFromRef(std::cout)) )
   : vb_(vb),
     fos_(fos)
  {
    bh_fos_ = Teuchos::getFancyOStream(Teuchos::rcpFromRef( myBHS_ ));
  };

  //! Destructor.
  virtual ~OutputManager() {};
  //@}
  
  //! @name Set/Get methods
  //@{ 

  //! Set the message output types for this manager.
  virtual void setVerbosity( int vb ) { vb_ = vb; }

  //! Get the message output types for this manager.
  virtual int getVerbosity( ) const { return vb_; }

  //! Set the formatted output stream object for this manager.
  virtual void setFancyOStream( const Teuchos::RCP<Teuchos::FancyOStream>& fos ) { fos_ = fos; }
  
  //! Get the formatted output stream object for this manager.
  virtual const Teuchos::RCP<Teuchos::FancyOStream>& getFancyOStream( ) const { return fos_; }

  //@}

  //! @name Output methods
  //@{ 

  //! Find out whether we need to print out information for this message type.
  /*! This method is used by the solver to determine whether computations are
      necessary for this message type.
  */
  virtual bool isVerbosity( MsgType type ) const;

  //! Send output to the output manager.
  virtual void print( MsgType type, const std::string output );

  //! Create a stream for outputting to.
  virtual Teuchos::FancyOStream &stream( MsgType type );

  //@}

  private:

  //! @name Undefined methods
  //@{ 

  //! Copy constructor.
  OutputManager( const OutputManager<ScalarType>& OM );

  //! Assignment operator.
  OutputManager<ScalarType>& operator=( const OutputManager<ScalarType>& OM );

  //@}

  protected:
  int vb_;
  Teuchos::RCP<Teuchos::FancyOStream> fos_, bh_fos_;
  Teuchos::oblackholestream myBHS_;
};

template<class ScalarType>
bool OutputManager<ScalarType>::isVerbosity( MsgType type ) const 
{
  if ( (type & vb_) == type ) {
    return true;
  }
  return false;
}

template<class ScalarType>
void OutputManager<ScalarType>::print( MsgType type, const std::string output ) 
{
  if ( (type & vb_) == type ) {
    *fos_ << output;
  }
}

template<class ScalarType>
Teuchos::FancyOStream & OutputManager<ScalarType>::stream( MsgType type ) {
  if ( (type & vb_) == type ) 
  {
    return *fos_;
  }
  return *bh_fos_;
}

} // end Anasazi namespace

#endif

// end of file AnasaziOutputManager.hpp
