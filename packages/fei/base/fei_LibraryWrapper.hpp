#ifndef _LibraryWrapper_hpp_
#define _LibraryWrapper_hpp_

/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include <fei_macros.hpp>
#include <fei_SharedPtr.hpp>

#include <fei_LinearSystemCore.hpp>
#include <fei_FiniteElementData.hpp>

class LibraryWrapper {
 public:
  LibraryWrapper(fei::SharedPtr<LinearSystemCore> lsc);
  LibraryWrapper(fei::SharedPtr<FiniteElementData> feData);
  virtual ~LibraryWrapper();

  bool haveLinearSystemCore() { return( haveLinearSystemCore_ ); }
  bool haveFiniteElementData(){ return( haveFiniteElementData_); }

  fei::SharedPtr<LinearSystemCore> getLinearSystemCore() { return( lsc_ ); }
  fei::SharedPtr<FiniteElementData> getFiniteElementData() { return( feData_ ); }

 private:
  bool haveLinearSystemCore_;
  bool haveFiniteElementData_;
  fei::SharedPtr<LinearSystemCore> lsc_;
  fei::SharedPtr<FiniteElementData> feData_;
};

#endif
