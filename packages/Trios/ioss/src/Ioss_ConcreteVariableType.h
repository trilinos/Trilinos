/*--------------------------------------------------------------------*/
/*    Copyright 2000 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#ifndef IOSS_Ioss_ConcreteVariableType_h
#define IOSS_Ioss_ConcreteVariableType_h

#include <Ioss_CodeTypes.h>
#include <string>

#include <Ioss_VariableType.h>

namespace Ioss {
  class StorageInitializer
    {
    public:
      StorageInitializer();
      // Assignment operator
      // Copy constructor
    };

#define MAKE_CLASS(X) \
class X : public VariableType {\
 public:\
  std::string label(int which, const char suffix_sep='_') const;	\
  static void factory();\
 protected:\
  X();\
 private:\
  X(const X&);\
}

  class Invalid_Storage : public VariableType {
  public:
    std::string label(int which, const char suffix_sep='_') const;
    std::string label_name(const std::string& base, int,
			      const char suffix_sep) const;
    int suffix_count() const {return 0;}
    static void factory();

  protected:
    Invalid_Storage();

  private:
    Invalid_Storage(const Invalid_Storage&);
  };

  class Scalar : public VariableType {
  public:
    std::string label(int which, const char suffix_sep='_') const;
    std::string label_name(const std::string& base, int,
			      const char suffix_sep) const;
    int suffix_count() const {return 0;}
    static void factory();

  protected:
    Scalar();

  private:
    Scalar(const Scalar&);
  };

  MAKE_CLASS(Vector_2D);
  MAKE_CLASS(Vector_3D);
  MAKE_CLASS(Quaternion_2D);
  MAKE_CLASS(Quaternion_3D);
  MAKE_CLASS(Full_Tensor_36);
  MAKE_CLASS(Full_Tensor_32);
  MAKE_CLASS(Full_Tensor_22);
  MAKE_CLASS(Full_Tensor_16);
  MAKE_CLASS(Full_Tensor_12);
  MAKE_CLASS(Sym_Tensor_33);
  MAKE_CLASS(Sym_Tensor_31);
  MAKE_CLASS(Sym_Tensor_21);
  MAKE_CLASS(Sym_Tensor_13);
  MAKE_CLASS(Sym_Tensor_11);
  MAKE_CLASS(Sym_Tensor_10);
  MAKE_CLASS(Asym_Tensor_03);
  MAKE_CLASS(Asym_Tensor_02);
  MAKE_CLASS(Asym_Tensor_01);
  MAKE_CLASS(Matrix_22);
  MAKE_CLASS(Matrix_33);
}
#endif
