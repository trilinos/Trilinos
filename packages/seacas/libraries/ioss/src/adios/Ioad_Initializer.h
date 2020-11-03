// Copyright(C) 1999-2020 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#ifndef IOSS_Ioad_Initializer_h
#define IOSS_Ioad_Initializer_h

namespace Ioad {
  /** \brief Initialization of the adios database parts of the Ioss library.
   *
   *  If any input or output type is adios, then an object of this type
   *  must be created before using any other functions or methods in the
   *  Ioss library except Ioss::Init::Initializer().
   */
  class Initializer
  {
  public:
    Initializer();
    ~Initializer();
    // Copy constructor
    // Assignment operator
  private:
    static int useCount;
  };
} // namespace Ioad
#endif
