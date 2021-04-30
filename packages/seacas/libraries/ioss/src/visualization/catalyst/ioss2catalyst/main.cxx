// Copyright(C) 1999-2021 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#include "IossApplication.h"

int main(int argc, char **argv) {
    IossApplication ioapp(argc, argv);
    ioapp.runApplication();
}
