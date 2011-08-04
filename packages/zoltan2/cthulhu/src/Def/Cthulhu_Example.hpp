#ifndef CTHULHU_EXAMPLES_HPP
#define CTHULHU_EXAMPLES_HPP

#include "Cthulhu_ConfigDefs.hpp"

// The purpose of this file is to allow a user to write program without worrying about templates parameters.
// For example, you can replace:
//
// typedef int LocalOrdinal;
// typedef int GlobalOrdinal;
//
// GlobalOrdinal nx = 4;
// GlobalOrdinal ny = 4;
// RCP<const Map<LocalOrdinal,GlobalOrdinal> > map = MapFactory<LocalOrdinal,GlobalOrdinal>::Build(cthulhuParameters.GetLib(), nx*ny, 0, comm);
//
// by:
//
// #include <Cthulhu_Example.hpp>
//
// GO nx = 4;
// GO ny = 4;
// RCP<const Map> map = MapFactory::Build(cthulhuParameters.GetLib(), nx*ny, 0, comm);
//
// The definition of template types can still be modified if needed.

// Cthulhu_DefaultTypes.hpp defines the types: Scalar, LocalOrdinal, GlobalOrdinal and Node.
// You can use your own definition of these types by defining them before including Cthulhu_Example.hpp. You just have to define also the macro CTHULHU_DEFAULTTYPES_HPP.
#include "Cthulhu_UseDefaultTypes.hpp"

// Cthulhu_UseShortNames.hpp get ride of template type.
// For example, instead of using the type 'Map<LocalOrdinal,GlobalOrdinal,Node>', you can simply use the short type name 'Map'.
#include "Cthulhu_UseShortNames.hpp"

#endif // CTHULHU_EXAMPLES_HPP
