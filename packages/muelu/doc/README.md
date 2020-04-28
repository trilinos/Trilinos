# Doxygen

`MueLu`'s documentation is built in the build tree:

1. Run your `cmake` configure line or script.
1. In the <buildDirectory>, type `make doc_muelu`.

The documentation will be placed in `<buildDirectory>/packages/muelu/doc/html`.

# User's Guide

The user's guide is in `muelu/doc/UsersGuide`.
To build the guide, simply run `make` in `muelu/doc/UsersGuide`, which should produce `mueluguide.pdf`.
