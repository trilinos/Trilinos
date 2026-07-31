# Doxygen

`MiniTensor`'s documentation is built in the build tree:

1. Run your `cmake` configure line or script.
1. In the <buildDirectory>, type `make doc_minitensor`.

The documentation will be placed in `<buildDirectory>/packages/minitensor/doc/html`.

The documentation can also be built without configuring, directly in
the source tree, by running `./build_docs` in this directory; the
HTML is then placed in `doc/html`.

## Publishing to the Trilinos website

The Trilinos website serves package documentation from the
[trilinos.github.io](https://github.com/trilinos/trilinos.github.io)
repository: the generated HTML for a package `<pkg>` lives under
`docs/<pkg>/` there and is served at
`https://trilinos.github.io/docs/<pkg>/index.html`.

Publishing is automated: the `documentation` workflow in
`trilinos.github.io` runs weekly (and on demand), executes every
package's `doc/build_docs` script from the Trilinos `develop` branch
via `doc/build_docs.pl`, copies each resulting `doc/html` into
`docs/<pkg>/`, and opens a pull request with the refreshed pages.
The `build_docs` script in this directory hooks MiniTensor into that
workflow; no manual publishing steps are required.
