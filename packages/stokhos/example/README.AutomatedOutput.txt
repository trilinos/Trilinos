--------------------------------------------------------------------------
Auto-generation of sample output and incorporation into HTML documentation
--------------------------------------------------------------------------

Here the steps needed to setup creating automated output, adding it to the
source tree, and then incorporating it into doxygen documentation.  To setup
this system for some package (called PACKAGE) follow the instructions belown
which are based on the example in stokhos:

1) Copy generate-sampe-output.pl.stub.in from stokhos/example to
PACKAGE/example and modify it to generate the desired output from your
package.  This will send output to files in the source tree for PACKAGE.

  cd PACKAGE/example;
  cp ../../stokhos/example/generate-sampe-output.pl.stub.in .

2) Add example/generate-sample-output.p.stub to PACKAGE/configure.ac to the
AC_CONFIG_FILES(...) M4 macro.  See stokhos/configure.ac for an example.

3) Add a make rule in PACKAGE/example/Makefile.am to create
generate-sample-output.pl from generate-sample-output.pl.stub and to make it
executable (see stokhos/example/Makefile.am as an example)

4) Bootstrap your PACKAGE

  cd PACKAGE
  ./bootstrap

5) Add an entry in PACKAGE/test/definition to run
../example/generate-sample-output.pl and give it the test category
DOCUMENTATION (see stokhos/test/definition for an example).

6) Update the makefiles in the build directory for your package and rebuild
the package (this will build PACKAGE/example/generate-sample-output.pl)

  cd TRILINOS_BUILD_DIR/packages/PACKAGE
  ./config.status
  make

7) Run test test harness for the package (usually 'make runtests-serial') and
define TRILINOS_TEST_CATEGORY=DOCUMENTATION.  This should run the
generate-sample-output.pl and therefore create output files in the source tree
for your package.

  cd TRILINOS_BUILD_DIR/packages/PACKAGE
  make runtests-serial TRILINOS_TEST_CATEGORY=DOCUMENTATION

8) In your Doxyfile, add directories where your output is written to the
variable EXAMPLE_PATH.  For example, in stokhos, the output is written to
stokhos/example and therefore the path ../example is set in the
EXAMPLE_PATH variable.

9) Put in calls to \verbinclude in your doxygen documentation to pull in the
generated output files.  See the file stokhos/src/Stokhos_Hello.h for an
example.

10) Build the doxygen documentation and verify that the generated output is
showing up in the HTML pages.

  cd stokhos/doc
  ./build_docs

11) Do a CVS add and commit for the output files that were generated.

Now you can automatically update sample output from any machine where your
package is built by simply running the test harness again and checking the new
files back into CVS.

  cd $TRILINOS_BUILD_DIR/packages/PACKAGE
  make runtests-serial TRILINOS_TEST_CATEGORY=DOCUMENTATION
  cd $TRILINOS_HOME/packages/PACKAGE
  cvs commit -m "Updating auto-generated sample output."

THE END
