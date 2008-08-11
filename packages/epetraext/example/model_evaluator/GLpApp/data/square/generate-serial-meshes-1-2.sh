#!/bin/sh
#
# This script generates two meshes for input to this driver program.
# In order to run this script you must have the program Triangle installed
# on your system.
#
# This script must be run from this directory to work currectly.
#
# This script should only be run when the square.poly file changes
# and not otherwise.
#

_OFILE=generate-serial-meshes-1-2.sh.out

echo > $_OFILE
echo "*** Generate the first mesh square.1.[poly|node|edge|ele] from this basic geometry square.poly" >> $_OFILE
echo >> $_OFILE
triangle -pe square.poly >> $_OFILE

echo >> $_OFILE
echo "*** Generate the formated input file square.1.000 from square.1 for input to NLPThyraEpetraAdvDiffReactModel.exe" >> $_OFILE
echo >> $_OFILE
../../from-triangle-to-serial-input-mesh.pl --input-base=square.1 >> $_OFILE

echo >> $_OFILE
echo "*** Generate a refinement of square.1.[poly|node|edge|ele] to square.2.[poly|node|edge|ele]" >> $_OFILE
echo >> $_OFILE
triangle -pera0.001 square.1 >> $_OFILE

echo >> $_OFILE
echo "*** Generate the formated input file square.2.000 from square.2 for input to NLPThyraEpetraAdvDiffReactModel.exe" >> $_OFILE
echo >> $_OFILE
../../from-triangle-to-serial-input-mesh.pl --input-base=square.2 >> $_OFILE

cat $_OFILE
