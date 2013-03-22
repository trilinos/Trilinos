
echo Running Tpetra/Doxygen ...
cd ../../tpetra/doc
sed -i 's/GENERATE_XML           = NO/GENERATE_XML           = YES/' Doxyfile
doxygen Doxyfile >/dev/null 2>/dev/null
git checkout Doxyfile
cd - >/dev/null

echo Generating Xpetra ...
python interfaces.py
python tpetra.py
python epetra.py

echo Removing trailing whitespace...
find ../src -iname "*.?pp" -exec sed -i 's/[ \t]*$//' {} \;

echo Removing tpetra/doc/xml ...
rm -rf ../../tpetra/doc/xml

# TMP
echo Warning\(TMP\): git checkout Xpetra_TpetraRowMatrix.hpp
git checkout ../src/RowMatrix/Xpetra_TpetraRowMatrix.hpp

echo done.