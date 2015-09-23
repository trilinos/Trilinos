#!/bin/sh

rm -f Doxyfile.input
rm -f index.dox

echo "/** \mainpage  Sierra Toolkit" > index.dox
echo " *" >> index.dox

for i in ../../stk_*/dox/Doxyfile.input; do 
  echo @INCLUDE = $i >> Doxyfile.input
done 

echo " * @section stk_products The Sierra Toolkit contains the following products:" >> index.dox
echo " *" >> index.dox
for i in ../../stk_*/dox/Doxyfile.input; do 
  gawk -- '/#BANNER/{print " * - "substr($0, 8)}' $i >> index.dox
  echo " *" >> index.dox
done 
echo " *" >> index.dox

echo " * @section stk_howto How to..." >> index.dox
echo " *" >> index.dox
for i in ../../stk_*/dox/Doxyfile.input; do 
  echo @INCLUDE = $i >> Doxyfile.input
  gawk -- '/#HOWTO/{print " * - "substr($0, 7)}' $i >> index.dox
  echo " *" >> index.dox
done 
echo " *" >> index.dox

echo " */" >> index.dox

doxygen

# Hello Carol
