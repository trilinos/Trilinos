cd ../..

# Generate short names
cd ./Headers/
./gen_UseShortNames.sh
cd ..

# Generate forward declarations
cd ./Utils/ForwardDeclaration/
./gen_fwdFiles.sh
cd ../..

# Generate explicit instantiation
cd ./Utils/ExplicitInstantiation/
./gen_cppFiles.sh
cd ../..

# Touch CMakeLists.txt for a new class
new=`git status -s ./Utils/ForwardDeclaration/ ./Utils/ExplicitInstantiation/ | grep "^??" | wc -l`
if [ "$new" != "0" ]; then
  echo "# touch CMakeLists.txt because a new file was created in Utils/ExplicitInstantiation of Utils/ForwardDeclaration" >> CMakeLists.txt
fi
