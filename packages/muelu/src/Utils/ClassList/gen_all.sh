cd ../..

cd ./Headers/
./gen_UseShortNames.sh
cd ..

cd ./Utils/ForwardDeclaration/
./gen_fwdFiles.sh
cd ../..

cd ./Utils/ExplicitInstantiation/
./gen_cppFiles.sh
cd ../..