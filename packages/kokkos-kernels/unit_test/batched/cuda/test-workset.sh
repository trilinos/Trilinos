exec=$1
testfile=$exec.workset.txt
rm -f $testfile

echo $exec > $testfile

for w in 1024 2048 4096 8192 16384 32768 65536 131072; do
    echo "./$exec -N $w  >> $testfile"
    ./$exec -N $w  >> $testfile
done
