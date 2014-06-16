##---------------------------------------------------------------------------##
## Install nemesis environment
## 
## install.sh <dir> where dir is the directory to install the shell script
##---------------------------------------------------------------------------##

# current directory
current=`dirname $0|pwd`

# environment
dirs="bibtex emacs fig latex templates"

# nemesis python scripts
pybin="nemesis-cc nemesis-class nemesis-fortran nemesis-header nemesis-python nemesis-package nemesis-note"
pybin="$pybin nemesis-test"
pysup="nemesis_files.py"
binfiles="nemesis-browse"

# find pythonpath
PYTHON=`which python`

# install directory
prefix=$1

remove='remove'
# if install directory exists remove it
if test -d "$prefix"; then
    echo ">>> Found $prefix, remove or exit [remove]"
    read remove
    if test "$remove" == 'exit'; then
	echo ">>> Exiting, to finish install, manually remove $prefix, or" 
	echo "    choose a different install location"
	exit 0
    else
	echo ">>> Removing $prefix"
	rm -rf $prefix
    fi
fi

# make the prefix directory
echo ">>> Making $prefix"
mkdir $prefix

# make directories and install files
for d in $dirs; do
    echo ">>> Installing $d in $prefix/$d"
    mkdir $prefix/$d
    cd $d
    files=`ls`
    for f in $files; do
	if test ! -d $f; then
	    cp $f $prefix/$d
	fi
    done
    cd $current
done

# install python binaries
echo ">>> Installing bin in $prefix/bin"
mkdir $prefix/bin
for b in $pybin; do
    sed 's!NEMESIS_PYTHON_DIR!'"$PYTHON"'!g' <bin/$b>$prefix/bin/$b
    chmod 755 $prefix/bin/$b
done

# install python helpers
for s in $pysup; do
    sed 's!NEMESIS_ENVIRONMENT_DIR!'"$prefix"'!g' <bin/$s>$prefix/bin/$s
done

# install other binaries
for f in $binfiles; do
    cp bin/$f $prefix/bin
    chmod 755 $prefix/bin/$f
done

echo
echo "========================================================================"
echo "Congratulations, you're almost finished the nemesis environment install!"
echo
echo "For complete operation you should set the following paths and "
echo "environment variables (depending upon your shell):"
echo 
echo " - Add $prefix/bin to PATH"
echo " - Add $prefix/latex to TEXINPUTS"
echo "   e.g. echo \$TEXINPUTS"
echo "        .:$prefix/latex:"
echo " - Add $prefix/bibtex to BSTINPUTS"
echo "   e.g. echo \$BSTINPUTS"
echo "        .:$prefix/bibtex:"
echo
echo "To use the nemesis GNU Emacs environment add the following to your"
echo ".emacs or .emacs.d/init.el file:"
echo 
echo " (setq nemesis-dir \"$prefix/emacs\")"
echo " (add-to-list 'load-path nemesis-dir)"
echo " (load-library \"nemesis\")"
echo "========================================================================"

##---------------------------------------------------------------------------##
## end of install.sh
##---------------------------------------------------------------------------##

