#!/bin/bash

#
# The purpose of this script is to create a compressed archive for
# PALMapper 
#

set -e

name=palmapper-0.5
tmp_dir=/tmp/$name
result_fn=/tmp/$name.tar.gz

# first we create a temporary directory to store the needed files:
mkdir $tmp_dir

rsync -a . $tmp_dir
cd $tmp_dir

# Clean palmapper
make clean

#Clean testcases
cd testcase
make clean
cd data
if [ ! -f "split_1m.000.bz2" -a -f "split_1m.000" ]
then
	bzip2 split_1m.000
fi
cd $tmp_dir

#clean doc
rm -f doc/palmapper-manual.{aux,log}
#Remove this script
rm -f makePackage.sh

#remove settings
rm -rf .settings
rm -f .cproject
rm -f .project

#remove galaxy test data
rm -rf galaxy/test_data

find . -name \*~ -exec rm {} \;

for i in `find . -name .svn`
do
    #echo $i
    rm -rf $i 
done


# create a gzipped tarball containing the PALMapper package
cd ..
tar -czf $result_fn $name

echo "Finished package-building result is under $result_fn"
