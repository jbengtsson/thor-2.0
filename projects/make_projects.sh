#!/bin/sh

dir=`pwd`

# Remove autoconf cashe.
rm -rf autom4te.cache
rm -rf aclocal.m4

./bootstrap
./configure --prefix=$dir/projects

make install
