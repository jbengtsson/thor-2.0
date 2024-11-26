#!/bin/sh

dir=`pwd`

# Remove autoconf cashe.
rm -rf autom4te.cache
rm -rf aclocal.m4

make distclean

./bootstrap
./configure --prefix=$dir/projects

make
