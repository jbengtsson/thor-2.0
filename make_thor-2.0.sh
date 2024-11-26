#!/bin/sh

dir=`pwd`

# Remove autoconf cashe.
rm -rf autom4te.cache
rm -rf aclocal.m4
rm -rf tho/lib/*

make distclean

mkdir -p config

# Configure libtool (for shared libraries).
#libtoolize

./bootstrap
./configure --prefix=$dir/thor

make install
