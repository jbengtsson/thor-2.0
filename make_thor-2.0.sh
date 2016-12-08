#!/bin/sh

dir=`pwd`

cd "$THOR_LIB"

rm -rf autom4te.cache
rm -rf aclocal.m4
rm -rf thor/bin/*
rm -rf thor/lib/*

make distclean

./bootstrap
./configure --prefix=$THOR_LIB/thor

make install
