#!/bin/bash
echo "start"
export OMP_NUM_THREADS=1
echo "Thread zahl=1" >gcc2.txt
echo "Thread zahl=1" >icc2.txt
./gcc 2 >> gcc2.txt
./gcc 2 >> gcc2.txt
./gcc 2 >> gcc2.txt

./icc 2 >> icc2.txt
./icc 2 >> icc2.txt
./icc 2 >> icc2.txt

export OMP_NUM_THREADS=4
echo "Thread zahl=4" >>gcc2.txt
echo "Thread zahl=4" >>icc2.txt
./gcc 2 >> gcc2.txt
./gcc 2 >> gcc2.txt
./gcc 2 >> gcc2.txt

./icc 2 >> icc2.txt
./icc 2 >> icc2.txt
./icc 2 >> icc2.txt

export OMP_NUM_THREADS=8
echo "Thread zahl=8" >>gcc2.txt
echo "Thread zahl=8" >>icc2.txt
./gcc 2 >> gcc2.txt
./gcc 2 >> gcc2.txt
./gcc 2 >> gcc2.txt

./icc 2 >> icc2.txt
./icc 2 >> icc2.txt
./icc 2 >> icc2.txt

export OMP_NUM_THREADS=16
echo "Thread zahl=16" >>gcc2.txt
echo "Thread zahl=16" >>icc2.txt
./gcc 2 >> gcc2.txt
./gcc 2 >> gcc2.txt
./gcc 2 >> gcc2.txt

./icc 2 >> icc2.txt
./icc 2 >> icc2.txt
./icc 2 >> icc2.txt
echo "start 4"
#4
export OMP_NUM_THREADS=1
echo "Thread zahl=1" >gcc4.txt
echo "Thread zahl=1" >icc4.txt
./gcc 4 >> gcc4.txt
./gcc 4 >> gcc4.txt
./gcc 4 >> gcc4.txt

./icc 4 >> icc4.txt
./icc 4 >> icc4.txt
./icc 4 >> icc4.txt

export OMP_NUM_THREADS=4
echo "Thread zahl=4" >>gcc4.txt
echo "Thread zahl=4" >>icc4.txt
./gcc 4 >> gcc4.txt
./gcc 4 >> gcc4.txt
./gcc 4 >> gcc4.txt

./icc 4 >> icc4.txt
./icc 4 >> icc4.txt
./icc 4 >> icc4.txt

export OMP_NUM_THREADS=8
echo "Thread zahl=8" >>gcc4.txt
echo "Thread zahl=8" >>icc4.txt
./gcc 4 >> gcc4.txt
./gcc 4 >> gcc4.txt
./gcc 4 >> gcc4.txt

./icc 4 >> icc4.txt
./icc 4 >> icc4.txt
./icc 4 >> icc4.txt

export OMP_NUM_THREADS=16
echo "Thread zahl=16" >>gcc4.txt
echo "Thread zahl=16" >>icc4.txt
./gcc 4 >> gcc4.txt
./gcc 4 >> gcc4.txt
./gcc 4 >> gcc4.txt

./icc 4 >> icc4.txt
./icc 4 >> icc4.txt
./icc 4 >> icc4.txt

echo "start32"
#32
export OMP_NUM_THREADS=1
echo "Thread zahl=1" >gcc32.txt
echo "Thread zahl=1" >icc32.txt
./gcc 32 >> gcc32.txt
./gcc 32 >> gcc32.txt
./gcc 32 >> gcc32.txt

./icc 32 >> icc32.txt
./icc 32 >> icc32.txt
./icc 32 >> icc32.txt

export OMP_NUM_THREADS=4
echo "Thread zahl=4" >>gcc32.txt
echo "Thread zahl=4" >>icc32.txt
./gcc 32 >> gcc32.txt
./gcc 32 >> gcc32.txt
./gcc 32 >> gcc32.txt

./icc 32 >> icc32.txt
./icc 32 >> icc32.txt
./icc 32 >> icc32.txt

export OMP_NUM_THREADS=8
echo "Thread zahl=8" >>gcc32.txt
echo "Thread zahl=8" >>icc32.txt
./gcc 32 >> gcc32.txt
./gcc 32 >> gcc32.txt
./gcc 32 >> gcc32.txt

./icc 32 >> icc32.txt
./icc 32 >> icc32.txt
./icc 32 >> icc32.txt

export OMP_NUM_THREADS=16
echo "Thread zahl=16" >>gcc32.txt
echo "Thread zahl=16" >>icc32.txt
./gcc 32 >> gcc32.txt
./gcc 32 >> gcc32.txt
./gcc 32 >> gcc32.txt

./icc 32 >> icc32.txt
./icc 32 >> icc32.txt
./icc 32 >> icc32.txt

echo "start128"
#128
export OMP_NUM_THREADS=1
echo "Thread zahl=1" >gcc128.txt
echo "Thread zahl=1" >icc128.txt
./gcc 128 >> gcc128.txt
./gcc 128 >> gcc128.txt
./gcc 128 >> gcc128.txt

./icc 128 >> icc128.txt
./icc 128 >> icc128.txt
./icc 128 >> icc128.txt

export OMP_NUM_THREADS=4
echo "Thread zahl=4" >>gcc128.txt
echo "Thread zahl=4" >>icc128.txt
./gcc 128 >> gcc128.txt
./gcc 128 >> gcc128.txt
./gcc 128 >> gcc128.txt

./icc 128 >> icc128.txt
./icc 128 >> icc128.txt
./icc 128 >> icc128.txt

export OMP_NUM_THREADS=8
echo "Thread zahl=8" >>gcc128.txt
echo "Thread zahl=8" >>icc128.txt
./gcc 128 >> gcc128.txt
./gcc 128 >> gcc128.txt
./gcc 128 >> gcc128.txt

./icc 128 >> icc128.txt
./icc 128 >> icc128.txt
./icc 128 >> icc128.txt

export OMP_NUM_THREADS=16
echo "Thread zahl=16" >>gcc128.txt
echo "Thread zahl=16" >>icc128.txt
./gcc 128 >> gcc128.txt
./gcc 128 >> gcc128.txt
./gcc 128 >> gcc128.txt

./icc 128 >> icc128.txt
./icc 128 >> icc128.txt
./icc 128 >> icc128.txt

echo "start256"
#256
export OMP_NUM_THREADS=1
echo "Thread zahl=1" >gcc256.txt
echo "Thread zahl=1" >icc256.txt
./gcc 256 >> gcc256.txt
./gcc 256 >> gcc256.txt


./icc 256 >> icc256.txt
./icc 256 >> icc256.txt


export OMP_NUM_THREADS=4
echo "Thread zahl=4" >>gcc256.txt
echo "Thread zahl=4" >>icc256.txt
./gcc 256 >> gcc256.txt
./gcc 256 >> gcc256.txt


./icc 256 >> icc256.txt
./icc 256 >> icc256.txt


export OMP_NUM_THREADS=8
echo "Thread zahl=8" >>gcc256.txt
echo "Thread zahl=8" >>icc256.txt
./gcc 256 >> gcc256.txt
./gcc 256 >> gcc256.txt


./icc 256 >> icc256.txt
./icc 256 >> icc256.txt


export OMP_NUM_THREADS=16
echo "Thread zahl=16" >>gcc256.txt
echo "Thread zahl=16" >>icc256.txt
./gcc 256 >> gcc256.txt
./gcc 256 >> gcc256.txt


./icc 256 >> icc256.txt
./icc 256 >> icc256.txt
