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

export OMP_NUM_THREADS=12
echo "Thread zahl=12" >>gcc2.txt
echo "Thread zahl=12" >>icc2.txt
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

export OMP_NUM_THREADS=12
echo "Thread zahl=12" >>gcc4.txt
echo "Thread zahl=12" >>icc4.txt
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

export OMP_NUM_THREADS=12
echo "Thread zahl=12" >>gcc32.txt
echo "Thread zahl=12" >>icc32.txt
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

echo "start7"
#32
export OMP_NUM_THREADS=1
echo "Thread zahl=1" >gcc7.txt
echo "Thread zahl=1" >icc7.txt
./gcc 7 >> gcc7.txt
./gcc 7 >> gcc7.txt
./gcc 7 >> gcc7.txt

./icc 7 >> icc7.txt
./icc 7 >> icc7.txt
./icc 7 >> icc7.txt
echo "plÃ¶p"
export OMP_NUM_THREADS=4
echo "Thread zahl=4" >>gcc7.txt
echo "Thread zahl=4" >>icc7.txt
./gcc 7 >> gcc7.txt
./gcc 7 >> gcc7.txt
./gcc 7 >> gcc7.txt

./icc 7 >> icc7.txt
./icc 7 >> icc7.txt
./icc 7 >> icc7.txt
echo "plÃ¶p"
export OMP_NUM_THREADS=8
echo "Thread zahl=8" >>gcc7.txt
echo "Thread zahl=8" >>icc7.txt
./gcc 7 >> gcc7.txt
./gcc 7 >> gcc7.txt
./gcc 7 >> gcc7.txt

./icc 7 >> icc7.txt
./icc 7 >> icc7.txt
./icc 7 >> icc7.txt
export OMP_NUM_THREADS=12
echo "Thread zahl=12" >>gcc7.txt
echo "Thread zahl=12" >>icc7.txt
./gcc 7 >> gcc7.txt
./gcc 7 >> gcc7.txt
./gcc 7 >> gcc7.txt

./icc 7 >> icc7.txt
./icc 7 >> icc7.txt
./icc 7 >> icc7.txt 
export OMP_NUM_THREADS=16
echo "Thread zahl=16" >>gcc7.txt
echo "Thread zahl=16" >>icc7.txt
./gcc 7 >> gcc7.txt
./gcc 7 >> gcc7.txt
./gcc 7 >> gcc7.txt

./icc 7 >> icc7.txt
./icc 7 >> icc7.txt
./icc 7 >> icc7.txt

