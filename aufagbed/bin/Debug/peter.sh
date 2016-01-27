#!/bin/bash
echo "start"

echo "Thread zahl=1" >gcc2.txt

./gcc 2 >> gcc2.txt
./gcc 2 >> gcc2.txt
./gcc 2 >> gcc2.txt

echo "start 4"
#4

echo "Thread zahl=1" >gcc4.txt

./gcc 4 >> gcc4.txt
./gcc 4 >> gcc4.txt
./gcc 4 >> gcc4.txt


echo "start32"
#32
export OMP_NUM_THREADS=1
echo "Thread zahl=1" >gcc32.txt

./gcc 32 >> gcc32.txt
./gcc 32 >> gcc32.txt
./gcc 32 >> gcc32.txt

echo "start128"
#128
export OMP_NUM_THREADS=1
echo "Thread zahl=1" >gcc128.txt

./gcc 128 >> gcc128.txt
./gcc 128 >> gcc128.txt
./gcc 128 >> gcc128.txt


echo "start1024"
#1024
export OMP_NUM_THREADS=1
echo "Thread zahl=1" >gcc1024.txt

./gcc 1024 >> gcc1024.txt
./gcc 1024 >> gcc1024.txt
./gcc 1024 >> gcc1024.txt
