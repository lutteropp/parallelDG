#!/bin/bash
echo "start"

export OMP_NUM_THREADS=1
echo "Thread zahl=1" >gcc65.txt
echo "Thread zahl=1"
./gcc 65 >> gcc65.txt

export OMP_NUM_THREADS=4
echo "Thread zahl=4" >>gcc65.txt
echo "Thread zahl=4"
./gcc 65 >> gcc65.txt
./gcc 65 >> gcc65.txt


export OMP_NUM_THREADS=8
echo "Thread zahl=8" >>gcc65.txt
echo "Thread zahl=8"
./gcc 65 >> gcc65.txt
./gcc 65 >> gcc65.txt


echo "289"
export OMP_NUM_THREADS=1
echo "Thread zahl=1" >gcc289.txt

./gcc 289 >> gcc289.txt
./gcc 289 >> gcc289.txt


export OMP_NUM_THREADS=4
echo "Thread zahl=4" >>gcc289.txt
echo "Thread zahl=4" 
./gcc 289 >> gcc289.txt
./gcc 289 >> gcc289.txt

echo "Thread zahl=8" 
export OMP_NUM_THREADS=8
echo "Thread zahl=8" >>gcc289.txt

./gcc 289 >> gcc289.txt
./gcc 289 >> gcc289.txt




echo "1089" 
export OMP_NUM_THREADS=1
echo "Thread zahl=1" >gcc1089.txt


./gcc 1089 >> gcc1089.txt

echo "Thread zahl=4" 

export OMP_NUM_THREADS=4
echo "Thread zahl=4" >>gcc1089.txt

./gcc 1089 >> gcc1089.txt

echo "Thread zahl=8" 
export OMP_NUM_THREADS=8
echo "Thread zahl=8" >>gcc1089.txt

./gcc 1089 >> gcc1089.txt




echo "3969" 
export OMP_NUM_THREADS=1
echo "Thread zahl=1" >gcc3969.txt

./gcc 3969 >> gcc3969.txt

echo "Thread zahl=4" 
export OMP_NUM_THREADS=4
echo "Thread zahl=4" >>gcc3969.txt

./gcc 3969 >> gcc3969.txt

echo "Thread zahl=8" 
export OMP_NUM_THREADS=8
echo "Thread zahl=8" >>gcc3969.txt

./gcc 3969 >> gcc3969.txt






