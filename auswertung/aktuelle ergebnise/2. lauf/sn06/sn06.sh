#!/bin/bash
echo "start"

export OMP_NUM_THREADS=1
echo "Thread zahl=1" >gcc33.txt
echo "Thread zahl=1"
./seriel 33 >> gcc33.txt
./gcc 33 >> gcc33.txt

./gcc 33 >> gcc33.txt
export OMP_NUM_THREADS=4
echo "Thread zahl=4" >>gcc33.txt
echo "Thread zahl=4"
./gcc 33 >> gcc33.txt
./gcc 33 >> gcc33.txt

export OMP_NUM_THREADS=8
echo "Thread zahl=8" >>gcc33.txt
echo "Thread zahl=8"
./gcc 33 >> gcc33.txt
./gcc 33 >> gcc33.txt

export OMP_NUM_THREADS=12
echo "Thread zahl=12" >>gcc33.txt
echo "Thread zahl=12"
./gcc 33 >> gcc33.txt
./gcc 33 >> gcc33.txt

export OMP_NUM_THREADS=1
echo "Thread zahl=1" >gcc65.txt
echo "Thread zahl=1"
./seriel 65 >> gcc65.txt
./gcc 65 >> gcc65.txt

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

export OMP_NUM_THREADS=12
echo "Thread zahl=12" >>gcc65.txt
echo "Thread zahl=12"
./gcc 65 >> gcc65.txt
./gcc 65 >> gcc65.txt

echo "129"
export OMP_NUM_THREADS=1
echo "Thread zahl=1" >gcc129.txt
./seriel 129 >> gcc129.txt
./gcc 129 >> gcc129.txt

./gcc 129 >> gcc129.txt
export OMP_NUM_THREADS=4
echo "Thread zahl=4" >>gcc129.txt
echo "Thread zahl=4" 
./gcc 129 >> gcc129.txt
./gcc 129 >> gcc129.txt
echo "Thread zahl=8" 
export OMP_NUM_THREADS=8
echo "Thread zahl=8" >>gcc129.txt
./gcc 129 >> gcc129.txt
./gcc 129 >> gcc129.txt

echo "Thread zahl=12" 
export OMP_NUM_THREADS=12
echo "Thread zahl=1129" >>gcc129.txt
./gcc 129 >> gcc129.txt
./gcc 129 >> gcc129.txt

echo "257" 
export OMP_NUM_THREADS=1
echo "Thread zahl=1" >gcc257.txt
./seriel 257 >> gcc257.txt
./gcc 257 >> gcc257.txt

echo "Thread zahl=4" 

export OMP_NUM_THREADS=4
echo "Thread zahl=4" >>gcc257.txt
./gcc 257 >> gcc257.txt


echo "Thread zahl=8" 
export OMP_NUM_THREADS=8
echo "Thread zahl=8" >>gcc257.txt
./gcc 257 >> gcc257.txt


echo "Thread zahl=12" 
export OMP_NUM_THREADS=12
echo "Thread zahl=1257" >>gcc257.txt
./gcc 257 >> gcc257.txt

echo "513" 
export OMP_NUM_THREADS=1
echo "Thread zahl=1" >gcc513.txt
./seriel 513 >> gcc513.txt
./gcc 513 >> gcc513.txt

echo "Thread zahl=4" 

export OMP_NUM_THREADS=4
echo "Thread zahl=4" >>gcc513.txt
./gcc 513 >> gcc513.txt


echo "Thread zahl=8" 
export OMP_NUM_THREADS=8
echo "Thread zahl=8" >>gcc513.txt
./gcc 513 >> gcc513.txt


echo "Thread zahl=12" 
export OMP_NUM_THREADS=12
echo "Thread zahl=1513" >>gcc513.txt
./gcc 513 >> gcc513.txt


echo "1025" 
export OMP_NUM_THREADS=1
echo "Thread zahl=1" >gcc1025.txt
./seriel 1025 >> gcc1025.txt
./gcc 1025 >> gcc1025.txt

echo "Thread zahl=4" 
export OMP_NUM_THREADS=4
echo "Thread zahl=4" >>gcc1025.txt

./gcc 1025 >> gcc1025.txt

echo "Thread zahl=8" 
export OMP_NUM_THREADS=8
echo "Thread zahl=8" >>gcc1025.txt

./gcc 1025 >> gcc1025.txt

echo "Thread zahl=12" 
export OMP_NUM_THREADS=12
echo "Thread zahl=11025" >>gcc1025.txt

./gcc 1025 >> gcc1025.txt




