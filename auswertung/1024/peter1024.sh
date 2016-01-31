#!/bin/bash
echo "start"






export OMP_NUM_THREADS=4
echo "Thread zahl=4" >>gcc1024.txt
echo "Thread zahl=4" >>icc1024.txt
./gcc 1024 >> gcc1024.txt
echo "start2"
./gcc 1024 >> gcc1024.txt
echo "start3"

./icc 1024 >> icc1024.txt
echo "start4"
./icc 1024 >> icc1024.txt
echo "start5"

export OMP_NUM_THREADS=8
echo "Thread zahl=8" >>gcc1024.txt
echo "Thread zahl=8" >>icc1024.txt
./gcc 1024 >> gcc1024.txt
echo "start6"
./gcc 1024 >> gcc1024.txt

echo "start7"
./icc 1024 >> icc1024.txt
echo "start8"
./icc 1024 >> icc1024.txt




