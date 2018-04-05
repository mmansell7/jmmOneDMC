#!/bin/bash

echo Sorting results...
sort -k1,1g -k2,2g Summary_tmp.txt | uniq > tmp.txt
cp tmp.txt Summary_sorted.txt
rm tmp.txt
