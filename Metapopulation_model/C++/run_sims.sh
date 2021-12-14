#!/bin/bash

max=640
for i in $(seq 1 $max)
do
    row=$(($max * $1 + $i + 1))
    ./a.out parameters.txt $row 2> /dev/null
done
