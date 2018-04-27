#!/bin/bash
arrayp=(2 3 4 6 9 12 18)
arrayt=(18 12 9 6 4 3 2)
len=${#arrayt[@]}
for ((i = 0; i<${len}; i++));
do
    echo ${arrayp[$i]} ${arrayt[$i]}
done