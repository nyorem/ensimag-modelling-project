#! /usr/bin/env bash

for mode in $(seq 0 4); do
    ./main 10 $mode
    poly="poly$mode".png
    angle="angle$mode".png

    echo "set terminal png size 400, 300 enhanced
    set output '$poly'
    plot 'Subdivision0.txt' with lines, 'Subdivision9.txt' with lines
    set output '$angle'
    plot 'AngleSubdivision9.txt' with lines" | gnuplot
done

