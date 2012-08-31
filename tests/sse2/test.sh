#!/bin/bash

n=6000
s=10

bins=( le2d_f le2d_d le2d_f_sse le2d_d_sse le2d_f_avx le2d_d_avx )

for i in ${bins[@]}
do
	echo "" > $i.log
	cp ../../src/$i ./
done

for ((c=16;c<7000;c+=64))
do
	for i in ${bins[@]}
	do
		echo $c $i
		./$i $c $c $((7000/$c)) >> $i.log
	done
done
