#!/bin/bash

n=7000
s0=16
s=64

bins=( le2d_base_f le2d_base_d le2d_f le2d_d le2d_f_sse le2d_d_sse le2d_f_avx le2d_d_avx )

for i in ${bins[@]}
do
	echo "" > $i.log
	cp ../../src/$i ./
done

for ((c=$s0;c<$n;c+=$s))
do
	for i in ${bins[@]}
	do
		echo $c $i
		./$i $c $c $(($n/$c)) >> $i.log
	done
done
