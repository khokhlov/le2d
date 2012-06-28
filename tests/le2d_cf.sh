#!/bin/bash

echo "" > le2d_cf.log

for ((c=4;c<7000;c+=64))
do
	echo $c
	./le2d_cf $c $c $((7000/$c)) $c >> le2d_cf.log
done

for cfs in 4096 4 16 32 64 128 256 512 1024 2048 
do
	echo "" > le2d_cf.$cfs.log

	for ((c=4;c<7000;c+=64))
	do
		if [ $c -ge $cfs ]
		then
			echo $c
			./le2d_cf $c $c $((7000/$c)) $cfs >> le2d_cf.$cfs.log
		fi
	done
done
