#!/bin/bash

echo "" > le2d_sse.log

for ((c=4;c<5000;c+=64))
do
	echo $c
	./le2d_sse $c $c $((7000/$c)) >> le2d_sse.log
done
