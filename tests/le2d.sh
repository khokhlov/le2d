#!/bin/bash

echo "" > le2d.log

for ((c=4;c<7000;c+=64))
do
	echo $c
	./le2d $c $c $((7000/$c)) >> le2d.log
done
