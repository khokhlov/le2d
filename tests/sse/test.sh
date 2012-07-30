#!/bin/bash

n=6000
s=10

for i in le2d le2d_soa le2d_sse le2d_
do
    echo $i
    double/$i $n $n $s > $i.d
    float/$i $n $n $s > $i.f
done

for i in le2d_cf
do
    echo $i
    double/$i $n $n $s $n > $i.d
    float/$i $n $n $s $n > $i.f
done
