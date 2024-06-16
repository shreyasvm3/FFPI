#!/bin/bash

for i in {2..4..1}
do 
#	j=$( echo "scale = 2; $i/10" | bc)
	mkdir -p $i/
	cd $i/
	cp ../jobrun.sh ../dvr.x ../dyn.x ../input .
	sed -i "s/input-mass/$i/g" input
        qsub jobrun.sh	
	cd ../
done
