#!/bin/bash

exp=$1

if [ $exp = "cluster" ];
then
	bash 5_cluster.sh
elif [ $exp = "AD" ];
then
	echo "waiting for join..."
else
	echo "please choose one experiments!"
fi
