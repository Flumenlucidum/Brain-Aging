#!/bin/bash

for ((i=1;i<=22;i++));
do
	nohup /media/leelabsg-storage0/jangho/plink2 --bfile MR_not_pruned_chr"$i"  --extract MR_list_pruned_chr"$i".prune.in  --make-bed  --out MR_pruned_chr"$i" > log_chr"$i" &
done
