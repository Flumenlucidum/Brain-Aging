#!/bin/bash

for ((i=1;i<=22;i++));
do
	nohup /media/leelabsg-storage0/jangho/plink2 --bfile MR_not_pruned_chr"$i"  --indep-pairwise 50 5 0.01    --out MR_list_pruned_chr"$i" > log_chr"$i" &
done
