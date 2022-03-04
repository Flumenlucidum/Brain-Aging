#!/bin/bash

for ((i=1;i<=22;i++));
do
	nohup /media/leelabsg-storage0/jangho/plink2 --bfile MR_pruned_chr"$i"  --keep /media/leelabsg-storage0/jangho/ukb_imp/samplelist_34129.txt  --make-bed  --out MR_pruned_chr"$i" > log_chr"$i" &

	nohup /media/leelabsg-storage0/jangho/plink2 --bfile MR_not_pruned_chr"$i"  --keep /media/leelabsg-storage0/jangho/ukb_imp/samplelist_34129.txt  --make-bed  --out MR_not_pruned_chr"$i" > log_not_chr"$i" &

done
