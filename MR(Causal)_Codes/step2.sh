#!/bin/bash

for ((i=1;i<=22;i++));
do
	nohup /media/leelabsg-storage0/jangho/plink2 --bed /media/leelabsg-storage0/DATA_1/UKBB/cal/ukb_cal_chr"$i"_v2.bed --bim /media/leelabsg-storage0/DATA_1/UKBB/cal/ukb_snp_chr"$i"_v2.bim --fam /media/leelabsg-storage0/DATA_1/UKBB/cal/ukb45227_cal_chr20_v2_s488264.fam --extract markers_list.txt  --make-bed   --out MR_not_pruned_chr"$i" > log_chr"$i" &
done
