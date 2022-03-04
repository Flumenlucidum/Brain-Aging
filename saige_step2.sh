#!/bin/bash
source /home/leelabsg/anaconda3/etc/profile.d/conda.sh
conda activate saige_new

for((i=1;i<=22;i++))
do
	nohup step2_SPAtests.R --bgenFile=/home/lee7801/DATA/UKBB/imp/ukb_imp_chr"$i"_v3.bgen \
	--bgenFileIndex=/home/lee7801/DATA/UKBB/imp/ukb_imp_chr"$i"_v3.bgen.bgi \
	--minMAF=0.0001 \
	--minMAC=1 \
	--sampleFile=/media/leelabsg_storage01/jangho/ukb_imp/saige_step1_sample_for_bgen.txt \
	--GMMATmodelFile=/media/leelabsg_storage01/jangho/ukb_imp/ukb_whole_saige_step1_result_LOCO.rda \
	--varianceRatioFile=/media/leelabsg_storage01/jangho/ukb_imp/ukb_whole_saige_step1_result_LOCO.varianceRatio.txt \
	--SAIGEOutputFile=/media/leelabsg_storage01/jangho/ukb_imp/ukb_whole_saige_step2_result_LOCO_chr${i}.txt \
	--numLinesOutput=2 \
	--chrom="$i" \
	--LOCO=TRUE \
	--IsOutputAFinCaseCtrl=TRUE > /media/leelabsg_storage01/jangho/ukb_imp/log_ukb_whole_saige_step2_LOCO_chr"$i" &
done
