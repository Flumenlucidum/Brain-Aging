nohup step1_fitNULLGLMM.R
--plinkFile=/media/leelabsg_storage01/jangho/ukb_imp/UKB_CAL_for_step1 \
--phenoFile=/media/leelabsg_storage01/jangho/ukb_imp/ukb_whole_pheno.txt \
--phenoCol=corrected_delta \
--covarColList=COVAR1,COVAR2,COVAR3,COVAR4,COVAR5,COVAR6,COVAR7,COVAR8,COVAR9,COVAR10,COVAR11,COVAR12,cv1,cv2,cv3,cv4 \   
--sampleIDColinphenoFile=IID \
--traitType=quantitative \
--invNormalize=FALSE \
--outputPrefix=/media/leelabsg_storage01/jangho/ukb_imp/ukb_whole_saige_step1_result_LOCO \
--nThreads=4 \
--LOCO=TRUE \
--tauInit=1,0 >log_0812_LOCOtrue  &


