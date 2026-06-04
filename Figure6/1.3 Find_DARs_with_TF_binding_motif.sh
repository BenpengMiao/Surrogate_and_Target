#Candidate.txt:
#MORE_AL_DAR_TCDD_Li_F_adt_Ctrl_Li_F_adt_RUVr_k3_filter  motif3  Fox
#MORE_AL_DAR_TCDD_Li_F_adt_Ctrl_Li_F_adt_RUVr_k3_filter  motif4  Thrb_Rxr_Nr2f2_Hnf4
#MORE_AL_DAR_TCDD_Li_F_adt_Ctrl_Li_F_adt_RUVr_k3_filter  motif5  Cebp_Nfil3_Hlf
#MORE_AL_DAR_TCDD_Li_M_adt_Ctrl_Li_M_adt_RUVr_k3_filter  motif2  Cebp_Nfil3_Hlf
#MORE_AL_DAR_TCDD_Li_M_adt_Ctrl_Li_M_adt_RUVr_k3_filter  motif3  Thrb_Rxr_Nr2f2_Hnf4
#MORE_AL_DAR_TCDD_Li_M_adt_Ctrl_Li_M_adt_RUVr_k3_filter  motif4  Fox


cat Candidate.txt|while read l
do

file=$(echo $l|cut -d " " -f 1)
motif=$(echo $l|cut -d " " -f 2)
gene=$(echo $l|cut -d " " -f 3)

echo $file
echo $motif
echo $gene

~/software/homer/bin/findMotifsGenome.pl ${file}.txt mm10 ${file} -size given -find ${file}/homerResults/${motif}.motif > outputfile.txt

mv outputfile.txt ${file}_${gene}.txt

cut -f 1 ${file}_${gene}.txt|sed "1d"|sort -u > ${file}_${gene}_peak.txt

done
