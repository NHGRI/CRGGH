#######################
##Prepare Instruments##
#######################

### Correct for multiple testing

# add header to mQTL output
sed -i '1i CpG CpG_chr CpG_start CpG_end strand n_cis dist snp snp_chr snp_start snp_end nom_pval r2 beta se best_hit' mQTLs_out.txt

# Use script to calculate FDR for large files from http://blog.mcbryan.co.uk/2013/01/false-discovery-rates-and-large-files.html
./PtoFDR mQTLs_out.txt 5 > mQTLs_out_fdr.txt

# Select FDR < 0.05 and grep list with CpGs of interest 
awk '$6 < 0.05' mQTLs_out_fdr.txt > mQTLs_out_fdr.05.txt
grep -f List_cpgs_n57.txt  mQTLs_out_fdr.05.txt >  mQTLs_out_fdr.05_57CpGs.txt

less -S mQTLs_out_fdr.05_57CpGs.txt | sort | uniq -c | awk '{print $2}' > List_cpgs_wmQTL_forclumping.txt #These are the potential instruments

#Make file per cpg - get cpg ids from List_cpgs_wmQTL.txt
while read -r cpg; do
    head -n1 mQTLs_out_fdr.05_57CpGs.txt > ${cpg}_mqtls_fdr0.05.txt
    awk -F'\t' -v cpg="$cpg" '$1 == cpg' mQTLs_out_fdr.05_57CpGs.txt >> ${cpg}_mqtls_fdr0.05.txt
done < List_cpgs_wmQTL_forclumping.txt



#### LD Clumping

# Convert files to plink format
plink2 --vcf Genotyping.vcf.gz --make-bed --double-id --out Genotyping.plink

# Perform clumping per CpG [ get cpgs from List_cpgs_wmQTL.txt]
while read -r cpg; do
plink \
    --bfile /Data/Genotyping.plink \
    --clump-p1 0.05 \
    --clump-r2 0.01 \
    --clump-kb 10000 \
    --clump /Data/${cpg}_mqtls_fdr0.05.txt \
    --clump-snp-field snp \
    --clump-field fdr \
    --out /Clumping/${cpg}_mqtls_fdr0.05.r01;
done < List_cpgs_wmQTL_forclumping.txt


# Make a file with all clumped mQTLs combined
head -n1 cpg1_mqtls_fdr0.05.r01.clumped.txt > allcpgs_mqtls_fdr0.05.r01txt
while read -r cpg; do
 awk 'NR >1 {print}' ${cpg}_mqtls_fdr0.05.r01.clumped.txt |  grep -v '^$' >> allcpgs_mqtls_fdr0.05.r01txt;
done < List_cpgs_wmQTL_forclumping.txt
