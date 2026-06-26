###5-3. Polygenic risk score (PRS)
plink --bfile ##name of bim/bed/fam file \
--score PGS002361.txt ##rsID_column ##effect_allele_column ##effect_weight_column header sum \
--no-sex \
--out PGS002361_PRS_sum