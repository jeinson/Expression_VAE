gtex_vcf="/gpfs/commons/datasets/controlled/GTEx/dbgap_restricted/data/gtex/exchange/GTEx_phs000424/exchange/analysis_releases/GTEx_Analysis_2017-06-05_v8/genotypes/WGS/variant_calls/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.vcf.gz"
muscle_expr="/gpfs/commons/datasets/controlled/GTEx/dbgap_restricted/data/gtex/exchange/GTEx_phs000424/exchange/analysis_releases/GTEx_Analysis_2017-06-05_v8/eqtl/GTEx_Analysis_v8_eQTL_expression_matrices/Muscle_Skeletal.v8.normalized_expression.bed.gz"
covs="/gpfs/commons/groups/lappalainen_lab/jeinson/random/ML_genomics/FAKE_Muscle_Skeletal_encodings.txt"

date=$(printf '%(%m.%d.%y)T\n' -1)

module load qtltools

# Run eQTL analysis using the fake encodings as covariates. 
out_dir=$date"_fake_eqtl_results"
mkdir $out_dir
for i in $(seq 1 22); do
chr="chr"$i
echo "\
module load qtltools; \
qtltools cis \
--vcf $gtex_vcf \
--bed $muscle_expr \
--cov $covs \
--permute 100 1000 \
--region $chr \
--out $out_dir/Muscle_Skeletal_eqtl_FAKE_covars.$chr.$date.txt \
&> 'logs/skeletal_muscle_eqtl_FAKE.$chr.$date.txt' \
" | qsub -cwd -V -l mem=15G -N 'eqtl_'$chr'_'$date -o log_files -e log_files
done
