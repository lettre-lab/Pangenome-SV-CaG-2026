ml load plink/2.00a5.8
ml load python/3.13.2
ml load scipy-stack

# Parameters
# $1 = chrom (ex. "chr21")
# $2 = $2

path=path/to/files/ # adjust this

# Step 1: Linear regression between phenotypes and standard approach variants from CaG

plink2 --bfile $path/$1/genorefq.500.maf0.01.$1 \
--glm no-x-sex zs \
--covar-variance-standardize \
--covar $path/covariates.txt \
--pheno $path/Matrix.gene_expression.INT_Norm.txt \
--pheno-name $2 \
--out $path/$1/genorefq.500.maf0.01.$1.25nov25.$2

# Step 1: Conditional analysis with python script using linear regression from step 1

# Create a pheno file from Matrix. ...

awk -v target="$2" '
NR==1 {
    for (i=1; i<=NF; i++)
        if ($i == target) keep=i
}
{
    print $1, $2, $keep
}' $path/Matrix.gene_expression.INT_Norm.txt > \
$path/$1/Matrix.gene_expression.INT_Norm.$2.txt

num="${1#chr}" # Removing the 'chr' from the chrom variable

python3 $path/cond.13nov25.py \
$path/$1/genorefq.500.maf0.01.$1.bed \
$path/$1/genorefq.500.maf0.01.$1.25nov25.$2.$2.glm.linear.zst \
$num \
$num \
$2 \
$path/$1/Matrix.gene_expression.INT_Norm.$2.txt \
$path/covariates.txt \
> $path/$1/conditional.$1.$2.txt

# Check if file is empty
if [ ! -s "$path/$1/conditional.$1.$2.txt" ]; then
  echo "file empty"
  plink2 --bfile $path/$1/hprc2.500.maf0.01.$1 \
  --covar $path/covariates.txt \ # this changes
  --covar-variance-standardize \
  --glm no-x-sex \
  --pheno $path/Matrix.gene_expression.INT_Norm.txt \
  --pheno-name $2 \
  --out $path/$1/hprc2.500.maf0.01.$1.$2.genorefq_conditions.25nov25
  
  head -n 1 $path/$1/hprc2.500.maf0.01.$1.$2.genorefq_conditions.25nov25.$2.glm.linear > $path/$1/tmp.$1.$2.header
  awk '$10 == "ADD"' $path/$1/hprc2.500.maf0.01.$1.$2.genorefq_conditions.25nov25.$2.glm.linear > $path/$1/tmp.$1.$2.ADD
  cat $path/$1/tmp.$1.$2.header $path/$1/tmp.$1.$2.ADD > $path/$1/hprc2.500.maf0.01.$1.$2.genorefq_conditions.25nov25.ADD.glm.linear
  
  # Clean up
  rm $path/$1/Matrix.gene_expression.INT_Norm.$2.txt
  rm $path/$1/genorefq.500.maf0.01.$1.25nov25.$2.log
  rm $path/$1/tmp.$1.$2.header
  rm $path/$1/tmp.$1.$2.ADD
  rm $path/$1/hprc2.500.maf0.01.$1.$2.genorefq_conditions.25nov25.$2.glm.linear
  
  
else

# Step 3: Extract lead variants that will be used as covariates

  awk '$19 <= 5e-8' $path/$1/conditional.$1.$2.txt | cut -f 3 > $path/$1/ids.$1.$2.list

  plink2 --bfile $path/$1/genorefq.500.maf0.01.$1 \
  --export A \
  --extract $path/$1/ids.$1.$2.list \
  --out $path/$1/genotypes.$1.$2

  # Merge with covariate files
  head -n 1 $path/$1/genotypes.$1.$2.raw | cut -f 7- > $path/$1/genotypes.$1.$2.header
  sed '1d' $path/$1/genotypes.$1.$2.raw > $path/$1/genotypes.$1.$2.raw2
  sort -k2,2 $path/$1/genotypes.$1.$2.raw2 > $path/$1/genotypes.$1.$2.sorted
  cut -f 7- $path/$1/genotypes.$1.$2.sorted > $path/$1/genotypes.$1.$2.tmp
  cat $path/$1/genotypes.$1.$2.header $path/$1/genotypes.$1.$2.tmp > $path/$1/genotypes.$1.$2.tmp2
  paste $path/covariates.txt $path/$1/genotypes.$1.$2.tmp2 > $path/$1/covariates.w_conditions.$1.$2.tsv

  # Etape 4: Linear regression with HPRC SVs
  plink2 --bfile $path/$1/hprc2.500.maf0.01.$1 \
  --covar $path/$1/covariates.w_conditions.$1.$2.tsv \
  --covar-variance-standardize \
  --glm no-x-sex \
  --pheno $path/Matrix.gene_expression.INT_Norm.txt \
  --pheno-name $2 \
  --out $path/$1/hprc2.500.maf0.01.$1.$2.genorefq_conditions.25nov25

  head -n 1 $path/$1/hprc2.500.maf0.01.$1.$2.genorefq_conditions.25nov25.$2.glm.linear > $path/$1/tmp.$1.$2.header
  awk '$10 == "ADD"' $path/$1/hprc2.500.maf0.01.$1.$2.genorefq_conditions.25nov25.$2.glm.linear > $path/$1/tmp.$1.$2.ADD
  cat $path/$1/tmp.$1.$2.header $path/$1/tmp.$1.$2.ADD > $path/$1/hprc2.500.maf0.01.$1.$2.genorefq_conditions.25nov25.ADD.glm.linear

  # Clean up
  rm $path/$1/genotypes.$1.$2.header
  rm $path/$1/genotypes.$1.$2.raw2
  rm $path/$1/genotypes.$1.$2.raw
  rm $path/$1/genotypes.$1.$2.sorted
  rm $path/$1/genotypes.$1.$2.tmp
  rm $path/$1/genotypes.$1.$2.tmp2
  rm $path/$1/Matrix.gene_expression.INT_Norm.$2.txt
  rm $path/$1/genorefq.500.maf0.01.$1.25nov25.$2.log
  rm $path/$1/covariates.w_conditions.$1.$2.tsv
  rm $path/$1/genotypes.$1.$2.log
  rm $path/$1/tmp.$1.$2.header
  rm $path/$1/tmp.$1.$2.ADD
  rm $path/$1/hprc2.500.maf0.01.$1.$2.genorefq_conditions.25nov25.$2.glm.linear
  
fi
