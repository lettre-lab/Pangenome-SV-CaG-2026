//parameters

params.chrom_files = "chr_genes.hg38.txt" // this has two columns with chr and gene_id


//processes

process run_phenotype {

  tag "${chrom}_${pheno}" // label for job
  
  cpus 1
  memory '22.GB'
  time '40.m'
  executor 'slurm'
  clusterOptions '--account=ACCOUNT_NAME'
  
  errorStrategy 'ignore'

  input:
  tuple val(chrom), val(pheno)
  
  output:
  
  script:
  """

  question1_script.sh ${chrom} ${pheno}
  
  """
}

//workflow script

workflow { // Orchestrates channels and processes

  chrom_files = Channel.fromPath(params.chrom_files) // Create a channel of chromosome files
  
  // Expand each file into (chrom, pheno) tuples
  pairs = chrom_files
    .splitText()
    .map { line ->
      def (chrom, pheno) = line.tokenize('\t')
      return [chrom, pheno]
    }
    
  // Run the process for each pair 
  run_phenotype(pairs)

}