import sys
from pyspark import SparkConf
from pyspark.sql import SparkSession
import pyspark.sql.functions as F
import pandas as pd

# Read the LD proxies ids
proxies = []
with open("ld_proxies.ids", "r") as proxies_file:
  line = proxies_file.readline()
  while line:
    proxies.append(line.strip())
    line = proxies_file.readline()

# directories/datasets downloaded from OpenTarget
evidencePath = "credible_set/"
studyPath = "study/"

# establish spark connection
spark = (
  SparkSession.builder
  .master('local[*]')
  .getOrCreate()
)

# read datasets
evd = spark.read.parquet(evidencePath)
stud = spark.read.parquet(studyPath)

proxies_lit = F.array([F.lit(p) for p in proxies])

# browse the schema
#evd.printSchema()

# select fields of interest from datasets
evdSelect = (evd
  .filter((F.col("variantId").isin(proxies)) | (F.size(F.array_intersect(F.col("locus.variantId"), proxies_lit)) > 0))
  .select("studyLocusId","studyId",F.col("variantId").alias("variantIdLead"),"chromosome", "position", F.col("beta").alias("betaLead"), 
  F.col("pValueExponent").alias("pValueExponentLead"),"locus.variantId", "locus.pValueExponent", "locus.beta")
)
#evdSelect.show()

studSelect = (stud
  .select("studyId", "traitFromSource", "diseaseIds", "studyType")
)

#studSelect.show()

# merge datasets (this is to have the reported trait with the credible set data)
evdWithTrait = evdSelect.join(studSelect, on="studyId", how="left")
evdWithTrait = (evdWithTrait
  .filter(F.col("studyType") == "gwas")
)
#evdWithTrait.show()


df_pandas = evdWithTrait.toPandas()
#print(df_pandas)
df_pandas_expanded = df_pandas.explode(['variantId', 'pValueExponent', 'beta'])
#print(df_pandas)
# Could add a filtering step, cause now the df has 5M lines

# Pretty print in a file
with open('output.txt', 'w') as f:
    f.write(df_pandas_expanded.to_csv(index=False, sep='\t'))







