ml load bedtools
  
# 53 unique SVs out of 163 intersect (at least 1 bp)
# $1 is output file (number of intersection for each permutation)
# $2 is prefix
# $3 is output file for bedtools intersection with -F 0.5

# Repeat this x times
for i in {1..1000}; do
  echo $i

# 1. Select 163 SVs randomly in the 139 706 that were tested and lifted over
# We match on SV type, so we need to pick 86 insertions and 77 del/complex
  shuf -n 86 hprc2.500.maf0.01.geno0.1.hg38.INS.bed > random.ins.$2.tmp
  
  # 10% of 105 014 pb is 10 501 pb
  # If the cummulative lenghts of the 77 picked del and complex is more than 115 515 pb or less than 94 513 pb, we reject the selection
  # Length is in the SV id
  while true; do
    shuf -n 77 hprc2.500.maf0.01.geno0.1.hg38.DEL.COMPLEX.bed > random.del.complex.$2.tmp
    
    cum_len=0
    while read line; do len=$(echo $line | cut -f 4 -d ' ' | cut -f 3 -d ':' | cut -f 4 -d '-'); cum_len=$(($cum_len + $len)); done < random.del.complex.$2.tmp
    
    if (( $cum_len >= 94513 && $cum_len <= 115515 )); then
      break
    fi
    done  

# 2. Intersection of the selected SVs
  bedtools intersect -a random.ins.$2.tmp -b ENCODE4_cCREs_cCREregistry.hg38.bed -wo > intersect.INS.$2.tmp
  bedtools intersect -a random.del.complex.$2.tmp -b ENCODE4_cCREs_cCREregistry.hg38.bed -wo > intersect.DEL.COMPLEX.$2.tmp

# 3. How many unique SVs are intersecting a cCRE?
  cat intersect.INS.$2.tmp intersect.DEL.COMPLEX.$2.tmp > intersect.$2.tmp
  cut -f 4 intersect.$2.tmp > svs.$2.tmp
  uniq svs.$2.tmp | wc -l >> $1
  
  # 3.5. How many unique SVs (not INS) overlap at least 50% of the cCRE?
  bedtools intersect -a random.del.complex.$2.tmp -b ENCODE4_cCREs_cCREregistry.hg38.bed -F 0.5 -wo > intersect.DEL.COMPLEX.$2.tmp
  cut -f 4 intersect.DEL.COMPLEX.$2.tmp > svs.$2.tmp
  uniq svs.$2.tmp | wc -l >> $3
  

  # Clean up
  rm random.ins.$2.tmp
  rm random.del.complex.$2.tmp
  rm intersect.INS.$2.tmp
  rm intersect.DEL.COMPLEX.$2.tmp
  rm intersect.$2.tmp
  rm svs.$2.tmp
  
done
