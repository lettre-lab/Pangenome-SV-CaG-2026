# Script for figures #

### Fig.1 ###
library(eulerr)
library(scales)
library(ggridges)
library(dplyr)
library(ggbreak)

# B #
fc <- read.delim("HPRC2.SVs.FC.afreq")
h <- read.delim("HPRC2.SVs.H.afreq")
m <- read.delim("HPRC2.SVs.M.afreq")
fc <- subset(fc, ALT_FREQS>0)
h <- subset(h, ALT_FREQS>0)
m <- subset(m, ALT_FREQS>0)

fit <- euler(list(
  QFC = fc$ID,
  QHA = h$ID,
  QMO = m$ID))

counts <- fit$original.values
total <- sum(counts)
labels <- paste0(comma(counts),"\n(",percent(counts / total, accuracy = 0.1),")")

B <- plot(
  fit,
  fills = c("#619cff", "#00ba38", "#f8766d"),
  quantities = list(labels = labels, cex=0.),
  edges = TRUE,
  labels = list(cex = 0.8))

# A #
fc <- read.delim("HPRC1.SVs.FC.afreq")
h <- read.delim("HPRC1.SVs.H.afreq")
m <- read.delim("HPRC1.SVs.M.afreq")
colnames(fc)[colnames(fc)=="ALT_FREQS"] <- "AF_FC"
colnames(h)[colnames(h)=="ALT_FREQS"] <- "AF_H"
colnames(m)[colnames(m)=="ALT_FREQS"] <- "AF_M"
af <- cbind(data.frame(fc$ID), data.frame(fc$AF_FC), data.frame(h$AF_H), data.frame(m$AF_M))
af$ID_new <- paste0(af$fc.ID, "_", sample(1:400000, nrow(af), replace = FALSE))
af <- af[,c(5,2,3,4)]

fc <- subset(af, fc.AF_FC>0)
h <- subset(af, h.AF_H>0)
m <- subset(af, m.AF_M>0)

fit <- euler(list(
  QFC = fc$ID_new,
  QHA = h$ID_new,
  QMO = m$ID_new))

counts <- fit$original.values
total <- sum(counts)
labels <- paste0(comma(counts),"\n(",percent(counts / total, accuracy = 0.1),")")

A <- plot(
  fit,
  fills = c("#619cff", "#00ba38", "#f8766d"),
  quantities = list(labels = labels, cex=0.8),
  edges = TRUE, labels = list(cex = 0.8))

# C #
fc <- read.delim("Standard.SVs.FC.afreq")
h <- read.delim("Standard.SVs.H.afreq")
m <- read.delim("Standard.SVs.M.afreq")
fc <- subset(fc, ALT_FREQS>0)
h <- subset(h, ALT_FREQS>0)
m <- subset(m, ALT_FREQS>0)

fit <- euler(list(
  QFC = fc$ID,
  QHA = h$ID,
  QMO = m$ID))

counts <- fit$original.values
total <- sum(counts) #25,474, the other SVs are only found in the other participants that are not part of the 3 groups
labels <- paste0(comma(counts),"\n(",percent(counts / total, accuracy = 0.1),")")

C <- plot(
  fit,
  fills = c("#619cff", "#00ba38", "#f8766d"),
  quantities = list(labels = labels, cex=0.8),
  edges = TRUE, labels = list(cex = 0.8))

# D #
dicho <- read.delim("CaG.Dicho_Traits.WGS.txt")
dicho <- dicho[,c(2,4,7)]

r1a <- read.delim("HPRC2.SVs.scount")
r1a <- merge(r1a, dicho, by.x="X.IID", by.y="code_111")
r1a$num_SV <- r1a$HET_REF_ALT_CT+r1a$HOM_ALT_CT
r1a <- subset(r1a, AncestralGroup=="FrenchCanada" | AncestralGroup=="Haiti" | AncestralGroup=="Morocco")

r1 <- read.delim("HPRC1.SVs.scount")
r1 <- merge(r1, dicho, by.x="X.IID", by.y="code_111")
r1$num_SV <- r1$HET_REF_ALT_CT+r1$HOM_ALT_CT
r1 <- subset(r1, AncestralGroup=="FrenchCanada" | AncestralGroup=="Haiti" | AncestralGroup=="Morocco")

s <- read.delim("Standard.SVs.scount")
ids <- read.delim("CaG.IDs.WGS.tsv")
ids <- ids[,c(1,3)]
ids <- merge(ids, dicho)
s <- merge(s, ids, by.x="X.IID", by.y="code_ns")
s$num_SV <- s$HET_REF_ALT_CT+s$HOM_ALT_CT
s <- subset(s, AncestralGroup=="FrenchCanada" | AncestralGroup=="Haiti" | AncestralGroup=="Morocco")

hp1 <- r1[,c("X.IID", "num_SV", "AncestralGroup")]
hp2 <- r1a[,c("X.IID", "num_SV", "AncestralGroup")]
st <- s[,c("code_111", "num_SV", "AncestralGroup")]
colnames(st)[colnames(st)=="code_111"] <- "X.IID"
hp1$Approach <- "HPRC1"
hp2$Approach <- "HPRC2"
st$Approach <- "Standard"
df <- bind_rows(hp1, hp2, st)
df$AncestralGroup <- ifelse(df$AncestralGroup=="FrenchCanada", "QFC",
                            ifelse(df$AncestralGroup=="Morocco", "QMO", "QHA"))
df$Approach <- factor(df$Approach, levels = c("Standard", "HPRC2", "HPRC1"))

D <- ggplot(df, aes(x = num_SV, y = Approach, fill = AncestralGroup)) +
  geom_density_ridges(alpha = 0.6) +
  scale_fill_manual(values = c("#619cff","#00ba38","#f8766d")) + 
  theme_ridges() +
  scale_x_break(c(3000, 24000)) + 
  scale_x_continuous(
    limits = c(1000, 40000),
    breaks = c(1000, 2000, 3000, 25000, 30000, 35000)) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.x.top = element_blank(), 
    axis.ticks.x.top = element_blank(),
    axis.title.x = element_text(hjust = 0.5, margin = margin(t = 10)),
    axis.title.y = element_text(hjust = 0.5, margin = margin(r = 10)),
    legend.position = "top", legend.justification = "center") +
  labs(x="Number of SVs per genome", y="Approach", fill="Ancestral Group")


### Fig. 2 ###
library(ggrepel)
library(ggplot2)

asso <- read.delim("Conditonal_Analysis.Association_Results.beta_for_minor_allele.tsv")
asso$CHROM <- mapply(function(id){return(as.numeric(str_split(id, "-")[[1]][1]))}, asso$ID)
asso$CHROM <- ifelse(is.na(asso$CHROM), 23, asso$CHROM)
asso$POS <- mapply(function(id){return(as.numeric(str_split(id, "-")[[1]][2]))}, asso$ID)

chr_lengths <- read.delim("chm13.lengths_chr.tsv")
genome_df <- chr_lengths %>% mutate(tot = lag(cumsum(as.numeric(chr_len)), default = 0))
df_plot <- asso %>% left_join(genome_df, by = c("CHROM" = "bin_chr")) %>%mutate(bp_cum = POS + tot)
axis_df <- genome_df %>%mutate(center = tot + chr_len / 2)

indep <- read.delim("unique_independent_SV_eQTLs.genes_paires.txt", header=FALSE)
colnames(indep) <- c("ID", "PHENO", "GENE")
top20 <- merge(df_plot, indep, by.x="ID", by.y="ID")
top20 <- subset(top20, PHENO.x==PHENO.y)
top20 <- top20[order(top20$P),]
top20 <- top20[c(1:20),]
colnames(top20)[colnames(top20)=="PHENO.x"] <- "PHENO"
colnames(top20)[colnames(top20)=="CHROM.x"] <- "CHROM"
df_plot <- df_plot %>% left_join(top20, by = c("ID", "PHENO"))
df_labels <- df_plot %>% filter(!is.na(GENE) & GENE != "")

df_plot$TYPE <- mapply(function(id){return(str_split(id, "-")[[1]][3])}, df_plot$ID)

A <- ggplot(df_plot, aes(x = bp_cum.x, y = -log10(P.x), shape=TYPE, color=TYPE)) +
  geom_vline(xintercept = axis_df$tot ,color = "grey90",linewidth = 0.8, lty="dotted") +
  geom_point(size = 3, alpha=0.6) + 
  scale_color_manual(
    values = c("INS"="#619cff", "DEL"="#00ba38", "COMPLEX"="#f8766d"), name = "SV type") +
  scale_shape_manual(values = c("INS" = 16,"COMPLEX" = 17,"DEL" = 15), name ="SV type") +
  scale_x_continuous(
    breaks = axis_df$center,
    labels = axis_df$bin_chr) +
  labs(x = "Chromosome", y = expression(-log[10](P))) +
  theme_bw() +
  theme(panel.grid = element_blank(), text = element_text(size = 13)) + 
  geom_hline(yintercept = -log10(5e-8), color="red", size=0.8, lty="dashed") +
  geom_text_repel(
    aes(label = GENE),
    size = 3.5,
    max.overlaps = Inf,
    show.legend = FALSE,
    box.padding = 0.8,
    point.padding = 0.5,
    force = 2,
    min.segment.length = 0,
    segment.color = "grey40",
    segment.size = 0.4, color="grey40")

# B and C #
b <- read.delim("hprc2.500.maf0.01.chr17.ENSG00000175711.genorefq_conditions.25nov25.ADD.glm.linear")
a <- read.delim("hprc2.500.maf0.01.chr17.ENSG00000175711.genorefq_conditions.25nov25.condition_DEL.ADD.glm.linear")
tmp_b <- subset(b, P<=5e-8)
tmp_a <- subset(a, ID %in% tmp_b$ID)

B <- ggplot(tmp_b, aes(x=POS, y=-log10(P))) + geom_point(size=3, alpha=0.5, color="blue") + 
  xlim(c(84071428-1000000, 84071428+1000000)) + 
  geom_hline(yintercept=-log10(5e-8), color="red", linetype="dashed") +
  geom_text_repel(data=subset(b, P<=5e-8), aes(label = ID),size = 3 ,max.overlaps = 50) +
  geom_vline(xintercept=84071428, color="darkgrey", linetype="dotted", size=0.7) +
  geom_point(data=subset(b, ID=="17-84059150-DEL-115"), aes(x=POS, y=-log10(P)), size=3, alpha=1, color="red") +
  geom_point(data=subset(b, P>5e-8), aes(x=POS, y=-log10(P)), size=3, alpha=0.5, color="darkgrey") +
  labs(x="Position on chr17") +
  theme_bw() + theme(panel.grid = element_blank(), text = element_text(size = 13)) + 
  scale_y_continuous(breaks = c(0, 2, 4, 6, 8, 10, 12), limits = c(0, 12.5))

C <- ggplot(data=subset(a, !(ID %in% tmp_b$ID)), aes(x=POS, y=-log10(P))) + geom_point(size=3, alpha=0.5, color="darkgrey") +
  geom_point(data=tmp_a, aes(x=POS, y=-log10(P)), size=3, alpha=0.5, color="blue") + 
  xlim(c(84071428-1000000, 84071428+1000000)) + 
  geom_hline(yintercept=-log10(5e-8), color="red", linetype="dashed") +
  geom_text_repel(data=subset(a, P<=5e-8 | ID %in% tmp_b$ID), aes(label = ID),size = 3,max.overlaps = 50) +
  geom_vline(xintercept=84071428, color="darkgrey", linetype="dotted", size=0.7) +
  labs(x="Position on chr17") +
  theme_bw() + theme(panel.grid = element_blank(), text = element_text(size = 13)) + 
  scale_y_continuous(breaks = c(0, 2, 4, 6, 8, 10, 12), limits = c(0, 12.5))

# D and E #
i <- read.table("10000.num_intersects", quote="\"", comment.char="")
D <- ggplot(i, aes(x=V1)) +geom_histogram(alpha=0.3, color="#9990AB", fill="#9990AB", bins=30) + 
  geom_vline(xintercept = 36, color="red", lty="dashed", size=0.75) +
  labs(y="Random sets count", x="Number of SVs intersecting a cCRE") + #theme_classic() +
  annotate("text", x = 36, y = Inf,label = "36", color = "red", fontface = "bold",vjust = 1.5, hjust = -0.5) +
  theme_bw() + theme(panel.grid = element_blank(), text = element_text(size = 13)) + ylim(c(0,2500))

i5 <- read.table("10000.num_intersects.F0.5", quote="\"", comment.char="")
E <- ggplot(i5, aes(x=V1)) +geom_histogram(alpha=0.3, color="#9990AB", fill="#9990AB",bins=30) + 
  geom_vline(xintercept = 11, color="red", lty="dashed", size=0.75) +
  labs(y="Random sets count", x="Number of SVs intersecting 50% of a cCRE") + #theme_classic() +
  annotate("text", x = 11, y = Inf,label = "11", color = "red", fontface = "bold",vjust = 1.5, hjust = -0.5) +
  theme_bw() + theme(panel.grid = element_blank(), text = element_text(size = 13)) + ylim(c(0,2500))


### Fig. 3 ###
# Was made using the IGV Web App (https://igv.org/app/)

### Fig. 4A ###
# Was made using the IGV Web App (https://igv.org/app/) and Clustal Omega (https://www.ebi.ac.uk/jdispatcher/msa/clustalo)

### Fig. 4B ###
indep <- read.table("unique_independent_sv_eqtls.26jun26.txt", quote="\"", comment.char="")
notmostcommon <- read.delim("multiallelic_SVeQTLs.not_most_common.TRF.tsv")
trf <- read.delim("TRF.outpt.tsv")
trf$period_size <- as.numeric(trf$period_size)
bigmat <- read.delim("Gene_exp.CaG.INT_Norm.txt")

convert_to_dosage <- function(gt){
  if (gt=="0/0") {return(0)}
  if (gt=="0/1" | gt=="1/0"){return(1)}
  if (gt=="1/1"){return(2)}
}

foo <- function(group, asso_allele, period_size, gt_file, covar_file, pheno, test) {
  # Tests
  #1. Cummulative Number of Repeats ~ GE
  #2. Cummulative Number of Repeats ~ GE (with samples with 2 INS alleles)
  #3. Number of INS alleles ~ GE 
  #4. Number of INS alleles ~ GE (without samples with associated allele)
  #5. Number of INS alleles ~ GE + Cumulative Number of Repeats
  #6. Just return the data (merged df)
  #7. Cummulative Number of Repeats ~ GE + Number of associated allele
  
  # Genotypes
  gt <- read.delim(gt_file)
  colnames(gt) <- mapply(function(id){return(str_split(id, "\\.")[[1]][3])},colnames(gt))
  colnames(gt)[colnames(gt)==1] <- "ID"
  gt_long <- gt %>% pivot_longer(cols = -ID, names_to = "sample", values_to = "gt")
  gt_long$dosage <- mapply(convert_to_dosage, gt_long$gt)
  # Remove samples if at least one GT is missing 
  miss <- subset(gt_long, gt==".")
  gt_long <- subset(gt_long, !(sample %in% miss$sample))
  
  # Calculate the cummulative length of repeats per sample
  tmp_trf <- subset(trf, period_size==period_size & group == group)
  
  result <- gt_long %>% left_join(tmp_trf, by = c("ID" = "sequence")) %>%
    mutate(total_copies_per_allele = as.numeric(dosage) * as.numeric(copy_number)) %>%
    group_by(sample) %>%
    summarise(cumulative_repeats = sum(total_copies_per_allele, na.rm = TRUE)) %>%
    pivot_wider(
      names_from = sample, 
      values_from = cumulative_repeats)
  result <- data.frame(t(result))
  result$IID <- rownames(result)
  result <- result[,c(2,1)]
  colnames(result) <- c("IID", "cum_num_repeats")
  
  # Calculate number of INS alleles per sample
  allele_counts <- gt_long %>%
    mutate(dosage = unlist(dosage)) %>%
    group_by(sample) %>%
    summarise(total_alleles = sum(dosage), .groups = "drop")
  colnames(allele_counts) <- c("IID", "total_alleles")
  
  # Gene Expression Data
  mat <- bigmat[,c("Sample", pheno)]
  mat <- subset(mat, Sample %in% result$IID)
  colnames(mat) <- c("IID", pheno)
  
  # Covariates
  covar <- read.delim(covar_file)
  covar <- covar[,-1]
  
  # Merge everything
  merged <- merge(result, allele_counts)
  merged <- merge(merged, mat)
  merged <- merge(merged, covar)
  merged <- merged[, colSums(!is.na(merged)) > 0]
  
  if (test == 1) {
    # Create formula
    covars <- setdiff(names(merged), c("IID", "cum_num_repeats", pheno, "total_alleles"))
    tmp <- paste0("cum_num_repeats ~ ", pheno, " +")
    form <- as.formula(paste(tmp, paste(covars, collapse = " + ")))
    
    # Linear regression
    model <- glm(data = merged, formula = form)                        
  }
  
  if (test == 2) {
    # Modify data
    merged2 <- subset(merged, total_alleles==2)
    
    # Create formula
    covars <- setdiff(names(merged2), c("IID", "cum_num_repeats", pheno, "total_alleles"))
    tmp <- paste0("cum_num_repeats ~ ", pheno, " +")
    form <- as.formula(paste(tmp, paste(covars, collapse = " + ")))
    model <- glm(data = merged2, formula = form)
    
    return(summary(model))
  }
  
  if (test == 3) {
    # Create formula
    covars <- setdiff(names(merged), c("IID", "cum_num_repeats", pheno, "total_alleles"))
    tmp <- paste0("total_alleles ~ ", pheno, " +")
    form <- as.formula(paste(tmp, paste(covars, collapse = " + ")))
    model <- glm(data = merged, formula = form)
    return(summary(model))
  }
  
  if (test == 4) {
    # Modify data
    remove <- subset(gt_long, ID==asso_allele)
    remove <- subset(remove, dosage>0)
    remove <- remove$sample
    merged3 <- subset(merged, !(IID %in% remove))
    
    # Create formula
    covars <- setdiff(names(merged3), c("IID", "cum_num_repeats", pheno, "total_alleles"))
    tmp <- paste0("total_alleles ~ ", pheno, " +")
    form <- as.formula(paste(tmp, paste(covars, collapse = " + ")))
    model <- glm(data = merged3, formula = form)
    
    return(summary(model))
  }
  
  if (test == 5) {
    # Create formula
    covars <- setdiff(names(merged), c("IID", "total_alleles", pheno))
    tmp <- paste0("total_alleles ~ ", pheno, " +")
    form <- as.formula(paste(tmp, paste(covars, collapse = " + ")))
    model <- glm(data = merged, formula = form)
    
    return(summary(model))
  }
  
  if (test==6) {
    return(list(merged = merged, gt_long = gt_long))
  }
}

gt <- read.delim("chr10.104182044.GT")
covar <- read.delim("covariates.w_conditions.chr10.ENSG00000166275.tsv")

# Associated allele: chr10-104182044-INS-1019-25 (2.6%)
# Most common allele: chr10-104182044-INS-1023-44 (26%)
chr10 <- foo("chr10-104182044", "chr10-104182044-INS-1019-25", 2,
             "chr10.104182044.GT",
             "covariates.w_conditions.chr10.ENSG00000166275.tsv",
             "ENSG00000166275", 6)
gt_long <- chr10$gt_long
merged <- chr10$merged
samples_asso <- subset(gt_long, ID=="chr10-104182044-INS-1019-25")
samples_asso <- subset(samples_asso, dosage>0)
samples_common <- subset(gt_long, ID=="chr10-104182044-INS-1023-44")
samples_common <- subset(samples_common, dosage>0)
samples_rest <- subset(gt_long, !(sample %in% samples_asso$sample & !(sample %in% samples_common$sample)))
merged2 <- merged
merged2$status <- ifelse(merged2$IID %in% samples_asso$sample, "SV-eQTL", 
                         ifelse(merged2$IID %in% samples_common$sample, "Most Common", 
                                ifelse(merged2$IID %in% samples_common$sample & merged2$IID %in% samples_asso$sample, "Both", 
                                       ifelse(merged2$IID %in% samples_rest$sample, "Control", NA))))
merged2 <- subset(merged2, !is.na(status))
merged2$status <- factor(merged2$status,
                         levels = c("Control", "Most Common", "SV-eQTL"))

B <- ggplot(merged2, aes(x=status, y=ENSG00000166275)) + geom_point(alpha=0.6) + 
  stat_summary(fun = mean,geom = "crossbar",width = 0.5,color = "red") +
  labs(x="", y= expression(italic("BORCS7")*" expression")) + theme_minimal()



### Supp. Fig. 1 ###
library(ggplot2)
library(scales)
a <- read.delim("data_supp_fig1.txt")

SF1 <- ggplot(a, aes(x=Approach, y=Count, color=SV_type, fill=SV_type)) + 
  geom_col(position = position_dodge(width = 0.9)) +
  geom_text(
    aes(label = comma(Count)),
    position = position_dodge(width = 0.9),
    vjust = -0.3,
    size = 4,
    color = "black") +
  scale_y_continuous(labels = comma, 
                     breaks = seq(0, 500000, by = 100000),
                     limits = c(0, 500000)) +
  labs(fill="SV Type", color="SV Type")


### Supp. Fig. 2 ###
library(dplyr)
library(tidyr)
library(ggplot2)

# A #
all <- read.delim("HPRC1.SVs.ALL.afreq")
fc <- read.delim("HPRC1.SVs.FC.afreq")
h <- read.delim("HPRC1.SVs.H.afreq")
m <- read.delim("HPRC1.SVs.M.afreq")
fc <- fc[,c(2,5)]
h <- h[,c(2,5)]
m <- m[,c(2,5)]
all <- all[,c(2,5)]
colnames(fc) <- c("ID", "AF_FC")
colnames(h) <- c("ID", "AF_H")
colnames(m) <- c("ID", "AF_M")
colnames(all) <- c("ID", "AF_ALL")

af <- merge(all, fc)
af <- merge(af, h)
af <- merge(af, m)

af_long <- af %>%
  pivot_longer(cols = starts_with("AF_"),names_to = "AncestryGroup",values_to = "AF") %>% mutate(AncestryGroup = sub("AF_", "", AncestryGroup))

af_long <- af_long %>%
  mutate(AF_Category = case_when(
      AF == 0              ~ "0%",
      AF > 0 & AF <= 0.001 ~ "]0%,0.1%]",
      AF > 0.001 & AF <= 0.01 ~ "]0.1%,1%]",
      AF > 0.01 & AF <= 0.05 ~ "]1%,5%]",
      AF > 0.05            ~ ">5%"))

af_long$AncestryGroup <- ifelse(af_long$AncestryGroup=="FC", "QFC", ifelse(af_long$AncestryGroup=="M", "QMO", ifelse(af_long$AncestryGroup=="H", "QHA", "ALL")))

plot_df <- af_long %>% count(AncestryGroup, AF_Category)

plot_df <- af_long %>%
  count(AncestryGroup, AF_Category) %>%
  group_by(AncestryGroup) %>%
  mutate(pct = n / sum(n),label = paste0(round(pct,2))) %>% ungroup()

A <- ggplot(plot_df,aes(x = AncestryGroup,y = n,fill = AF_Category)) +
  geom_col(position = "fill") +
  geom_text(aes(label = label),position = position_fill(vjust = 0.5),color = "white",size = 3.5) +
  labs(x = "Ancestral Group", y = "Proportion of SVs",fill = "Allele Frequencies")

# B #
all <- read.delim("HPRC2.SVs.ALL.afreq")
fc <- read.delim("HPRC2.SVs.FC.afreq")
h <- read.delim("HPRC2.SVs.H.afreq")
m <- read.delim("HPRC2.SVs.M.afreq")
fc <- fc[,c(2,6)]
h <- h[,c(2,6)]
m <- m[,c(2,6)]
all <- all[,c(2,6)]
colnames(fc) <- c("ID", "AF_FC")
colnames(h) <- c("ID", "AF_H")
colnames(m) <- c("ID", "AF_M")
colnames(all) <- c("ID", "AF_ALL")

af <- merge(all, fc)
af <- merge(af, h)
af <- merge(af, m)

af_long <- af %>% 
  pivot_longer(cols = starts_with("AF_"),names_to = "AncestryGroup",values_to = "AF") %>%
  mutate(AncestryGroup = sub("AF_", "", AncestryGroup))

af_long <- af_long %>%mutate(AF_Category = case_when(
      AF == 0              ~ "0%",
      AF > 0 & AF <= 0.001 ~ "]0%,0.1%]",
      AF > 0.001 & AF <= 0.01 ~ "]0.1%,1%]",
      AF > 0.01 & AF <= 0.05 ~ "]1%,5%]",
      AF > 0.05            ~ ">5%"))

af_long$AncestryGroup <- ifelse(af_long$AncestryGroup=="FC", "QFC",ifelse(af_long$AncestryGroup=="M", "QMO",ifelse(af_long$AncestryGroup=="H", "QHA", "ALL")))

plot_df <- af_long %>% count(AncestryGroup, AF_Category)

plot_df <- af_long %>% count(AncestryGroup, AF_Category) %>% group_by(AncestryGroup) %>%
  mutate(pct = n / sum(n),label = paste0(round(pct,2))) %>% ungroup()

B <- ggplot(plot_df,aes(x = AncestryGroup,y = n,fill = AF_Category)) +
  geom_col(position = "fill") +
  geom_text(aes(label = label),position = position_fill(vjust = 0.5),color = "white", size = 3.5) +
  labs( x = "Ancestral Group",y = "Proportion of SVs",fill = "Allele Frequencies")

# C #
all <- read.delim("Standard.SVs.ALL.afreq")
fc <- read.delim("Standard.SVs.FC.afreq")
h <- read.delim("Standard.SVs.H.afreq")
m <- read.delim("Standard.SVs.M.afreq")
all <- all[,c(2,5)]
fc <- fc[,c(2,5)]
h <- h[,c(2,5)]
m <- m[,c(2,5)]
colnames(fc) <- c("ID", "AF_FC")
colnames(h) <- c("ID", "AF_H")
colnames(m) <- c("ID", "AF_M")
colnames(all) <- c("ID", "AF_ALL")

af <- merge(all, fc)
af <- merge(af, h)
af <- merge(af, m)

af_long <- af %>% 
  pivot_longer(cols = starts_with("AF_"),names_to = "AncestryGroup",values_to = "AF") %>%
  mutate(AncestryGroup = sub("AF_", "", AncestryGroup))

af_long <- af_long %>% mutate(AF_Category = case_when(
      AF == 0              ~ "0%",
      AF > 0 & AF <= 0.001 ~ "]0%,0.1%]",
      AF > 0.001 & AF <= 0.01 ~ "]0.1%,1%]",
      AF > 0.01 & AF <= 0.05 ~ "]1%,5%]",
      AF > 0.05            ~ ">5%"))

af_long$AncestryGroup <- ifelse(af_long$AncestryGroup=="FC", "QFC", ifelse(af_long$AncestryGroup=="M", "QMO", ifelse(af_long$AncestryGroup=="H", "QHA", "ALL")))

plot_df <- af_long %>% count(AncestryGroup, AF_Category)

plot_df <- af_long %>% 
  count(AncestryGroup, AF_Category) %>%
  group_by(AncestryGroup) %>%
  mutate(pct = n / sum(n),label = paste0(round(pct,2))) %>%ungroup()

C <- ggplot(plot_df, aes(x = AncestryGroup,y = n,fill = AF_Category)) +
  geom_col(position = "fill") +
  geom_text(aes(label = label), position = position_fill(vjust = 0.5),color = "white", size = 3.5) +
  labs(x = "Ancestral Group",y = "Proportion of SVs",fill = "Allele Frequencies")

### Supp. Fig. 3 ###
library(arrow)
library(stringr)
library(scattermore)
library(dplyr)
library(ggplot2)

combine_parquet <- function(path, prefix, filter){
  # path is the path where the parquet files are
  # prefix is what is before '.cis_qtl_pairs.chr1.parquet'
  
  output_df <- data.frame()
  
  for (i in c(1:22)) {
    file_name <- paste0(path, prefix, ".cis_qtl_pairs.chr", as.character(i), ".parquet")
    df <- read_parquet(file_name, as_data_frame = TRUE)
    
    if (i == 1){
      output_df <- df
    }
    else {
      output_df <- rbind(output_df, df)
    }
  }
  
  # Add chrX
  file_name <- paste0(path, prefix, ".cis_qtl_pairs.chrX.parquet")
  df <- read_parquet(file_name, as_data_frame = TRUE)
  
  output_df <- rbind(output_df, df)
  
  return(output_df)
}

# A (CHM13)
hp <- combine_parquet(path="/path/to/files",
                      prefix="hprc2.13nov25", filter=TRUE)
hp$CHROM <- mapply(function(id){return(str_split(id, "-")[[1]][1])}, hp$variant_id)
hp$CHROM <- gsub("chr","", hp$CHROM)
hp$CHROM <- ifelse(hp$CHROM=="X", 23, as.numeric(hp$CHROM))
hp$POS <- mapply(function(id){return(as.numeric(str_split(id, "-")[[1]][2]))}, hp$variant_id)

chr_lengths <- read.delim("chm13.lengths_chr.tsv")
genome_df <- chr_lengths %>% mutate(tot = lag(cumsum(as.numeric(chr_len)), default = 0))

df_plot <- hp %>% left_join(genome_df, by = c("CHROM" = "bin_chr")) %>% mutate(bp_cum = POS + tot)

axis_df <- genome_df %>% mutate(center = tot + chr_len / 2)

A <- ggplot(df_plot, aes(x = bp_cum, y = -log10(pval_nominal))) +
  geom_scattermore(pointsize = 6, alpha = 1, interpolate = TRUE, pixels = c(2000, 2000), aes(color = factor(CHROM)),) + 
  scale_color_manual(
    values = rep(c("black", "grey50"), length.out = length(unique(df_plot$CHROM)))) + 
  scale_x_continuous(
    breaks = axis_df$center,
    labels = axis_df$bin_chr) +
  labs(x = "Chromosome", y = expression(-log[10](P))) +
  theme_bw() + theme(legend.position = "none",panel.grid.major.x = element_blank(),panel.grid.minor.x = element_blank()) + 
  geom_hline(yintercept=-log10(5e-8), color="red", linewidth=0.8, lty="dashed")

# B (GRCh38)
st <- combine_parquet(path="C:/path/to/files",
                      prefix="standard.13nov25", filter=TRUE)

extract_chr <- function(variant_id){
  return(str_split(str_split(variant_id,":")[[1]][1], "_")[[1]][2])
}
extract_bp <- function(variant_id){
  return(as.numeric(str_split(str_split(variant_id,":")[[1]][2], "-")[[1]][1]))
}

st$CHROM <- mapply(extract_chr, st$variant_id)
st$POS <- mapply(extract_bp, st$variant_id)
st$CHROM <- ifelse(st$CHROM=='X', 23, st$CHROM)

chr_lengths <- read.delim("hg38.chr_sizes.txt", header=FALSE)
colnames(chr_lengths) <- c("bin_chr", "chr_len")
chr_lengths$bin_chr <- gsub("chr", "", chr_lengths$bin_chr)
genome_df <- chr_lengths %>% mutate(tot = lag(cumsum(as.numeric(chr_len)), default = 0))

df_plot <- st %>% left_join(genome_df, by = c("CHROM" = "bin_chr")) %>% mutate(bp_cum = POS + tot)

axis_df <- genome_df %>% mutate(center = tot + chr_len / 2)

B <- ggplot(df_plot, aes(x = bp_cum, y = -log10(pval_nominal))) +
  geom_scattermore(pointsize = 6, alpha = 1, interpolate = TRUE, pixels = c(2000, 2000), aes(color = factor(CHROM)),) + 
  scale_color_manual(
    values = rep(c("black", "grey50"), length.out = length(unique(df_plot$CHROM)))) + 
  scale_x_continuous(
    breaks = axis_df$center,
    labels = axis_df$bin_chr) +
  labs(x = "Chromosome", y = expression(-log[10](P))) +
  theme_bw() + theme(legend.position = "none",panel.grid.major.x = element_blank(),panel.grid.minor.x = element_blank()) + 
  geom_hline(yintercept=-log10(5e-8), color="red", linewidth=0.8, lty="dashed")


### Supp. Fig. 4 ###
library(ggplot2)
library(arrow)
library(dplyr)
indep <- read.delim("unique_independent_SV_eQTLs.genes_paires.txt", header=FALSE)
colnames(indep) <- c("ID", "PHENO" ,"GENE")

# Extract pvalues from before cond. analysis
combine_parquet <- function(path, prefix, filter){
  # path is the path where the parquet files are
  # prefix is what is before '.cis_qtl_pairs.chr1.parquet'
  
  output_df <- data.frame()
  
  for (i in c(1:22)) {
    file_name <- paste0(path, prefix, ".cis_qtl_pairs.chr", as.character(i), ".parquet")
    df <- read_parquet(file_name, as_data_frame = TRUE)
    df$variant_id <- gsub("chr", "", df$variant_id)
    df <- subset(df, variant_id %in% indep$ID)
    df <- subset(df, phenotype_id %in% indep$PHENO)
    
    if (i == 1){
      output_df <- df
    }
    else {
      output_df <- rbind(output_df, df)
    }
  }
  
  # Add chrX
  file_name <- paste0(path, prefix, ".cis_qtl_pairs.chrX.parquet")
  df <- read_parquet(file_name, as_data_frame = TRUE)
  df <- subset(df, variant_id %in% indep$ID)
  df <- subset(df, phenotype_id %in% indep$PHENO)
  
  output_df <- rbind(output_df, df)
  
  return(output_df)
}

cis <- combine_parquet(path="path/to/files",
                       prefix="hprc2.13nov25", filter=TRUE)
x <- read_parquet("hprc2.13nov25.cis_qtl_pairs.chrX.parquet", as_data_frame=TRUE)
x$variant_id <- gsub("chr", "", x$variant_id)
x <- subset(x, phenotype_id == "ENSG00000147164" & variant_id == "X-69505643-INS-63-5")
cis <- rbind(cis, x)

# Extract p-values after conditional analysis
asso <- read.delim("Conditonal_Analysis.Association_Results.tsv")
asso <- asso[,c("ID", "Pheno", "BETA", "P")]
colnames(asso)[colnames(asso)=="Pheno"] <- "PHENO"
asso <- subset(asso, ID %in% indep$ID)

# Merge before and after
merged <- merge(asso, cis, by.x="ID", by.y="variant_id")
merged <- subset(merged, PHENO == phenotype_id)
merged <- merged %>% semi_join(indep, by = c("ID", "PHENO"))
merged$Significance <- ifelse(merged$pval_nominal<=5e-8, "Yes", "No")

ggplot(merged, aes(y=-log10(P), x=-log10(pval_nominal), color=Significance)) + geom_point(alpha=0.5, size=2.5) +
  labs(y="-log10(P) after conditional analysis", x="-log10(P) before conditional analysis", color="Significant before\nconditional analysis") +
  geom_vline(xintercept = -log10(5e-8), color="red", lty="dotted") + 
  geom_hline(yintercept = -log10(5e-8), color="red", lty="dotted")


### Supp. Fig. 5 ###
library(ggplot2)
library(dplyr)

tested <- read.delim("HPRC2.SVs.tested_in_eQTL.afreq")
eqtl <- read.delim("HPRC2.SVs.indep_SVeQTL.afreq")

df_plot <- bind_rows(data.frame(AF = tested, Group = "SV-eQTLs"),data.frame(AF = eqtl, Group = "Tested SVs"))

SF6 <- ggplot(df_plot, aes(x = AF.ALT_FREQS, fill = Group, color = Group)) +
  geom_density(alpha = 0.5, linewidth = 1) +
  facet_wrap(~Group, ncol = 1, scales = "free_y") + 
  scale_fill_manual(values = c("red", "blue")) +
  scale_color_manual(values = c("red", "blue")) +
  labs(x = "Frequency of ALT allele",y = "Relative density", color="", fill="") +
  theme(legend.position = "none") + theme_minimal() +
  theme(axis.line = element_line(color = "black", linewidth = 0.5),
        axis.ticks = element_line(color = "black", linewidth = 0.5))

# Supp. Fig. 6 #
library(dplyr)
library(ggplot2)
library(ggbreak)
library(tidyr)

bins <- read.delim("test_eqtls.GQ.bins", header=FALSE)
colnames(bins) <- c("GQ", "count")
qtls <- read.delim("unique_indep_sv_eqtls.GQ.bins", header=FALSE)
colnames(qtls) <- c("GQ", "count")  

df_combined <- bind_rows(bins %>% mutate(Source = "Tested SVs"),qtls %>% mutate(Source = "SV-eQTLs")) %>%
  group_by(Source) %>% mutate(prop = count / sum(count)) %>% ungroup()

# A #
A <- ggplot(df_combined, aes(x = GQ, y = prop, fill = Source)) +
  geom_col(position = "identity", alpha = 0.5) +
  labs(x = "Genotype Quality (GQ)", y = "Proportion", fill="") +
  scale_x_break(c(200,9990)) +
  scale_x_continuous(breaks=c(0,50,100,150,200,10000)) + 
  theme_minimal() +
  theme(axis.text.x.top = element_blank(), axis.ticks.x.top = element_blank(),
        axis.line.x.top = element_blank(),
        axis.line = element_line(color = "black", linewidth = 0.5),
        axis.ticks = element_line(color = "black", linewidth = 0.5)) +
  scale_fill_manual(values=c("red", "blue")) +
  scale_y_break(c(0.085,0.16)) +
  scale_y_continuous(breaks=c(0,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.16,0.17))

# B #
eqtl_summary <- qtls %>% summarise(success = sum(as.numeric(count[GQ >= 30])),fail= sum(as.numeric(count[GQ <= 30])))
bins_summary <- bins %>% summarise(success = sum(as.numeric(count[GQ >= 30])),fail= sum(as.numeric(count[GQ <= 30])))
successes <- c(eqtl_summary$success, bins_summary$success)
totals <- c(eqtl_summary$success + eqtl_summary$fail,bins_summary$success + bins_summary$fail)

plot_prop_data <- data.frame(
  Dataset = c("SV-eQTLs", "Tested SVs"),
  `GQ >= 30` = c(eqtl_summary$success, bins_summary$success),
  `GQ < 30` = c(eqtl_summary$fail, bins_summary$fail),
  check.names = FALSE) %>%
  pivot_longer(cols = c(`GQ >= 30`, `GQ < 30`), names_to = "GQ", values_to = "Count") %>%
  group_by(Dataset) %>%
  mutate(Proportion = Count / sum(as.numeric(Count)))

B <- ggplot(plot_prop_data, aes(x = Dataset, y = Proportion, fill = GQ)) +
  geom_bar(stat = "identity", position = "fill", width = 0.6) +
  scale_fill_manual(values = c("GQ >= 30" = "#2ca02c", "GQ < 30" = "#d62728")) +
  labs(x = "", y = "Proportion",fill = "Genotype Quality (GQ)") +
  theme(legend.position = "right") + theme_minimal() +
  theme(axis.line = element_line(color = "black", linewidth = 0.5),
        axis.ticks = element_line(color = "black", linewidth = 0.5)) +
  geom_text(
    aes(label = round(Proportion, 2)),
    position = position_stack(vjust = 0.5),
    color = "white",
    fontface = "bold",
    size = 4)

# Supp. Fig. 7 #
library(ggplot2)

ccres <- read.delim("Conditonal_Analysis.ENCODE_4_cCREs.tsv")
ccres <- ccres[,c(1,10)]

del <- read.delim("intersect_indep_sveqtls.ccres.DEL_COMPLEX.txt", header=FALSE)
del$id <- mapply(function(x){return(str_split(x, ":")[[1]][3])}, del$V4)
del$type <- mapply(function(x){return(str_split(x, "-")[[1]][3])}, del$id)
del$len_ccre <- del$V8-del$V7
del$per_ccre <- del$V10/del$len_ccre
del <- merge(del, ccres, by.x="V9", by.y="id_cCRE")
del$pos_chm13 <- mapply(function(id){return(str_split(id, "-")[[1]][2])},del$id)
query <- read.delim("HPRC2.SVs.tested_in_eQTL.query", header=FALSE)
del <- merge(del, query, by.x="id", by.y="V3")
del$len_sv <- nchar(del$V4.y)
del$len_sv <- ifelse(del$type=="DEL", del$len_sv-1, del$len_sv)
del$id2 <- paste0(del$V1.x, "-",del$pos_chm13, "-", del$type, "-",del$len_sv)

del <- del %>%
  mutate(chr = str_extract(id2, "chr[^-]+") |> str_remove("^chr"), pos = as.numeric(str_extract(id2, "-\\d+") |> str_remove("-"))) %>%
  mutate(chr = factor(chr, levels = c(as.character(1:22), "X", "Y", "M"))) %>%arrange(chr, pos)

del$id2 <- factor(del$id2, levels = unique(del$id2))
del<- del[!(rownames(del) %in% as.character(85:88)), ]
del<- del[!(rownames(del) %in% as.character(90:165)), ]
del<- del[!(rownames(del) %in% as.character(82:83)), ]
del<- del[!(rownames(del) %in% as.character(168:170)), ]

SF7 <-ggplot(del, aes(x = per_ccre, y = V9, fill = class_cCRE)) +
  geom_col(width = 0.8, color = "black") + 
  coord_cartesian(xlim = c(0, 1), clip = "off") +
  facet_grid(rows = vars(id2), scales = "free_y", space = "free_y", switch = "y") +
  labs(x = "Proportion of cCRE overlapped", y = "SV", fill = "cCRE class") +
  theme_grey() +
  theme(
    axis.text.y = element_blank(),        
    axis.ticks.y = element_blank(),       
    strip.placement = "outside",          
    strip.background = element_blank(),   
    strip.text.y.left = element_text(angle = 0, hjust = 1, size = 9), 
    panel.spacing = unit(0.1, "lines"),   
    panel.grid.major.y = element_blank(),
    strip.clip = "off",
    plot.margin = margin(t = 10, r = 10, b = 10, l = 150, unit = "pt"))

### Supp. Fig 8. ###
indep <- read.table("unique_independent_sv_eqtls.26jun26.txt", quote="\"", comment.char="")
asso <- read.delim("Conditonal_Analysis.Association_Results.beta_for_minor_allele.tsv")
asso <- asso[,-c(10,11)]
asso <- subset(asso, asso$ID %in% indep$V1)
vars <- unique(asso$ID) 
extract_multi <- function(ids) {
  output <- c()
  for (i in ids) {
    i_ <- strsplit(i, "-")
    if (!(is.na(i_[[1]][5]))){
      output <- c(output, i)
    }
  }
  return(output)
}

multi <- extract_multi(vars)
multi_asso <- subset(asso, ID %in% multi) 
multi_asso$group <- paste(multi_asso$chr, multi_asso$sv_start, sep="-")

extract_group <- function(id) {
  id_ <- strsplit(id, '-')
  output <- paste(id_[[1]][1], id_[[1]][2], sep='-')
  return(output)
}

extract_pos <- function(id){
  id_ <- strsplit(id, '-')
  return(id_[[1]][2])
}

af_all <- read.delim("HPRC2.SVs.multiallelic.AFREQ")
af_all$group <- mapply(extract_group, af_all$ID)
af_all$POS <- mapply(extract_pos, af_all$ID)
af_all$ID <- gsub("chr", "", af_all$ID) 
af_all <- subset(af_all, group %in% multi_asso$group)

bim_500 <- read.delim("HPRC2.SVs.tested_in_eQTL.bim", header=FALSE)
af_all_in_500 <- subset(af_all, ID %in% bim_500$V2)
af_all_asso <- subset(af_all, ID %in% asso$ID)

af_all_asso <- af_all_asso %>%
  mutate(chr = str_extract(group, "(?<=chr)[^-]+"),chr = factor(chr,levels = c(as.character(1:22), "X", "Y", "MT")),
    pos = as.integer(str_extract(group, "(?<=-)\\d+"))) %>%
  arrange(chr, pos) %>%
  mutate(group = factor(group, levels = unique(group)))

SF8 <- ggplot(af_all, aes(x=ALT_FREQS, y=group))+geom_point(alpha=0.7, size=1, color="darkgrey") +
  geom_point(data=af_all_in_500, aes(ALT_FREQS, y=group), color="black", alpha=1, size=1) +
  geom_point(data=af_all_asso, aes(x=ALT_FREQS, y=group), color="blue", alpha=1, size=1) +
  geom_vline(xintercept=0.01, color="red") +  theme(axis.text.y = element_text(size = 5)) +
  labs(y="Loci", x="AF of ALT alleles")


### Supp. Fig. 9 ###
library(stringi)
library(stringr)
library(ggplot2)
library(dplyr)
library(tidyr)

ins <- read.delim("TRF_INS.tsv")
del <- read.delim("TRF_DEL.tsv")
comref <- read.delim("TRF_COMPLEX_REF.tsv")
comalt<- read.delim("TRF_COMPLEX_ALT.tsv")
comref$allele <- "REF"
comalt$allele <- "ALT"
ins$allele <- NA
del$allele <- NA

tr <- rbind(ins, del, comref, comalt)
tr <- tr[,c(2:8,9)]
tr$type <- mapply(function(id){return(str_split(id, "-")[[1]][3])}, tr$sequence)
tr$group <- mapply(function(id){return(paste0(str_split(id, "-")[[1]][1],"-", str_split(id, "-")[[1]][2]))},tr$sequence)
tr <- tr[,c(10,1:9)]

associated <- read.table("Conditonal_Analysis.Association_Results.beta_for_minor_allele.tsv", quote="\"", comment.char="")
associated$V1 <- paste0("chr",associated$V1)
associated_tr <- subset(tr, sequence %in% associated$V1)

tested <- read.table("HPRC2.SVs.tested_in_eQTL.ids", quote="\"", comment.char="")
tested$V1 <- paste0("chr", tested$V1)
trt <- subset(tr, tr$sequence %in% tested$V1) 
trt <- trt[,c(10,1:9)]
trt <- trt[order(trt$group, trt$period_size),] 
trt <- subset(trt, group %in% multi_asso$group)

extract_group <- function(id) {
  id_ <- strsplit(id, '-')
  output <- paste(id_[[1]][1], id_[[1]][2], sep='-')
  return(output)
}

extract_pos <- function(id){
  id_ <- strsplit(id, '-')
  return(id_[[1]][2])
}

gt <- read.delim("hr10.104182044.GT")
colnames(gt) <- mapply(function(id){return(str_split(id, "\\.")[[1]][3])},colnames(gt))
colnames(gt)[colnames(gt)==1] <- "ID"
gt_long <- gt %>% pivot_longer(cols = -ID, names_to = "sample", values_to = "gt")

miss <- subset(gt_long, gt==".")
gt_long <- subset(gt_long, !(sample %in% miss$sample)) # 2,171 samples remaining

convert_to_dosage <- function(gt){
  if (gt=="0/0") {return(0)}
  if (gt=="0/1" | gt=="1/0"){return(1)}
  if (gt=="1/1"){return(2)}
}

gt_long$dosage <- mapply(convert_to_dosage, gt_long$gt)

result <- gt_long %>% left_join(tr, by = c("ID" = "sequence")) %>%
  mutate(total_copies_per_allele = as.numeric(dosage) * as.numeric(copy_number)) %>%
  group_by(sample) %>%
  summarise(cumulative_repeats = sum(total_copies_per_allele, na.rm = TRUE)) %>%
  pivot_wider(
    names_from = sample, 
    values_from = cumulative_repeats)
result <- data.frame(t(result))
result$IID <- rownames(result)
result <- result[,c(2,1)]
colnames(result) <- c("IID", "cum_num_repeats")
 
bigmat <- read.delim("Gene_exp.CaG.INT_Norm.txt")
mat <- bigmat[,c("Sample", "ENSG00000166275")]
mat <- subset(mat, Sample %in% result$IID)
colnames(mat) <- c("IID", "ENSG00000166275")
covar <- read.delim("covariates.w_conditions.chr10.ENSG00000166275.tsv")
rownames(covar) <- covar[,1]
covar <- covar[,-1]

# A #
merged <- merge(result, mat)
merged <- merge(merged, covar)

covars <- setdiff(names(merged), c("IID", "cum_num_repeats", "ENSG00000166275"))
form <- as.formula(paste("ENSG00000166275 ~ cum_num_repeats +", paste(covars, collapse = " + ")))

model <- glm(data = merged, formula = form)
summary(model)

A <- subset(gt_long, ID=="chr10-104182044-INS-1019-25" & dosage > 0)
merged_s <- subset(merged, IID %in% s$sample) # 27 samples with at least one associated allele
ggplot(merged, aes(x=cum_num_repeats, y=ENSG00000166275)) + geom_point(alpha=0.5, size=2) +
  geom_smooth(method = "lm", se = FALSE, color = "red") + 
  labs(x="Cumulative number of repeats", y = expression(italic("BORCS7") ~ "expression")) +
  geom_point(data=merged_s, aes(x=cum_num_repeats, y=ENSG00000166275), color="blue", alpha=0.9, size=2) +
  theme_bw() + scale_x_continuous(breaks = seq(0, 60, by = 5))

# B #
gt_long2 <- gt_long %>% group_by(sample) %>%
  filter(sum(dosage, na.rm = TRUE) == 2) %>% ungroup()

result2 <- gt_long2 %>% left_join(tr, by = c("ID" = "sequence")) %>%
  mutate(total_copies_per_allele = as.numeric(dosage) * as.numeric(copy_number)) %>%
  group_by(sample) %>%
  summarise(cumulative_repeats = sum(total_copies_per_allele, na.rm = TRUE)) %>%
  pivot_wider(
    names_from = sample, 
    values_from = cumulative_repeats)
result2 <- data.frame(t(result2))
result2$IID <- rownames(result2)
result2 <- result2[,c(2,1)]
colnames(result2) <- c("IID", "cum_num_repeats")
merged2 <- merge(result2, mat)
merged2 <- merge(merged2, covar)
covars <- setdiff(names(merged2), c("IID", "cum_num_repeats", "ENSG00000166275"))
form <- as.formula(paste("ENSG00000166275 ~ cum_num_repeats +", paste(covars, collapse = " + ")))
model <- glm(data = merged2, formula = form)
summary(model)

B <- subset(gt_long2, ID=="chr10-104182044-INS-1019-25" & dosage > 0)
merged_s <- subset(merged2, IID %in% s$sample) # 21 samples
ggplot(merged2, aes(x=cum_num_repeats, y=ENSG00000166275)) + geom_point(alpha=0.5, size=2) +
  geom_smooth(method = "lm", se = FALSE, color = "red") + 
  labs(x="Cumulative number of repeats", y = expression(italic("BORCS7") ~ "expression")) +
  geom_point(data=merged_s, aes(x=cum_num_repeats, y=ENSG00000166275), color="blue", alpha=0.9, size=2) +
  theme_bw() + scale_x_continuous(breaks = seq(30, 60, by = 1))

# C #
allele_counts <- gt_long %>%
  mutate(dosage = unlist(dosage)) %>%
  group_by(sample) %>%
  summarise(total_alleles = sum(dosage), .groups = "drop")

merged3 <- merge(allele_counts, mat, by.x="sample", by.y="IID")
merged3 <- merge(merged3, covar, by.x="sample", by.y="IID")
covars <- setdiff(names(merged3), c("sample", "total_alleles", "ENSG00000166275"))

form <- as.formula(paste("ENSG00000166275 ~ total_alleles +", paste(covars, collapse = " + ")))
model <- glm(data = merged3, formula = form)
summary(model)

C <- subset(gt_long, ID=="chr10-104182044-INS-1019-25" & dosage > 0)
merged_s <- subset(merged3, sample %in% s$sample) # 21 samples
ggplot(merged3, aes(x=total_alleles, y=ENSG00000166275)) + geom_point(alpha=0.5, size=2) +
  geom_smooth(method = "lm", se = FALSE, color = "red") + 
  labs(x="Number of alternate alleles", y = expression(italic("BORCS7") ~ "expression")) +
  geom_point(data=merged_s, aes(x=total_alleles, y=ENSG00000166275), color="blue", alpha=0.9, size=2) +
  theme_bw() + scale_x_continuous(breaks = seq(0, 2, by = 1))
