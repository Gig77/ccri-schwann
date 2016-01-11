# read input data
m <- read.delim("/mnt/projects/schwann/results/anduril/execute/qcReport-samplesClusterHeatmap/vst.csv", check.names = F, stringsAsFactors = F)
ann.var <- read.delim("/mnt/projects/schwann/data/ensembl.homo_sapiens_75_37.geneAnnotations.tsv", check.names = F, stringsAsFactors = F)
ann.sample <- read.delim("/mnt/projects/schwann/data/samples.csv", check.names = F, stringsAsFactors = F)
groups <- read.delim("/mnt/projects/schwann/results/anduril/execute/sampleGroups/table.csv", check.names = F, stringsAsFactors = F)

# add group memberships as columns to sample annotation
for(g in groups$ID) {
  ann.sample[,g] <- NA
  members <- unlist(strsplit(groups$Members[groups$ID==g], ","))
  ann.sample[,g] <- ifelse(ann.sample$Alias %in% members, "yes", "no")
}

# annotate expression matrix
sample.names <- colnames(m)[-1]
ann.var.names <- c("Ensembl", "Description", "Band", "HGNC", "Chr", "Start", "End")
names(ann.var) <- ann.var.names
m.ann <- merge(m, ann.var[,ann.var.names], all.x = T)
m.ann <- m.ann[,c(ann.var.names, sample.names)]

# combine sample annotation with expression matrix into single data frame
ann.sample <- ann.sample[match(sample.names, ann.sample$Alias),]
ann.sample.transposed <- as.data.frame(t(ann.sample), stringsAsFactors = F)
colnames(ann.sample.transposed) <- ann.sample$Alias
ann.sample.transposed$Attribute <- rownames(ann.sample.transposed)

qlucore <- merge(m.ann, ann.sample.transposed, by.x=c("Ensembl", sample.names), by.y=c("Attribute", sample.names), all=T)
qlucore <- qlucore[order(qlucore$Description, na.last=F),c(ann.var.names, sample.names)]
qlucore <- rbind(qlucore[qlucore$Ensembl %in% ann.sample.transposed$Attribute,], c(ann.var.names, rep(NA, length(sample.names))), qlucore[!qlucore$Ensembl %in% ann.sample.transposed$Attribute,])
qlucore <- cbind(qlucore[,ann.var.names], data.frame(Annotation=NA), qlucore[,sample.names])
qlucore$Annotation <- as.character(ifelse(qlucore$Ensembl %in% ann.sample.transposed$Attribute, qlucore$Ensembl, NA))
qlucore$Ensembl <- ifelse(qlucore$Ensembl %in% ann.sample.transposed$Attribute, NA, qlucore$Ensembl)

# put annotation with sample id and alias as first row
qlucore <- qlucore[order(qlucore$Annotation %in% c("ID", "Alias"), decreasing = T),]

write.table(qlucore, "/mnt/projects/schwann/results/qlucore/schwann.qlucore.expression-matrix.tsv", quote=F, sep="\t", col.names=F, row.names=F)
