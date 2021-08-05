setwd("D:\\DATA\\pheweb_UKBB")
# install.packages("RSQLite")
#install.packages('Rcpp')
#install.packages("limma")

library(data.table)
dat<-fread(file = "reformat_all_withPhenoDetails.txt" )

head(dat)
length(unique(dat$pheno_code))

# Generate pheno table 
first_element = aggregate(dat, by=list(dat$pheno_code), FUN=first)
pheno_info = first_element[,2:6]
pheno_info$group_num = 1
pheno_info$group_num[pheno_info$pheno_group=="self-report"]=2
pheno_info$group_num[pheno_info$pheno_group=="quantitative"]=3
color_all=c("blue", "green", "red")
Pinfo = data.frame(id = 1:nrow(pheno_info),
                   phecode=substr(pheno_info$pheno_code, 2,1000),
                   phenostring = pheno_info$pheno_string,
                   category = pheno_info$pheno_group,
                   num_cases = pheno_info$Ncase, 
                   num_controls = pheno_info$Ncontrol,
                   sex = "Both",
                   color=color_all[pheno_info$group_num],
                   OrgCode=pheno_info$pheno_code,
                   groupnum =  pheno_info$group_num)

# Obtain gene info 
# Use both refseq and previously constructed DB since some of alias are not found in refseq..
gene_names_infile = data.frame(name = unique(dat$Gene))
dim(gene_names_infile)

refseq<-fread(file = "refseq_hg38_curated.txt" )
first_element_refseq = aggregate(refseq, by=list(refseq$name2), FUN=first)
gene_names_refseq = data.frame(name = first_element_refseq$name2, 
                               chrom = substr(first_element_refseq$chrom, 4,1000),
                               startpos    = first_element_refseq$txStart, 
                               endpos = first_element_refseq$txEnd)   


library(DBI)
library(RSQLite)
DBFile = "./assoc.db"
mydb <- dbConnect(SQLite(), DBFile)

res <- dbSendQuery(mydb, "SELECT * FROM gene")
gene_names_indb = dbFetch(res)
gene_names_indb$startpos=0
gene_names_indb$endpos=0
for(i in 1:nrow(gene_names_indb)){
  gene_id = gene_names_indb$id[i]
  query = sprintf("SELECT startpos, endpos FROM assoc where gene_id=%d LIMIT 1", gene_id)
  res1 <- dbSendQuery(mydb, query)
  out_res1 <- dbFetch(res1)
  gene_names_indb$startpos[i] = out_res1$startpos
  gene_names_indb$endpos[i] = out_res1$endpos
  
}
gene_names_indb = gene_names_indb[,-1]
out1 = setdiff(gene_names_infile$name, gene_names_refseq$name)
out2 = setdiff(out1, gene_names_indb$name)

gene_names_ref<-gene_names_refseq
gene_names_ref<-rbind(gene_names_ref, gene_names_indb[gene_names_indb$name %in% out1,])

Ginfo = merge(gene_names_infile, gene_names_ref, by="name")
Ginfo$id = 1:nrow(Ginfo)

# remove _GL000256v2_alt ...
chr = strsplit(Ginfo$chrom,"-")
chr_new = NULL
for(i in 1:length(chr)){
  chr_new[i] = chr[[i]][1]
}
Ginfo$chrom = chr_new
sapply(Pinfo, class)
dbDisconnect(mydb)


########################################
#
#
# Write into sqlite...
#
#
########################################

DBFile_New = "./assoc_new.db"
mydb_new <- dbConnect(SQLite(), DBFile_New)

dbExecute(mydb_new, 'drop table pheno')
dbExecute(mydb_new, 'drop table gene')
dbExecute(mydb_new, 'drop table assoc')
dbExecute(mydb_new, 'drop table assoc_group')

dbExecute(mydb_new, 
  'create table pheno (id INT PRIMARY KEY, 
  phecode TEXT, phenostring TEXT, category TEXT, 
  num_cases INT, num_controls INT, 
  sex TEXT, color TEXT)'
)

dbExecute(mydb_new, 
          'create table gene (
          id INT PRIMARY KEY, 
          name TEXT, chrom TEXT)'
)

dbExecute(mydb_new, 
          'create table assoc (
          id INT PRIMARY KEY, pheno_id INT, gene_id INT, 
          pval REAL, num_rare INT, startpos INT, endpos INT, 
          FOREIGN KEY(pheno_id) REFERENCES pheno(id), 
          FOREIGN KEY(gene_id) REFERENCES gene(id))'
)

dbExecute(mydb_new, 
          'create table assoc_group (
          id INT PRIMARY KEY, assoc_id INT, 
          description TEXT,
          pval REAL, 
          FOREIGN KEY(assoc_id) REFERENCES assoc(id))'
)

dbExecute(mydb_new,
          'create INDEX idx_assoc_id on assoc_group(assoc_id)'
)


dbWriteTable(mydb_new, "pheno", Pinfo[,1:8], overwrite=FALSE, append=TRUE)

Ginfo_db<-Ginfo[,c(5,1,2)]
dbWriteTable(mydb_new, "gene", Ginfo_db, overwrite=FALSE, append=TRUE)

#check 
#res <- dbSendQuery(mydb_new, "SELECT * FROM pheno")
#dbFetch(res)


##########################################################
# Write in assoc 
#######################################
# P-value check...
IDX_group_p<-8:16
N_group<-length(IDX_group_p)
Group_Name<-substr(colnames(dat)[IDX_group_p], 9, 1000)

gene_id_next = 1
gene_group_id_next=1

for(i in 1:nrow(Pinfo)){
  
  if(i==1){
    gene_id_next = 1
    gene_group_id_next=1  
  }
  code = Pinfo$OrgCode[i]
  dat1 = dat[dat$pheno_code == code]
  dat2 = merge(dat1, Ginfo, by.x="Gene", by.y="name")
  
  #startpos = sprintf("%s:%d", dat2$chrom, dat2$startpos)
  #endpos = sprintf("%s:%d", dat2$chrom, dat2$endpos)
  
  n1<-nrow(dat2)
  #id INT PRIMARY KEY, pheno_id INT, gene_id INT, pval REAL, num_rare INT, startpos INT, endpos INT, 
  assoc_df = data.frame(id = gene_id_next:(gene_id_next+ n1 -1),
                         pheno_id = Pinfo$id[i], gene_id = dat2$id,
                         pval = dat2$P_value,
                         startpos = dat2$startpos, endpos = dat2$endpos)
  

  MAT<-NULL
  for(j in 1:N_group){
    MAT<-cbind(MAT, dat2[[IDX_group_p[j]]])
  }
  Pval_vec<-as.vector(t(MAT))
  #id INT PRIMARY KEY, assoc_id INT, description TEXT, pval REAL, 
  n1_group = n1*N_group
  assoc_group_df = data.frame(id = gene_group_id_next:(gene_group_id_next+ n1_group-1), 
                               assoc_id = rep(assoc_df$id, each=N_group), 
                               description = rep(Group_Name, n1),
                               pval = Pval_vec)
  
  gene_id_next= gene_id_next+ n1 
  gene_group_id_next =  gene_group_id_next+ n1_group
  
  # Insert to DB 
  dbWriteTable(mydb_new, "assoc", assoc_df, overwrite=FALSE, append=TRUE)
  dbWriteTable(mydb_new, "assoc_group", assoc_group_df, overwrite=FALSE, append=TRUE)
  
}
dbDisconnect(mydb_new)


###########################################
# Check 

DBFile_New = "./assoc_new.db"
mydb_new <- dbConnect(SQLite(), DBFile_New)
dbListTables(mydb_new)

res <- dbSendQuery(mydb_new, "SELECT count(*) FROM pheno")
dbFetch(res)

res <- dbSendQuery(mydb_new, "SELECT count(*) FROM gene")
dbFetch(res)

res <- dbSendQuery(mydb_new, "SELECT count(*) FROM assoc_group")
dbFetch(res)

res <- dbSendQuery(mydb_new, "SELECT count(*) FROM assoc")
dbFetch(res)
dbDisconnect(mydb_new)
