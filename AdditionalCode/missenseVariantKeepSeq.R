####Wait to embed sequence vectors for Python instead####

library(VariantAnnotation)
library(RSQLite)
library(tidyverse) 
library(dbplyr)
library(magrittr)
library(stringr)
library(biomaRt)
library(Biostrings)
library(dplyr)
library(data.table)


applyVariant <- function(gene1, gene2) {
  
  geneList = c(gene1, gene2)
  db <- DBI::dbConnect(dbDriver("SQLite"), "/sdata/Simons/SPARK/DATA/SPARK_Freeze_Three_update_20181105/vcfdb/spark.freeze3.db") 
  cleanGeneVar <- function(geneList,db) {
    
    
    myList <- list()
    for (gene in geneList) {
      
      t <- gene
      #t <- 'ACOT7'
      #t <- 'CHD8'
      var_impact <- tbl(db, "variant_impact")
      
      gene_Var <-  var_impact %>%
        filter(symbol == t) %>%  ## filter() operates on rows
        collect()
      
      gene_useful <- gene_Var[c(3,4,6,7,11,14,15,16,29,32,36,49)]
      
      step_1 <- separate(gene_useful, amino_acids, into = c("Old_aa", "New_aa"), sep = "\\/")
      
      step_2 <- separate(step_1, protein_position, into = c("protein_pos", "seq_length"), sep = "\\/")
      
      step_3 <- separate(step_2, codons, into = c("Old_codon", "New_codon"), sep = "\\/")
      
      step_4 <- step_3[!is.na(step_3$New_aa),]
      
      step_5 <- subset(step_4, nchar(as.character(New_codon)) == 3)
      
      step_6 <- step_5 %>% 
        filter(str_detect(ensp, 'ENSP'))
      
      tmp <- step_6[!duplicated(step_6),]
      
      myList[[t]] <- tmp
    }
    full <- do.call(rbind, myList)
    ##DBI::dbDisconnect(db)
  }
  geneVariants <- cleanGeneVar(geneList,db)
  #Ignore warning#
  DBI::dbDisconnect(db)
  rownames(geneVariants) <- c()
  
  #Run to grab uniport ids based on list of uniparc ids of gene variants
  uniparcList <- geneVariants$uniparc
  grabUniprotID <- function(uniparcList) {
    
    id_list <- uniparcList
    listMarts()
    mart = useMart('ENSEMBL_MART_ENSEMBL','hsapiens_gene_ensembl')
    dataset = useDataset('hsapiens_gene_ensembl', mart)
    
    #note not all are human
    uniprot_id <- getBM(attributes = c('uniprotswissprot','uniparc'),
                        filter='uniparc',
                        values= id_list,
                        mart=mart)
  }
  
  variantIDs <- grabUniprotID(uniparcList)
  
  #match uniparc in geneVaraints to uniswiss in variantIDs
  matchProt <- function(variantIDs, geneVariants) {
    protDF <- dplyr::inner_join(variantIDs, geneVariants, by = "uniparc")
    names(protDF)[1] <- "Protein"
    protDF <- protDF[!(is.na(protDF$Protein) | protDF$Protein==""), ]
  }
  
  protein_var <- matchProt(variantIDs, geneVariants)
  
  
  #match sequence data from uniprot to respective protein id in variant data
  matchSeq <- function(protein_var, seq) {
    data <- data.frame("Protein" = (1:561176), "Sequence" = (1:561176), stringsAsFactors = FALSE)
    
    data$Protein <- seq %>% filter(str_detect(V1, ">")) %>% unlist() %>% as.vector()
    
    data$Sequence <- seq %>% filter(str_detect(V1,">", negate = TRUE)) %>% unlist() %>% as.vector()
    
    step1 <- separate(data, Protein, into = c("Drop1", "Protein", "Rest"), sep = "\\|")
    
    uprt <- step1[c(2,4)]
    
    seqDF <- dplyr::inner_join(uprt, protein_var, by = "Protein")
    
  }
  
  seq <- read.table("uniprot_sprot_flat.txt", sep="\t", header=FALSE,
                    na.strings=".", stringsAsFactors=FALSE,
                    quote="", fill=FALSE)
  
  fullProt <- matchSeq(protein_var, seq)
  
  
  fullProt$protein_pos <- as.numeric(fullProt$protein_pos)
  
  uni <- fullProt %>% distinct()
  
  fullProt$seq_matrix <-NA
  fullProt$variantSeq_matrix <- NA
  
  ##########################################################################################################################################
  ############Cleaned up to here, work on below to get into function(s) for the whole pipeline/allow for list of genes######################
  ##########################################################################################################################################
  
  #gene1 = "ACOT7"
  #gene2 = "CHD8"
  
  ################keep one instance of gene2###################
  gene1df <- fullProt[fullProt$symbol==gene1,]
  gene2df <- fullProt[fullProt$symbol==gene2,]
  one <- gene2df[1,]
  
  fullProt <- dplyr::bind_rows(gene1df,one)
  
  
  ######need to figure out why assignning same protein sequence###################
  for(v_1 in fullProt$hgvsp[fullProt$symbol==gene1]) {
    
    #gene2 = "CHD8"
    #iter 1:
    #v_1 = 'ENSP00000402532.1:p.Asn227Ser'
    
    #iter 2:
    #v_1 = 'ENSP00000402532.1:p.Lys225Asn'
    s_1 = fullProt$Sequence[fullProt$symbol==gene1 & fullProt$hgvsp==v_1] #ACOT7, just change this position,pick one vairant
    s_2 = fullProt$Sequence[fullProt$symbol==gene2] #CHD8
    print(s_1)
    print(s_2)
    
    s1 <- BStringSet(s_1)
    s2 <- BStringSet(s_2)
    
    x1 = as.matrix(stackStrings(s1,from=1,to=2500,shift=pmax(0,2500-width(s1))/2))
    
    x2 = as.matrix(stackStrings(s2,from=1,to=2500,shift=pmax(0,2500-width(s2))/2))
    
    ### getting the coordinates of the changes
    tmp = x1!=" " 
    st = apply(tmp,1,function(x) which(x)[1]) - 1
    pos = fullProt$protein_pos[fullProt$symbol==gene1 & fullProt$hgvsp==v_1]
    print(pos)
    pos = as.numeric(pos)
    idx = cbind(1:length(st),(st+pos)) ## 'pos' is the relative position of the AA change within the protein sequence
    
    xa = x1
    xa[idx] = fullProt$New_aa[fullProt$symbol==gene1 & fullProt$hgvsp==v_1] ## 'AAto' contains the AA that the missense variant changes *to*
    
    # load("AA_to_int_map.Rdata")
    # x1 = matrix(map[x1],nrow(x1),ncol(x1),dimnames=list(rownames(x1),1:2500)) ## matrix (integer) of protein seqs of interactor 1
    # x2 = matrix(map[x2],nrow(x2),ncol(x2),dimnames=list(rownames(x2),1:2500)) ## matrix (integer) of protein seqs of interactor 2
    # xa = matrix(map[xa],nrow(xa),ncol(xa),dimnames=list(rownames(xa),1:2500)) ## x1, with variant spiked in
    # 
    # 
    # x_1 = list(x1)
    # x_2 = list(x2)
    # x_a = list(xa)
   
    
    fullProt$seq_matrix[fullProt$symbol==gene1 & fullProt$hgvsp==v_1] = toString(x1)
    fullProt$seq_matrix[fullProt$symbol==gene2] = toString(x2)
    fullProt$variantSeq_matrix[fullProt$symbol==gene1 & fullProt$hgvsp==v_1] = toString(xa)
    fullProt$variantSeq_matrix[fullProt$symbol==gene2] = toString(x2)
    
  }
  
  new_variant_data <- fullProt[!duplicated(fullProt$variantSeq_matrix), ]
  
  
  gene_tab <- read.csv(file="/sdata/Simons/SPARK/DATA/SPARK_Freeze_Three_update_20181105/SPARK.30K.mastertable.20190208.csv", sep=",", header = TRUE)
  
  gene_Sibling <-  gene_tab[gene_tab$asd == 1 & gene_tab$role =='Sibling', 'spid'] %>% as.character()
  gene_Proband <-  gene_tab[gene_tab$asd == 2 & gene_tab$role == 'Proband', 'spid'] %>% as.character()
  
  var_of_interest = unique(new_variant_data$variant_id)
  
  library(dplyr)
  
  db <- DBI::dbConnect(dbDriver("SQLite"), "/sdata/Simons/SPARK/DATA/SPARK_Freeze_Three_update_20181105/vcfdb/spark.freeze3.db") 
  
  genos <- tbl(db, "variant_info") %>%
    dplyr::filter(variant_id %in% var_of_interest) %>%
    dplyr::select(geno) %>%
    collect() %>%
    unlist() %>% 
    data.table::tstrsplit(.,'/Dedicated/jmichaelson-') %>% 
    .[[2]] %>% 
    paste('/',.,sep='') %>% 
    as.data.frame() %>% 
    setNames(.,nm='geno') %>% 
    mutate(geno = map(as.character(geno), read_rds))
  
  
  #####from here get variant_id one by one and count how many times it shows up
  geno_sib = genos %>% 
    unnest('geno') %>% 
    dplyr::filter(sample %in% gene_Sibling) %>%
    dplyr::filter(gt != 0) %>%
    dplyr::select(variant_id) %>%
    unique()
  
  geno_prob = genos %>% 
    unnest('geno') %>% 
    dplyr::filter(sample %in% gene_Proband) %>%
    dplyr::filter(gt != 0) %>%
    dplyr::select(variant_id) %>%
    unique()
  
  foo = geno_prob
  geno_prob = geno_prob[!geno_prob %in% geno_sib]
  geno_sib = geno_sib[!geno_sib %in% foo]
  
  #this was wrong because it assume 'both' status when some variant_ids were not present in proband or sibling 
  new_variant_data$source = 'both'
  new_variant_data$source[new_variant_data$variant_id %in% geno_prob$variant_id] = 'proband'
  new_variant_data$source[new_variant_data$variant_id %in% geno_sib$variant_id] = 'sib'
  
  
  return(new_variant_data)
  DBI::dbDisconnect(db)
}

gene1 = "ACOT7"
gene2 = "CHD8"
python_ready_variants <- applyVariant(gene1, gene2)


library(stringr)
python_ready_variants$seq_matrix <- str_replace_all(python_ready_variants$seq_matrix, ",", "")
python_ready_variants$seq_matrix <- str_replace_all(python_ready_variants$seq_matrix, " ", "")

python_ready_variants$variantSeq_matrix <- str_replace_all(python_ready_variants$variantSeq_matrix, ",", "")
python_ready_variants$variantSeq_matrix <- str_replace_all(python_ready_variants$variantSeq_matrix, " ", "")

sapply(python_ready_variants, function(x) length(unique(x)))


write.csv(python_ready_variants, file="testSeqs.csv")





