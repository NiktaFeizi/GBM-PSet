#GBM PSets:

# Loading packages
library(openxlsx)
library(stringr)
library(Biobase)
library(GEOquery)#For unpakcing a tar file
library(hglue20hsensgcdf)
library(affy)
library(oligo)
library(affycoretools)
library(pd.hta.2.0)
library(RUVnormalize)
library(RUVnormalizeData)
library(PharmacoGx)
library(abind)
library(data.table)

screen2<-read.delim("Data/Screen2-drugData.txt", stringsAsFactors = FALSE)#Drug_dose_scr2 response data
screen3<-read.delim("Data/Screen3-drugData.txt", stringsAsFactors = FALSE)#Drug_dose_scr3 response data
cell<-read.xlsx("Data/mmc2-cell_lines.xlsx",rowNames = TRUE)#Cell_line names
load("/Data/Ensembl.v99.annotation.RData")

#Functions
######################################################### drug_name_correction ################################################
#Drug_names correction function
#Function to find the synonym drug names by not considering the salts and unwanted puncutations

drug_name_correction<-function(to_be_corrected_table, to_be_corrected_col){
  
  badchars <- "[\xb5]|[]|[ ,]|[;]|[:]|[-]|[+]|[*]|[%]|[$]|[#]|[{]|[}]|[[]|[]]|[|]|[\\^]|[/]|[\\]|[.]|[_]|[ ]|[(]|[)]|[8]"
  salts <- c("malate","sulfate","dihydrochloride","hydrochloride","citrate", "acid","ethanolamine",
             "2hcl", "oxalate" , "sodiumsalt" ,"bromide" , "monohydrate" ,"isethionate","sodium",
             "furoate","hcl","ca","dimaleate", "oxalate","dihydrate","maleate","fumarate","lactate","di")
  
  to_be_corrected_table$clean_drugID<-gsub("\\s*\\([^\\)]+\\)","",to_be_corrected_table[,to_be_corrected_col],ignore.case=T)
  to_be_corrected_table$clean_drugID<-gsub(badchars,"", to_be_corrected_table$clean_drugID, ignore.case=T)
  for (salt in salts){ 
    to_be_corrected_table$clean_drugID<- gsub(salt,"", to_be_corrected_table$clean_drugID, ignore.case=T)
    to_be_corrected_table$clean_drugID<- tolower(to_be_corrected_table$clean_drugID)
  }
  return(to_be_corrected_table)
}

######################################################### fdata_builder ######################################################
#Function to annotate the genes when creating feature data for a PSet. If annotation for a gene is not avialble 
#it will not be removed but its annotation will be presented by "NA".

fdata_builder<-function(annotation, assay,ID_column="gene_name"){
  feat_empty<-data.frame(matrix(NA, nrow = nrow(assay), ncol = ncol(annotation)))
  colnames(feat_empty)<-colnames(annotation)
  feat_empty[,ID_column]<-rownames(assay)
  
  annotated<-annotation[annotation[,ID_column] %in% rownames(assay),] #Subseting the genes from features_gene that belong to cnv
  feat<-merge(x=feat_empty, y=annotated, all=TRUE) #we we keep all the genes even the ones that don't have annotations.
  feat<-feat[!duplicated(feat[,ID_column]),]
  
  #Reformatting feature data to match the PSet requirements:
  feat<-feat[match(rownames(assay), feat[,ID_column]),]
  rownames(feat)<-feat[,ID_column]
  feat$Symbol<-feat$gene_name
  return(feat)}

######################################################### ph_data_builder #################################################
#Function to create the pheno data from cell-line object (given it includes the meta data)

ph_data_builder<- function(annotation,assay){
  
  phen<-as.data.frame(colnames(assay))
  phen$temp<-sub("U","",phen[,1])
  phen$temp<-substring(phen$temp, 1, 4)
  annotation$temp<-substring(rownames(annotation), 2, 5)
  phen<-merge(phen, annotation, by="temp" , all.x=TRUE)
  
  phen$batchid <- NA
  phen$cellid<-phen$`colnames(assay)`
  phen$Replicate <- substring(phen$cellid,6)
  phen$Replicate[phen$Replicate ==""]<- NA
  rownames(phen)<-phen$`colnames(assay)`
  
  phen<-phen[,c(15,18,16,17, 9:13)]
  return(phen)}

######################################################### eSetToSE ######################################################
# A function converting ExpressionSet to SummarizedExperiment

eSetToSE <- function(eSet , annot_name) {
  
  BiocGenerics::annotation(eSet) <- annot_name
  stopifnot(all(rownames(fData(eSet)) == rownames(exprs(eSet))))
  stopifnot(all(rownames(pData(eSet)) == colnames(exprs(eSet))))
  
  # Build summarized experiment from eSet
  SE <- SummarizedExperiment::SummarizedExperiment(
    assays=SimpleList(as.list(Biobase::assayData(eSet))),
    # Switch rearrange columns so that IDs are first, probes second
    rowData=S4Vectors::DataFrame(Biobase::fData(eSet),
                                 rownames=rownames(Biobase::fData(eSet))),
    colData=S4Vectors::DataFrame(Biobase::pData(eSet),
                                 rownames=rownames(Biobase::pData(eSet))),
    metadata=list("experimentData" = eSet@experimentData, 
                  "annotation" = Biobase::annotation(eSet), 
                  "protocolData" = Biobase::protocolData(eSet)))
  # Extract names from expression set                  
  SummarizedExperiment::assayNames(SE) <- Biobase::assayDataElementNames(eSet)
  
  stopifnot(all(rownames(colData(SE)) == rownames(pData(eSet))))
  stopifnot(all(rownames(rowData(SE)) == rownames(fData(eSet))))
  return(SE)
}

########################################### GSM_maping ########################################## 
# Cell_Maps_GEO.csv file is created using data from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE152160
# The file includes GSM id of the samples mapped to their corresponding cell ids
GSM_map<-read.csv("Data/Cell_Maps_GEO.csv", stringsAsFactors = FALSE, header = FALSE)

# Creating replicate column including replicate markers (e.g. "A") 
GSM_map$Replicate <- gsub(".*cells", "", GSM_map$V2)
GSM_map$Replicate [GSM_map$Replicate =="" | GSM_map$Replicate =="human astrocytes"]<-NA
GSM_map$Replicate <- gsub(" ", "", GSM_map$Replicate)


# Creating Patient-id column consistant with cell names in the main paper
GSM_map$Patient_id <-gsub("cells.*", "", GSM_map$V2)
GSM_map$Patient_id <- gsub(" ", "", GSM_map$Patient_id)
GSM_map$Patient_id [GSM_map$Patient_id =="humanastrocytes"]<-"human_astrocytes"
#GSM_map$V1 [GSM_map$Patient_id =="human_astrocytes"]<-"human_astrocytes"

# Creating unique cell-ids based on patient ids and replicates
GSM_map$cellid<-sub("MG cells", "", GSM_map$V2)
GSM_map$cellid<-sub(" ", "_", GSM_map$cellid)

# GSM ids are only used in the expression data
colnames(GSM_map)[colnames(GSM_map)=="V1"]<- "array_id"

GSM_map<-GSM_map[, c("Patient_id","Replicate","cellid", "array_id")]

########################################### Molecular Profiles ##########################################################
#################### expression data #####################
#Creating eset from raw expression data 

#untar("Data/Raw_expression/GSE152160_RAW.tar", exdir="GSE152160_RAW")#Unpack the CEL files
#cels<-list.files("GSE152160_RAW/", pattern = "CEL")
#sapply(paste("GSE152160_RAW", cels, sep="/"), gunzip)
#install.packages("/Data/Raw_expression/Source_hglue20hsensgcdf_25.0.0.tar.gz", repos = NULL, type = "source")

####Gene-level expression
cdf <- "hglue20hsensgcdf"
cels <- list.celfiles("Data/GSE152160_RAW", full.names = TRUE)#Raw_expression Folder includes 145 CEL files
expr_cel <- justRMA(filenames = cels, verbose = TRUE, cdfname = cdf)

# Assay data expression 
assay_exp<-as.data.frame(exprs(expr_cel))
rownames(assay_exp)<-sub("_at","",rownames(assay_exp))
colnames(assay_exp)<-sub("_PA.*","",colnames(assay_exp))

# Feature data expression
feat_exp<-fdata_builder(annotation=features_gene, assay=assay_exp,ID_column="gene_id")#features_gene from Ensembel.v99.annotation.RData
feat_exp$BEST<-NA

# Pheno data expression
cell$Patient_id<-rownames(cell)
phen_exp<-merge(GSM_map , cell, by="Patient_id" , all.x=TRUE)
phen_exp<-phen_exp[, c(1:4 , 11:15)]
phen_exp$batchid <- NA
rownames(phen_exp)<-phen_exp$array_id


########## Normalizing based on negative control genes ##########
ctrl<-read.delim("Data/HK_genes_2013.txt",stringsAsFactors = FALSE, header = FALSE) # Negative control genes from https://www.tau.ac.il/~elieis/HKG/.
ctrl$gene_name<-sub(" ","",ctrl$V1)
ctrl<-merge(ctrl, feat_exp[,c("gene_name","gene_id")], by="gene_name")
ctrl_ind<-which(rownames(assay_exp) %in% sub(" ","",ctrl$gene_id))

# Identifying cells that might contain unwanted sources of variation
XX<-as.numeric(ifelse(is.na(phen_exp$Subtype),0,1)) # Cells that are not among the original 101 GBM cells
XX<- XX - mean(XX)
XX <- cbind(XX/(sqrt(sum(XX^2))))

# Center gene expressions
YY <- t(scale(assay_exp, scale=FALSE)) 

## Prepare control samples
# A table that has as many columns as the largest set of replicates for one sample. Each
# row corresponds to a set of replicates of the same sample and gives the row indices of the
# replicates in the gene expression matrix, padded with -1 entries.
# See https://www.bioconductor.org/packages/release/bioc/vignettes/RUVnormalize/inst/doc/RUVnormalize.pdf

rep_cIdx <- matrix(-1,145,2) # 145 is number of cell lines and 2 is max number of replicates a cell has
added_pat<-c()
count <-1#row number of rep_cIdx

for (i in 1:nrow(rep_cIdx)){
  pat = phen_exp$Patient_id[i]
  if (pat %in% added_pat)
    next
  GSM=phen_exp$array_id[phen_exp$Patient_id==pat]
  if(length(GSM)==2){
    rep_cIdx[count,1]=which(rownames(YY)==GSM[1])
    rep_cIdx[count,2]=which(rownames(YY)==GSM[2])
  }
  else(rep_cIdx[count,1]=which(rownames(YY)==GSM))
  
  added_pat=c(added_pat,pat)  
  count=count+1
}

## Replicate-based
# Set k to the number of samples / 4 or to the number of replicates, if the latter is smaller than the former. 
# See https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4896498/
Res <- naiveReplicateRUV(YY, ctrl_ind, rep_cIdx, k=17)

# Rebuilding assay_exp based on normalized data
assay_exp<-data.frame(t(Res$cY))

# Protocol data expression
protocol_exp<-as.data.frame(rep("Affymetrix HTA 2.0 array",ncol(assay_exp)),row.names=colnames(assay_exp))
colnames(protocol_exp)<-"Array"
protocol_exp$URL_raw_data<-"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE152160" #Link to access raw expression data
protocol_exp$Folder_raw_data<-"GSE152160_RAW.tar" #Folder including the raw data in the above link
protocol_exp$Negative_control_genes<-"https://www.tau.ac.il/~elieis/HKG/"

  #Creating ExpressionSet 
assay_exp<-assay_exp[,rownames(phen_exp)] #rearranging the colnames so it is similar to pheno data
protocol_exp<-protocol_exp[rownames(phen_exp),]#rearranging the rownames so it is similar to pheno data
exp_eSet<- ExpressionSet(assayData = as.matrix(assay_exp), phenoData = AnnotatedDataFrame(phen_exp), 
                         featureData = AnnotatedDataFrame(feat_exp),
                         protocolData=AnnotatedDataFrame(protocol_exp)) 

####Probe-level expression
raw_data<- oligo::read.celfiles(cels) # If required run this: memory.limit(size=76000) , memory.limit()
norm_data<-rma(raw_data,target="core") # Perform RMA normalization
norm_main<- getMainProbes(norm_data, level = "core")#Remove the control transcripts (Controls do not match to any genes)
norm_main_annot<- annotateEset(norm_main, pd.hta.2.0)#Annotating the eset on transcript level : "https://support.bioconductor.org/p/89308/"

# Assay data expression-probe level
assay_exp_probe<-as.data.frame(exprs(norm_main_annot))
colnames(assay_exp_probe)<-sub("_PA.*","",colnames(assay_exp_probe))

# Feature data expression-probe level
feat_exp_probe<-fData(norm_main_annot)
feat_exp_probe$BEST<-NA

# Pheno data expression-probe level is same as Pheno data expression

# Protocol data expression-probe level
protocol_exp_probe<-as.data.frame(rep("Affymetrix HTA 2.0 array",ncol(assay_exp_probe)),row.names=colnames(assay_exp_probe))
colnames(protocol_exp_probe)<-"Array"
protocol_exp_probe$URL_raw_data<-"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE152160" #Link to access raw expression data
protocol_exp_probe$annotation<-"pd.hta.2.0"

#Creating ExpressionSet 
assay_exp_probe<-assay_exp_probe[,rownames(phen_exp)]#rearranging the colnames so it is similar to pheno data
protocol_exp_probe<-protocol_exp_probe[rownames(phen_exp),]#rearranging the rownames so it is similar to pheno data
exp_eSet_probe<- ExpressionSet(assayData = as.matrix(assay_exp_probe), phenoData = AnnotatedDataFrame(phen_exp), 
                               featureData = AnnotatedDataFrame(feat_exp_probe),
                               protocolData=AnnotatedDataFrame(protocol_exp_probe)) 

#################### CNV data #####################

assay_cnv<- read.delim("Data/HGCC_DNA_copy_number_gene_level.txt", header=T, sep="\t", stringsAsFactors = FALSE)
rownames(assay_cnv)<-assay_cnv$Row
assay_cnv<-assay_cnv[,-1]

feat_cnv<-fdata_builder(annotation=features_gene, assay=assay_cnv,ID_column="gene_name")
phen_cnv<-ph_data_builder(annotation=cell,assay=assay_cnv)
phen_cnv$Replicate[phen_cnv$cellid=="U10000" |phen_cnv$cellid=="U10001" ]<-NA

#Creating ExpressionSet 
assay_cnv<-assay_cnv[,rownames(phen_cnv)]#rearranging the colnames so it is similar to pheno data
cnv_eSet<- ExpressionSet(assayData = as.matrix(assay_cnv), phenoData = AnnotatedDataFrame(phen_cnv), featureData = AnnotatedDataFrame(feat_cnv)) 

#################### Mutation data ######################
#Assay data
mutation<- read.delim("Data/HGCC_WES_mutations_variants.txt", header=T, sep="\t", stringsAsFactors = FALSE,na.strings=c("", " ", "NA","NaN"))

mutation$CELL_LINE<- paste("U", mutation $CELL_LINE, sep="")
assay_mut<-data.frame(matrix("wt", nrow = length(unique(mutation$Gene.refGene[!is.na(mutation$Gene.refGene)])), ncol = length(unique(mutation$CELL_LINE))))
rownames(assay_mut)<-unique(mutation$Gene.refGene[!is.na(mutation$Gene.refGene)])
colnames(assay_mut)<-unique(mutation$CELL_LINE)

#Labeling wild type cells as "wt"
for (i in 1:nrow(assay_mut)){
  if(i%%100 ==0 ){print(i)}
  gene<-rownames(assay_mut)[i]
  
  indS<-which(mutation$Gene.refGene == gene)
  for (ind in indS){
    cell <- mutation$CELL_LINE [ind]
    func <- mutation$ExonicFunc.refGene [ind]
    
    if (!is.na(func)){
      if(assay_mut[gene, cell]=="wt"){assay_mut[gene, cell]=func}
      else if(grepl(func, assay_mut[gene, cell], perl=TRUE)) {next} #If the mutation effect is already reported for a gene the we skip
      else{assay_mut[gene, cell] = paste(assay_mut[gene,cell],func, sep="///")
      }} }}

assay_mut<-assay_mut[,-2]

feat_mutation<-fdata_builder(annotation=features_gene, assay=assay_mut)
phen_mutation<-ph_data_builder(annotation=cell,assay=assay_mut)
phen_mutation$Patient_id[phen_mutation$cellid=="U3067"]<-"U3067MG"# Based on "GSM_map" file

#Creating ExpressionSet 
assay_mut<-assay_mut[,rownames(phen_mutation)]#rearranging the colnames so it is similar to pheno data
mutation_eSet<- ExpressionSet(assayData = as.matrix(assay_mut), phenoData = AnnotatedDataFrame(phen_mutation), featureData = AnnotatedDataFrame(feat_mutation)) 

#################### Methylation data ####################
#Assay data
assay_methyl<- read.delim("Data/HGCC_DNA_methylation.txt", header=T, sep="\t",stringsAsFactors = FALSE)
illum<-read.csv("Data/MethylationEPIC_v-1-0_B2.csv",skip=7, stringsAsFactors = FALSE) #The first 7 rows are not informative 

#Removing SNP probes and sex chromosome probes
assay_methyl<-assay_methyl[-c(grep("rs",rownames(assay_methyl))),]
assay_methyl<-assay_methyl[-c(grep("ch.X",rownames(assay_methyl))),]

#Removing cross-reactive and polymorphic probes
polymorph_probes<-read.delim("Data/EPIC_ polymorphic_probes.txt", stringsAsFactors = FALSE)
cross_probes<-read.delim("Data/cross-hybridising_CpG-targeting probes.txt", header=F , stringsAsFactors = FALSE)
assay_methyl<-assay_methyl[rownames(assay_methyl) %in% polymorph_probes$IlmnID == FALSE,]
assay_methyl<-assay_methyl[rownames(assay_methyl) %in% cross_probes$V1 == FALSE,]
colnames(assay_methyl)[colnames(assay_methyl) == "hAstro"]<-"human_astrocytes"

#Assay data gene-level
assay_methyl_gene<-read.delim("Data/methMat.txt", sep="\t", stringsAsFactors = FALSE)
colnames(assay_methyl_gene)[colnames(assay_methyl_gene) == "hAstro"]<-"human_astrocytes"

#Feature data
feat_methyl<-fdata_builder(annotation=illum, assay=assay_methyl,ID_column="IlmnID")
feat_methyl_gene<-fdata_builder(annotation=features_gene, assay=assay_methyl_gene)

#Pheno data
phen_methyl<-ph_data_builder(annotation=cell, assay=assay_methyl)
phen_methyl$Patient_id[phen_methyl$cellid=="human_astrocytes"]<-"human_astrocytes"
phen_methyl$Replicate[phen_methyl$Replicate=="_astrocytes"]<-NA

#Pheno data methylation gene level is same as Pheno data probe level

#Protocol data
protocol_methyl<-as.data.frame(rep("Infinium® MethylationEPIC BeadChip Kit",ncol(assay_methyl)),row.names=colnames(assay_methyl))
colnames(protocol_methyl)<-"Array"
protocol_methyl$Provider<-"Illumina"
protocol_methyl$URL_GpG<-"http://portal.hgcc.se/" #Link to access methylation data provided by the authors
protocol_methyl$File_GpG<-"HGCC_DNA_methylation.txt" #File including the data in the above link
protocol_methyl$URL_annotation<-"https://support.illumina.com/array/array_kits/infinium-methylationepic-beadchip-kit/downloads.html"#Link to access annotations by illumina
protocol_methyl$File_annotation<-"Infinium MethylationEPIC v1.0 B5 Manifest File (CSV Format)"#File including annotations in the above link

protocol_methyl_gene<-as.data.frame(rep("Infinium® MethylationEPIC BeadChip Kit",ncol(assay_methyl_gene)),row.names=colnames(assay_methyl_gene))
colnames(protocol_methyl_gene)<-"Array"
protocol_methyl_gene$Provider<-"Illumina"

#Creating ExpressionSet 
assay_methyl<-assay_methyl[,rownames(phen_methyl)]#rearranging the colnames so it is similar to pheno data
protocol_methyl<-protocol_methyl[rownames(phen_methyl),]#rearranging the rownames so it is similar to pheno data
methyl_eSet<- ExpressionSet(assayData = as.matrix(assay_methyl), phenoData = AnnotatedDataFrame(phen_methyl), 
                            featureData = AnnotatedDataFrame(feat_methyl), protocolData=AnnotatedDataFrame(protocol_methyl)) 

assay_methyl_gene<-assay_methyl_gene[,rownames(phen_methyl)]#rearranging the colnames so it is similar to pheno data
protocol_methyl_gene<-protocol_methyl_gene[rownames(phen_methyl),]#rearranging the rownames so it is similar to pheno data
methyl_gene_eSet<- ExpressionSet(assayData = as.matrix(assay_methyl_gene), 
                                 phenoData = AnnotatedDataFrame(phen_methyl), 
                                 featureData = AnnotatedDataFrame(feat_methyl_gene),
                                 protocolData=AnnotatedDataFrame(protocol_methyl_gene)) 

#################### SE objects ####################

#Checks are included in the eSetToSE function
expression_SE<- eSetToSE(exp_eSet,annot_name="rna")
expression_probe_SE<- eSetToSE(exp_eSet_probe,annot_name="rna_probe")
cnv_SE <- eSetToSE(cnv_eSet,annot_name="cnv")
mutation_SE <- eSetToSE(mutation_eSet,annot_name="mut")
methyl_gene_SE <- eSetToSE(methyl_gene_eSet,annot_name="methyl_gene")
methyl_SE <- eSetToSE(methyl_eSet,annot_name="methyl_probe")

########################################### sensitivity data ###########################################
#Published AUC info
drug_cell<-read.delim("Data/HGCC_drug_response_AUC.txt", header=T, sep="\t", stringsAsFactors = FALSE)

#It is required to match the drugnames from published AUCs with the drugnames from sensitivity raw data 
#As we need to compare the published-AUCs with the computed-AUCs
#Synonym correction:
drug_cell["Azacytidine-5","X"] <- "5-azacytidine"
drug_cell["Azaguanine-8","X"] <- "8-Azaguanine"
drug_cell["Azelastine HCl","X"] <- "azelastine"
drug_cell["Balsalazide Sodium","X"] <- "balsalazide"
drug_cell["Bepridil hydrochloride ","X"] <- "bepridil"
drug_cell["Bortezomib  ","X"] <- "bortezomib [velcade ]"
drug_cell["Bucladesine sodium salt","X"] <- "bucladesine"
drug_cell["Cabazitaxel  ","X"] <- "cabazitaxel"
drug_cell["Carboplatin  ","X"] <- "carboplatin"
drug_cell["Cefazolin sodium salt","X"] <- "cefazolin"
drug_cell["Ciclopirox ethanolamine","X"] <- "ciclopirox"
drug_cell["Clofarabine  ","X"] <- "clofarabine"
drug_cell["Clomiphene citrate (Z.E) ","X"] <- "clomiphene citrate (z.e)"
drug_cell["Crizotinib  ","X"] <- "crizotinib"
drug_cell["Dichlorbenzamil","X"] <- "dichlorbenzamide" #to be duble checked
drug_cell["Doxorubicin","X"] <- "doxorubicin hydrochloride"
drug_cell["Everolimus  ","X"] <- "everolimus"
drug_cell["Lomerizine 2HCl","X"] <- "Lomerizine dihydrochloride"
drug_cell["Mometasone furoate","X"] <- "mometasone"
drug_cell["Nilotinib  ","X"] <- "nilotinib"
drug_cell["Pentamidine isethionate","X"] <- "pentamidine"
drug_cell["Pitavastatin Ca","X"] <- "pitavastatin"
drug_cell["Propantheline bromide","X"] <- "propantheline"
drug_cell["Sunitinib Malate  ","X"] <- "sunitinib malate"
drug_cell["Tenovin-6 derivat A","X"] <- "tenovin-6 derivative 39" #to be duble checked
drug_cell["Tenovin-6 derivat B","X"] <- "tenovin-6 derivative 50" #to be duble checked
drug_cell["Terreic acid","X"] <- "terreic-acid-(-)"
drug_cell["Tetrahydrozoline hydrochloride","X"] <- "tetrahydrozoline"
drug_cell["Theophylline monohydrate","X"] <- "theophylline"
drug_cell["Vemurafenib  ","X"] <- "vemurafenib"
drug_cell["Vinblastine Sulfate  ","X"] <- "vinblastine sulfate"
drug_cell["Ziprasidone  Hydrochloride","X"] <- "ziprasidone hydrochloride"
drug_cell["(R)-Duloxetine hydrochloride","X"] <- "Duloxetine"
drug_cell["(S)-propranolol hydrochloride","X"] <- "propranolol-(S)"
drug_cell["5-Nonyloxytryptamine oxalate","X"] <- "5-nonyloxytryptamine"
drug_cell["AMA404","X"] <- "AM404"
drug_cell["CGS 12066B dimaleate","X"] <- "cgs 12066b"
drug_cell["Dasatinib  ","X"] <- "dasatinib"
drug_cell["Vinorelbinetartrate","X"] <- "vinorelbine tartrate"
drug_cell["Mitoxantrone dihydrochloride","X"] <- "lomerizine" #based on molecular formula (C22H30Cl2N4O6)
drug_cell["Vakuinol-1","X"] <- "GLN-1001"#based on molecular formula (ClC1=CC=C(C2=NC(C=CC=C3)=C3C(C(O)C4NCCCC4)=C2)C=C1), The correct dictation is "Vacquinol-1"
drug_cell["Mitoxantrone dihydrochloride","X"] <-"mitoxantrone"
drug_cell["prm-116","X"] <- "tv-001"
drug_cell["prm-122","X"] <- "tv-002"
drug_cell["PRM-123","X"] <- "tv-003"

#Transposing from wide to long
drug_cell_long<- melt(setDT(drug_cell), id.vars = "X", variable.name = "Pat")
drug_cell_long$EXP_details<-paste(drug_cell_long$X, drug_cell_long$Pat , sep ="_")

######################## Screen2 ########################
screen2<-screen2[rowSums(is.na(screen2[,3:ncol(screen2)])) <ncol(screen2)-2, ] #Keeping the rows that are not all NaNs.(2 represents the two first columns which are drugs and pat names)
screen2$Drug<-sub("  ","", screen2$Drug)

#Multiplying the viabilites by 100
screen2[,3:ncol(screen2)]<-apply(screen2[,3:ncol(screen2)],2,function(x){x*100})

screen2$EXP_details<-paste(screen2$Drug, screen2$Pat , sep ="_")
colnames(screen2) <- c("drugid","cellid", colnames(screen2)[3:length(colnames(screen2))])
rownames(screen2)<-screen2$EXP_details

sen_info_scr2<-screen2[,c(1,2,14)]
sen_info_scr2$Duration<-"3 Days"

###sensitivity_raw
dose_scr2<-as.data.frame(matrix(nrow=nrow(screen2)))
i=1
for (col in colnames(screen2)[3:13]){
  micro_mol<-as.numeric(substr(col, start = 6, stop = nchar(col)))
  name<-paste("doses",i,sep="") #doses need to be named as (doses1, doses2,..)
  dose_scr2[,name]<-micro_mol
  i=i+1}

dose_scr2<-dose_scr2[-1]
rownames(dose_scr2)<-rownames(screen2)
viability_scr2<-screen2[,3:13] #Only keeping the doses 
colnames(viability_scr2)<-colnames(dose_scr2) 
sen_raw_scr2<-abind(dose_scr2, viability_scr2, along = 3) #Creating a 3d array
dimnames(sen_raw_scr2)[[3]] <- c("Dose","Viability")

###Sensitivity_profile
AAC_IC50_scr2<-PharmacoGx:::.calculateFromRaw(sen_raw_scr2, cap = NA, nthread = 1, family = c("normal","Cauchy"), scale = 0.07, n=1 )
sen_profile_scr2<-cbind(as.data.frame(AAC_IC50_scr2["AUC"]), as.data.frame(AAC_IC50_scr2["IC50"]))
colnames(sen_profile_scr2)<-c("Computed_AAC","Computed_IC50")

sen_profile_scr2$Published_AUC<-drug_cell_long$value[match(rownames(sen_profile_scr2),drug_cell_long$EXP_details)]
sen_profile_scr2$Published_AUC<-sen_profile_scr2$Published_AUC*100

######################## Screen3 ########################
screen3<-screen3[rowSums(is.na(screen3[,3:ncol(screen3)])) <ncol(screen3)-2, ] #Keeping the rows that are not all NaNs.(2 represents the two first columns which are drugs and pat names)
screen3$Drug<-sub("  ","", screen3$Drug)
screen3$Pat[screen3$Pat == "hAst"] <-"human_astrocytes"

#Multiplying the viabilites by 100
screen3[,3:ncol(screen3)]<-apply(screen3[,3:ncol(screen3)],2,function(x){x*100})

screen3$EXP_details<-paste(screen3$Drug, screen3$Pat , sep ="_")
colnames(screen3) <- c("drugid","cellid", colnames(screen3)[3:length(colnames(screen3))])
rownames(screen3)<-screen3$EXP_details

sen_info_scr3<-screen3[,c(1,2,14)]
sen_info_scr3$Duration<-"3 Days"

###sensitivity_raw
dose_scr3<-as.data.frame(matrix(nrow=nrow(screen3)))
i=1
for (col in colnames(screen3)[3:13]){
  micro_mol<-as.numeric(substr(col, start = 6, stop = nchar(col)))
  name<-paste("doses",i,sep="") #dose_scr3s need to be named as (dose_scr3s1, dose_scr3s2,..)
  dose_scr3[,name]<-micro_mol 
  i=i+1}

dose_scr3<-dose_scr3[-1]
rownames(dose_scr3)<-rownames(screen3)
viability_scr3<-screen3[,3:13] #Only keeping the dose_scr3s 
colnames(viability_scr3)<-colnames(dose_scr3) 
sen_raw_scr3<-abind(dose_scr3, viability_scr3, along = 3) #Creating a 3d array
dimnames(sen_raw_scr3)[[3]] <- c("Dose","Viability")

###Sensitivity_profile

AAC_IC50_scr3<-PharmacoGx:::.calculateFromRaw(sen_raw_scr3, cap = NA, nthread = 1, family = c("normal","Cauchy"), scale = 0.07, n=1 )
sen_profile_scr3<-cbind(as.data.frame(AAC_IC50_scr3["AUC"]), as.data.frame(AAC_IC50_scr3["IC50"]))
colnames(sen_profile_scr3)<-c("Computed_AAC","Computed_IC50")

sen_profile_scr3$Published_AUC<-drug_cell_long$value[match(rownames(sen_profile_scr3),drug_cell_long$EXP_details)]
sen_profile_scr3$Published_AUC<-sen_profile_scr3$Published_AUC*100

#Removing the published AUC values for the cell-drug pairs mutual between screen2 and screen3
#The published AUCs have a higher correlation with calculated AUCs from screen2
#Since the original paper has not specified what screen do the published AUCs belong to
#It is more conservative to report the published AUCs for screen2 

Mutual_pairs<-merge(sen_profile_scr2,sen_profile_scr3, by="row.names")
#Correlation between Published_AUC and Computed_AAC from screen2
cor(Mutual_pairs$Computed_AAC.x , Mutual_pairs$Published_AUC.x , use= "complete.obs")#-0.9051596
#Correlation between Published_AUC and Computed_AAC from screen3
cor(Mutual_pairs$Computed_AAC.y , Mutual_pairs$Published_AUC.y , use= "complete.obs")#-0.7312538

sen_profile_scr3$Published_AUC[rownames(sen_profile_scr3) %in% Mutual_pairs$Row.names]<-NA

########################################### Cell object #######################################################
#Gathering cell lines from all the experiments
phen_exp<-phen_exp[,colnames(phen_cnv)] # Reordeing the columns of phen_exp for rbind
cell_obj<-as.data.frame(unique(rbindlist( list(phen_exp,phen_cnv,phen_mutation,phen_methyl))))
rownames(cell_obj)<-cell_obj$cellid

######Screen2
#Adding cell-lines from screen2 obj
scr2_cell_obj<-cell_obj
diff_scr2<-unique(screen2$cellid[screen2$cellid %in% scr2_cell_obj$cellid ==FALSE])

#For i in diff_scr2, adds a new row 
for(i in 1:length(diff_scr2)){
  cid=substring(diff_scr2[i],1,5) 
  if (cid %in% scr2_cell_obj$cellid){
    scr2_cell_obj[nrow(scr2_cell_obj)+1,]<- scr2_cell_obj[scr2_cell_obj$cellid==cid,][1,]
    scr2_cell_obj[nrow(scr2_cell_obj),2]<- substring(diff_scr2[i],6) #Replicate name
    scr2_cell_obj[nrow(scr2_cell_obj),4]<- diff_scr2[i] #cell id
  }
  else{scr2_cell_obj[nrow(scr2_cell_obj)+1,]<- c(NA, NA,NA, diff_scr2[i], rep(NA,length(colnames(scr2_cell_obj))-4))}
}

rownames(scr2_cell_obj)<-scr2_cell_obj$cellid

######Screen3
#Adding cell-lines from screen3 obj
scr3_cell_obj<-cell_obj
diff_scr3<-unique(screen3$cellid[screen3$cellid %in% scr3_cell_obj$cellid ==FALSE])

#For i in diff_scr3, adds a new row 
for(i in 1:length(diff_scr3)){
  cid=substring(diff_scr3[i],1,5) 
  if (cid %in% scr3_cell_obj$cellid){
    scr3_cell_obj[nrow(scr3_cell_obj)+1,]<- scr3_cell_obj[scr3_cell_obj$cellid==cid,][1,]
    scr3_cell_obj[nrow(scr3_cell_obj),2]<- sub("_","",substring(diff_scr3[i],6)) #Replicate name
    scr3_cell_obj[nrow(scr3_cell_obj),4]<- diff_scr3[i] #cell id
  }
  else{scr3_cell_obj[nrow(scr3_cell_obj)+1,]<- c(NA, NA,NA, diff_scr3[i], rep(NA,length(colnames(scr3_cell_obj))-4))}
}

rownames(scr3_cell_obj)<-scr3_cell_obj$cellid

########################################### Drug object #######################################################

drugs<- read.xlsx("Data/mmc3-drugs.xlsx",rowNames = TRUE , startRow = 2)

drugs[duplicated(drugs[ , "Compound.name"]),] #Checking duplications in drug names: 3 duplications found

#Replacing duplicated drug names in "Compound.name" column with the correct drug names:
drugs$Compound.name[drugs$Molecular.Formula=="C26H27ClN2O"]<-"LOFEPRAMINE"
drugs$Compound.name[drugs$`Chemical.Abstracts.Service.(CAS).code`=="101477-54-7"]<-"Lomerizine 2hcl"
drugs$Compound.name [drugs$ Compound.name == "lomerizine"] <-"Mitoxantrone dihydrochloride" 

#To make sure that the drug-object contains ALL the drug names used in screen2 and screen3 
#the drug-names in drug-object are mapped to drug-names from screen2 and screen3.

drugs$Compound.name[drugs$Compound.name == "bortezomib [velcade ]"]<-"Bortezomib  "
drugs$Compound.name[drugs$Compound.name == "5-azacytidine"]<-"Azacytidine-5"
drugs$Compound.name[drugs$Compound.name == "TV-001"]<-"prm-116"
drugs$Compound.name[drugs$Compound.name == "TV-002"]<-"prm-122"
drugs$Compound.name[drugs$Compound.name == "TV-003"]<-"PRM-123"
drugs$Compound.name[drugs$Compound.name == "tenovin-6 derivative 39"] <- "Tenovin-6 derivat A"#To be double checked
drugs$Compound.name[drugs$Compound.name == "tenovin-6 derivative 50"] <- "Tenovin-6 derivat B"#To be double checked 
drugs$Compound.name[drugs$Compound.name == "GLN-1001"]<-"Vakuinol-1" # Based on molecular formula (ClC1=CC=C(C2=NC(C=CC=C3)=C3C(C(O)C4NCCCC4)=C2)C=C1), The correct dictation is "Vacquinol-1"
drugs$Compound.name[drugs$Compound.name == "AM404"]<-"AMA404"
drugs$Compound.name[drugs$Compound.name == "dichlorbenzamide"]<-"Dichlorbenzamil"#To be double checked 

drugs<-drug_name_correction(to_be_corrected_table=drugs , to_be_corrected_col="Compound.name") 

###################### Screen2
screen2<-drug_name_correction(to_be_corrected_table= screen2, to_be_corrected_col="drugid")
scr2_drugs<-merge(screen2[!duplicated(screen2$drugid),c("drugid", "clean_drugID")],drugs, by="clean_drugID",all.x = TRUE)
rownames(scr2_drugs)<-scr2_drugs$drugid
scr2_drugs<-scr2_drugs[,4:ncol(scr2_drugs)]


###################### Screen3
screen3<-drug_name_correction(to_be_corrected_table= screen3, to_be_corrected_col="drugid")
scr3_drugs<-merge(screen3[!duplicated(screen3$drugid),c("drugid", "clean_drugID")],drugs, by="clean_drugID",all.x = TRUE)
rownames(scr3_drugs)<-scr3_drugs$drugid
scr3_drugs<-scr3_drugs[,4:ncol(scr3_drugs)]


########################################### Curation ##########################################################

########### Cell dataframe ###########
scr2_cur_cell <- data.frame(cellid=rownames(scr2_cell_obj),
                            unique.cellid=rownames(scr2_cell_obj),
                            row.names = rownames(scr2_cell_obj))

scr3_cur_cell <- data.frame(cellid=rownames(scr3_cell_obj),
                            unique.cellid=rownames(scr3_cell_obj),
                            row.names = rownames(scr3_cell_obj))


###### Updating "Drugs_with_ids" files based on GBM data ######

####Function for extracting data using Pubchem API#
#from "https://github.com/sisiranair/Data-curation/blob/master/R/GetPubChemDetails.R"#

#to retrieve Pubchem names
# library("httr")
# library("stringi")
# library("XML")
# 
# getPubChemDetails <- function(drugList){
#   require("httr")
#   require("stringi")
#   require("XML")
#   baseURL <- "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/DRUG_NAME/property/Inchikey,CanonicalSMILES/XML"
#   baseMatrix <- matrix("", nrow = length(drugList), ncol = 4)
#   
#   for(i in 1:length(drugList)){
#     if(i %% 5 == 0){
#       #Wait a second before next execution as PubChem restricts the requests at 5 requests per second
#       Sys.sleep(1)}
#     currentDrugName <- drugList[i]
#     requestURL <- gsub("DRUG_NAME",currentDrugName, baseURL)
#     #fetching drug data from PubChem using HTTP GET request
#     requestURL <- URLencode(requestURL)
#     response <- GET(requestURL)
#     responseText <- content(response,"text")
#     responseXML <- XML::xmlParse(responseText)
#     pubChemList <- XML::xmlToList(responseXML)
#     baseMatrix[i, 1] <- drugList[i]
#   
#     if(!is.null(pubChemList$Properties$CID[1])){
#       baseMatrix[i, 2] <- pubChemList$Properties$CID[1]
#     }else{
#       print(paste(drugList[i]," ","fetched null values"))
#     }
#     if(!is.null(pubChemList$Properties$CanonicalSMILES[1])){
#       baseMatrix[i, 3] <- pubChemList$Properties$CanonicalSMILES[1]
#     }else{
#       print(paste(drugList[i]," ","fetched null values"))
#     }
#     if(!is.null(pubChemList$Properties$InChIKey[1])){
#       baseMatrix[i, 4] <- pubChemList$Properties$InChIKey[1]
#     }else{
#       print(paste(drugList[i]," ","fetched null values"))
#     }
#   }
#   
#   drugDataFrame <- as.data.frame(baseMatrix)
#   colnames(drugDataFrame) <- c("Drug", "cid","smiles","inchikey")
#   return(drugDataFrame)
# }

#Adding GBM drug ids to "drugs_with_ids" file:
# drugs_with_ids <- read.csv("Data/drugs_with_ids.csv", na.strings=c("", " ", "NA"), stringsAsFactors = FALSE, row.names = 1)
# screen2<-read.delim("Data/Screen2-drugData.txt", stringsAsFactors = FALSE)
# screen3<-read.delim("Data/Screen3-drugData.txt", stringsAsFactors = FALSE)
# 
# drugs_with_ids <-drug_name_correction(to_be_corrected_table= drugs_with_ids, to_be_corrected_col="unique.drugid")
# screen2<-drug_name_correction(to_be_corrected_table= screen2, to_be_corrected_col="Drug")
# screen3<-drug_name_correction(to_be_corrected_table= screen3, to_be_corrected_col="Drug")
# 
# both_screens<-unique(rbind(screen2[,c("Drug", "clean_drugID")] , screen3[,c("Drug", "clean_drugID")]))
# both_screens<-both_screens[!duplicated(both_screens$clean_drugID),]
# 
# drugs_with_ids_GBM<-merge(both_screens,drugs_with_ids, by="clean_drugID", all = TRUE)
# colnames(drugs_with_ids_GBM)<-c("clean_drugID", "GBM.drugid", colnames(drugs_with_ids_GBM)[3:ncol(drugs_with_ids_GBM)])
# 
# #Function for replacing the GBM.ids with their synonym unique ids
# replacement<-function(data=drugs_with_ids_GBM, unique, GBM){
#   i= which(data$unique.drugid == unique)
#   j= which(data$GBM.drugid == GBM)
#   print(i)
#   print(j)
#   data$GBM.drugid[i]= GBM
#   data<-data[-j,]
#   return(data)
# }
# 
# 
# drugs_with_ids_GBM <-replacement(unique="Amodiaquine", GBM="Amodiaquin dihydrochloride dihydrate")
# drugs_with_ids_GBM <-replacement(unique="Meclizine" , GBM="Meclozine dihydrochloride")
# drugs_with_ids_GBM <-replacement(unique="Cephalotaxine" , GBM="Homoharringtonine")
# drugs_with_ids_GBM <-replacement(unique="Carboplatinum" , GBM="Carboplatin  ")
# drugs_with_ids_GBM <-replacement(unique="Tyrphostin AG 1024" , GBM="Tyrphostin AG 1295")
# drugs_with_ids_GBM <-replacement(unique="Seliciclib" , GBM="Roscovitine")
# drugs_with_ids_GBM <-replacement(unique="Z-Leu-leu-leu-al" , GBM="MG-132")
# drugs_with_ids_GBM <-replacement(unique="Withaferin A" , GBM="withaferin" )#to be double checked
# drugs_with_ids_GBM <-replacement(unique="Ribavirin" , GBM="M15")
# drugs_with_ids_GBM <-replacement(unique="Vincaleukoblastine" , GBM= "Vinblastine Sulfate  ")
# drugs_with_ids_GBM <-replacement(unique="7H-Purine-6-thiol, 2-amino-" , GBM= "Thioguanine")
# drugs_with_ids_GBM <-replacement(unique="Mithramycin", GBM= "Plicamycin")
# drugs_with_ids_GBM <-replacement(unique="1H-Purine-2,6-dione, 3,7-dihydro-1,3-dimethyl-" , GBM="Theophylline monohydrate")
# drugs_with_ids_GBM <-replacement(unique="Azacitidine" , GBM= "Azacytidine-5")
# drugs_with_ids_GBM <-replacement(unique="Mestranol" , GBM= "Ethynylestradiol 3-methyl ether")
# drugs_with_ids_GBM <-replacement(unique="Dibutyryl cAMP" , GBM= "Bucladesine sodium salt")
# drugs_with_ids_GBM <-replacement(unique="Ramifenazone" , GBM= "Isopyrin hydrochloride")
# drugs_with_ids_GBM <-replacement(unique="Vinburnine" , GBM= "Eburnamonine (-)")
# drugs_with_ids_GBM <-replacement(unique="Hydrocortisone" , GBM="Cortisol acetate")
# drugs_with_ids_GBM <-replacement(unique="Vanoxerine" , GBM="GBR 12909 dihydrochloride")
# drugs_with_ids_GBM <-replacement(unique="Taractan" , GBM="Chlorprothixene hydrochloride")
# drugs_with_ids_GBM <-replacement(unique="6-Thioguanosine" , GBM= "Thioguanosine")


#Matching GBM.ids with unique ids based on "cid"
# unmatched_uniques<-drugs_with_ids_GBM$GBM.drugid[which(!is.na(drugs_with_ids_GBM$GBM.drugid) & is.na(drugs_with_ids_GBM$unique.drugid))]
# 
# Pubchem_details<-getPubChemDetails(unmatched_uniques)
# colnames(Pubchem_details)<-c("GBM.drugid","cid","smiles","inchikey")
# 
# drugs_with_ids_GBM<-merge(drugs_with_ids_GBM,Pubchem_details, by=c("GBM.drugid","cid","smiles","inchikey"), all = TRUE)
# drugs_with_ids_GBM<-drugs_with_ids_GBM[!duplicated(drugs_with_ids_GBM[, c("GBM.drugid" , "unique.drugid")]),]
# drugs_with_ids_GBM$unique.drugid[which(!is.na(drugs_with_ids_GBM$GBM.drugid) & is.na(drugs_with_ids_GBM$unique.drugid))]<-drugs_with_ids_GBM$GBM.drugid[which(!is.na(drugs_with_ids_GBM$GBM.drugid) & is.na(drugs_with_ids_GBM$unique.drugid))]
# 
# #Rearranging the columns of drugs_with_ids_GBM dataframe
# drugs_with_ids_GBM<-drugs_with_ids_GBM[, c("unique.drugid","CGP.drugid","CCLE.drugid","GSK.drugid","NCI60.drugid","GRAY.drugid",
#                                            "GNE.drugid","CTRP.drugid","CMAP.drugid","CTRPv2.drugid","gCSI.drugid","UHNBreast.drugid","MERCK.drugid",
#                                            "GDSC1000.drugid","NCIALM.drugid","YALE_TNBC.drugid","GBM.drugid","smiles","inchikey","cid","GDSC2019.drugid",
#                                            "FIMM.drugid","FDA","PharmacoDB.uid","PharmacoDB.uid")]
# 
# drugs_with_ids_GBM$GBM.drugid<-sub("  ","",drugs_with_ids_GBM$GBM.drugid)
# 
# write.csv(drugs_with_ids_GBM, "Data/Generated_files/Drugs_with_ids.csv", row.names = FALSE)


########### Drug dataframe ###########
drug_with_ids<-read.delim("Data/Generated_files/Drugs_with_ids.csv" , stringsAsFactors = FALSE, sep = ",")
scr2_cur_drug<-drug_with_ids[!is.na(drug_with_ids$GBM.drugid) & drug_with_ids$GBM.drugid %in% screen2$drugid, c("unique.drugid", "GBM.drugid")]
scr3_cur_drug<-drug_with_ids[!is.na(drug_with_ids$GBM.drugid) & drug_with_ids$GBM.drugid %in% screen3$drugid, c("unique.drugid", "GBM.drugid")]

#"Vinorelbinetartrate" is differently spelled in GBM_scr2 (with space) and GBM_scr3 (without space)
#To avoid defining a unique id for a same drug it is manually added to scr3_cur_drug through the below line.
scr3_cur_drug[nrow(scr3_cur_drug) + 1 ,] <- c("Vinorelbine tartrate", "Vinorelbinetartrate")

########################################### Creating PSet ###########################################
################Screen2

GBM_scr2_PSet<- PharmacoGx::PharmacoSet("GBM_scr2_PSet",
                                        molecularProfiles <- list( "rna" = expression_SE, "rna_probe" = expression_probe_SE,
                                                                   "cnv"= cnv_SE, "mut" = mutation_SE, "methyl_probe" = methyl_SE,
                                                                   "methyl_gene" = methyl_gene_SE),
                                        cell <- scr2_cell_obj,
                                        drug <- scr2_drugs,
                                        sensitivityInfo <- sen_info_scr2,
                                        sensitivityRaw <-  sen_raw_scr2,
                                        sensitivityProfiles <- sen_profile_scr2,
                                        curationDrug=scr2_cur_drug,
                                        curationCell=scr2_cur_cell,
                                        datasetType <- "sensitivity",
                                        verify = TRUE)

################Screen3

GBM_scr3_PSet<- PharmacoGx::PharmacoSet("GBM_scr3_PSet",
                                        molecularProfiles <- list( "rna" = expression_SE, "rna_probe" = expression_probe_SE,
                                                                   "cnv"= cnv_SE, "mut" = mutation_SE,"methyl_probe" = methyl_SE, 
                                                                   "methyl_gene" = methyl_gene_SE),
                                        cell <- scr3_cell_obj,
                                        drug <- scr3_drugs,
                                        sensitivityInfo <- sen_info_scr3,
                                        sensitivityRaw <-  sen_raw_scr3,
                                        sensitivityProfiles <- sen_profile_scr3,
                                        curationDrug=scr3_cur_drug,
                                        curationCell=scr3_cur_cell,
                                        datasetType <- "sensitivity",
                                        verify = TRUE)

