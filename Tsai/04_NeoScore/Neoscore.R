args = commandArgs(trailingOnly=TRUE)
if (length(args)<8){
    stop("Eight arguments must be provided: NetMHCpan file, NetMHCstabpan file, 
        RNA quant file, wildtype file, transcript names file, full output file, restricted output file, 
        and binding neoantigen output file")
} 

#install.packages("stringr",repos = "http://cran.us.r-project.org")
library(stringr)

################# Define input files #################
netMHC_file <- args[1]
netMHCstab_file <- args[2]
quants_file <- args[3]
wildtype_file <- args[4]
transcript_names <- args[5]
output_file <- args[6]
output_file_2 <- args[7]
output_file_3 <- args[8]

################# Read in input files #################
netMHC <-as.data.frame(read.table(netMHC_file, header=TRUE,skip=1))
netMHCstab <-as.data.frame(read.table(netMHCstab_file, header=TRUE,skip=1))
quants <-as.data.frame(read.table(quants_file, header=TRUE))
quants$Name=gsub("\\..*$","",quants$Name, fixed=FALSE)
wildtype <- as.data.frame(read.table(wildtype_file, header=FALSE))
transcripts <- as.data.frame(read.table(transcript_names, header=FALSE))

################# Format NetMHC file #################
netMHC <- unique(netMHC)
netMHC$ID=gsub("MT_","",netMHC$ID, fixed=TRUE)
netMHC$ID=gsub("_E.*","",netMHC$ID)

################# Ensure that neoantigen is not present in wildtpe #################
wt <- as.matrix(read.table(wildtype_file, header = FALSE))
wt <- unique(wt)
# Define a function that splits the wildtype peptides into every possible 9mer for comparison
combinated_letters <- function(string, n) {
    length_ <- str_length(string)
    str_sub(string, seq(1, length_ + 1 - n), seq(n, length_))
}
# Split wildtype peptides into every possible 9mer and create a matrix to compare to
wt_expanded <- matrix(NA, nrow=1)
colnames(wt_expanded) <- "peptides"
for (g in 1:length(wt)){
    temp = combinated_letters(wt[g], 9)
    temp <- as.data.frame(temp)
    colnames(temp) <- "peptides"
    wt_expanded <- rbind(wt_expanded, temp)
}
# Take out peptides that exist in wildtype library
print(length(netMHC[which(netMHC$Peptide %in% wt_expanded$peptides),1]))
if(length(netMHC[which(netMHC$Peptide %in% wt_expanded$peptides),1])!=0){
    netMHC <- netMHC[-which(netMHC$Peptide %in% wt_expanded$peptides),]
}
################# Format NetMHCstabpan file #################  
netMHCstab <- netMHCstab[match(netMHC$Peptide, netMHCstab$Peptide),]
  identical(netMHC$Peptide, netMHCstab$Peptide)
  total_MHC <- cbind(netMHC, netMHCstab)
  total_MHC <- total_MHC[,-4]

################# Separate out the binding and stability for each of the alleles #################
a=1
h=length(netMHC[1,])+4
for (j in seq(6,(length(netMHC[1,])-3),5)){
    assign(paste("total_MHC_",a, sep=""), total_MHC[, c(2,3,j,h)])
    a=a+1
    h=h+3
}
a=a-1
print(a)
# Select the minimum Kd and keep only that allele, simultaneously assign expression
total_MHC_summarized <- matrix(NA, nrow=length(total_MHC_1[,1]), ncol=6)
mhc_record <- matrix(NA, nrow=length(total_MHC_1[,1]), ncol=1)
total_MHC_summarized <- as.data.frame(total_MHC_summarized)
for (k in 1:length(total_MHC_1[,1])){
    # Find minimum allele
    if (a==1){
        total_MHC_summarized[k,] <- total_MHC_1[k,]
        mhc_record[k]=1
    }
    if (a==2){
        x= min(c(total_MHC_1$nM[k],total_MHC_2$nM[k]))
        if (total_MHC_1$nM[k]==x){total_MHC_summarized[k,] <- total_MHC_1[k,]
        mhc_record[k]=1}
        if (total_MHC_2$nM[k]==x){total_MHC_summarized[k,] <- total_MHC_2[k,]
        mhc_record[k]=2}
    }
    if (a==3){
        x= min(c(total_MHC_1$nM[k],total_MHC_2$nM[k],total_MHC_3$nM[k]))
        if (total_MHC_1$nM[k]==x){total_MHC_summarized[k,] <- total_MHC_1[k,]
        mhc_record[k]=1}
        if (total_MHC_2$nM[k]==x){total_MHC_summarized[k,] <- total_MHC_2[k,]
        mhc_record[k]=2}
        if (total_MHC_3$nM[k]==x){total_MHC_summarized[k,] <- total_MHC_3[k,]
        mhc_record[k]=3}
    }
    if (a==4){
        x= min(c(total_MHC_1$nM[k],total_MHC_2$nM[k],total_MHC_3$nM[k],total_MHC_4$nM[k]))
        if (total_MHC_1$nM[k]==x){total_MHC_summarized[k,] <- total_MHC_1[k,]
        mhc_record[k]=1}
        if (total_MHC_2$nM[k]==x){total_MHC_summarized[k,] <- total_MHC_2[k,]
        mhc_record[k]=2}
        if (total_MHC_3$nM[k]==x){total_MHC_summarized[k,] <- total_MHC_3[k,]
        mhc_record[k]=3}
        if (total_MHC_4$nM[k]==x){total_MHC_summarized[k,] <- total_MHC_4[k,]
        mhc_record[k]=4}
    }
    if (a==5){
        x= min(c(total_MHC_1$nM[k],total_MHC_2$nM[k],total_MHC_3$nM[k],total_MHC_4$nM[k],total_MHC_5$nM[k]))
        if (total_MHC_1$nM[k]==x){total_MHC_summarized[k,] <- total_MHC_1[k,]
        mhc_record[k]=1}
        if (total_MHC_2$nM[k]==x){total_MHC_summarized[k,] <- total_MHC_2[k,]
        mhc_record[k]=2}
        if (total_MHC_3$nM[k]==x){total_MHC_summarized[k,] <- total_MHC_3[k,]
        mhc_record[k]=3}
        if (total_MHC_4$nM[k]==x){total_MHC_summarized[k,] <- total_MHC_4[k,]
        mhc_record[k]=4}
        if (total_MHC_5$nM[k]==x){total_MHC_summarized[k,] <- total_MHC_5[k,]
        mhc_record[k]=5}
    }
    if (a==6){
        x= min(c(total_MHC_1$nM[k],total_MHC_2$nM[k],total_MHC_3$nM[k],total_MHC_4$nM[k],total_MHC_5$nM[k],total_MHC_6$nM[k]))
        if (total_MHC_1$nM[k]==x){total_MHC_summarized[k,] <- total_MHC_1[k,]
        mhc_record[k]=1}
        if (total_MHC_2$nM[k]==x){total_MHC_summarized[k,] <- total_MHC_2[k,]
        mhc_record[k]=2}
        if (total_MHC_3$nM[k]==x){total_MHC_summarized[k,] <- total_MHC_3[k,]
        mhc_record[k]=3}
        if (total_MHC_4$nM[k]==x){total_MHC_summarized[k,] <- total_MHC_4[k,]
        mhc_record[k]=4}
        if (total_MHC_5$nM[k]==x){total_MHC_summarized[k,] <- total_MHC_5[k,]
        mhc_record[k]=5}
        if (total_MHC_6$nM[k]==x){total_MHC_summarized[k,] <- total_MHC_6[k,]
        mhc_record[k]=6}
    }
    if (a==7){
        x= min(c(total_MHC_1$nM[k],total_MHC_2$nM[k],total_MHC_3$nM[k],total_MHC_4$nM[k],total_MHC_5$nM[k],total_MHC_6$nM[k],
            total_MHC_7$nM[k]))
        if (total_MHC_1$nM[k]==x){total_MHC_summarized[k,] <- total_MHC_1[k,]
        mhc_record[k]=1}
        if (total_MHC_2$nM[k]==x){total_MHC_summarized[k,] <- total_MHC_2[k,]
        mhc_record[k]=2}
        if (total_MHC_3$nM[k]==x){total_MHC_summarized[k,] <- total_MHC_3[k,]
        mhc_record[k]=3}
        if (total_MHC_4$nM[k]==x){total_MHC_summarized[k,] <- total_MHC_4[k,]
        mhc_record[k]=4}
        if (total_MHC_5$nM[k]==x){total_MHC_summarized[k,] <- total_MHC_5[k,]
        mhc_record[k]=5}
        if (total_MHC_6$nM[k]==x){total_MHC_summarized[k,] <- total_MHC_6[k,]
        mhc_record[k]=6}
        if (total_MHC_7$nM[k]==x){total_MHC_summarized[k,] <- total_MHC_7[k,]
        mhc_record[k]=7}
    }
    if (a==8){
        x= min(c(total_MHC_1$nM[k],total_MHC_2$nM[k],total_MHC_3$nM[k],total_MHC_4$nM[k],total_MHC_5$nM[k],total_MHC_6$nM[k],
            total_MHC_7$nM[k],total_MHC_8$nM[k]))
        if (total_MHC_1$nM[k]==x){total_MHC_summarized[k,] <- total_MHC_1[k,]
        mhc_record[k]=1}
        if (total_MHC_2$nM[k]==x){total_MHC_summarized[k,] <- total_MHC_2[k,]
        mhc_record[k]=2}
        if (total_MHC_3$nM[k]==x){total_MHC_summarized[k,] <- total_MHC_3[k,]
        mhc_record[k]=3}
        if (total_MHC_4$nM[k]==x){total_MHC_summarized[k,] <- total_MHC_4[k,]
        mhc_record[k]=4}
        if (total_MHC_5$nM[k]==x){total_MHC_summarized[k,] <- total_MHC_5[k,]
        mhc_record[k]=5}
        if (total_MHC_6$nM[k]==x){total_MHC_summarized[k,] <- total_MHC_6[k,]
        mhc_record[k]=6}
        if (total_MHC_7$nM[k]==x){total_MHC_summarized[k,] <- total_MHC_7[k,]
        mhc_record[k]=7}
        if (total_MHC_8$nM[k]==x){total_MHC_summarized[k,] <- total_MHC_8[k,]
        mhc_record[k]=8}
    }
    if (a==9){
        x= min(c(total_MHC_1$nM[k],total_MHC_2$nM[k],total_MHC_3$nM[k],total_MHC_4$nM[k],total_MHC_5$nM[k],total_MHC_6$nM[k],
            total_MHC_7$nM[k],total_MHC_8$nM[k],total_MHC_9$nM[k]))
        if (total_MHC_1$nM[k]==x){total_MHC_summarized[k,] <- total_MHC_1[k,]
        mhc_record[k]=1}
        if (total_MHC_2$nM[k]==x){total_MHC_summarized[k,] <- total_MHC_2[k,]
        mhc_record[k]=2}
        if (total_MHC_3$nM[k]==x){total_MHC_summarized[k,] <- total_MHC_3[k,]
        mhc_record[k]=3}
        if (total_MHC_4$nM[k]==x){total_MHC_summarized[k,] <- total_MHC_4[k,]
        mhc_record[k]=4}
        if (total_MHC_5$nM[k]==x){total_MHC_summarized[k,] <- total_MHC_5[k,]
        mhc_record[k]=5}
        if (total_MHC_6$nM[k]==x){total_MHC_summarized[k,] <- total_MHC_6[k,]
        mhc_record[k]=6}
        if (total_MHC_7$nM[k]==x){total_MHC_summarized[k,] <- total_MHC_7[k,]
        mhc_record[k]=7}
        if (total_MHC_8$nM[k]==x){total_MHC_summarized[k,] <- total_MHC_8[k,]
        mhc_record[k]=8}
        if (total_MHC_9$nM[k]==x){total_MHC_summarized[k,] <- total_MHC_9[k,]
        mhc_record[k]=9}
    }
    
    total_MHC_summarized <- as.data.frame(total_MHC_summarized)
    # Find expression value
    colnames(total_MHC_summarized) <- c("Peptide", "Gene", "Kd", "Stability", "Gene1", "Expression")
    transcript_name <- transcripts[grep(total_MHC_summarized$Peptide[k], as.character(transcripts$V2), value=FALSE),]
    if (length(unique(transcript_name$V2))==1){
        total_MHC_summarized$Gene1[k] <- unique(transcript_name$V2)
    }
    if (length(unique(transcript_name$V2))>1){
        total_MHC_summarized$Gene1[k] <- paste(unique(transcript_name$V2), collapse="_ ")
    }
    expression <- quants[which(quants$Name %in%  transcript_name$V1),4]
    print(quants[which(quants$Name %in% transcript_name$V1),])
    if (length(expression) !=0){
        expression <- na.omit(expression)
        total_MHC_summarized$Expression[k] <- sum(as.numeric(expression))
    }
}
total_MHC_summarized <- cbind(total_MHC_summarized, mhc_record)
colnames(total_MHC_summarized) <- c("Peptide", "Gene", "Kd", "Stability", "Peptide1", "Expression", "HLA")
#total_MHC_summarized <- total_MHC_summarized[-(which(is.na(total_MHC_summarized$Gene1))),]
immunogenicity <- (-1.7951005 + 1.2868273*log(as.numeric(total_MHC_summarized$Expression)+0.1,10)
    -1.8477746*log(as.numeric(total_MHC_summarized$Kd),10)
    +0.7698449*log(as.numeric(total_MHC_summarized$Stability)+0.01,10))
neoantigen_candidates <- cbind(total_MHC_summarized, immunogenicity)
neoantigen_candidates <- as.data.frame(neoantigen_candidates)
colnames(neoantigen_candidates) <- c("Peptide", "Gene", "Kd", "Stability", "Gene1", "Expression", "HLA", "Immunogenicity_Score")
neoantigen_candidates <- neoantigen_candidates[order(neoantigen_candidates$Immunogenicity_Score, decreasing=TRUE),]
neoantigen_candidates_restricted <- neoantigen_candidates[which(neoantigen_candidates$Immunogenicity_Score > -2.4777),]
binding_neoantigens <- neoantigen_candidates[which(neoantigen_candidates$Kd < 500),]

write.csv(neoantigen_candidates, file=output_file, quote=FALSE, row.names=FALSE)
write.csv(neoantigen_candidates_restricted, file=output_file_2, quote=FALSE, row.names=FALSE)
write.csv(binding_neoantigens, file=output_file_3, quote=FALSE, row.names=FALSE)
