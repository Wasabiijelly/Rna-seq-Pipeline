#####################
## project: Hye yoon's RNA-seq Pipeline
## 
#####################

# Set Working Env
setRepositories(ind=1:8)

# Set library
library(stringr)

############## Get Argument ##############
Args <- commandArgs(trailingOnly = T)

# Args <- paste0('time Rscript khyPipeline.R',' -Fa http://ftp.ensembl.org/pub/release-105/fasta/gallus_gallus/dna/",Spc,".GRCg6a.dna.toplevel.fa.gz',
#                ' -Gtf http://ftp.ensembl.org/pub/release-105/gtf/gallus_gallus/Gallus_gallus.GRCg6a.105.gtf.gz',
#                ' -Spc Gallus_gallus -Input /disk4/ -Base /disk4/bikhy -Trim /program/Trimmomatic/trimmomatic-0.39.jar ',
#                 '-End PE -t 8 -Adapt /program/Trimmomatic/adapters/TruSeq3-PE.fa -QC /program/FastQC/fastqc ',
#                 '-HISAT /program/HISAT2/hisat2 -Samtool /program/samtools/bin/samtools ',
#                 '-Featurecount /program/subread/bin/featureCounts')
# Args <- unlist(strsplit(Args,' '))

if(!which(Args=='-Fa')){
  print('no Fa !')
}else{
  referenceLink <- Args[which(Args=='-Fa')+1]
}

if(!which(Args=='-Gtf')){
  print('no Gtf !')
}else{
  GeneLink <- Args[which(Args=='-Gtf')+1]
}

if(!which(Args=='-Spc')){
  print('no Spc !')
}else{
  species <- Args[which(Args=='-Spc')+1]
}

if(!which(Args=='-Input')){
  print('no Input !')
}else{
  inputfileLoc <- Args[which(Args=='-Input')+1]
}

if(!which(Args=='-Base')){
  print('no Base !')
}else{
  baseLoc <- Args[which(Args=='-Base')+1]
}

if(!which(Args=='-Trim')){
  print('no Trim !')
}else{
  TrimPrgmLoc <- Args[which(Args=='-Trim')+1]
}

if(!which(Args=='-End')){
  print('no End !')
}else{
  Endtype <- Args[which(Args=='-End')+1]
}

if(!which(Args=='-t')){
  print('no t !')
}else{
  thread <- Args[which(Args=='-t')+1]
}

if(!which(Args=='-Adapt')){
  print('no Adapt !')
}else{
  AdapterLoc <- Args[which(Args=='-Adapt')+1]
}

if(!which(Args=='-QC')){
  print('no QC !')
}else{
  QCPrgmLoc <- Args[which(Args=='-QC')+1]
}

if(!which(Args=='-HISAT')){
  print('no HISAT !')
}else{
  HISATPrgmLoc <- Args[which(Args=='-HISAT')+1]
}

if(!which(Args=='-Samtool')){
  print('no Samtool !')
}else{
  SamtoolPrgmLoc <- Args[which(Args=='-Samtool')+1]
}

if(!which(Args=='-Featurecount')){
  print('no Featurecount !')
}else{
  FeatureCountPrgmLoc <- Args[which(Args=='-Featurecount')+1]
}

####################### Set Working Directory #######################
setwd(baseLoc)

####################### Download reference file & gene annotation file #######################
EnsemblLoc <- paste0(baseLoc,'/1.Reference_annotation/')
command <- paste0('mkdir ',EnsemblLoc)
#system(command)

# Get reference genome file
setwd(EnsemblLoc)
command <- paste0("wget ",referenceLink," &")
#system(command)

# Get gene annotation file
command <- paste0("wget ",GeneLink," &")
#system(command)
setwd(baseLoc)
######################## Trimmomatic #########################
fastqName <- list.files(inputfileLoc,pattern="\\.fastq\\.gz")
sampleName <- substr(fastqName[1],1,str_locate(fastqName[1],'_')-1)

TrimOutLoc <- paste0(baseLoc,'/2.Trimming_Trimmomatic/')
command <- paste0('mkdir ',TrimOutLoc)
#system(command)

inputfile1 <- paste0(inputfileLoc,fastqName[1])
inputfile2 <- paste0(inputfileLoc,fastqName[2])
outputfile1 <- paste0(TrimOutLoc,sampleName,'_forward_paired.fastq.gz')
outputfile2 <- paste0(TrimOutLoc,sampleName,'_forward_unpaired.fastq.gz')
outputfile3 <- paste0(TrimOutLoc,sampleName,'_reverse_paired.fastq.gz')
outputfile4 <- paste0(TrimOutLoc,sampleName,'_reverse_unpaired.fastq.gz')

Trimcommand <- paste0("java -jar ",TrimPrgmLoc," ",Endtype," -threads ",thread," -phred33 ", 
                      inputfile1,' ',inputfile2,' ',
       outputfile1,' ', outputfile2,' ',outputfile3,' ',outputfile4,' ',
       "ILLUMINACLIP:",AdapterLoc,":2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36")
#system(Trimcommand)

######################## Quality Control by fastQC ########################
QCOutLoc <- paste0(baseLoc,'/3.QC_fastQC/')
command <- paste0('mkdir ',QCOutLoc)
#system(command)

seqfile1 <- paste0(TrimOutLoc,list.files(TrimOutLoc,pattern = 'forward_paired'))
seqfile2 <- paste0(TrimOutLoc,list.files(TrimOutLoc,pattern = 'reverse_paired'))

QCcommand <- paste0(QCPrgmLoc,' -o ',QCOutLoc,' -t ',thread," ",seqfile1," ",seqfile2 )
#system(QCcommand)

######################## HISAT & SAM TO BAM########################
HisatOutLoc <- paste0(baseLoc,'/4.HISAT/')
command <- paste0('mkdir ',HisatOutLoc)
system(command)
BamOutLoc <- paste0(baseLoc,'/5.BAM/')
command <- paste0('mkdir ',BamOutLoc)
system(command)

referenceFile <- paste0(EnsemblLoc, list.files(EnsemblLoc, pattern = '\\.fa\\.gz'))
system(paste0('gzip -d ',referenceFile))
referenceFile <- paste0(EnsemblLoc, list.files(EnsemblLoc, pattern = '\\.fa'))

#command <- paste0(HISATPrgmLoc,'-build -p 16 ',referenceFile,' ',paste0(HisatOutLoc,'genome'))
#system(command)

genome <- paste0(HisatOutLoc,'genome')
HISATcommand <- paste0(HISATPrgmLoc," -p ",thread," -x ",genome," -1 ",seqfile1," -2 ",seqfile2,
                       " | ",SamtoolPrgmLoc," sort -@ ",thread," -o ",paste0(BamOutLoc,sampleName,".bam"))

#2> ",Hisat2_id[d],".log | ",PutYouSamtools_Pro_Dir," sort -@ ",PutYourThreadsNumber," -o ",PutYourSamtools_Sort_Dir,"/",Hisat2_id[d],".Sorted.bam &")
#system(HISATcommand)

bamFile <- paste0(BamOutLoc, list.files(BamOutLoc, pattern = '\\.bam'))
sortBamcommand <- paste0(SamtoolPrgmLoc,' sort ',bamFile,' -o ',paste0(BamOutLoc,sampleName,"_sorted.bam"))

#system(sortBamcommand)

######################## Feature Count ########################
FeatureCountOutLoc <- paste0(baseLoc,'/6.FeatureCount/')
command <- paste0('mkdir ',FeatureCountOutLoc)
system(command)

FeatureCountFile <- paste0(FeatureCountOutLoc,sampleName,'_counts.txt')
sortedBamFile <- paste0(BamOutLoc,list.files(BamOutLoc,pattern = '_sorted\\.bam'))

gtfFile <- paste0(EnsemblLoc,list.files(EnsemblLoc,pattern = '\\.gtf.\\gz'))
system(paste0('gzip -d ',gtfFile)) 
gtfFile <- paste0(EnsemblLoc,list.files(EnsemblLoc,pattern = '\\.gtf'))

FeatureCountcommand <- paste0(FeatureCountPrgmLoc,' -T ',thread,' -p -a ',gtfFile,' -t exon -g gene_id -s 0 -o ',FeatureCountFile,' ', sortedBamFile)
system(FeatureCountcommand)







