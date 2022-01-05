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

 # Args <- paste0('time Rscript khyPipeline.R',' --Fa http://ftp.ensembl.org/pub/release-105/fasta/gallus_gallus/dna/gallus_gallus.GRCg6a.dna.toplevel.fa.gz',
 #                ' --Gtf http://ftp.ensembl.org/pub/release-105/gtf/gallus_gallus/Gallus_gallus.GRCg6a.105.gtf.gz',
 #                ' --Spc Gallus_gallus -Input /disk4/bikhy/oneSampleTest -Base /disk4/bikhy -Trim /program/Trimmomatic/trimmomatic-0.39.jar ',
 #                 '-End PE --t 8 -Adapt /program/Trimmomatic/adapters/TruSeq3-PE.fa -QC /program/FastQC/fastqc ',
 #                 '-HISAT /program/HISAT2/hisat2 -Samtool /program/samtools/bin/samtools ',
 #                 '-Featurecount /program/subread/bin/featureCounts')
 #Args <- unlist(strsplit(Args,' '))
 
Args <- 'time Rscript khyPipeline.R --Fa 1.Reference_annotation --Gtf 1.Reference_annotation -Input /disk4/bikhy/oneSampleTest -Base /disk4/bikhy -Trim /program/Trimmomatic/trimmomatic-0.39.jar -End PE --t 8 -Adapt /program/Trimmomatic/adapters/TruSeq3-PE.fa -QC /program/FastQC/fastqc -HISAT /program/HISAT2/hisat2 -Samtool /program/samtools/bin/samtools -Featurecount /program/subread/bin/featureCounts'
Args <- unlist(strsplit(Args,' '))

if(!'--Fa' %in% Args){
  print('no Fa !')
}else{
  FaFileLoc <- Args[which(Args=='--Fa')+1]
}

if(!'--Gtf' %in% Args){
  print('no Gtf !')
}else{
  GtfFileLoc <- Args[which(Args=='--Gtf')+1]
}

if('--EblFa' %in% Args){
  referenceLink <- Args[which(Args=='--EblFa')+1]
}else{
  referenceLink <- NULL
}

if('--EblGtf' %in% Args){
  GeneLink <- Args[which(Args=='--EblGtf')+1]
}else{
  GeneLink <- NULL
}

if('--Spc' %in% Args){
  species <- Args[which(Args=='--Spc')+1]
}else{
  species <- NULL
}

if(!'-Input' %in% Args){
  print('no Input !')
}else{
  inputfileLoc <- Args[which(Args=='-Input')+1]
}

if(!'-Base' %in% Args){
  print('no Base !')
}else{
  baseLoc <- Args[which(Args=='-Base')+1]
}

if(!'-Trim' %in% Args){
  print('no Trim !')
}else{
  TrimPrgmLoc <- Args[which(Args=='-Trim')+1]
}

if(!'-End' %in% Args){
  print('no End !')
}else{
  Endtype <- Args[which(Args=='-End')+1]
}

if(!'--t' %in% Args){
  thread=NULL
  CPUCMD <- paste0("top -n 1 -b| grep -i cpu\\(s\\)| awk \'{print $8}\'")
  Current_CPU_remain <- system(CPUCMD,intern=T)
  Current_CPU_remain_n <- as.numeric(Current_CPU_remain)
}else{
  thread <- Args[which(Args=='--t')+1]
}

if(!'-Adapt' %in% Args){
  print('no Adapt !')
}else{
  AdapterLoc <- Args[which(Args=='-Adapt')+1]
}

if(!'-QC' %in% Args){
  print('no QC !')
}else{
  QCPrgmLoc <- Args[which(Args=='-QC')+1]
}

if(!'-HISAT' %in% Args){
  print('no HISAT !')
}else{
  HISATPrgmLoc <- Args[which(Args=='-HISAT')+1]
}

if(!'-Samtool' %in% Args){
  print('no Samtool !')
}else{
  SamtoolPrgmLoc <- Args[which(Args=='-Samtool')+1]
}

if(!'-Featurecount' %in% Args){
  print('no Featurecount !')
}else{
  FeatureCountPrgmLoc <- Args[which(Args=='-Featurecount')+1]
}

####################### Set Working Directory #######################
setwd(baseLoc)

####################### Define Fuction #######################
fastqName <- list.files(inputfileLoc,pattern="\\.fastq\\.gz")
sampleNum <- length(fastqName)%/%2

setThread <- function(thread){
  if(!is.null(thread)){
    return(thread)
  }
  
  if(Current_CPU_remain_n < 50){
    thread <- 16%/%sampleNum
    print("Threads number is changed : 16")
  }else{
    thread <- 32%/%sampleNum
    print("Threads number is unchanged : 32")
  }
  
  return(thread)
}

####################### Parallel #######################
#MAXthread <- detectCores() # 쓰레드 수 확인
MAXthread <- sampleNum*2 # 12
library(doParallel) # 백엔드 생성
cl=makeCluster(MAXthread %/% as.numeric(thread))
registerDoParallel(cl)

library(foreach)

####################### 1.Download reference file & gene annotation file #######################
RefAnnLoc <- paste0(baseLoc,'/1.Reference_annotation/')

if(!is.null(referenceLink) | !is.null(GeneLink)){
  command <- paste0('mkdir ',RefAnnLoc)
  system(command)
 }
 
if(!is.null(referenceLink)){
  setwd(RefAnnLoc)
  command <- paste0("wget ",referenceLink," &")
  system(command)
    
  referenceFile <- paste0(RefAnnLoc, list.files(RefAnnLoc, pattern = '\\.fa\\.gz'))
  system(paste0('gzip -d ',referenceFile))
    
  FaFileLoc <- paste0(RefAnnLoc, list.files(RefAnnLoc,pattern = '\\.fa'))
  setwd(baseLoc)
}
  
if(!is.null(GeneLink)){
  setwd(RefAnnLoc)
  command <- paste0("wget ",GeneLink," &")
  system(command)
    
  gtfFile <- paste0(RefAnnLoc,list.files(RefAnnLoc,pattern = '\\.gtf.\\gz'))
  system(paste0('gzip -d ',gtfFile))
    
  GtfFileLoc <- paste0(RefAnnLoc, list.files(RefAnnLoc,pattern = '\\.gtf'))
  setwd(baseLoc)
}

setwd(baseLoc)
######################## 2.Trimmomatic #########################
TrimOutLoc <- paste0(baseLoc,'/2.Trimming_Trimmomatic/')

sampleName <- c()
for(i in 1:sampleNum){
  temp <- substr(fastqName[2*i],1,str_locate(fastqName[2*i],'_')-1)
  sampleName <- append(sampleName,temp)
}


if(dir.exists(TrimOutLoc) && (length(list.files(TrimOutLoc))!=0)){   #dir 이름 & 내용물 있음
  print("!!Trimming exist!!")
  }else{
    
    if(!dir.exists(TrimOutLoc)){
      command <- paste0('mkdir ',TrimOutLoc)
      system(command)
    }
    
    print("-------------START TRIMMING-------------")
    
    thread <- setThread(thread)
    
    Trimcommand <- c()
    for(i in 1:sampleNum){
      inputfile1 <- paste0(inputfileLoc,'/',fastqName[2*i-1])
      inputfile2 <- paste0(inputfileLoc,'/',fastqName[2*i])
      
      outputfile1 <- paste0(TrimOutLoc,sampleName[i],'_forward_paired.fastq.gz')
      outputfile2 <- paste0(TrimOutLoc,sampleName[i],'_forward_unpaired.fastq.gz')
      outputfile3 <- paste0(TrimOutLoc,sampleName[i],'_reverse_paired.fastq.gz')
      outputfile4 <- paste0(TrimOutLoc,sampleName[i],'_reverse_unpaired.fastq.gz')
      
      temp <- paste0("java -jar ",TrimPrgmLoc," ",Endtype," -threads ",thread," -phred33 ", 
                     inputfile1,' ',inputfile2,' ',
                     outputfile1,' ', outputfile2,' ',outputfile3,' ',outputfile4,' ',
                     "ILLUMINACLIP:",AdapterLoc,":2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36")
      Trimcommand <- append(Trimcommand, temp)
    }
    
    foreach(i = 1:sampleNum) %dopar% {
      system(Trimcommand[i])
    }   
  }

######################## 3.Quality Control by fastQC ########################
QCOutLoc <- paste0(baseLoc,'/3.QC_fastQC/')

seqfile1 <- paste0(TrimOutLoc,list.files(TrimOutLoc,pattern = 'forward_paired'))
seqfile2 <- paste0(TrimOutLoc,list.files(TrimOutLoc,pattern = 'reverse_paired'))

if(dir.exists(QCOutLoc)&& (length(list.files(QCOutLoc))!=0)){
  print("!!fastQC exist!!")
}else{
  if(!dir.exists(QCOutLoc)){
    command <- paste0('mkdir ',QCOutLoc)
    system(command)
  }
  print("-------------START FASTQC-------------")

  thread <- setThread(thread)

  QCcommand <- c()
  for(i in 1:sampleNum){
    temp <- paste0(QCPrgmLoc,' -o ',QCOutLoc,' -t ',thread," ",seqfile1[i]," ",seqfile2[i] )
    QCcommand <- append(QCcommand, temp)
  }
  foreach(i = 1:sampleNum) %dopar% {
    system(QCcommand[i])
  }
}

######################## 4.HISAT & 5.SAM TO BAM########################
referenceFile <- FaFileLoc

HisatOutLoc <- paste0(baseLoc,'/4.HISAT/')

if(dir.exists(HisatOutLoc) && (length(list.files(HisatOutLoc))!=0)){
  print("!!HISAT exist!!")
}else{
  if(!dir.exists(HisatOutLoc)){
    command <- paste0('mkdir ',HisatOutLoc)
    system(command)
  }
  print("-------------START HISAT-------------")

  thread <- setThread(thread)

  Refcommand <- paste0(HISATPrgmLoc,'-build -p 16 ',referenceFile,' ',paste0(HisatOutLoc,'genome'))
  system(Refcommand)
}

BamOutLoc <- paste0(baseLoc,'/5.BAM/')

if(dir.exists(BamOutLoc)&& (length(list.files(BamOutLoc))!=0)){
  print("!!BAM exist!!")
}else{
  if(!dir.exists(BamOutLoc)){
    command <- paste0('mkdir ',BamOutLoc)
    system(command)
  }
  print("-------------START BAM-------------")

  genome <- paste0(HisatOutLoc,'genome')

  HISATcommand <- c()
  for(i in 1:sampleNum){
    temp <- paste0(HISATPrgmLoc," -p ",thread," -x ",genome," -1 ",seqfile1[i]," -2 ",seqfile2[i],
                   " | ",SamtoolPrgmLoc," sort -@ ",thread," -o ",paste0(BamOutLoc,sampleName[i],".bam"))

    HISATcommand <- append(HISATcommand,temp)
  }

  foreach(i= 1:sampleNum) %dopar% {
    system(HISATcommand[i])
  }

  print("-------------START SORTED BAM-------------")

  bamFile <- paste0(BamOutLoc, list.files(BamOutLoc, pattern = '\\.bam'))
  
  sortBamcommand <- c()
  for(i in 1:sampleNum){
    temp <- paste0(SamtoolPrgmLoc,' sort ',bamFile[i],' -o ',paste0(BamOutLoc,sampleName[i],"_sorted.bam"))
    sortBamcommand <- append(sortBamcommand, temp)
  }

  foreach(i= 1:sampleNum) %dopar% {
    system(sortBamcommand[i])
  }

}

######################## 6.Feature Count ########################
thread <- setThread(thread)
FeatureCountOutLoc <- paste0(baseLoc,'/6.FeatureCount/')

if(dir.exists(FeatureCountOutLoc) && (length(list.files(FeatureCountOutLoc))!=0)){
  print("!!FeatureCount exist!!")
}else{
  if(!dir.exists(FeatureCountOutLoc)){
    command <- paste0('mkdir ',FeatureCountOutLoc)
    system(command)
  }

  print("-------------START FEARUTECOUNT-------------")

  FeatureCountFile <- paste0(FeatureCountOutLoc,sampleName,'_counts.txt')
  sortedBamFile <- paste0(BamOutLoc,list.files(BamOutLoc,pattern = '_sorted\\.bam'))

  gtfFile <- GtfFileLoc

  FeatureCountcommand <- c()

  for(i in 1:sampleNum){
    temp <- paste0(FeatureCountPrgmLoc,' -T ',thread,' -p -a ',gtfFile,' -t exon -g gene_id -s 0 -o ',FeatureCountFile[i],' ', sortedBamFile[i])
    FeatureCountcommand <- append(FeatureCountcommand,temp)
  }

  foreach(i= 1:sampleNum) %dopar% {
    system(FeatureCountcommand[i])
  }

}









