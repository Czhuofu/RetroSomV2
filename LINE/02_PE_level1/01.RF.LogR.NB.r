#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
set.seed(1)
datpath<-args[1]
retro<-args[2]
strand<-args[3]
sub<-args[4]
masterpath<-args[5]
set.seed(1)

library('randomForest')
library('e1071')
library('glmnet')

### read in the RF, LogR and NB models ###
filename <- paste(as.character(masterpath),"/LINE/02_PE_level1/RFI_1.rds", sep="")
rf1 <- readRDS(filename)
filename <- paste(as.character(masterpath),"/LINE/02_PE_level1/RFI_2.rds", sep="")
rf2 <- readRDS(filename)
filename <- paste(as.character(masterpath),"/LINE/02_PE_level1/RFI_3.rds", sep="")
rf3 <- readRDS(filename)
filename <- paste(as.character(masterpath),"/LINE/02_PE_level1/RFI_4.rds", sep="")
rf4 <- readRDS(filename)
filename <- paste(as.character(masterpath),"/LINE/02_PE_level1/RFI_5.rds", sep="")
rf5 <- readRDS(filename)
filename <- paste(as.character(masterpath),"/LINE/02_PE_level1/RFI_6.rds", sep="")
rf6 <- readRDS(filename)
filename <- paste(as.character(masterpath),"/LINE/02_PE_level1/RFI_7.rds", sep="")
rf7 <- readRDS(filename)
filename <- paste(as.character(masterpath),"/LINE/02_PE_level1/RFI_8.rds", sep="")
rf8 <- readRDS(filename)
filename <- paste(as.character(masterpath),"/LINE/02_PE_level1/LogRI.rds", sep="")
logr <- readRDS(filename)
filename <- paste(as.character(masterpath),"/LINE/02_PE_level1/NBI.rds", sep="")
nb <- readRDS(filename)

### read in the test data
filename <- paste(as.character(datpath),"/",as.character(sub),"/retro_v",as.character(retro),"_",as.character(strand),"/",as.character(sub),".pe.LINE.matrix", sep="")
clone <- read.table(filename, sep="\t", header=T)
subclone <- subset(clone, (clone[,7] == 1 | clone[,8] == 1) & ref == 0)
subclone$gap <- subclone$gap/(abs(subclone$refpos1 - subclone$refpos2))

### make predictions ###
### RF models ###
    test1 <- subset(subclone, !is.na(ACAG))
    nx1 <- test1[,c("seg", "map", "anchor_align", "anchor_len", "anchor_mm", "anchor_XS", "anchor_AS", "anchor_mapQ", "depth","map_len","short","end3","end5","direction", "refpos", "dist", "c1", "upstream", "ACAG")]
    pred.RF1 <- predict(rf1,newdata=nx1,type='prob')[,2]
    pred1 <- data.frame(cbind(as.character(test1$chr), test1$cord1, test1$cord2, as.character(test1$read), pred.RF1))
    filename <- paste(as.character(datpath),"/",as.character(sub),"/retro_v",as.character(retro),"_",as.character(strand),"/LINE/",as.character(sub),".pe.pred.G1.txt", sep="")
    write(t(pred1), file=filename, ncol=5, sep="\t")

    test2 <- subset(subclone, !is.na(TAG))
    nx2 <- test2[,c("seg", "map", "anchor_align", "anchor_len", "anchor_mm", "anchor_XS", "anchor_AS", "anchor_mapQ", "depth","map_len","short","end3","end5","direction", "refpos", "dist", "c1", "upstream", "TAG")]
    pred.RF2 <- predict(rf2,newdata=nx2,type='prob')[,2]
    pred2 <- data.frame(cbind(as.character(test2$chr), test2$cord1, test2$cord2, as.character(test2$read), pred.RF2))
    filename <- paste(as.character(datpath),"/",as.character(sub),"/retro_v",as.character(retro),"_",as.character(strand),"/LINE/",as.character(sub),".pe.pred.G2.txt", sep="")
    write(t(pred2), file=filename, ncol=5, sep="\t")

    test3 <- subset(subclone, is.na(TAG) & is.na(ACAG) & (ORF==1))
    nx3 <- test3[,c("seg", "map", "anchor_align", "anchor_len", "anchor_mm", "anchor_XS", "anchor_AS", "anchor_mapQ","depth","map_len","short","end3","end5","direction", "refpos", "dist", "c1", "upstream", "Xscore")]
    pred.RF3 <- predict(rf3,newdata=nx3,type='prob')[,2]
    pred3 <- data.frame(cbind(as.character(test3$chr), test3$cord1, test3$cord2, as.character(test3$read), pred.RF3))
    filename <- paste(as.character(datpath),"/",as.character(sub),"/retro_v",as.character(retro),"_",as.character(strand),"/LINE/",as.character(sub),".pe.pred.G3.txt", sep="")
    write(t(pred3), file=filename, ncol=5, sep="\t")

    test4 <- subset(subclone, is.na(TAG) & is.na(ACAG) & (ORF==0) & (refpos < 6015))
    nx4 <- test4[,c("seg", "map", "anchor_align", "anchor_len", "anchor_mm", "anchor_XS", "anchor_AS", "anchor_mapQ", "depth","map_len","short","end3","end5","direction", "refpos", "dist", "c1", "upstream")]
    pred.RF4 <- predict(rf4,newdata=nx4,type='prob')[,2]
    pred4 <- data.frame(cbind(as.character(test4$chr), test4$cord1, test4$cord2, as.character(test4$read), pred.RF4))
    filename <- paste(as.character(datpath),"/",as.character(sub),"/retro_v",as.character(retro),"_",as.character(strand),"/LINE/",as.character(sub),".pe.pred.G4.txt", sep="")
    write(t(pred4), file=filename, ncol=5, sep="\t")

    test5 <- subset(subclone, is.na(TAG) & is.na(ACAG) & !is.na(GC5389))
    nx5 <- test5[,c("seg", "map", "anchor_align", "anchor_len", "anchor_mm", "anchor_XS", "anchor_AS", "anchor_mapQ", "depth","map_len","short","end3","end5","direction", "refpos", "dist", "c1", "upstream", "GC5389")]
    pred.RF5 <- predict(rf5,newdata=nx5,type='prob')[,2]
    pred5 <- data.frame(cbind(as.character(test5$chr), test5$cord1, test5$cord2, as.character(test5$read), pred.RF5))
    filename <- paste(as.character(datpath),"/",as.character(sub),"/retro_v",as.character(retro),"_",as.character(strand),"/LINE/",as.character(sub),".pe.pred.G5.txt", sep="")
    write(t(pred5), file=filename, ncol=5, sep="\t")

    test6 <- subset(subclone, is.na(TAG) & is.na(ACAG) & !is.na(G5533))
    nx6 <- test6[,c("seg", "map", "anchor_align", "anchor_len", "anchor_mm", "anchor_XS", "anchor_AS", "anchor_mapQ","depth","map_len","short","end3","end5","direction", "refpos", "dist", "c1", "upstream", "G5533", "C5533")]
    pred.RF6 <- predict(rf6,newdata=nx6,type='prob')[,2]
    pred6 <- data.frame(cbind(as.character(test6$chr), test6$cord1, test6$cord2, as.character(test6$read), pred.RF6))
    filename <- paste(as.character(datpath),"/",as.character(sub),"/retro_v",as.character(retro),"_",as.character(strand),"/LINE/",as.character(sub),".pe.pred.G6.txt", sep="")
    write(t(pred6), file=filename, ncol=5, sep="\t")

    test7 <- subset(subclone, is.na(TAG) & is.na(ACAG) & !is.na(AT5710))
    nx7 <- test7[,c("seg", "map", "anchor_align", "anchor_len", "anchor_mm", "anchor_XS", "anchor_AS", "anchor_mapQ", "depth","map_len","short","end3","end5","direction", "refpos", "dist", "c1", "upstream", "AT5710")]
    pred.RF7 <- predict(rf7,newdata=nx7,type='prob')[,2]
    pred7 <- data.frame(cbind(as.character(test7$chr), test7$cord1, test7$cord2, as.character(test7$read), pred.RF7))
    filename <- paste(as.character(datpath),"/",as.character(sub),"/retro_v",as.character(retro),"_",as.character(strand),"/LINE/",as.character(sub),".pe.pred.G7.txt", sep="")
    write(t(pred7), file=filename, ncol=5, sep="\t")

    test8 <- subset(subclone, anchor_split==1)
    nx8 <- test8[,c("seg", "map", "anchor_align", "anchor_len", "anchor_mm", "anchor_direct", "anchor_insert", "anchor_seg", "anchor_map", "anchor_XS", "anchor_AS", "anchor_mapQ", "depth","map_len","short","end3","end5","direction", "refpos", "dist", "c1", "upstream")]
    pred.RF8 <- predict(rf8,newdata=nx8,type='prob')[,2]
    pred8 <- data.frame(cbind(as.character(test8$chr), test8$cord1, test8$cord2, as.character(test8$read), pred.RF8))
    filename <- paste(as.character(datpath),"/",as.character(sub),"/retro_v",as.character(retro),"_",as.character(strand),"/LINE/",as.character(sub),".pe.pred.G8.txt", sep="")
    write(t(pred8), file=filename, ncol=5, sep="\t")

    test9 <-subclone
    nx9 <- test9[,c("seg", "map", "anchor_align", "anchor_len", "anchor_mm", "anchor_XS", "anchor_AS", "anchor_mapQ", "depth","map_len","short","end3","end5","direction", "refpos", "dist", "c1", "upstream")]
    pred.RF9 <- predict(nb,newdata=nx9,type="raw")[,2]
    pred9 <- data.frame(cbind(as.character(test9$chr), test9$cord1, test9$cord2, as.character(test9$read), pred.RF9))
    filename <- paste(as.character(datpath),"/",as.character(sub),"/retro_v",as.character(retro),"_",as.character(strand),"/LINE/",as.character(sub),".pe.pred.NB.txt", sep="")
    write(t(pred9), file=filename, ncol=5, sep="\t")

    pred.RF10 <- predict(logr,newdata=nx9,type="response")
    pred10 <- data.frame(cbind(as.character(test9$chr), test9$cord1, test9$cord2, as.character(test9$read), pred.RF10))
    filename <- paste(as.character(datpath),"/",as.character(sub),"/retro_v",as.character(retro),"_",as.character(strand),"/LINE/",as.character(sub),".pe.pred.logR.txt", sep="")
    write(t(pred10), file=filename, ncol=5, sep="\t")

