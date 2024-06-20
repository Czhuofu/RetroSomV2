#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)
set.seed(1)
datpath<-args[1]
retro<-args[2]
strand<-args[3]
sub<-args[4]
masterpath<-args[5]

library('randomForest')
library('e1071')
library('glmnet')

### read in the RF, LogR and NB models ###
filename <- paste(as.character(masterpath),"/ALU/02_PE_level1/RFIII_1.rds", sep="")
rf1 <- readRDS(filename)
filename <- paste(as.character(masterpath),"/ALU/02_PE_level1/RFIII_2.rds", sep="")
rf2 <- readRDS(filename)
filename <- paste(as.character(masterpath),"/ALU/02_PE_level1/LogRIII.rds", sep="")
logr <- readRDS(filename)
filename <- paste(as.character(masterpath),"/ALU/02_PE_level1/NBIII.rds", sep="")
nb <- readRDS(filename)

    filename <- paste(as.character(datpath),"/",as.character(sub),"/retro_v",as.character(retro),"_",as.character(strand),"/",as.character(sub),".pe.ALU.matrix", sep="")
    clone1 <- read.table(filename, sep="\t", header=T)
    filename <- paste(as.character(datpath),"/",as.character(sub),"/retro_v",as.character(retro),"_0/",as.character(sub),".pe.ALU.matrix", sep="")
    clone2 <- read.table(filename, sep="\t", header=T)
    clone <- rbind(clone1, clone2)
    subclone <- subset(clone, (clone[,7] == 1 | clone[,8] == 1) & ref == 0)
    subclone$gap <- subclone$gap/(abs(subclone$refpos1 - subclone$refpos2))
    subclone$frag <- abs(subclone$refpos1 - subclone$refpos2) / 151

    test1 <- subset(subclone, anchor_split==1)
    nx1 <- test1[,c("seg", "map", "anchor_align", "anchor_len", "anchor_mm", "anchor_direct", "anchor_insert", "anchor_seg", "anchor_map", "anchor_XS", "anchor_AS", "anchor_mapQ", "depth","map_len","end3", "end5", "trans2","direction", "refpos", "dist", "c1", "upstream", "frag")]
    pred.RF1 <- predict(rf1,newdata=nx1,type='prob')[,2]
    pred1 <- data.frame(cbind(as.character(test1$chr), test1$cord1, test1$cord2, as.character(test1$read), pred.RF1))
    filename <- paste(as.character(datpath),"/",as.character(sub),"/retro_v",as.character(retro),"_",as.character(strand),"/ALU/",as.character(sub),".pe.pred.A1.txt", sep="")
    write(t(pred1), file=filename, ncol=5, sep="\t")

    test2 <- subset(subclone, anchor_split==0)
    nx2 <- test2[,c("seg", "map", "anchor_align", "anchor_len", "anchor_mm", "anchor_XS", "anchor_AS", "anchor_mapQ", "depth","map_len","end3", "end5", "trans2","direction", "refpos", "dist", "c1", "upstream", "frag")]
    pred.RF2 <- predict(rf2,newdata=nx2,type='prob')[,2]
    pred2 <- data.frame(cbind(as.character(test2$chr), test2$cord1, test2$cord2, as.character(test2$read), pred.RF2))
    filename <- paste(as.character(datpath),"/",as.character(sub),"/retro_v",as.character(retro),"_",as.character(strand),"/ALU/",as.character(sub),".pe.pred.A2.txt", sep="")
    write(t(pred2), file=filename, ncol=5, sep="\t")

    test3 <-subclone
    nx3 <- test3[,c("seg", "map", "anchor_align", "anchor_len", "anchor_mm", "anchor_XS", "anchor_AS", "anchor_mapQ", "depth","map_len","end3","end5","trans2", "direction", "refpos", "dist", "c1", "upstream", "frag")]
    pred.RF3 <- predict(logr,newdata=nx3,type="response")
    pred3 <- data.frame(cbind(as.character(test3$chr), test3$cord1, test3$cord2, as.character(test3$read), pred.RF3))
    filename <- paste(as.character(datpath),"/",as.character(sub),"/retro_v",as.character(retro),"_",as.character(strand),"/ALU/",as.character(sub),".pe.pred.logR.txt", sep="")
    write(t(pred3), file=filename, ncol=5, sep="\t")

    pred.RF4 <- predict(nb,newdata=nx3,type="raw")[,2]
    pred4 <- data.frame(cbind(as.character(test3$chr), test3$cord1, test3$cord2, as.character(test3$read), pred.RF4))
    filename <- paste(as.character(datpath),"/",as.character(sub),"/retro_v",as.character(retro),"_",as.character(strand),"/ALU/",as.character(sub),".pe.pred.NB.txt", sep="")
    write(t(pred4), file=filename, ncol=5, sep="\t")
