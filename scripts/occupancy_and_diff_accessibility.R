library(DiffBind)

### Set the working directory
setwd("/path/to/mapping/diffbind/")

### Load the peak files from MACS2 to DiffBind
SC_vs_hybrid<-dba.peakset(NULL,
                      	peaks="../SRR10261591_peak_peaks.xls",
                      	peak.caller="macs", sampID="SC1",condition = "Parent", replicate=1, bamReads = "../SRR10261591_SC_subset_dedup.bam")
SC_vs_hybrid<-dba.peakset(SC_vs_hybrid,
                      	peaks="../SRR10261592_peak_peaks.xls",
                      	peak.caller="macs", sampID="SC2",condition = "Parent", replicate=2, bamReads = "../SRR10261592_SC_subset_dedup.bam")
SC_vs_hybrid<-dba.peakset(SC_vs_hybrid,
                      	peaks="../SRR10261593_peak_peaks.xls",
                      	peak.caller="macs", sampID="SC3",condition = "Parent", replicate=3, bamReads = "../SRR10261593_SC_subset_dedup.bam")


SC_vs_hybrid<-dba.peakset(SC_vs_hybrid,
                      	peaks="../SRR10261594_peak_peaks.xls",
                      	peak.caller="macs", sampID="hybrid1",condition = "Hybrid", replicate=1, bamReads = "../SRR10261594_SC_subset_dedup.bam")
SC_vs_hybrid<-dba.peakset(SC_vs_hybrid,
                      	peaks="../SRR10261596_peak_peaks.xls",
                      	peak.caller="macs", sampID="hybrid3",condition = "Hybrid", replicate=2, bamReads = "../SRR10261596_SC_subset_dedup.bam")


### Occupancy analysis

SC_vs_hybrid_occupancy<-dba.peakset(SC_vs_hybrid, consensus = DBA_CONDITION,minOverlap = 0.33)
 

png('./occupancy_consensus_peaks_SC_vs_hybrid.png', units="in", width=7, heigh=7, res=800)
dba.plotVenn(SC_vs_hybrid_occupancy,SC_vs_hybrid_occupancy$masks$Consensus, main = "Overlap of consensus peaks of SC and hybrid")
dev.off()

#### Affinity analysis - differential accessibility analysis
counts<-dba.count(SC_vs_hybrid, bParallel = TRUE, score=DBA_SCORE_READS)
SC_vs_hybrid_counts<-dba.contrast(counts, categories=DBA_CONDITION,minMembers = 2)


SC_vs_hybrid_analyzed<-dba.analyze(SC_vs_hybrid_counts)
 
SC_vs_hybrid_DE_peaks <- dba.report(SC_vs_hybrid_analyzed)
SC_vs_hybrid_DE_peaks_DF<-as.data.frame(SC_vs_hybrid_DE_peaks)
# Filtering
write.table(SC_vs_hybrid_DE_peaks_DF[abs(SC_vs_hybrid_DE_peaks_DF$Fold)>1 & SC_vs_hybrid_DE_peaks_DF$FDR<0.05,], "./diff_regions_hybrid_vs_parent.txt", sep = "\t", quote = F)
