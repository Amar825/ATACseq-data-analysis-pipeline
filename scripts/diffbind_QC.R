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
                       	peaks="../SRR10261595_peak_peaks.xls",
                       	peak.caller="macs", sampID="hybrid2",condition = "Hybrid", replicate=2, bamReads = "../SRR10261595_SC_subset_dedup.bam")
SC_vs_hybrid<-dba.peakset(SC_vs_hybrid,
                      	peaks="../SRR10261596_peak_peaks.xls",
                      	peak.caller="macs", sampID="hybrid3",condition = "Hybrid", replicate=3, bamReads = "../SRR10261596_SC_subset_dedup.bam")



### Plot Venn diagrams of replicates
png('./SC_replicates_atacseq.png', units="in", width=7, heigh=7, res=800)
dba.plotVenn(SC_vs_hybrid,SC_vs_hybrid$masks$Parent, main = "Open chromatic region overlaps in S. cerevisiae replicates")
dev.off()

png('./hybrid_SC_replicates_atacseq.png', units="in", width=7, heigh=7, res=800)
dba.plotVenn(SC_vs_hybrid,SC_vs_hybrid$masks$Hybrid, main = "Open chromatic region overlaps in the hybrid replicates")
dev.off()

#### Plot heatmap of correlations
SC_vs_hybrid_counts<-dba.count(SC_vs_hybrid, bParallel = TRUE, score=DBA_SCORE_READS)

png('./heatmap.png', units="in", width=7, heigh=7, res=800)
plot(SC_vs_hybrid_counts, main="Correlation plot of studied samples")
dev.off()
