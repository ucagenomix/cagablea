
#####################################################
#generate .bed files with the coordinates of apmlified fragments:

primers <- read.table("C:\\Users\\IPMC\\Desktop\\nCoV-2019.bed", header = F, sep = "\t")
names(primers) <- c("Ref", "first", "last", "L_or_R", "pair_N", "plusminus")

#divide in two pools:
P1 <- primers[primers$pair_N == "nCoV-2019_1",]
P2 <- primers[primers$pair_N == "nCoV-2019_2",]

#coordinates of left and right primers:
P1_left <-  P1[c(grep ("LEFT", P1$L_or_R)), ]
P1_right <- P1[c(grep ("RIGHT", P1$L_or_R)), ]

P2_left <-  P2[c(grep ("LEFT", P2$L_or_R)), ]
P2_right <- P2[c(grep ("RIGHT", P2$L_or_R)), ]

#coordinates of amplified fragments:
#pool 1:
L_P1 <- P1_left$last + 1 #left coordinate
R_P1 <- P1_right$first - 1 #right coordinate

#pool 2:
L_P2 <- P2_left$last + 1 #left coordinate
R_P2 <- P2_right$first - 1 #right coordinate

#bed file for Pool1:
Pool1 <- as.data.frame(cbind
                       (Ref = P1_right$Ref,
                        first = L_P1,
                        last = R_P1))
Pool1$first <- as.integer(Pool1$first)
Pool1$last <- as.integer(Pool1$last)
Pool1$Ref <- as.character(P1_right$Ref)

write.table(Pool1, "pool1.bed", sep = "\t", col.names = F, row.names = F)

#bed file for Pool2:
Pool2 <- as.data.frame(cbind
                       (Ref = P2_right$Ref,
                         first = L_P2,
                         last = R_P2))
Pool2$first <- as.integer(Pool2$first)
Pool2$last <- as.integer(Pool2$last)
write.table(Pool2, "pool2.bed", sep = "\t", col.names = F, row.names = F)

#visualize coverage

plot_depth <- function (mypath) {
  depth <- read.table(mypath, header = F, sep = "\t")
  
  #1) The number of features in A that overlapped the B interval.
  #2) The number of bases in B that had non-zero coverage.
  #3) The length of the entry in B.
  #4) The fraction of bases in B that had non-zero coverage.
  
  names(depth) <- c("Ref", "first", "last", "L_or_R", "pair_N", "plusminus", "feat_overlap","nbp_nonzero_cov", "length", "nonzero_fraction")
  
  #some plots:
  library(ggplot2)
  #ggplot(depth, aes(first, feat_overlap)) + geom_point() + facet_grid(~ pair_N)
  ggplot(depth, aes(x = pair_N, y= feat_overlap)) + geom_boxplot()
}

mypath <- ("C:\\Users\\IPMC\\Desktop\\depth.txt")

#plot_depth("C:\\Users\\IPMC\\Desktop\\depth.txt")
ggplot(depth, aes(x = first, y= feat_overlap)) + geom_point()+facet_grid(~ pair_N)
