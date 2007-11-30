library(dtw);

ldist<-matrix(1,nrow=6,ncol=6);  # Matrix of ones
ldist[2,]<-0; ldist[,5]<-0;      # Mark a clear path of zeroes
ldist[2,5]<-.01;		 # Forcely cut the corner

ds<-dtw(ldist);			 # DTW with user-supplied local cost matrix
da<-dtw(ldist,step="a");	 # Also compute the asymmetric


ds$distance
ds$index1
ds$index2

da$index1
da$index2
da$distance
