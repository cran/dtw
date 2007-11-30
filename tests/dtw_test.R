library(dtw);

### Synthetic example: check indexes, distance
ldist<-matrix(1,nrow=6,ncol=6);  # Matrix of ones
ldist[2,]<-0; ldist[,5]<-0;      # Mark a clear path of zeroes
ldist[2,5]<-.01;		 # Forcely cut the corner

ds<-dtw(ldist);			 # DTW with user-supplied local cost matrix
ds$distance			 # 2
ds$index1			 # 1 2 2 2 2 3 4 5 6 6
ds$index2			 # 1 1 2 3 4 5 5 5 5 6

da<-dtw(ldist,step="a");	 # Also compute the asymmetric
da$distance			 # 2
da$index1			 # 1 2 3 4 5 6
da$index2			 # 1 3 5 5 5 6

###  Synthetic example: verify native output
ds<- globalCostMatrix(ldist) 
dsn<- globalCostNative(ldist)
all.equal(ds,dsn)		 # TRUE


###  Sine/cosine example: verify native output
### there may be a random chance of failing due to rounding errors
idx<-seq(0,6.28,len=100);
query<-sin(idx)+runif(100)/10;	
template<-cos(idx)
ldist<-outer(query,template,FUN=euclideanSquared)
ds<- globalCostMatrix(ldist)
dsn<- globalCostNative(ldist)
all.equal(ds,dsn)		# TRUE

