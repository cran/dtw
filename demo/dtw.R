
## $Id:$
## Show off some capabilities, on the same sine/cosine
## alignment problem

library(dtw);
idx<-seq(0,6.28,len=100);

## A noisy sine wave as query
query<-sin(idx)+runif(100)/10;
template<-cos(idx);

## A cosine is for template; sin and cos are offset by 25 samples
dtw(query,template,keep=TRUE)->alignment;

## Display the mapping with reference
dtw(query,template,keep=TRUE)->alignment;
plot(alignment$index1,alignment$index2);
lines(1:100-25,col="red")

## We can also make a nice three-way plot
## Beware of the template's y axis, may be confusing
dtwPlot(alignment,x=query,y=template,type="threeway");


## A profile of the cumulative distance matrix
## similar to: dtwPlot(alignment,type="density");

image(alignment$costMatrix,col=terrain.colors(100),x=1:100,y=1:100,
	xlab="Query (noisy sine)",ylab="Template (cosine)");
contour(alignment$costMatrix,x=1:100,y=1:100,add=TRUE);

lines(alignment$index1,alignment$index2,col="red",lwd=2);


## Do the same with asymmetric step (takes ~5 s)
dtw(query,template,keep=TRUE,step=asymmetric)->ita;
dtwPlot(ita,type="density",main="Sine and cosine, asymmetric step");


## Windowing functions (global constraints) can be applied and plotted
dtwWindow.plot(itakuraWindow, main="So-called Itakura parallelogram window");


## Symmetric step with global parallelogram-shaped constraint
## Note how long (>2 steps) horizontal stretches are allowed within the window.
dtw(query,template,keep=TRUE,window=itakuraWindow)->ita;
dtwPlot(ita,type="density",main="Symmetric step with Itakura parallelogram window");


## Asymmetric step with slope constraint
## The local constraint: three sides of the parallelogram are seen
dtw(query,template,keep=TRUE,step=asymmetricItakura)->ita;
dtwPlot(ita,type="density",main="Slope-limited asymmetric step (Itakura)");

