#Generate clustered spatial point process
#
library(spatstat)
library(maptools)


#define variables
#variables that control the size and strength of clusters
val.at.center=1
effect.range=10
background=0.001

#variables that control the number of points and spatial dimensions
Pointnum=100
Xmin=-50
Xmax=50
Ymin=-50
Ymax=50

#define the center locations and set up the distance vector
centers=matrix(data=c(-25,-25,25,25,-25,25,-25,25),nrow=4,ncol=2)
dist=matrix(nrow=4,ncol=1)

#define outputs
output.X=matrix(nrow = Pointnum, ncol = 1)
output.Y=matrix(nrow = Pointnum, ncol = 1)

#precalcs - calculate the slope of the clustering effect
slope=(val.at.center-background)/effect.range

#set a counter
outcounter=0

#main for loop
for (i in 1:100000){
  #generate a random candidate point
  x.candidate=runif(1, min=Xmin, max=Xmax)
  y.candidate=runif(1, min=Ymin, max=Ymax)
  
  #calculate the distance between the candidate point and the nearest cluster center
  for (j in 1:4){
    dist[j]=sqrt((x.candidate-centers[j,1])^2+(y.candidate-centers[j,2])^2)
  }
  min.dist=min(dist)
  
  #calculate the probability of retaining the candidate point
  if(min.dist<effect.range){
    prob=val.at.center-slope*min.dist
  }
  else
    prob=background
  
  #roll the dice to see if you keep the candidate point
  testval=runif(1,min=0,max=1)  
  if (testval<prob){
    outcounter=outcounter+1
    keep=1
    output.X[outcounter]=x.candidate
    output.Y[outcounter]=y.candidate
  }
  
  #if you've reached your target number of points, break from the for loop
  if(outcounter==Pointnum){
    break
  }
}

#plot the points and the cluster centers
plot(output.X,output.Y)
points(centers,type="p", col="red", pch=21, bg="red")

#create a point pattern object for analysis using the spatstat library
output_ppp = ppp(output.X, output.Y, c(Xmin,Xmax), c(Ymin,Ymax))
plot(output_ppp)

#create a CSR object for analysis 

lambda.u <- function(x,y){1000 * x^2 * y^2 + 1000}
output_csr <- rpoispp(lambda = lambda.u, win = owin(c(0,1), c(0,1)))
plot(output_csr)


####Q1.1####

#The coordinates in the center matrix control the spatial extent of clusters, you can change where the centers
#are located by manipulating the coordinates (e.g. changing everything to 50 will space it out)

####Q1.2####
#effect.range and background have major controls over the strength of clustering 
#The Y and X min/max have some control over the strength of clustering 

####Q1.3####
#run these two plot functions 
plot(output_ppp)
plot(output_csr)

####Q2####
quadrat.test(output_ppp, alternative = "clustered")
quadrat.test(output_csr, alternative = "regular")
#PPP quadrat test X2 = 291.5, df = 24, p-value <2.2e-16
#CSR quadrat test X2 = 46.633, df = 24, p-value 0.9963
#P value close to 1 means this does not differ from CSR
#(which hopefully makes sense, because it was coded as a true CSR)

####Q3####

#Sample example
X <- runifpoint(100)
K <- Kest(X)
plot(K)

pest <- Kest(output_ppp)
plot(pest)
cest <- Kest(output_csr)
plot(cest)

plot(envelope(output_ppp))
plot(envelope(output_csr))
#For the point pattern we generated, at small scales the points are randomly distributed, 
#however as we increase scale it becomes obvious there is a pattern to the points - they're clustered 
