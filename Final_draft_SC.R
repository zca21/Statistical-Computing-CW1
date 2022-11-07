## Question 1: Wild Life Population Density

#a

#Initializing the environment
setwd("~/Desktop/Statistical Computing/CW/Statistical-Computing-CW1")
eco.df <- read.table("ecoStudy.txt",header = T)

#Turning habitat variable into factor to allow for plotting and analysis of that variable
eco.df$habitat <- factor(eco.df$habitat,levels=c("A","B","C"),labels=c("A","B","C"))
attach(eco.df)

#i
round(min(density),3)
#minimum density is 0

#Selecting 1st quantile from the density column of eco-study dataset
round(quantile(density)[[2]],3)
#1st quantile of density is 1.413 to (3 d.p.)

round(median(density),3)
#Median of density is 2.744 (3 d.p.)

round(mean(density),3)
#Mean of density is 3.103 (3 d.p.)

#Selecting 3rd quantile
round(quantile(density)[[4]],3)
#3rd quantile of density is 4.193 (3 d.p.)

round(max(density),3)
#maximum value of density is 20.497 (3 d.p.)

#ii
table(habitat)
#For habitat A, B and C there are 100, 125 and 75 observations in each respectively

#b

#Creating new column in eco-study dataframe which is log of density+1
eco.df$logDensity <- log(density + 1)
attach(eco.df)

#Saving plot created below to pdf file called Boxplot.pdf into current directory (set at start up)
pdf(file="Boxplot.pdf")

#Creating graphic of side by side box plots
par(mfrow=c(2,1))

boxplot(density~habitat,
        main="Boxplot of density in each habitat",
        xlab="Habitat type",
        ylab="Density",
        col=c("red","blue","green"))

boxplot(logDensity~habitat,
        main="Boxplot of log density in each habitat",
        xlab="Habitat type",
        ylab="Log density",
        col=c("red","blue","green"),
        sub="Question 1: part b",
        cex.sub = 0.7)

dev.off()

#c

#i
#Creating vectors to store the observed densities of each habitat type
density.A <- density[habitat=="A"]
density.B <- density[habitat=="B"]
density.C <- density[habitat=="C"]

#calculating the denominator of F statistic
F.stat.denom <- (sum((density.A-mean(density.A))^2)+
                   sum((density.B-mean(density.B))^2)+
                   sum((density.C-mean(density.C))^2))/(300-3)

#calculating the numerator of F statistic
F.stat.numerator <- ((100*((mean(density.A)-mean(density))^2))+
                       (125*((mean(density.B)-mean(density))^2))+
                       (75*((mean(density.C)-mean(density))^2)))/2

#Calculating the F statistic
F.stat <- F.stat.numerator/F.stat.denom
#F statistic is 10.5 to 3 s.f.

#ii
#Finding p-value, as pf is probability that x from a RV X~F_{2}{297} is less than F.stat value the probability to be equal or larger is 1-pf
1-pf(F.stat,2,297)
#p-value is 0.0000400 to 3 s.f.

#d
#First, I create a function to calculate the numerator and denominator for one habitat type
F.stat <- function(vec=NULL){
  
  numerator<-(length(vec)*(mean(vec)-mean(density))^2)/2
  denominator<-(sum((vec-mean(vec))^2))/(300-3)
  
  #Storing in a list so can extract both values from the function output
  output <- list("num"=numerator,"den"=denominator)
  
  return(list("num"=numerator,"den"=denominator))
}

#Applying the function to each habitat type separately using the tapply function
tapply.out.Fstat<-tapply(density,habitat,FUN=F.stat)

#Calculating the F statistic by extracting and summing the numerator and denominator parts of the tapply output
tapply.F.stat<-sum(c(tapply.out.Fstat$A$num,tapply.out.Fstat$B$num,tapply.out.Fstat$C$num))/
  sum(c(tapply.out.Fstat$A$den,tapply.out.Fstat$B$den,tapply.out.Fstat$C$den))
#The answer is 10.5 (3 s.f.) as before

#e
#Function to calculate the residuals for one habitat
residual.func<- function(x){
  
  residuals <- x-mean(x)
  
  return(residuals)
}

#Calculating residuals for each habitat using tapply
residuals.list <- tapply(density, habitat, residual.func)
#Extracting residuals values from list and putting into a vector
residuals<-c(residuals.list$A,residuals.list$B,residuals.list$C)

pdf(file="QQ plot.pdf")
#Creating QQ plot
sorted.residuals<-sort(residuals)
sorted.normal.dist<-sort(rnorm(300))

plot(x=sorted.normal.dist,y=sorted.residuals,
     xlab="Theoretical quantiles",
     ylab="Sample quantiles",
     main="QQ plot of residuals",
     sub="Question 1: part e",
     cex.sub = 0.7)
#Adding in reference line
qqline(residuals, col = "red", lwd = 2)

#Viewing the QQ plot the tails diverge from the line (especially the upper tail which diviates significantly) 
#and the shape has a slight U. This gives evidence that the normality assumption has not been satisfied

dev.off()

detach(eco.df)
#----------------------------------------------------------------------------------------------------# 
## Question 2: Mixture Models

#a
calcLogMix2Gamma <- function(x,alpha,beta,p){
  
  #checking input argument are valid 
  #(note check numeric before using setequal so dont non-numeric input error)
  if(is.numeric(x)==F || setequal(x,abs(x))==F){
    stop("x is not a numeric vector of only positive values")
    
  }else if(length(alpha)!=2 || is.numeric(alpha)==F || setequal(alpha,abs(alpha))==F){
    stop("alpha needs to a numeric vector of length 2 of positive values")
    
  }else if(length(beta)!=2 || is.numeric(beta)==F || setequal(beta,abs(beta))==F){
    stop("beta needs to a numeric vector of length 2 of positive values")
    
  }else if(is.numeric(p)==F || length(p)!=1 || p>1 || p<0){
    stop("p must be a single numeric value between 0 and 1 inclusive")
  }
  
  #calculating the pdf for each value of x using the provided formula
  prob.dist<-p*(x^(alpha[1]-1)*exp(-beta[1]*x)*beta[1]^alpha[1])/gamma(alpha[1])+
    (1-p)*(x^(alpha[2]-1)*exp(-beta[2]*x)*beta[2]^alpha[2])/gamma(alpha[2])
  
  #calculating the log of the joint probability distribution
  log.prob.dist <- log(prod(prob.dist))
  
  return(log.prob.dist)
}

#c
#Creating arguments for calcLogMix2Gamma function
x <- c(0.06, 2.32, 4.81, 0.02, 2.33, 2.18, 0.83, 2.45, 2.10, 3.27)
alpha <- c(0.5, 6)
beta <- c(0.5, 2.5)
p <- 0.35

#Using calcLogMixGamma function with required inputs and rounding to 3 d.p.
round(calcLogMix2Gamma(x,alpha,beta,p),3)

#The log of the joint pdf is -13.928 (to 3 d.p.) using the parameters and observed values above
#-------------------------------------------------------------------------------------------------#
## Question 3: Riemann Integration
#a
smoothCurve <- function(x,alpha,beta,p,q,k){
  #To calculate g(x) I split up the term to allow for easier manual 
  #checking of values and reduce errors due to missing a bracket
  
  #calculating numerator and denominator in log term of g(x) function
  numerator <- (gamma(x+beta)^k)*(p^alpha)*(q^beta)
  denominator <- (gamma(alpha)^k)*(gamma(beta)^k)
  
  #combining parts of g(x) to calculate value of g(x)
  g <- sin((1/(2*alpha))*log(numerator/denominator))
  
  return(g)
}

#b
calcMidRiemannLoop <- function(xVec,alpha,beta,p,q,k){
  #setting up vector to store areas of each (x_{i-1},x_{i}) segment
  interval.area<-rep(0,length(xVec)-1)
  
  #loop to calculate the area in each interval
  #As we are told that the x intercept is not crossed in any interval we do not need to worry about 
  #having to split the area into 2 and calculating the area separately for part above and below x-axis
  
  for (i in 2:length(xVec)){
    #Calculating x input for g(x) function
    x <- (xVec[i]+xVec[i-1])/2
    
    #using equation 6 with absolute value as area cannot be negative
    interval.area[i-1] <- abs((xVec[i]-xVec[i-1])*smoothCurve(x,alpha,beta,p,q,k))
  }
  
  #Calculating total area by summing areas within each segment
  approx.area <- sum(interval.area)
  
  return(approx.area)
}

#c
calcMidRiemann <- function(xVec,alpha,beta,p,q,k){
  #Checking xVec argument
  if (is.numeric(xVec)==F || length(xVec)<2 || setequal(xVec,sort(xVec)==F)){
    stop("xVec argument must be numeric, of length greater than 1 with values in
         increasing order")
    
  }else if (length(alpha)>1 || is.numeric(alpha)==F || alpha<=0){
    stop("alpha must be a single positive numeric value")
    
  }else if (length(beta)>1 || is.numeric(beta)==F || beta<=0){
    stop("beta must be a single positive numeric value")
    
  }else if (length(p)>1 || is.numeric(p)==F){
    stop("p must be a single numeric value")
    
  }else if (length(k)>1 || is.numeric(k)==F){
    stop("k must be a single numeric value")
    
  }else if (length(q)>1 || is.numeric(q)==F){
    stop("q must be a single numeric value")
  }
  #Creating vector containing the intervals which we need to calculate the area of
  interval.region<-xVec[-1]-xVec[1:length(xVec)-1]
  
  #values of x to input into g(x) for each of the corresponding intervals
  g.func.input <-(xVec[-1]+xVec[1:length(xVec)-1])/2
  
  #calculating each intervals approximate area using equation (6) without summation
  interval.area <- abs(interval.region*smoothCurve(g.func.input,alpha,beta,p,q,k))
  
  #summing individual areas to find total area 
  approx.area <- sum(interval.area)
  
  return(approx.area)
}
#e
#Creating vector of intervals as specified
xVec <- seq(from=2,to=8.5,by=0.01)

round(calcMidRiemann(xVec,2.1,0.5,3,6,2),3)
#The approximate area using middle Riemann sum is 4.737 to 3 d.p.

#f
calcMidRiemannAreas <- function(xSeqList,alpha,beta,p,q,k){
  #Using sapply to apply calcMidRiemann function to each list element and return each calculated area as an element of a vector
  vector.of.areas <- sapply(xSeqList,calcMidRiemann,alpha=alpha,beta=beta,p=p,q=q,k=k)
  
  #Summing areas calculated for each list element to find total area
  total.area<-sum(vector.of.areas)
  
  return(total.area)
}

#g
#creating list of sequences as defined in question
xSeqList <- list(seq(from=3.5,to=8,by=0.01),seq(from=2,to=4.3,by=0.1),seq(from=1.07,to=9.012,by=0.001))

round(calcMidRiemannAreas(xSeqList,1.9,0.75,2.15,7,1.65),3)
#Total area is 10.847 to 3 d.p.

#-------------------------------------------------------------------------------------------------#
## Question 4: Hidden Markov Model
#a
weatherSeqProb <- function(weatherSeq,trProbs,initProbs){
  
  #Setting up vector to store each transition probability in
  transition.prob<-rep(0,length(weatherSeq))
  
  #Finding initial probability of chosen sequence
  if (weatherSeq[1]=="s"){
    transition.prob[1] <- initProbs[1]
  }else if (weatherSeq[1]=="c"){
    transition.prob[1]  <- initProbs[2]
  }else if (weatherSeq[1]=="r"){
    transition.prob[1]  <- initProbs[3]
  }
  #Assigning numeric values to weather sequence letters to be able to call indexes of matrix
  weatherSeq.numeric <- ifelse(weatherSeq=="s",1,ifelse(weatherSeq=="c",2,3))
  
  for (i in 2:length(weatherSeq)){
    #trProbs[i,j] is the probability to go from state i to state j and
    #i-1 in weatherSeq is last state and i is current state, thus can use
    #wetherSeq.numeric to access the transitional probability
    transition.prob[i] <- trProbs[weatherSeq.numeric[i-1],weatherSeq.numeric[i]]
  }
  
  #Calculating log total probability
  total.prob.log <- log(prod(transition.prob))
  
  return(total.prob.log)
}

#b
weatherSeq <- c("c","s","c","r","s","s")
trProbs <- matrix(c(0.5,0.4,0.1,0.33,0.35,0.32,0.3,0.3,0.4),nrow=3,ncol=3,byrow = T)
initProbs <- c(0.45,0.25,0.3)

round(weatherSeqProb(weatherSeq,trProbs,initProbs),3)
#log probability of weather sequence is -6.448 (3 d.p.)

#c
weatherColourProbs <- function(colourSeq,emitProbs,weatherSeq,trProbs,initProbs){
  #We can rearrange the probabilities 
  
  #calculate probabilities for weather sequence as before
  total.weather.prob <- exp(weatherSeqProb(weatherSeq,trProbs,initProbs))
  
  #Creating numeric vectors of weather and jacket sequence
  weatherSeq.numeric <- ifelse(weatherSeq=="s",1,ifelse(weatherSeq=="c",2,3))
  colourSeq.numeric <- ifelse(colourSeq=="B",1,2)
  
  #Creating vector to store the probabilities of selecting the ith jacket colour  
  #given the ith weather condition in the sequence
  colour.trans.prob <- rep(0,length(colourSeq))
  for (i in 1:length(colourSeq)){
    #Storing prob of selecting the jacket in the ith sequence position for the ith weather condition in the sequence
    colour.trans.prob[i] <- emitProbs[weatherSeq.numeric[i],colourSeq.numeric[i]]  
  }
  
  #calculating total probability of observing sequence of weather and jacket colour
  total.prob <- total.weather.prob*prod(colour.trans.prob)
  
  log.total.prob <- log(total.prob)
  return(log.total.prob)
}

#d
#Setting up arguments to input into function
weatherSeq <- c("r","s","c","r","c","r","s","s")
colourSeq <- c("B","W","W","B","B","W","W","W")
initProbs <- c(0.35,0.45,0.2)
trProbs <- matrix(c(0.55,0.25,0.2,0.25,0.35,0.4,0.2,0.15,0.65),nrow=3,ncol=3,byrow = T)
emitProbs <- matrix(c(0.2,0.8,0.55,0.45,0.9,0.1),nrow=3,ncol=2,byrow = T)

round(weatherColourProbs(colourSeq,emitProbs,weatherSeq,trProbs,initProbs),3)
#Answer is -15.121 (3 d.p.)
