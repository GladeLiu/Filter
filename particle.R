set.seed(2);
#simulate the system
#initial state value
t<-100; #final time
f_state<-function(x,k){
   output<-x/2+25*x/(1+x^2)+8*cos(1.2*k);
}
f_observe<-function(x){
   output<-x^2/20;
}
x<-rep(rnorm(1,0,sqrt(5)),t)
y<-f_observe(x[1]);
for(i in 2:t){
   x[i]<-f_state(x[i-1],i-1)+rnorm(1,0,sqrt(10));
   y[i]<-f_observe(x[i])+rnorm(1,0,1);
}

#Particle filter
library(pracma);
ParticleFilter<-function(f,h,pdf_v,Q,P0,M,y,algorithm=1){
 Q<-as.matrix(Q);
 P0<-as.matrix(P0);
 M<-as.matrix(M);
 n<-size(P0,2);
 x<-sqrt(P0)%*%randn(n,M); #Initialize particles
 tf<-size(y,2);
 xhat<-rep(0,tf);
 if(algorithm==1){
   for(t in 1:tf){
      e<-repmat(y[t],1,M)-f_observe(x); #Calculate weights
      w<-pdf_v(e); #The likelihood function
      w<-w/sum(w); #Normalize importance weights
      xhat[t]<-sum(repmat(w,n,1)*x,2);
      ind<-resampling(w); #Resampling
      x<-x[,ind]; #The new particles
      x<-f_state(x,t)+sqrt(Q)%*%randn(n,M); #Time update
   }
 }else{
   for(t in 1:tf){
     x<-f_state(x,t)+sqrt(Q)%*%randn(n,M); #Time update
     e<-repmat(y[t],1,M)-f_observe(x); #Calculate weights
     w<-pdf_v(e); #The likelihood function
     w<-w/sum(w); #Normalize importance weights
     xhat[t]<-sum(repmat(w,n,1)*x,2);
     ind<-resampling(w); #Resampling
     x<-x[,ind]; #The new particles
   }
 }
 xhat;
}

#resampling
resampling<-function(w){
   wc<-cumsum(w);
   M<-length(w);
   u<-(seq(0,M-1)+runif(1,0,1))/M;
   i<-matrix(0,1,M); 
   k<-1;
   for(j in 1:M){
      while(wc[k]<u[j]){
        k<-k+1;
      }
    i[j]<-k;
   }
   i;
}

pdf_v<-function(x){
  1/(2*pi*1)^(1/2)*exp(-(x^2)/(2*1))
}

xTrue<-x;
xhat<-PF(f_state,f_observe,pdf_v,Q=10,P0=5,M=1000,y,algorithm=2);
draw<-data.frame(t=c(1:100,1:100),data=c(as.vector(xhat),xTrue),
      type=c(rep("filter",100),rep("truth",100) ));
library(ggplot2);
qplot(t,data,data=draw,colour=type,geom="line",)+geom_line(size=1)
