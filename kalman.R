#Kalman Filter
Kalman_filter<-function(Data,H,R,M,Q){
  if(class(Data)=="matrix"){
    size<-dim(Data);  
  }else{
    stop("Please ensure the class of input!");
  }  
  m<-size[1];T<-size[2]; 
  x_hat<-Data;
  P<-array(rep(1,T*m^2),dim=c(m,m,T));
  P_up<-P_pred<-P;
  for(i in 1:(T-1)){
    #Prediction phase
    x_hat[,i+1]<-M%*%x_hat[,i];        
    # State estimate predict
    P[,,i+1]= M%*%P[,,i]%*%t(M)+Q;    
    # Tracking error covariance predict
    P_pred[,,i+1]<-P[,,i+1];    
    #Kalman gain
    K=P[,,i+1]%*%t(H)%*%solve(t(H)%*%P[,,i+1]%*%H+R);
    #Updata step
    x_hat[,i+1]<-x_hat[,i+1]+K%*%(Data[,i+1]-H%*%x_hat[,i+1]);            
    # State estimate update
    P[,,i+1]<-P[,,i+1]-K%*%H%*%P[,,i+1];                       
    # Tracking error covariance update
    P_up[,,i+1]<-P[,,i+1];   
 }
 output<-t(x_hat);
}

#a simple example
R=0.1;H=1;Q=0.5;M=0.7;T=100;
Xt<-arima.sim(list(order=c(1,0,0),ar=M,sd=sqrt(Q)),n=T);
Xt<-as.vector(Xt);
Data<-H*Xt+rnorm(100,0,R);
Data<-matrix(Data,nrow=1);
output<-Kalman_filter(Data,H,R,M,Q);
output;
draw<-data.frame(t=c(1:100,1:100,1:100),data=c(as.vector(output),
      as.vector(Data),Xt),
      type=c(rep("filter",100),rep("data",100),rep("truth",100) ));
library(ggplot2);
qplot(t,data,data=draw,colour=type,geom="line",)+geom_line(size=1)
