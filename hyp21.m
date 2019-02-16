clc
clear all
close all
 
ds = datastore('heart_DD.csv','TreatAsMissing','NA',.....
    'MissingValue',0,'ReadSize',250);
T = read(ds);

m= length(T{:,1});
Alpha=0.01;
lambda=100;
crossVal=ceil(m*0.8);

U=T{1:crossVal,1:13};

Y=T{1:crossVal,14};
 
X=[ones(crossVal,1) U U.^2 exp(-U)];

n=length(X(1,:)); 
 for w=2:n
    if max(abs(X(:,w)))~=0
    X(:,w)=(X(:,w)-mean((X(:,w))))./std(X(:,w));
    end
 end
 
 theta=zeros(n,1);
 
 h=1./(1+exp(-X*theta));  %sigmoid function
 
 k=1;
 J(k)=-(1/crossVal)*sum(Y.*log(h)+(1-Y).*log(1-h))+(lambda/(2*crossVal))*sum((theta).^2);  %cost function
 
 grad=zeros(size(theta,1),1);     %gradient vector
  
 for i=1:size(grad)
     grad(i)=(1/crossVal)*sum((h-Y)'*X(:,i));
 end
 

R=1;
while R==1
Alpha=Alpha*1;
theta=theta-(Alpha/crossVal)*X'*(h-Y);
h=1./(1+exp(-X*theta));  %sigmoid function
k=k+1

J(k)=(-1/crossVal)*sum(Y.*log(h)+(1-Y).*log(1-h))+(lambda/(2*crossVal))*sum((theta).^2);
if J(k-1)-J(k) <0 
    break
end 
q=(J(k-1)-J(k))./J(k-1);
if q <.00001
    R=0;
end
end

plot (J)

%%%%%%%%%%%%%%%%%%%%%% the first 10%
crossValnew=ceil(m*0.1);

U1=T{crossVal+1:crossVal+crossValnew,1:13};
Y1=T{crossVal+1:crossVal+crossValnew,14};
 
X1=[ones(crossValnew,1) U1 U1.^2 exp(-U1)];

for w=2:n
    if max(abs(X1(:,w)))~=0
    X1(:,w)=(X1(:,w)-mean((X1(:,w))))./std(X1(:,w));
    end
 end
 
theta1=theta; 
h1=1./(1+exp(-X1*theta1));  %sigmoid function
 
k1=1;
JJ(k1)=-(1/crossValnew)*sum(Y1.*log(h1)+(1-Y1).*log(1-h1))+(lambda/(2*crossValnew))*sum((theta1).^2);  %cost function
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% the second 10%

crossValnew2=m-(crossVal+crossValnew);

U2=T{crossVal+crossValnew+1:end,1:13};
Y2=T{crossVal+crossValnew+1:end,14};
X2=[ones(crossValnew2,1) U2 U2.^2 exp(-U2)];

for w=2:n
    if max(abs(X2(:,w)))~=0
    X2(:,w)=(X2(:,w)-mean((X2(:,w))))./std(X2(:,w));
    end
 end
 
theta2=theta; 
h2=1./(1+exp(-X2*theta2));  %sigmoid function
 
k2=1;
JJJ(k2)=-(1/crossValnew2)*sum(Y2.*log(h2)+(1-Y2).*log(1-h2))+(lambda/(2*crossValnew2))*sum((theta2).^2); %cost function
 

