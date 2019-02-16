clear all
ds = datastore('house_data_complete.csv','TreatAsMissing','NA',.....
'MissingValue',0,'ReadSize',25000);
T = read(ds);
size(T);
Alpha=0.01;
m=length(T{:,1});

crossVal=ceil(m*0.6);

U0=T{1:crossVal,2};
U=T{1:crossVal,4:19};
U1=T{1:crossVal,20:21};

X=[ones(crossVal,1) U U1 U.^2 U.^3 U.^4];

n=length(X(1,:));

%Form of normalization using mean
for w=2:n
    if max(abs(X(:,w)))~=0
    X(:,w)=(X(:,w)-mean((X(:,w))))./std(X(:,w));
    end
end

Y=T{1:crossVal,3}/mean(T{1:crossVal,3});
Theta=zeros(n,1);
k=1;

E(k)=(1/(2*crossVal))*sum((X*Theta-Y).^2);

R=1;
while R==1
Alpha=Alpha*1;
Theta=Theta-(Alpha/crossVal)*X'*(X*Theta-Y);
k=k+1;
E(k)=(1/(2*crossVal))*sum((X*Theta-Y).^2);
if E(k-1)-E(k)<0
    break
end 
q=(E(k-1)-E(k))./E(k-1);
if q <.0001;
    R=0;
end
end

figure(1)
plot(E)


%%%%%%% first 20% of data
crossValnew=ceil(m*0.2);
U0new=T{crossVal+1:crossVal+crossValnew,2};
Unew=T{crossVal+1:crossVal+crossValnew,4:19};
U1new=T{crossVal+1:crossVal+crossValnew,20:21};

U2new=Unew.^2;

Xnew=[ones(crossValnew,1) Unew U1new Unew.^2 Unew.^3 Unew.^4];
ThetaNew=Theta;
Ynew=T{crossVal+1:crossVal+crossValnew,3}/mean(T{crossVal+1:crossVal+crossValnew,3});
for w=2:n
    if max(abs(Xnew(:,w)))~=0
    Xnew(:,w)=(Xnew(:,w)-mean((Xnew(:,w))))./std(Xnew(:,w));
    end
end
k1=1;
Enew(k1)=(1/(2*crossValnew))*sum((Xnew*ThetaNew-Ynew).^2);


%%%%%%% second 20% of data
crossValnew2=m-(crossVal+crossValnew);
U0new2=T{crossVal+crossValnew+1:end,2};
Unew2=T{crossVal+crossValnew+1:end,4:19};
U1new2=T{crossVal+crossValnew+1:end,20:21};

U2new2=Unew2.^2;

Xnew2=[ones(crossValnew2,1) Unew2 U1new2 Unew2.^2 Unew2.^3 Unew2.^4];
ThetaNew2=Theta;
Ynew2=T{crossVal+crossValnew+1:end,3}/mean(T{crossVal+crossValnew+1:end,3});
for w=2:n
    if max(abs(Xnew2(:,w)))~=0
    Xnew2(:,w)=(Xnew2(:,w)-mean((Xnew2(:,w))))./std(Xnew2(:,w));
    end
end
k2=1;
Enew2(k2)=(1/(2*crossValnew2))*sum((Xnew2*ThetaNew2-Ynew2).^2);

