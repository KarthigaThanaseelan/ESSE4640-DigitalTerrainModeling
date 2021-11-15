%% ESSE4640 DTM Lab 3 
clc
clear all
%% Read Data
% read file
data = xlsread('lab3esse4640');
%%Initial (Given) Coordinates
xi = data(:,1);
yi = data(:,2);
zi = data(:,3);


%% Generate 2nd and 3rd order trend surfaces
% Initial coeffiecnt values 
n = 92;
x = sum(xi);
y = sum(yi);
z = sum(zi);
% Estimating the corresponding polynomial coefficient
X = xi;
Y = yi;
Z = zi;
xy = xi.* yi;
xz = xi.* zi;
yz = yi.* zi;
for i = 2 : 8
    X(:,i) = x.^i;
    Y(:,i) = y.^i;
end
for i = 1 : 8
    SumX(i) = sum(X(:,i));
    SumY(i) = sum(Y(:,i));
end
xySum = sum(xy);
xzSum = sum(xz);
yzSum = sum(yz);

% Estimate using the normalized coordinates
Xmx = max(xi); Xmn = min(xi);
Ymx = max(yi); Ymn = min(yi);
Zmx = max(zi); Zmn = min(zi);

Xn = (xi-Xmn)/(Xmx-Xmn);
Yn = (yi-Xmn)/(Ymx-Ymn);
SumZn(:,1) = (zi-Xmn)/(Zmx-Zmn);
xyn = Xn.* Yn;
xzn = Xn.* SumZn(:,1);
yzn = Yn.* SumZn(:,1);
for i = 2 : 8
    Xn(:,i) = Xn(i,1).^i;
    Yn(:,i) = Yn(i,1).^i;
end
for i = 1 : 8
    SumXn(i) = sum(Xn(:,i));
    SumYn(i) = sum(Yn(:,i));
end
xynSum = sum(xyn);
xznSum = sum(xzn);
yznSum = sum(yzn);
SumZn = sum(SumZn(:,1));

% Estimate using the centroid coordinates
Xcent = xi/n;
Ycent = yi/n;
Zcent = zi/n;
SumXcent = sum(Xcent);
SumYcent = sum(Ycent);

for i = 2 : 8
    Xcent(:,i) = x.^i;
    Ycent(:,i) = y.^i;
end

xycent = xi.* yi;
xzcent = xi.* zi;
yzcent = yi.* zi;

for i = 1 : 8
    SumXcent(i) = sum(Xcent(:,i));
    SumYcent(i) = sum(Ycent(:,i));
end
xycentSum = sum(xycent);
xzcentSum = sum(xzcent);
yzcentSum = sum(yzcent);
SumZcent = sum(zi);

%% Normalized First Order Polynomial Equation
% Design Matrix
An = [n   SumXn(:,1)      SumYn(:,1);
     SumXn(:,1) SumXn(2) xynSum;
     SumYn(:,1)  xynSum  SumYn(2)];
Zon = [SumZn xznSum yznSum];
con = An/Zon
condn = abs (max(con)/min(con))
 
%% Centroid First Order Polynomial Equation
% Design Matrix
Acent = [n     SumXcent(:,1)       SumYcent(:,1);
         SumXcent(:,1)  SumXcent(2) xycentSum;
         SumYcent(:,1)  xycentSum  SumYcent(2)];
Zocent = [SumZcent xzcentSum yzcentSum];
cocent = Acent/Zocent
condcent = abs (max(cocent/min(cocent)));
 
%% Second Order Polynomial Equation
% Design Matrix
A2 = [n         x           y           SumX(2)         xySum           SumY(2);
      x         SumX(2)     xySum       SumX(3)         SumX(2)*y       SumY(2)*x;
      y         xySum       SumY(2)     SumX(2)*y       SumY(2)*x       SumY(3);
   SumX(2)      SumX(3)     SumX(2)*y   SumX(4)         SumX(3)*y       SumX(2)*SumY(2);
   xySum        SumX(2)*y   SumY(2)*x   SumX(3)*y       SumX(2)*SumY(2) SumY(3)*x;
   SumY(2)      SumY(2)*x   SumY(3)     SumX(2)*SumY(2) SumY(3)*x       SumY(4)];

Zo2 = [z z*x z*y SumX(2)*z  xySum*z SumY(2)*z];
co2 = A2/Zo2; 
cond2 = abs (max(co2)/min(co2));

%% Normalized Second Order Polynomial Equation
% Design Matrix
An2 = [n    SumXn(:,1)    SumYn(:,1)       SumXn(2)     xynSum      SumYn(2);
      SumXn(:,1)   SumXn(2) xynSum    SumXn(3)    SumXn(2)*SumYn(:,1)    SumY(2)*SumXn(:,1);
      SumYn(:,1)    xynSum  SumYn(2) SumXn(2)*SumYn(:,1)   SumY(2)*SumXn(:,1)    SumYn(3);
   SumXn(2) SumXn(3) SumXn(2)*SumYn(:,1) SumXn(4) SumXn(3)*SumYn(:,1)  SumXn(2)*SumYn(2);
   xynSum   SumXn(2)*SumYn(:,1)  SumYn(2)*SumXn(:,1) SumXn(3)*SumYn(:,1) SumXn(2)*SumYn(2) SumYn(3)*SumXn(:,1);
   SumYn(2) SumYn(2)*SumXn(:,1)  SumYn(3) SumXn(2)*SumYn(2) SumYn(3)*SumXn(:,1) SumYn(4)];
Zon2 = [SumZn(:,1) SumZn(:,1)*SumXn(:,1) SumZn(:,1)*SumYn(:,1) SumXn(2)*SumZn(:,1)  xynSum*SumZn(:,1) SumYn(2)*SumZn(:,1)];
con2 = An2/Zon2; 
condn2 = abs (max(con2)/min(con2));
 
%% Centroid Second   Order Polynomial Equation
% Design Matrix
Acent2 = [n             SumXcent(:,1)               SumYcent(:,1)               SumXcent(2)                 xycentSum                   SumYcent(2);
         SumXcent(:,1)  SumXcent(2)                 xycentSum                   SumXcent(3)                 SumXcent(2)*SumYcent(:,1)   SumYcent(2)*SumXn(:,1);
         SumYcent(:,1)  xycentSum                   SumYcent(2)                 SumXcent(2)*SumYcent(:,1)   SumYcent(2)*SumXcent(:,1)   SumYcent(3);
         SumXcent(2)    SumXcent(3)                 SumXcent(2)*SumYcent(:,1)   SumXcent(4)                 SumXcent(3)*SumYcent(:,1)   SumXcent(2)*SumYcent(2);
         xycentSum      SumXcent(2)*SumYcent(:,1)   SumYcent(2)*SumXcent(:,1)   SumXcent(3)*SumYcent(:,1)   SumXcent(2)*SumYcent(2)     SumYcent(3)*SumXcent(:,1);
         SumYcent(2)    SumYcent(2)*SumXcent(:,1)   SumYcent(3)                 SumXcent(2)*SumYcent(2)     SumYcent(3)*SumXcent(:,1)   SumYcent(4)];

Zocent2 = [SumZcent(:,1) SumZcent(:,1)*SumXcent(:,1) SumZcent(:,1)*SumYcent(:,1) SumXcent(2)*SumZcent(:,1)  xycentSum*SumZcent(:,1) SumYcent(2)*SumZcent(:,1)];
cocent2 = Acent2/Zocent2
condcent2 = abs (max(cocent2/min(cocent2)));
 

%% Third Order Polynomial Equation
%Design Matrix
A3 = [ n           x               y               SumX(2)              xySum          SumY(2)          SumX(3)          SumX(2)*y        x*SumY(2)        SumY(3);  
        x           SumX(2)          xySum          SumX(3)              SumX(2)*y        x*SumY(2)        SumX(4)          SumX(3)*y        SumX(2)*SumY(2)   x*SumY(3) ;
        y           xySum          SumY(2)          SumX(2)*y            x*SumY(2)        SumY(3)          SumX(3)*y        SumX(2)*SumY(2)   x*SumY(3)        SumY(4) ;
        SumX(2)      SumX(3)          SumX(2)*y        SumX(4)              SumX(3)*y        SumX(2)*SumY(2)   SumX(5)          SumX(4)*y        SumX(3)*SumY(2)   SumX(2)*SumY(3) ;
        xySum      SumX(2)*y        x*SumY(2)        SumX(3)*y            SumX(2)*SumY(2)   x*SumY(3)        SumX(4)*y        SumX(3)*SumY(2)   SumX(2)*SumY(3)   x*SumY(4);
        SumY(2)      x*SumY(2)        SumY(3)          SumX(2)*SumY(2)       x*SumY(3)        SumY(4)          SumX(3)*SumY(2)   SumX(2)*SumY(3)   x*SumY(4)        SumY(5);
        SumX(3)      SumX(4)          SumX(3)*y        SumX(5)              SumX(4)*y        SumX(3)*SumY(2)   SumX(6)          SumX(5)*y        SumX(4)*SumY(2)   SumX(3)*SumY(3);
        SumX(2)*y    SumX(3)*SumY(2)   SumX(2)*SumY(2)   SumX(4)*y            SumX(3)*SumY(2)   SumX(2)*SumY(3)   SumX(5)*y        SumX(4)*SumY(2)   SumX(3)*SumY(3)   SumX(2)*SumY(4);
        x*SumY(2)    SumX(2)*SumY(2)   x*SumY(3)        SumX(3)*SumY(2)       SumX(2)*SumY(3)   x*SumY(4)        SumX(4)*SumY(2)   SumX(3)*SumY(3)   SumX(2)*SumY(4)   x*SumY(5);
        SumY(3)      x*SumY(3)        SumY(4)          SumX(2)*SumY(3)       x*SumY(4)        SumY(5)          SumX(3)*SumY(3)   SumX(2)*SumY(4)   x*SumY(5)        SumY(6)];
Zo3 = [z z*x y*z SumX(2)*z xySum*z SumY(2)*z SumX(3)*z SumX(2)*y*z x*SumY(2)*z SumY(3)*z];    
co3 = A3/Zo3;
cond3 = abs(max(co3)/min(co3));

%% Normalized Third Order Polynomial Equation
%Design Matrix
An3 = [ n           SumXn(:,1)               SumYn(:,1)               SumXn(2)              xynSum          SumYn(2)          SumXn(3)          SumXn(2)*SumYn(:,1)        SumXn(:,1)*SumYn(2)        SumYn(3);  
        SumXn(:,1)           SumXn(2)          xynSum          SumXn(3)              SumXn(2)*SumYn(:,1)        SumXn(:,1)*SumYn(2)        SumXn(4)          SumXn(3)*SumYn(:,1)        SumXn(2)*SumYn(2)   SumXn(:,1)*SumYn(3) ;
        SumYn(:,1)           xynSum          SumYn(2)          SumXn(2)*SumYn(:,1)            SumXn(:,1)*SumYn(2)        SumYn(3)          SumXn(3)*SumYn(:,1)        SumXn(2)*SumYn(2)   SumXn(:,1)*SumYn(3)        SumYn(4) ;
        SumXn(2)      SumXn(3)          SumXn(2)*SumYn(:,1)        SumXn(4)              SumXn(3)*SumYn(:,1)        SumXn(2)*SumYn(2)   SumXn(5)          SumXn(4)*SumYn(:,1)        SumXn(3)*SumYn(2)   SumXn(2)*SumYn(3) ;
        xynSum      SumXn(2)*SumYn(:,1)        SumXn(:,1)*SumYn(2)        SumXn(3)*SumYn(:,1)            SumXn(2)*SumYn(2)   SumXn(:,1)*SumYn(3)        SumXn(4)*SumYn(:,1)        SumXn(3)*SumYn(2)   SumXn(2)*SumYn(3)   SumXn(:,1)*SumYn(4);
        SumYn(2)      SumXn(:,1)*SumYn(2)        SumYn(3)          SumXn(2)*SumYn(2)       SumXn(:,1)*SumYn(3)        SumYn(4)          SumXn(3)*SumYn(2)   SumXn(2)*SumYn(3)   SumXn(:,1)*SumYn(4)        SumYn(5);
        SumXn(3)      SumXn(4)          SumXn(3)*SumYn(:,1)        SumXn(5)              SumXn(4)*SumYn(:,1)        SumXn(3)*SumYn(2)   SumXn(6)          SumXn(5)*SumYn(:,1)        SumXn(4)*SumYn(2)   SumXn(3)*SumYn(3);
        SumXn(2)*SumYn(:,1)    SumXn(3)*SumYn(2)   SumXn(2)*SumYn(2)   SumXn(4)*SumYn(:,1)            SumXn(3)*SumYn(2)   SumXn(2)*SumYn(3)   SumXn(5)*SumYn(:,1)        SumXn(4)*SumYn(2)   SumXn(3)*SumYn(3)   SumXn(2)*SumYn(4);
        SumXn(:,1)*SumYn(2)    SumXn(2)*SumYn(2)   SumXn(:,1)*SumYn(3)        SumXn(3)*SumYn(2)       SumXn(2)*SumYn(3)   SumXn(:,1)*SumYn(4)        SumXn(4)*SumYn(2)   SumXn(3)*SumYn(3)   SumXn(2)*SumYn(4)   SumXn(:,1)*SumYn(5);
        SumYn(3)      SumXn(:,1)*SumYn(3)        SumYn(4)          SumXn(2)*SumYn(3)       SumXn(:,1)*SumYn(4)        SumYn(5)          SumXn(3)*SumYn(3)   SumXn(2)*SumYn(4)   SumXn(:,1)*SumYn(5)        SumYn(6)];
Zon3 = [SumZn(:,1) SumZn(:,1)*SumXn(:,1) SumYn(:,1)*SumZn(:,1) SumXn(2)*SumZn(:,1) xynSum*SumZn(:,1) SumYn(2)*SumZn(:,1) SumXn(3)*SumZn(:,1) SumXn(2)*SumYn(:,1)*SumZn(:,1) SumXn(:,1)*SumYn(2)*SumZn(:,1) SumYn(3)*SumZn(:,1)];    
con3 = An3/Zon3;
condn3 = abs(max(con3)/min(con3));

%% Centroid Third Order Polynomial Equation
%Design Matrix
Acent3 = [ n           SumXcent(:,1)               SumYcent(:,1)               SumXcent(2)              xycentSum          SumYcent(2)          SumXcent(3)          SumXcent(2)*SumYcent(:,1)        SumXcent(:,1)*SumYcent(2)        SumYcent(3);  
        SumXcent(:,1)           SumXcent(2)          xycentSum          SumXcent(3)              SumXcent(2)*SumYcent(:,1)        SumXcent(:,1)*SumYcent(2)        SumXcent(4)          SumXcent(3)*SumYcent(:,1)        SumXcent(2)*SumYcent(2)   SumXcent(:,1)*SumYcent(3) ;
        SumYcent(:,1)           xycentSum          SumYcent(2)          SumXcent(2)*SumYcent(:,1)            SumXcent(:,1)*SumYcent(2)        SumYcent(3)          SumXcent(3)*SumYcent(:,1)        SumXcent(2)*SumYcent(2)   SumXcent(:,1)*SumYcent(3)        SumYcent(4) ;
        SumXcent(2)      SumXcent(3)          SumXcent(2)*SumYcent(:,1)        SumXcent(4)              SumXcent(3)*SumYcent(:,1)        SumXcent(2)*SumYcent(2)   SumXcent(5)          SumXcent(4)*SumYcent(:,1)        SumXcent(3)*SumYcent(2)   SumXcent(2)*SumYcent(3) ;
        xycentSum      SumXcent(2)*SumYcent(:,1)        SumXcent(:,1)*SumYcent(2)        SumXcent(3)*SumYcent(:,1)            SumXcent(2)*SumYcent(2)   SumXcent(:,1)*SumYcent(3)        SumXcent(4)*SumYcent(:,1)        SumXcent(3)*SumYcent(2)   SumXcent(2)*SumYcent(3)   SumXcent(:,1)*SumYcent(4);
        SumYcent(2)      SumXcent(:,1)*SumYcent(2)        SumYcent(3)          SumXcent(2)*SumYcent(2)       SumXcent(:,1)*SumYcent(3)        SumYcent(4)          SumXcent(3)*SumYcent(2)   SumXcent(2)*SumYcent(3)   SumXcent(:,1)*SumYcent(4)        SumYcent(5);
        SumXcent(3)      SumXcent(4)          SumXcent(3)*SumYcent(:,1)        SumXcent(5)              SumXcent(4)*SumYcent(:,1)        SumXcent(3)*SumYcent(2)   SumXcent(6)          SumXcent(5)*SumYcent(:,1)        SumXcent(4)*SumYcent(2)   SumXcent(3)*SumYcent(3);
        SumXcent(2)*SumYcent(:,1)    SumXcent(3)*SumYcent(2)   SumXcent(2)*SumYcent(2)   SumXcent(4)*SumYcent(:,1)            SumXcent(3)*SumYcent(2)   SumXcent(2)*SumYcent(3)   SumXcent(5)*SumYcent(:,1)        SumXcent(4)*SumYcent(2)   SumXcent(3)*SumYcent(3)   SumXcent(2)*SumYcent(4);
        SumXcent(:,1)*SumYcent(2)    SumXcent(2)*SumYcent(2)   SumXcent(:,1)*SumYcent(3)        SumXcent(3)*SumYcent(2)       SumXcent(2)*SumYcent(3)   SumXcent(:,1)*SumYcent(4)        SumXcent(4)*SumYcent(2)   SumXcent(3)*SumYcent(3)   SumXcent(2)*SumYcent(4)   SumXcent(:,1)*SumYcent(5);
        SumYcent(3)      SumXcent(:,1)*SumYcent(3)        SumYcent(4)          SumXcent(2)*SumYcent(3)       SumXcent(:,1)*SumYcent(4)        SumYcent(5)          SumXcent(3)*SumYcent(3)   SumXcent(2)*SumYcent(4)   SumXcent(:,1)*SumYcent(5)        SumYcent(6)];
Zocent3 = [SumZcent(:,1) SumZcent(:,1)*SumXcent(:,1) SumYcent(:,1)*SumZcent(:,1) SumXcent(2)*SumZcent(:,1) xycentSum*SumZcent(:,1) SumYcent(2)*SumZcent(:,1) SumXcent(3)*SumZcent(:,1) SumXcent(2)*SumYcent(:,1)*SumZcent(:,1) SumXcent(:,1)*SumYcent(2)*SumZcent(:,1) SumYcent(3)*SumZcent(:,1)];    
cocent3 = Acent3/Zocent3;
condcent3 = abs(max(cocent3)/min(cocent3));


%Second order Design Matrix
for i=1:1:92
    A2(i,1) = (1);
    A2(i,2) = xi(i);
    A2(i,3) = yi(i);
    A2(i,4)= xi(i)^2;
    A2(i,5) = xi(i)*yi(i);
    A2(i,6)= yi(i)^2;
end
[X2, L2, R2, apost2] = LnrPrmtrcLSA(zi,A2)

%Third Order Design Matrix
for i=1:1:92
    A3(i,1) = (1);
    A3(i,2) = xi(i);
    A3(i,3) = yi(i);
    A3(i,4)= xi(i)^2;
    A3(i,5) = xi(i)*yi(i);
    A3(i,6)= yi(i)^2;
    A3(i,7)= xi(i)^3;
    A3(i,8)= xi(i)^2*yi(i);
    A3(i,9)= yi(i)^2*xi(i);
    A3(i,10)= yi(i)^3;
end
[X3, L3, R3, apost3] = LnrPrmtrcLSA(zi,A3)
%% Estimate Accuracy
apost2; apost3;

%% Elevations and Residuals


%%Second order
for i=1:1:92
    A2(i,1) = (1);
    A2(i,2) = xi(i);
    A2(i,3) = yi(i);
    A2(i,4)= xi(i)^2;
    A2(i,5) = xi(i)*yi(i);
    A2(i,6)= yi(i)^2;
end
p1 = inv(transpose(A2)*A2);
p2 = transpose(A2);
Xhat2 = A2*p1*p2*zi;


%Second Order
Mean_R2 = mean(R2)
Max_R2 = max(R2)
Min_R2 = min(R2)
STD_R2 = std (R2)
RMSE_R2 = sqrt(mean(((zi-Xhat2).^2)))

%%Third Order
for i=1:1:92
    A3(i,1) = (1);
    A3(i,2) = xi(i);
    A3(i,3) = yi(i);
    A3(i,4)= xi(i)^2;
    A3(i,5) = xi(i)*yi(i);
    A3(i,6)= yi(i)^2;
    A3(i,7)= xi(i)^3;
    A3(i,8)= xi(i)^2*yi(i);
    A3(i,9)= yi(i)^2*xi(i);
    A3(i,10)= yi(i)^3;
end
Xhat3 = A3*inv(transpose(A3)*A3)*transpose(A3)*zi;

%Third Order
Mean_R3 = mean(R3)
Max_R3 = max(R3)
Min_R3 = min(R3)
STD_R3 = std (R3)
RMSE_R3 = sqrt(mean(((zi-Xhat3).^2)))

% Worthiness
%%Second Order
avg_Z = mean(zi);
for i=1:1:92
SST_or2(i) = sum((zi(i)-avg_Z)^2);
SSR_or2(i) = sum ((Xhat2(i)-avg_Z)^2);
SSD_or2(i) = sum((zi(i)-Xhat2(i))^2);
end
SST_o2 = sum(SST_or2)
SSR_o2 = sum(SSR_or2)
SSD_o2 = sum(SSD_or2)
%Goodness-of-Fit
R2 = SSR_o2/SST_o2
nn = 92
m = 6
MST_o2 = SST_o2 / nn
MSR_o2 = SSR_o2 / m
MSD_o2 = SSD_o2 / (n-m)
% Fstat
F2_test = MSR_o2/MSD_o2
%%Third Order
for i=1:1:92
SST_or3(i) = sum((zi(i)-avg_Z)^2);
SSR_or3(i) = sum ((Xhat3(i)-avg_Z)^2);
SSD_or3(i) = sum((zi(i)-Xhat3(i))^2);
end
SST_o3 = sum(SST_or3)
SSR_o3 = sum(SSR_or3)
SSD_o3 = sum(SSD_or3)
% of Goodness-of-Fit
R3 = SSR_o3/SST_o3
n_3 = 92
m_3 = 10
MST_o3 = SST_o3 / n_3
MSR_o3 = SSR_o3 / m_3
MSD_o3 = SSD_o3 / (n_3-m_3)
% Fstat
F3_test = MSR_o3/MSD_o3

%% Elevation 
point = ['A','B','C'];
x_elevLoc = [7000, 6350, 8800];
y_elevLoc = [3500, 5500, 2400];
for i = 1:3
    point(i)
    elev_2(i) = cocent2(1)+(cocent2(2)*x_elevLoc(i)) + (cocent2(3)*y_elevLoc(i)) + (cocent2(4)*(x_elevLoc(i)^2)) + (cocent2(5)*x_elevLoc(i)*y_elevLoc(i)) + (cocent2(6)*(y_elevLoc(i)^2))
    elev_3(i) = cocent3(1)+(cocent3(2)*x_elevLoc(i)) +(cocent3(3)*y_elevLoc(i)) + (cocent3(4)*(x_elevLoc(i)^2)) +(cocent3(5)*x_elevLoc(i)*y_elevLoc(i)) + (cocent3(6)*(y_elevLoc(i)^2))+(cocent3(7)*(x_elevLoc(i)^3))+(cocent3(8)*(x_elevLoc(i)^2)*y_elevLoc(i)) + (cocent3(9)*x_elevLoc(i)*(y_elevLoc(i)^2))+(cocent3(10)*(y_elevLoc(i)^3))
end
   

%% Visualize 1
figure(1)  
PlanimetricPlot = scatter(xi,yi,zi);
hold on 
title('Plainmetric Distribution of the Given Elevation Points') 
xlabel('X (m)') 
ylabel('Y (m)')

%% Visualize 2
figure;
stem3(xi, yi, zi)
grid on
xv = linspace(min(xi), max(xi), 20);
yv = linspace(min(yi), max(yi), 20);
[X,Y] = meshgrid(xv, yv);
Z = griddata(xi,yi,zi,X,Y);
%figure(2)
sc = surfc(X, Y, Z);
grid on
set(gca, 'ZLim',[0 100])
shading interp
zlim([200,300])
title('Original Trend Surface') 
xlabel('X (m)') 
ylabel('Y (m)') 
zlabel('Z (m)') 
%% Visualize 3
figure;
stem3(xi, yi, Xhat2)
grid on
xv = linspace(min(xi), max(xi), 20);
yv = linspace(min(yi), max(yi), 20);
[X,Y] = meshgrid(xv, yv);
Z = griddata(xi,yi,Xhat2,X,Y);
sc = surfc(X, Y, Z);
grid on
set(gca, 'ZLim',[0 100])
shading interp
zlim([200,300])
title('Trend Surface: 2nd order') 
xlabel('X (m)') 
ylabel('Y (m)') 
zlabel('Z (m)') 
%% Visualize 4
figure;
stem3(xi, yi, Xhat3)
grid on
xv = linspace(min(xi), max(xi), 20);
yv = linspace(min(yi), max(yi), 20);
[X,Y] = meshgrid(xv, yv);
Z = griddata(xi,yi,Xhat3,X,Y);
sc = surfc(X, Y, Z);
grid on
set(gca, 'ZLim',[0 100])
zlim([200,300])
title('Trend Surface: 3rd order') 
xlabel('X (m)') 
ylabel('Y (m)') 
zlabel('Z (m)') 

%% Adjustment
function [X, L, R,apost] = LnrPrmtrcLSA(l,A)
    P = diag(ones (length(l),1));
    X =  inv(A'*P*A)*A'*inv(P)*l;
    L =  A*X;
    R =  l - L;
    dof =  length(L) - length(X);
    apost = (R'*P*R)/dof;
end
