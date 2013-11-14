% Composites Project 2B
% Austin Stoker
%

close all
clearvars
clc

readFromFile=false;

if readFromFile==false
    
% Material Write out


% Write the material properties
            %E1 , E2   ,  G12, v12, alpha1 ,alpha2 , v13  ,v23 ,alpha3    
%Mat_Types=[155E9,12.1E9,4.4E9,.248,-.018E-6,24.3E-6,.248,.458,24.3E-6];
Mat_Types=[19.2E6,1.56E6,.82E6,.24,-.43E-6,13.6E-6,.24,.59,13.6E-6];
MatHeader =    'E1 , E2   ,  G12, v12, alpha1 ,alpha2 , v13  ,v23 ,alpha3'; 
dlmwrite('materials.txt',MatHeader,'');
dlmwrite('materials.txt',Mat_Types,'-append');


% Geometry Write out
r_in=1.5;

Angles=[26]; %

t=.005;

Materials=ones(1,numel(Angles)); % [1,1,1,1,1,1,1,1];
Thicknesses= t*ones(1,numel(Angles)); %[t,t,t,t,t,t,t,t];

Geometery(1,:)=Angles;
Geometery(2,:)=Materials;
Geometery(3,:)=Thicknesses;

GeoHeader ='The rows are: Angles (degrees), Material numbers, thicknesses';
dlmwrite('LayerGeometery.dat',GeoHeader,'');
dlmwrite('LayerGeometery.dat',Geometery,'-append');

PipeDimHeader='Inner Radius';
dlmwrite('PipeDimensions.dat',PipeDimHeader,'');
dlmwrite('PipeDimensions.dat',r_in,'-append');


% Loads Write out
Pin=1;  %pressure inside the cylinder
Pout=2; % Pressure on the outside
dT=0;   %change in Temperature

GivenPx=true;
Px_ex=1;

GivenTx=true;
Tx_Gamx=1;

LoadsTemp=[Pin,Pout,dT,Px_ex,Tx_Gamx,GivenPx,GivenTx];
LoadHeader='Pin, Pout, dT, Px/ex, Tx/Gamx, GivenPx, GivenTx';
dlmwrite('Loads.txt',LoadHeader,'');
dlmwrite('Loads.txt',LoadsTemp,'-append');

end

%% read in
% Geometery read in
Geometery =  dlmread('LayerGeometery.dat',',',1,0);
Angles= Geometery(1,:);
Materials=Geometery(2,:);
Thicknesses=Geometery(3,:);
NL=size(Angles,2);
rsteps=20;
rList=-t*NL/2+t/rsteps:t/rsteps:t*NL/2-t/rsteps;

% Pipe Geometery Read in
r_in= dlmread('PipeDimensions.dat',',',1,0);
r(1)=r_in;
for i=2:NL+1
    r(i)=r(i-1)+Thicknesses(i-1);
end

% Material Read in
Mat_Types = dlmread('materials.txt',',',1,0);

% Loads Read in
LoadsTemp=dlmread('Loads.txt',',',1,0);
Loads=LoadsTemp(1:5);
GivenPx=LoadsTemp(6);
GivenTx=LoadsTemp(7);
Pin=Loads(1);
Pout=Loads(2);
dT=Loads(3);
if(GivenPx)
    Px=Loads(4);
else
    Epsx=Loads(4);
end

if(GivenTx)
    Tx=Loads(5);
else
    Gamx=Loads(5);
end



%% Calculations
%Get the C_bar matrices
C_bar=zeros(6,6,NL);
for i=1:NL
    C_bar=get_Cbar3D(Mat_Types(Materials(i),:),Angles(i));
end

% Get Gamma, Omega, Psi etc

%%
Gam=zeros(1,NL);
Omega=zeros(1,NL);
Psi=zeros(1,NL);
Sigma=zeros(1,NL);
lambda=zeros(1,NL);
al_1=zeros(1,NL);
al_2=zeros(1,NL);
al_x=zeros(1,NL);
al_th=zeros(1,NL);
al_r=zeros(1,NL);
al_1=zeros(1,NL);
al_x_th=zeros(1,NL);

for i=1:NL
    al_1(i)=Mat_Types(Materials(i),5);
    al_2(i)=Mat_Types(Materials(i),6);
    al_r(i)=Mat_Types(Materials(i),9);
    th(i)=Angles(i);
    al_x(i)=al_1(i)*cosd(th(i))^2+al_2(i)*sind(th(i))^2;
    al_th(i)=al_1(i)*sind(th(i))^2+al_2(i)*cosd(th(i))^2;
    al_x_th(i)=2*sind(th(i))*cosd(th(i))*(al_1(i)-al_2(i));
    
    if(Angles(i)==0||Angles(i)==90)
        Gam(i)=0;
        Omega(i)=0;
        Psi(i)=0;
        Sigma(i)=0;
        lambda(i)=0;
    else
        Gam(i)=(C_bar(1,2,i)-C_bar(1,3,i))/(C_bar(3,3,i)-C_bar(2,2,i));
        Omega(i)=(C_bar(2,6,i)-2*C_bar(3,6,i))/(4*C_bar(3,3,i)-C_bar(2,2,i));
        Sigma(i)=(C_bar(1,3,i)-C_bar(1,2,i))*al_x(i)+...
            (C_bar(2,3,i)-C_bar(2,2,i))*al_th(i) +...
            (C_bar(3,3,i)-C_bar(3,2,i))*al_r(i)+...
            (C_bar(6,3,i)-C_bar(6,2,i))*al_x_th(i);
        Psi(i)=Sigma(i)/(C_bar(3,3,i)-C_bar(2,2,i));
        lambda(i)=sqrt(C_bar(2,2,i)/C_bar(3,3,i));
    end
end

%% Get the K matrix

K=zeros(2*NL+2,2*NL+2);
F=zeros(2*NL+2,1);

K(1,1)=(C_bar(2,3,1)+lambda(1)*C_bar(3,3,1))*r(1)^(lambda(1)-1);
K(1,2)=(C_bar(2,3,1)-lambda(1)*C_bar(3,3,1))*r(1)^(-lambda(1)-1);

if(GivenPx)
    K(1,2*NL+1)=C_bar(1,3,1)+(C_bar(2,3,1)+C_bar(3,3,1))*Gam(1);
else
    K(1,2*NL+1)=0;
end

if(GivenTx)
    K(1,2*NL+2)=(C_bar(3,6,1)+(C_bar(2,3,1)+2*C_bar(3,3,1))*Omega(1))*r(1);
else
    K(1,2*NL+2)=0;
end

bob=0; %TODO fix the c*alpha sum
F(1)=-Pin-((C_bar(2,3,1)+C_bar(3,3,1))*Psi(1)-bob)*dT;
if(~GivenPX)
    F(1)=F(1)-(C_bar(1,3,1)+(C_bar(2,3,1)+C_bar(3,3,1))*Gam(1))*Epsx;
end
if(~GivenTx)
    F(1)=F(1)-(C_bar(3,6,1)+(C_bar(2,3,1)+2*C_bar(3,3,1))*Omega(1))*Gamx*r(1);
end

 
for k=1:NL-1
%even rows
    K(2*k,2*k-1)=r(k+1)^lambda(k);
    K(2*k,2*k)=r(k+1)^-lambda(k);
    K(2*k,2*k+1)=-r(k+1)^lambda(k+1);
    K(2*k,2*k+2)=-r(k+1)^-lambda(k+1);
    
    
    F(2*k)=r(k+1)*(Psi(k+1)-Psi(k))*dT;
    
    if(GivenPx)
        K(2*k,2*NL+1)=(Gam(k)-Gam(k+1))*r(k+1);
    else
        K(2*k,2*NL+1)=0;
        F(2*k)=F(2*k)-(Gam(k)-Gam(k+1))*r(k+1)*Epsx;
    end
    
    if(GivenTx)
        K(2*k,2*NL+2)=(Omega(k)-Omega(k+1))*r(k+1)^2;
    else
        K(2*k,2*NL+2)=0;
        F(2*k)=F(2*k)-(Omega(k)-Omega(k+1))*r(k+1)^2*Gamx;
    end
    
% odd rows
    K(2*k+1,2*k-1)=r(k+1)^(lambda(k)-1)*C_bar(2,3,k)+lambda(k)*C_bar(3,3,k));
    K(2*k+1,2*k)=r(k+1)^(-lambda(k)-1)*C_bar(2,3,k)-lambda(k)*C_bar(3,3,k));
    
    F(2*k+1)=((C_bar(2,3,k+1)+C_bar(3,3,k+1))*Phi(k+1)-bob(k+1)-((C_bar(2,3,k)+C_bar(3,3,k))*Phi(k)-bob(k)))*dT;
    
    
    if(GivenPx)
    K(2*k+1,2*NL+1)=C_bar(1,3,k)+(C_bar(2,3,k)+C_bar(3,3,k))*Gam(k)-(C_bar(1,3,k+1)+(C_bar(2,3,k+1)+C_bar(3,3,k+1))*Gam(k+1));
    else
    K(2*k+1,2*NL+1)=0;
    F(2*k+1)=F(2*k+1)-(C_bar(1,3,k)+(C_bar(2,3,k)+C_bar(3,3,k))*Gam(k)-(C_bar(1,3,k+1)+(C_bar(2,3,k+1)+C_bar(3,3,k+1))*Gam(k+1)))*Epsx;
    
    end
    
    if(GivenTx)
    K(2*k+1,2*NL+2)=r(k+1)*(C_bar(3,6,k)+(C_bar(2,3,k)+2*C_bar(3,3,k))*Omega(k)-(C_bar(3,6,k+1)+(C_bar(2,3,k+1)+2*C_bar(3,3,k+1))*Omega(k+1)));
    else
    K(2*k+1,2*NL+2)=0;
    F(2*k+1)=F(2*k+1)-r(k+1)*(C_bar(3,6,k)+(C_bar(2,3,k)+2*C_bar(3,3,k))*Omega(k)-(C_bar(3,6,k+1)+(C_bar(2,3,k+1)+2*C_bar(3,3,k+1))*Omega(k+1)))*Gamx;
    end
    
end


K(2*NL,2*NL-1)=(C_bar(2,3,NL)+lambda(NL)*C_bar(3,3,NL))*r(NL+1)^(lambda(NL)-1);
K(2*NL,2*NL)=(C_bar(2,3,NL)-lambda(NL)*C_bar(3,3,NL))*r(NL+1)^(-lambda(NL)-1);



