% Composites Project 2B
% Austin Stoker
%

readFromFileOnly=true;
useOldMod=false;

if readFromFileOnly==false
    
% Material Write out


% Write the material properties
            %E1 , E2   ,  G12, v12, alpha1 ,alpha2 , v13  ,v23 ,alpha3    
%Mat_Types=[155E9,12.1E9,4.4E9,.248,-.018E-6,24.3E-6,.248,.458,24.3E-6];
Mat_Types=[19.2E6,1.56E6,.82E6,.24,-.43E-6,13.6E-6,.24,.59,13.6E-6];
MatHeader =    'E1 , E2   ,  G12, v12, alpha1 ,alpha2 , v13  ,v23 ,alpha3'; 
dlmwrite('materials.txt',MatHeader,'');
dlmwrite('materials.txt',Mat_Types,'-append');


% Geometry Write out
r_in=30;

Angles=[0,90,90,0]; %

t=.025;

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
Pin=0;  %pressure inside the cylinder
Pout=0; % Pressure on the outside
dT=0;   %change in Temperature

GivenPx=true;
Px_ex=0;

GivenTx=false;
Tx_Gamx=.001;

LoadsTemp=[Pin,Pout,dT,Px_ex,Tx_Gamx,GivenPx,GivenTx];
LoadHeader='Pin, Pout, dT, Px/ex, Tx/Gamx, GivenPx, GivenTx';
dlmwrite('Loads.txt',LoadHeader,'');
dlmwrite('Loads.txt',LoadsTemp,'-append');

MaterialsFile='materials.txt';
LayerGeoFile='LayerGeometery.dat';
LoadFile='Loads.txt';
PipeFile='PipeDimensions.dat';

end



%% read in
% Geometery read in
Geometery =  dlmread(LayerGeoFile,',',1,0);
Angles= Geometery(1,:);
Materials=Geometery(2,:);
Thicknesses=Geometery(3,:);
NL=size(Angles,2);


% Pipe Geometery Read in
r_in= dlmread(PipeFile,',',1,0);
rsteps=2; %points in mid layer.
r=zeros(NL+1,1);
r(1)=r_in;
rList=zeros((rsteps+1)*NL+1,1);
rLayerList=ones((rsteps+1)*NL+1,1);
rList(1)=r_in;
for i=1:NL 
    for j=1:rsteps+1
        k=(i-1)*(rsteps+1)+j+1;
        if j==rsteps+1
            rList(k)=rList(k-1);
        else
            rList(k)=rList(k-1)+Thicknesses(i)/rsteps;
        end
            rLayerList(k-1)=i;
    end
end


for i=2:NL+1
    r(i)=r(i-1)+Thicknesses(i-1);
end

% Material Read in
Mat_Types = dlmread(MaterialsFile,',',1,0);

% Loads Read in
LoadsTemp=dlmread(LoadFile,',',1,0);
Loads=LoadsTemp(1:5);
GivenPx=LoadsTemp(6);
GivenTx=LoadsTemp(7);
Pin=Loads(1);
Pout=Loads(2);
dT=Loads(3);
if(GivenPx)
    Px=Loads(4);
    Epsx=0;
else
    Epsx=Loads(4);
    Px=0;
end

if(GivenTx)
    Tx=Loads(5);
    Gamx=0;
else
    Gamx=Loads(5);
    Tx=0;
end



%% Calculations
%Get the C_bar matrices
C_bar=zeros(6,6,NL);
for i=1:NL
    C_bar(:,:,i)=get_Cbar3D(Mat_Types(Materials(i),:),Angles(i));
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

al_on=zeros(3,NL);
al_off=zeros(4,NL);

for i=1:NL
    al_1(i)=Mat_Types(Materials(i),5);
    al_2(i)=Mat_Types(Materials(i),6);
    al_r(i)=Mat_Types(Materials(i),9);
    th=Angles(i);
    al_x(i)=al_1(i)*cosd(th)^2+al_2(i)*sind(th)^2;
    al_th(i)=al_1(i)*sind(th)^2+al_2(i)*cosd(th)^2;
    al_x_th(i)=2*sind(th)*cosd(th)*(al_1(i)-al_2(i));
    al_on(:,i)=[al_1(i);al_2(i);al_r(i)];
    al_off(:,i)=[al_x(i);al_th(i);al_r(i);al_x_th(i)];

    if(Angles(i)==0)
        Gam(i)=0;
        Omega(i)=0;
        Psi(i)=0;
        Sigma(i)=0;
    else
        Gam(i)=(C_bar(1,2,i)-C_bar(1,3,i))/(C_bar(3,3,i)-C_bar(2,2,i));
        Omega(i)=(C_bar(2,6,i)-2*C_bar(3,6,i))/(4*C_bar(3,3,i)-C_bar(2,2,i));
        Sigma(i)=(C_bar(1,3,i)-C_bar(1,2,i))*al_x(i)+...
            (C_bar(2,3,i)-C_bar(2,2,i))*al_th(i) +...
            (C_bar(3,3,i)-C_bar(3,2,i))*al_r(i)+...
            (C_bar(6,3,i)-C_bar(6,2,i))*al_x_th(i);
        Psi(i)=Sigma(i)/(C_bar(3,3,i)-C_bar(2,2,i));
    end
    lambda(i)=sqrt(C_bar(2,2,i)/C_bar(3,3,i));
end

%% Get the K matrix

K=zeros(2*NL+2,2*NL+2);
F=zeros(2*NL+2,1);

% first row k

K(1,1)=(C_bar(2,3,1)+lambda(1)*C_bar(3,3,1))*r(1)^(lambda(1)-1);
K(1,2)=(C_bar(2,3,1)-lambda(1)*C_bar(3,3,1))*r(1)^(-lambda(1)-1);

%if(GivenPx)
K(1,2*NL+1)=C_bar(1,3,1)+(C_bar(2,3,1)+C_bar(3,3,1))*Gam(1);

%if(GivenTx)
K(1,2*NL+2)=(C_bar(3,6,1)+(C_bar(2,3,1)+2*C_bar(3,3,1))*Omega(1))*r(1);


bob=0; %TODO fix the c*alpha sum
F(1)=-Pin-((C_bar(2,3,1)+C_bar(3,3,1))*Psi(1)-bob)*dT;

 
for k=1:NL-1
%even rows 2k
    K(2*k,2*k-1)=r(k+1)^lambda(k);
    K(2*k,2*k)=r(k+1)^-lambda(k);
    K(2*k,2*k+1)=-r(k+1)^lambda(k+1);
    K(2*k,2*k+2)=-r(k+1)^-lambda(k+1);
    
    
    F(2*k)=r(k+1)*(Psi(k+1)-Psi(k))*dT;
    
    %if(GivenPx)
    K(2*k,2*NL+1)=(Gam(k)-Gam(k+1))*r(k+1);

    
    %if(GivenTx)
    K(2*k,2*NL+2)=(Omega(k)-Omega(k+1))*r(k+1)^2;

    
% odd rows 2k+1
    K(2*k+1,2*k-1)=r(k+1)^(lambda(k)-1)*(C_bar(2,3,k)+lambda(k)*C_bar(3,3,k));
    K(2*k+1,2*k)=r(k+1)^(-lambda(k)-1)*(C_bar(2,3,k)-lambda(k)*C_bar(3,3,k));
    
    K(2*k+1,2*k+1)=-r(k+1)^(lambda(k+1)-1)*(C_bar(2,3,k+1)+lambda(k+1)*C_bar(3,3,k+1));
    K(2*k+1,2*k+2)=-r(k+1)^(-lambda(k+1)-1)*(C_bar(2,3,k+1)-lambda(k+1)*C_bar(3,3,k+1));
    
    bob(k)=0;
    bob(k+1)=0;
    F(2*k+1)=((C_bar(2,3,k+1)+C_bar(3,3,k+1))*Psi(k+1)-bob(k+1)-((C_bar(2,3,k)+C_bar(3,3,k))*Psi(k)-bob(k)))*dT;
    
    
    %if(GivenPx)
    K(2*k+1,2*NL+1)=C_bar(1,3,k)+(C_bar(2,3,k)+C_bar(3,3,k))*Gam(k)-(C_bar(1,3,k+1)+(C_bar(2,3,k+1)+C_bar(3,3,k+1))*Gam(k+1));
    
    %if(GivenTx)
    K(2*k+1,2*NL+2)=r(k+1)*(C_bar(3,6,k)+(C_bar(2,3,k)+2*C_bar(3,3,k))*Omega(k)-(C_bar(3,6,k+1)+(C_bar(2,3,k+1)+2*C_bar(3,3,k+1))*Omega(k+1)));
end

for k=1:NL
    
    % 2N row
    K(2*NL,2*NL-1)=(C_bar(2,3,NL)+lambda(NL)*C_bar(3,3,NL))*r(NL+1)^(lambda(NL)-1);
    K(2*NL,2*NL)=(C_bar(2,3,NL)-lambda(NL)*C_bar(3,3,NL))*r(NL+1)^(-lambda(NL)-1);
    K(2*NL,2*NL+1)=C_bar(1,3,NL)+(C_bar(2,3,NL)+C_bar(3,3,NL))*Gam(NL);
    K(2*NL,2*NL+2)=(C_bar(3,6,NL)+(C_bar(2,3,NL)+2*C_bar(3,3,NL))*Omega(NL))*r(NL+1);
    
    F(2*NL)=-Pout-((C_bar(3,6,NL)+C_bar(3,3,NL))*Psi(NL)-bob(k))*dT;
    
    % 2N+1 row
    if(abs(lambda(k)+1)>1E-10)
        K(2*NL+1,2*k-1)=2*pi*(C_bar(1,2,k)+lambda(k)*C_bar(1,3,k))*(r(k+1)^(lambda(k)+1)-r(k)^(lambda(k)+1))/(lambda(k)+1);
    else
        K(2*NL+1,2*k-1)=0;
    end
    
    if(abs(1-lambda(k))>1E-10)
        K(2*NL+1,2*k)=2*pi*(C_bar(1,2,k)-lambda(k)*C_bar(1,3,k))*(r(k+1)^(1-lambda(k))-r(k)^(1-lambda(k)))/(1-lambda(k));
    else
        K(2*NL+1,2*k)=0;
    end
    
    
    % 2N+2 row
    if(abs(lambda(k)+2)>1E-10)
        K(2*NL+2,2*k-1)=2*pi*(C_bar(2,6,k)+lambda(k)*C_bar(3,6,k))*(r(k+1)^(lambda(k)+2)-r(k)^(lambda(k)+2))/(lambda(k)+2);
    else
        K(2*NL+2,2*k-1)=0;
    end
    
    if(abs(2-lambda(k))>1E-10)
        K(2*NL+2,2*k)=2*pi*(C_bar(2,6,k)-lambda(k)*C_bar(3,6,k))*(r(k+1)^(2-lambda(k))-r(k)^(2-lambda(k)))/(2-lambda(k));
    else
        K(2*NL+2,2*k)=0;
    end
end

%2N+1 and 2n+2 ends
K(2*NL+1,2*NL+1)=0;
K(2*NL+1,2*NL+2)=0;
if(GivenPx)
    F(2*NL+1)=Px;
else
    F(2*NL+1)=0;
end

if(GivenTx)
    F(2*NL+2)=Tx;
else
    F(2*NL+2)=0;
end
    
for k=1:NL
    K(2*NL+1,2*NL+1)=K(2*NL+1,2*NL+1)+...
        pi*(C_bar(1,1,k)+(C_bar(1,3,k)+C_bar(1,2,k))*Gam(k))*(r(k+1)^2-r(k)^2);

    K(2*NL+1,2*NL+2)=K(2*NL+1,2*NL+2)+...
        2*pi/3*(C_bar(1,6,k)+(C_bar(1,2,k)+2*C_bar(1,3,k))*Omega(k))*(r(k+1)^3-r(k)^3);
    
    F(2*NL+1)=F(2*NL+1)-...
        dT*pi*((C_bar(1,2,k)+C_bar(1,3,k))*Psi(k)-bob(k))*(r(k+1)^2-r(k)^2);


    
    K(2*NL+2,2*NL+1)=K(2*NL+2,2*NL+1)+...
        2*pi/3*(C_bar(1,6,k)+(C_bar(2,6,k)+C_bar(3,6,k))*Gam(k))*(r(k+1)^3-r(k)^3);

    K(2*NL+2,2*NL+2)=K(2*NL+2,2*NL+2)+...
        pi/2*(C_bar(6,6,k)+(C_bar(2,6,k)+2*C_bar(3,6,k))*Omega(k))*(r(k+1)^4-r(k)^4);

    F(2*NL+2)=F(2*NL+2)-...
        dT*2*pi/3*((C_bar(2,6,k)+C_bar(3,6,k))*Psi(k)-bob(k))*(r(k+1)^3-r(k)^3);
        
end




% Modify K and F
if useOldMod
    K_mod=K;
    F_mod=F;
    if ~GivenPx
        for k=1:NL-1
            K_mod(1,2*NL+1)=0;
            F_mod(1)=F_mod(1)-(C_bar(1,3,1)+(C_bar(2,3,1)+C_bar(3,3,1))*Gam(1))*Epsx;

            K_mod(2*k,2*NL+1)=0;
            F_mod(2*k)=F_mod(2*k)-(Gam(k)-Gam(k+1))*r(k+1)*Epsx;

            K_mod(2*k+1,2*NL+1)=0;
            F_mod(2*k+1)=F_mod(2*k+1)-(C_bar(1,3,k)+(C_bar(2,3,k)+C_bar(3,3,k))*Gam(k)-(C_bar(1,3,k+1)+(C_bar(2,3,k+1)+C_bar(3,3,k+1))*Gam(k+1)))*Epsx;

            K_mod(2*NL,2*NL+1)=0;
            F_mod(2*NL)=F_mod(2*NL)-(C_bar(1,3,NL)+(C_bar(2,3,NL)+C_bar(3,3,NL))*Gam(NL))*Epsx;

            K_mod(2*NL+1,2*NL+1)=-1;

            K_mod(2*NL+2,2*NL+1)=0;
        end
        for k=1:NL
            F_mod(2*NL+1)=F_mod(2*NL+1)-...
                    pi*(Epsx*(C_bar(1,1,k)+(C_bar(1,3,k)+C_bar(1,2,k))*Gam(k))+dT*((C_bar(1,2,k)+C_bar(1,3,k))*Psi(k)-bob(k)))*(r(k+1)^2-r(k)^2);
            F_mod(2*NL+1)=F_mod(2*NL+1)-...
                    Epsx*2*pi/3*(C_bar(1,6,k)+(C_bar(2,6,k)+C_bar(3,6,k))*Gam(k))*(r(k+1)^3-r(k)^3);  

        end

    end

    if ~GivenTx
        for k=1:NL-1
            K_mod(1,2*NL+2)=0;
            F_mod(1)=F_mod(1)-(C_bar(3,6,1)+(C_bar(2,3,1)+2*C_bar(3,3,1))*Omega(1))*Gamx*r(1);

            K_mod(2*k,2*NL+2)=0;
            F_mod(2*k)=F_mod(2*k)-(Omega(k)-Omega(k+1))*r(k+1)^2*Gamx;

            K_mod(2*k+1,2*NL+2)=0;
            F_mod(2*k+1)=F_mod(2*k+1)-r(k+1)*(C_bar(3,6,k)+(C_bar(2,3,k)+2*C_bar(3,3,k))*Omega(k)-(C_bar(3,6,k+1)+(C_bar(2,3,k+1)+2*C_bar(3,3,k+1))*Omega(k+1)))*Gamx;

            K_mod(2*NL,2*NL+2)=0;
            F_mod(2*NL)=F_mod(2*NL)-(C_bar(3,6,NL)+(C_bar(2,3,NL)+2*C_bar(3,3,NL))*Omega(NL))*Gamx*r(NL+1);

            K_mod(2*NL+1,2*NL+2)=0;

            K_mod(2*NL+2,2*NL+2)=-1;
        end
        for k=1:NL
            F_mod(2*NL+1)=F_mod(2*NL+1)-Gamx*2*pi/3*(C_bar(1,6,k)+(C_bar(1,2,k)+2*C_bar(1,3,k))*Omega(k))*(r(k+1)^3-r(k)^3);
            F_mod(2*NL+2)=F_mod(2*NL+2)-...
                Gamx*pi/2*(C_bar(6,6,k)+(C_bar(2,6,k)+2*C_bar(3,6,k))*Omega(k))*(r(k+1)^4-r(k)^4)-...
                dT*2*pi/3*((C_bar(2,6,k)+C_bar(3,6,k))*Psi(k)-bob(k))*(r(k+1)^3-r(k)^3);

        end
    end
else
    K_mod=K;
    F_mod=F;
    if ~GivenPx
        F_mod=F_mod-K_mod(:,2*NL+1)*Epsx;
        K_mod(:,2*NL+1)=0;
        K_mod(2*NL+1,2*NL+1)=-1;
    end
    
    if ~GivenTx
        F_mod=F_mod-K_mod(:,2*NL+2)*Gamx;
        K_mod(:,2*NL+2)=0;
        K_mod(2*NL+2,2*NL+2)=-1;
    end
end

%% Use K and F mod to find Strains and Stresses

R=K_mod^-1*F_mod;

if(GivenPx)
    Epsx=R(2*NL+1);
else
    Px=R(2*NL+1);
end

if(GivenTx)
    Gamx=R(2*NL+2);
else
    Tx=R(2*NL+2);
end

A1=zeros(NL,1);
A2=A1;
C_reduced=zeros(4,4,NL);
strains=zeros(4,(rsteps+1)*NL);
w=zeros((rsteps+1)*NL,1);
stresses=strains;

for i=1:NL
    A1(i)=R(2*i-1);
    A2(i)=R(2*i);
    C_reduced(1:3,1:3,i)=C_bar(1:3,1:3,i);
    C_reduced(1:3,4,i)=C_bar(1:3,6,i);
    C_reduced(4,1:3,i)=C_bar(6,1:3,i);
    C_reduced(4,4,i)=C_bar(6,6,i);
    
end

for i=1:(rsteps+1)*NL
    k=rLayerList(i);
    r_=rList(i);
    strains(1,i)=Epsx;
    strains(2,i)=A1(k)*r_^(lambda(k)-1)+A2(k)*r_^(-lambda(k)-1)+Gam(k)*Epsx+Omega(k)*Gamx*r_+Psi(k)*dT; %Eps_Th
    strains(3,i)=A1(k)*lambda(k)*r_^(lambda(k)-1)-A2(k)*lambda(k)*r_^(-lambda(k)-1)+Gam(k)*Epsx+2*Omega(k)*Gamx*r_+Psi(k)*dT; %Eps_r
    strains(4,i)=Gamx*r_;%gamma_x_th
    w(i)=A1(k)*r_^(lambda(k))+A2(k)*r_^(-lambda(k))+Gam(k)*Epsx*r_+Omega(k)*Gamx*r_^2+Psi(k)*r_*dT; %w

    stresses(:,i)=C_reduced(:,:,k)*(strains(:,i)-dT*al_off(:,k));
end


