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
            %E1 , E2   ,  G12, v12, alpha1 ,alpha2 ,    ,v23 ,alpha3    
Mat_Types=[155E9,12.1E9,4.4E9,.248,-.018E-6,24.3E-6,.248,.458,24.3E-6];
MatHeader =    'E1 , E2   ,  G12, v12, alpha1 ,alpha2 ,    ,v23 ,alpha3'; 
dlmwrite('materials.txt',MatHeader,'');
dlmwrite('materials.txt',Mat_Types,'-append');


% Geometry Write out
r_in=.05;
Angles=[0,90];%

t=.00015;

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
numLayers=size(Angles,2);
rsteps=20;
rList=-t*numLayers/2+t/rsteps:t/rsteps:t*numLayers/2-t/rsteps;

% Pipe Geometery Read in

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
C_bar=zeros(4,4,numLayers);
for i=1:numLayers
    C3D=get_Cbar3D(Mat_Types(Materials(i),:),Angles(i));
    C_bar(1:3,1:3,i)=C3D(1:3,1:3);
    C_bar(1:3,4,i)=C3D(1:3,6);
    C_bar(4,1:3,i)=C3D(6,1:3);
    C_bar(4,4,i)=C3D(6,6);
end


    




