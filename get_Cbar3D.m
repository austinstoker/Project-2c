function [ Cbar3D ] = get_Cbar3D( MatProp, Angle)
Trad = (Angle*pi)/180;
m = cos(Trad);
n = sin(Trad);

C=get_C3D(MatProp);

T1 = [m^2 n^2 0 0 0 2*m*n; n^2 m^2 0 0 0 -2*n*m; 0 0 1 0 0 0; 0 0 0 m -n 0; 0 0 0 n m 0; -m*n m*n 0 0 0 m^2-n^2];
T2 = [m^2 n^2 0 0 0 m*n; n^2 m^2 0 0 0 -n*m; 0 0 1 0 0 0; 0 0 0 m -n 0; 0 0 0 n m 0; -2*m*n 2*m*n 0 0 0 m^2-n^2];

CbTemp = T1^-1*C*T2;
Cbar3D = CbTemp;
end

