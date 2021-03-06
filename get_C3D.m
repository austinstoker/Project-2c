function [ C3D ] = get_C3D(MatProp)
E1=MatProp(1);
E2=MatProp(2);
G12=MatProp(3);
v12=MatProp(4);
v13=MatProp(7);
v23=MatProp(8);

E3=E2;
v21=v12*E2/E1;
v31=v21;
v32=v23*E3/E2;
nu=[0,v12,v13;v21,0,v23;v31,v32,0];

v=v12*v21+v23*v32+v31*v13+2*v21*v32*v13;

G13=G12;
G23=E2/(2*(1+v23));
G32=G23;


den=(1-v);
C11=(1-v23*v32)*E1/den;
C12=(v12+v32*v13)*E2/den;
C13=(v13+v12*v23)*E3/den;

C22=(1-v13*v31)*E2/den;
C23=(v23+v21*v13)*E3/den;
C33=(1-v12*v21)*E3/den;

C44=G23;
C55=G13;
C66=G12;

C3D=[C11 C12 C13  0  0  0
     C12 C22 C23  0  0  0
     C13 C23 C33  0  0  0
     0    0   0  C44 0  0
     0    0   0   0 C55 0
     0    0   0   0  0 C66];
end

