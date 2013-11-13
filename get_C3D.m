function [ C3D ] = get_C3D(MatProp)
E1=MatProp(1);
E2=MatProp(2);
G12=MatProp(3);
v12=MatProp(4);
v13=MatProp(7);
v23=MatProp(8);

v21=v12;
v31=v13;
v32=v23;

v=v12*v21+v23*v32+v31*v13+2*v21*v32*v13;

G13=G12;
G23=E2/(2*(1+v23));
E3=E2;

C11=(1-v21*v32)*E1/(1-v);
C12=(v12-v32*v13)*E2/(1-v);
C13=(v13-v12*v23)*E3/(1-v);

C22=(1-v13*v31)*E2/(1-v);
C23=(v23-v21*v13)*E3/(1-v);
C33=(1-v12*v21)*E3/(1-v);

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

