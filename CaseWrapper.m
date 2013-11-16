close all
clearvars
clc

myCase=4;


switch myCase
    case 1
        Results=zeros(12,6);
        MaterialsFile='VerificationMaterials.txt';
        PipeFile='VerificationPipe.dat';
        for i1=1:4
            LoadFile=strcat('VerificationLoads',num2str(i1),'.txt');
            for j1=1:6
                LayerGeoFile=strcat('VerificationLayerGeo',num2str(j1),'.dat');
                Project_2B_Laminated_Cylinder_Axisymm
                if i1==1
                        Results(i1*3-2,j1)=Gamx;
                        Results(i1*3-1,j1)=w(1);
                        Results(i1*3,j1)=Px;
                elseif i1==2
                        Results(i1*3-2,j1)=Epsx;
                        Results(i1*3-1,j1)=w(1);
                        Results(i1*3,j1)=Tx;
                elseif i1== 3
                        Results(i1*3-2,j1)=Gamx;
                        Results(i1*3-1,j1)=w(1);
                        Results(i1*3,j1)=Epsx;
                elseif i1== 4
                        Results(i1*3-2,j1)=Epsx;
                        Results(i1*3-1,j1)=Gamx;
                        Results(i1*3,j1)=w(1);
                end
                        
            end
        end
        Results=Results
    case 2
        ResultsB=zeros(3,5);
        MaterialsFile='VerificationMaterials.txt';
        PipeFile='VerificationPipe.dat';
        LoadFile=strcat('VerificationLoadsB.txt');
        for j1=1:3
            LayerGeoFile=strcat('VerificationLayerGeo',num2str(j1+3),'.dat');
            Project_2B_Laminated_Cylinder_Axisymm
            ResultsB(j1,1)=Px;
            ResultsB(j1,2)=Tx;
            ResultsB(j1,3)=Epsx;
            ResultsB(j1,4)=Gamx;
            ResultsB(j1,5)=w(1);

        end
        ResultsB=ResultsB
        
    case 3 % Strains and Displacements
        NP=80;
        Ang=-90:180/(NP-1):90;
        res=zeros(NP,3);
        MaterialsFile='StrainDispMaterials.txt';
        PipeFile='StrainDispPipe.dat';
        LoadFile=strcat('StrainDispLoads.txt');
        
        t=.025;
        Materials=ones(1,4); % [1,1,1,1,1,1,1,1];
        Thicknesses= t*ones(1,4); %[t,t,t,t,t,t,t,t];
        Geometery(2,:)=Materials;
        Geometery(3,:)=Thicknesses;
        for j1=1:NP
            Angles=[Ang(j1),-Ang(j1),-Ang(j1),Ang(j1)]; %

            Geometery(1,:)=Angles;
            GeoHeader ='The rows are: Angles (degrees), Material numbers, thicknesses';
            dlmwrite('StrainDisplGeo.dat',GeoHeader,'');
            dlmwrite('StrainDisplGeo.dat',Geometery,'-append');
            LayerGeoFile=strcat('StrainDisplGeo.dat');
            
            Project_2B_Laminated_Cylinder_Axisymm
            res(j1,1)=w(1);
            res(j1,2)=Epsx;
            res(j1,3)=Gamx;
        end 
        figure(1)
        plot(Ang,res(:,1))
        ylabel('w_i_n')
        xlabel('\theta (deg)')
        title('w_i_n of a [\theta,-\theta,-\theta,\theta]')
        
        figure(2)
        plot(Ang,res(:,2))
        ylabel('\epsilon_x')
        xlabel('\theta (deg)')
        title('\epsilon_x of a [\theta,-\theta,-\theta,\theta]')
        
        figure(3)
        plot(Ang,res(:,3))
        ylabel('\gamma_x')
        xlabel('\theta (deg)')
        title('\gamma_x of a [\theta,-\theta,-\theta,\theta]')
    case 4
        use_many_r_steps=true;
        rsteps=8;
        
        MaterialsFile='StressMaterials.txt';
        PipeFile='StressPipe.dat';
        LoadFile=strcat('StressLoads.txt');
        LayerGeoFile=strcat('StressGeo.dat');

        Project_2B_Laminated_Cylinder_Axisymm
        
        NR=numel(rList);
        
        figure(1)
        plot(stresses(1,1:NR-1)',rList(1:end-1))
        ylabel('r')
        xlabel('\sigma_x')
        title('\sigma_x of a [26,-26,-26,26]')
        
        figure(2)
        plot(stresses(2,1:NR-1)',rList(1:end-1))
        ylabel('r')
        xlabel('\sigma_\theta')
        title('\sigma_\theta of a [26,-26,-26,26]')
        
        figure(3)
        plot(stresses(3,1:NR-1)',rList(1:end-1))
        ylabel('r')
        xlabel('\sigma_r')
        title('\sigma_r of a [26,-26,-26,26]')
        
        figure(4)
        plot(stresses(4,1:NR-1)',rList(1:end-1))
        ylabel('r')
        xlabel('\tau_x_\theta')
        title('\tau_x_\theta of a [26,-26,-26,26]')
        
        use_many_r_steps=false;
    case 5
        MaterialsFile='SmearedMaterials.txt';
        PipeFile='SmearedPipe.dat';
        LoadFile=strcat('SmearedLoads.txt');
        LayerGeoFile=strcat('SmearedGeo.dat');

        Project_2B_Laminated_Cylinder_Axisymm
        
        
end
