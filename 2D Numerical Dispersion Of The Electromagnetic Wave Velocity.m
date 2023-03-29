clc; clear; clf;

for num=1:2                   
    
 if (num==1)                    %first part with dx=dy
  lo=1;                         %wave length
  fi=0:pi/100:pi/2;             %calculate angle aray from 0 to pi/2 with step pi/100
  Nl=[5 10 20];
  deltax=lo./Nl;                %calculate dx, step of x axes , dx=dy
  
 K=zeros(length(fi),3);         %array for every final K, after Newton - Raphson 
 Vp=zeros(length(fi),3);        %array for wave velocity 
 normalize=zeros(length(fi),3); %array for normalization wave velocity 
  for j=1:3
    
   c=1;
   dt=(deltax(j)*deltax(j))/2*c; %calculate dt, for stabillity system  
   omega=(2*pi*c)/lo;            %calculate omega
   S=0.5;                        %factor stability,  S=cdt/dx  

   k=2*pi;                       %start value for Newton Raphson Method
   N=3;                          %times for calculate k for Newton Raphson
   
                                 %calculate array A,B and C from 4.16b in A. Taflove and S.C Hagness,
                                 %Computational Electrodynamics. The
                                 %Finite-Difference Time-Domain Method, 2nd edition, 
   A=deltax(j)*cos(fi)/2;        %Artech House Inc., Norwood,2000 
   B=deltax(j)*sin(fi)/2;
   C=(1/S^2)*(sin(pi*S/Nl(j)))^2;
    
   %Newton Raphson method 4.16a in A. Taflove and S.C Hagness,
   %Computational Electrodynamics. The
   %Finite-Difference Time-Domain Method, 2nd edition, 
   %Artech House Inc., Norwood,2000 
   %UNIFORM GRID

    for i=1:N
      k=k - ((sin(A.*k)).^2+(sin(B.*k)).^2-C)./(A.*sin(2.*A.*k)+B.*sin(2.*B.*k));
        
    end

  K(:,j)=k;
  Vp(:,j)=(2*pi*c)./K(:,j);
  normalize(:,j)=Vp(:,j)./c;

  end
%figure_1 Uniform Grid
  
 fidegree=fi.*(180/pi);

 figure(1)
 plot(fidegree,normalize(:,1),'-')
 hold on
 plot(fidegree,normalize(:,2),'-')
 plot(fidegree,normalize(:,3),'-')
 title("Numerical DispersionIn 2D FTDT Uniform Grid")
 xlabel("Wave Angle, fi (degrees)")
 ylabel("Normalized Phase Velocity Vp/c")
 legend({'lo/5','lo/10','lo/20'},'Location','east')
 grid on

end 

if (num==2)   %case 2 => dy=dx/2
    
 deltay=deltax./2;
 
 S=0.5;
  for j=1:3   %calculate k with Newton Raphson from the eq.4.20
              %Stephen D. Gedney University of Kentucky
              %Introduction to the Finite-Difference Time-Domain (FDTD)
              %Method for Electromagnetics
              % UN-UNIFORM GRID
         
  
     for ang=1:length(fi)
       k=2*pi;  
       A=(1/deltax(j)^2)*sin(k*cos(fi(ang))*deltax(j)/2)^2;
       B=(1/deltay(j)^2)*sin(k*sin(fi(ang))*deltay(j)/2)^2;
       C=(1/(S*deltay(j))^2)*(sin(pi/lo*S*deltay(j)))^2;
       D=(cos(fi(ang))/(2*deltax(j)))*sin(k*cos(fi(ang))*deltax(j));
       E=(sin(fi(ang))/(2*deltay(j)))*sin(k*sin(fi(ang))*deltay(j));
       
        for i=1:N
           k=k - (A+B-C)/(D+E);
        end
 
    K(ang,j)=k;

     end 

   Vp(:,j)=omega./K(:,j);
   normalize(:,j)=Vp(:,j)./c;

  end
  
  %figure_2 Non-Uniform Grid

  figure(2)
  plot(fidegree,normalize(:,1),'-')
  hold on
  plot(fidegree,normalize(:,2),'-')
  plot(fidegree,normalize(:,3),'-')
  title("Numerical Dispersion In 2D FTDT Non-Uniform Grid")
  xlabel("Wave Angle, fi (degrees)")
  ylabel("Normalized Phase Velocity Vp/c")
  legend({'lo/5','lo/10','lo/20'},'Location','southeast')
  grid on
  
 end

end    
    
    

