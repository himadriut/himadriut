# Matlab code for solving stochastic partial differential equations

nx=3. ;      %Size of the box-x .  this code is only for a cubic box.

ny=3. ;      %Size of the box-y

nz=3. ;     %Size of the box-y

nt=700;        %total_time step

delx = 0.1 ;   %minimal length_scale or grid size

dely = 0.1 ;     %minimal length_scale or grid size

delz =0.1 ;      %minimal length_scale or grid size

delt = 0.01 ; %time step

nx=round(nx/delx);
ny=round(ny/dely);
nz=round(nz/delz);

Diff_coeff = 0.0000054 ;     
k_a = 0.1 ;
k_b = 1 ;    % birth-rate . 10 times bigger than death rate
sigma =0.1;        %shoart range interaction scale
adhesion=3.0;%0.0000001 ;         %due to cell cell 

ff = delt*Diff_coeff/(delx*delx); %dimensionless number
%fileID = fopen('diffusion_data.txt','w');

 

phi_n=zeros(nx,ny,nz);      % cc
psi_n=zeros(nx,ny,nz);   %tp

X=(1:nx)*delx;

Y=(1:ny)*delx;

Z=(1:nz)*delz;

for i=1:nx

  for j=1:ny
      for k2=1:nz
             sigma1=0.005;
             
             phi_n(i,j,k2) = exp(-(i-nx/2)*(i-nx/2)*delx^2/sigma1^2)*exp(-(j-ny/2)*(j-ny/2)*dely^2/sigma1^2)*exp(-(k2-nz/2)*(k2-nz/2)*delz^2/sigma1^2)/10000;
             psi_n(i,j,k2) =(exp(-(i-nx/2)*(i-nx/2)*delx^2/sigma1^2)*exp(-(j-ny/2)*(j-ny/2)*dely^2/sigma1^2)*exp(-(k2-nz/2)*(k2-nz/2)*delz^2/sigma1^2))/10000000;  %tp density assumed to be small
             %phi_n(i,j,k2) =8.5204e-15*rand;%  exp(-(i-nx/2)*(i-nx/2)*delx^2/sigma^2)*exp(-(j-ny/2)*(j-ny/2)*dely^2/sigma^2)*exp(-(k2-nz/2)*(k2-nz/2)*delz^2/sigma^2)/10000;
             %psi_n(i,j,k2) =8.5204e-17*rand;% (exp(-(i-nx/2)*(i-nx/2)*delx^2/sigma^2)*exp(-(j-ny/2)*(j-ny/2)*dely^2/sigma^2)*exp(-(k2-nz/2)*(k2-nz/2)*delz^2/sigma^2))/100000000;  %tp density assumed to be small
      end
  end

end

%norm_phi = sum(sum(sum(phi_n)));
%phi_n = phi_n/norm_phi;

%norm_psi = sum(sum(sum(psi_n)));
%psi_n = psi_n/norm_psi;



mean_phi=mean(mean(mean(phi_n)));
mean_psi=mean(mean(mean(psi_n)));
%b_noise =sqrt(k_b*(mean_phi)+k_a*(mean_phi)*(mean_phi)/2);

eta=zeros(nx,ny,nz);
eta_phi=zeros(nx,ny,nz);
b_noise=zeros(nx,ny,nz);


phi=zeros(nx,ny,nz);
total_phi=zeros(nx,ny,nz,nt);

psi=zeros(nx,ny,nz);
total_psi=zeros(nx,ny,nz,nt);

for k1= 1:nt                  %beginning the dynamics

    k1
        phi(i,j)

  for j=1:nx
      for k=1:nz

    phi_n(1,j,k) = phi_n(nx,j,k);
    psi_n(1,j,k) = psi_n(nx,j,k);
      end

  end
  for j=1:ny
      for i=1:nx

    phi_n(i,j,1)=phi_n(i,j,nz);
    psi_n(i,j,1)=psi_n(i,j,nz);
      end

  end
  for k=1:nz
      for i=1:nx

    phi_n(i,1,k)=phi_n(i,ny,k);
    psi_n(i,1,k)=psi_n(i,ny,k);
      end

  end
  
  for i=1:nx
    for j=1:ny
        for k=1:nz
            eta_phi(i,j,k)=0.000000001*sqrt(2.0*Diff_coeff*delt*(phi_n(i,j,k)))*randn; %noise for CC
            eta(i,j,k)=0.0000000001*sqrt(2.0*Diff_coeff*delt*(psi_n(i,j,k)))*randn;  %noise for TP
            b_noise(i,j,k)=sqrt(k_b*phi_n(i,j,k)+k_a/2*(phi_n(i,j,k))^2)*randn*sqrt(delt);
        end
    end
 end
    potential=zeros(nx,ny,nz);
    gx=zeros(nx,ny,nz);
    gz=zeros(nx,ny,nz);
    gzz=zeros(nx,ny,nz);
    fzz=zeros(nx,ny,nz);
    fz=zeros(nx,ny,nz);
    gy=zeros(nx,ny,nz);
    %pot=zeros(nx,ny,nz);
    fx=zeros(nx,ny,nz);
    fy=zeros(nx,ny,nz);
    %pot1=zeros(nx,ny,nz);
    %eta=zeros(nx,ny,nz);
    %pot2=zeros(nx,ny,nz);
    %pot3=zeros(nx,ny,nz);
    
     for i=1:nx
        %i

        for j=1:ny
            
            
            for k=1:nz                    %(updates the density at each site. we need to worry about the short range interactions)
                
                                                                                            
                
                            for k2=i-round((sigma/delx)):i+round((sigma/delx))
               
                                 for l=j-round((sigma/dely)):j+round((sigma/dely))
                   
                    
                                        for l1=k-round((sigma/delz)):k+round((sigma/delz))
                                            
                                               potential(i,j,k)=(adhesion/(pi*sigma*sigma)^(3/2))*exp(-(i-k2)*(i-k2)*delx*delx/(2*sigma*sigma))* exp(-(j-l)*(j-l)*dely*dely/(2*sigma*sigma))* exp(-(k-l1)*(k-l1)*delz*delz/(2*sigma*sigma));
                                               
                                                 if(k2>=1 && k2<=nx && l>=1 && l<=ny && l1>=1 && l1 <=nz)    %1,1,1
                    
                                                         gx(i,j,k)=gx(i,j,k)-(0.0*psi_n(k2,l,l1)+phi_n(k2,l,l1))*potential(i,j,k)*(i-k2)*delx/(sigma*sigma);   
                                                         gy(i,j,k)=gy(i,j,k)-(0.0*psi_n(k2,l,l1)+phi_n(k2,l,l1))*potential(i,j,k)*(j-l)*dely/(sigma*sigma);
                                                         gz(i,j,k)=gz(i,j,k)-(0.0*psi_n(k2,l,l1)+phi_n(k2,l,l1))*potential(i,j,k)*(k-l1)*delz/(sigma*sigma);
                                                         fx(i,j,k)=fx(i,j,k)-(psi_n(k2,l,l1)+0.0*phi_n(k2,l,l1))*(i-k2)*potential(i,j,k)*delx/(sigma*sigma);    
                                                         fy(i,j,k)=fy(i,j,k)-(psi_n(k2,l,l1)+0.0*phi_n(k2,l,l1))*(j-l)*potential(i,j,k)*dely/(sigma*sigma);
                                                         fz(i,j,k)=fz(i,j,k)-(psi_n(k2,l,l1)+0.0*phi_n(k2,l,l1))*(k-l1)*potential(i,j,k)*delz/(sigma*sigma); 
                                                         %pot(i,j,k)=pot(i,j,k)+ psi_n(k2,l,l1)*potential(i,j,k)*((i-k2)*(i-k2))*delx*delx/(sigma*sigma*sigma*sigma);
                                                         %gzz(i,j,k)=gzz(i,j,k)+ psi_n(k2,l,l1)*potential(i,j,k)*((k-l1)*(k-l1))*delz*delz/(sigma*sigma*sigma*sigma);
                                                         %pot1(i,j,k)=pot1(i,j,k)+ psi_n(k2,l,l1)*potential(i,j,k)*((j-l)*(j-l))*dely*dely/(sigma*sigma*sigma*sigma);
                                                         %pot2(i,j,k)=pot2(i,j,k)+ phi_n(k2,l,l1)*potential(i,j,k)*((i-k2)*(i-k2))*delx*delx/(sigma*sigma*sigma*sigma);
                                                         %fzz(i,j,k)=fzz(i,j,k)+ phi_n(k2,l,l1)*potential(i,j,k)*((k-l1)*(k-l1))*delz*delz/(sigma*sigma*sigma*sigma);
                                                         %pot3(i,j,k)=pot3(i,j,k)+ phi_n(k2,l,l1)*potential(i,j,k)*((j-l)*(j-l))*dely*dely/(sigma*sigma*sigma*sigma);
                                                         %elseif((i-k)>(sigma/delx) && (j-l)>(sigma/dely))
                                                         % potential(i,j)=0;
                                                 end
                                                   
                                                  if(k2>=1 && k2<=nx && l>=1 && l<=ny && l1<1) %1,1,2
                                                         l1_copy=l1;
                                                         l1=l1+nz;
                                                         gx(i,j,k)=gx(i,j,k)-(0.0*psi_n(k2,l,l1)+phi_n(k2,l,l1))*potential(i,j,k)*(i-k2)*delx/(sigma*sigma);   
                                                         gy(i,j,k)=gy(i,j,k)-(0.0*psi_n(k2,l,l1)+phi_n(k2,l,l1))*potential(i,j,k)*(j-l)*dely/(sigma*sigma);
                                                         gz(i,j,k)=gz(i,j,k)-(0.0*psi_n(k2,l,l1)+phi_n(k2,l,l1))*potential(i,j,k)*(k-l1)*delz/(sigma*sigma);
                                                         fx(i,j,k)=fx(i,j,k)-(psi_n(k2,l,l1)+0.0*phi_n(k2,l,l1))*(i-k2)*potential(i,j,k)*delx/(sigma*sigma);    
                                                         fy(i,j,k)=fy(i,j,k)-(psi_n(k2,l,l1)+0.0*phi_n(k2,l,l1))*(j-l)*potential(i,j,k)*dely/(sigma*sigma);
                                                         fz(i,j,k)=fz(i,j,k)-(psi_n(k2,l,l1)+0.0*phi_n(k2,l,l1))*(k-l1)*potential(i,j,k)*delz/(sigma*sigma);
                                                         l1=l1_copy;
                                                      
                                                  end
                                                  
                                                  if(k2>=1 && k2<=nx && l>=1 && l<=ny && l1>nz) %1,1,3
                                                         l1_copy=l1;
                                                         l1=l1-nz;
                                                         gx(i,j,k)=gx(i,j,k)-(0.0*psi_n(k2,l,l1)+phi_n(k2,l,l1))*potential(i,j,k)*(i-k2)*delx/(sigma*sigma);   
                                                         gy(i,j,k)=gy(i,j,k)-(0.0*psi_n(k2,l,l1)+phi_n(k2,l,l1))*potential(i,j,k)*(j-l)*dely/(sigma*sigma);
                                                         gz(i,j,k)=gz(i,j,k)-(0.0*psi_n(k2,l,l1)+phi_n(k2,l,l1))*potential(i,j,k)*(k-l1)*delz/(sigma*sigma);
                                                         fx(i,j,k)=fx(i,j,k)-(psi_n(k2,l,l1)+0.0*phi_n(k2,l,l1))*(i-k2)*potential(i,j,k)*delx/(sigma*sigma);    
                                                         fy(i,j,k)=fy(i,j,k)-(psi_n(k2,l,l1)+0.0*phi_n(k2,l,l1))*(j-l)*potential(i,j,k)*dely/(sigma*sigma);
                                                         fz(i,j,k)=fz(i,j,k)-(psi_n(k2,l,l1)+0.0*phi_n(k2,l,l1))*(k-l1)*potential(i,j,k)*delz/(sigma*sigma);
                                                         l1=l1_copy;
                                                  end
                                                  
                                                  if(k2>=1 && k2<=nx && l<1 && l1>=1 && l1 <=nz)    %1,2,1
                                                         l_copy=l;
                                                         l=ny+l;
                                                         gx(i,j,k)=gx(i,j,k)-(0.0*psi_n(k2,l,l1)+phi_n(k2,l,l1))*potential(i,j,k)*(i-k2)*delx/(sigma*sigma);   
                                                         gy(i,j,k)=gy(i,j,k)-(0.0*psi_n(k2,l,l1)+phi_n(k2,l,l1))*potential(i,j,k)*(j-l)*dely/(sigma*sigma);
                                                         gz(i,j,k)=gz(i,j,k)-(0.0*psi_n(k2,l,l1)+phi_n(k2,l,l1))*potential(i,j,k)*(k-l1)*delz/(sigma*sigma);
                                                         fx(i,j,k)=fx(i,j,k)-(psi_n(k2,l,l1)+0.0*phi_n(k2,l,l1))*(i-k2)*potential(i,j,k)*delx/(sigma*sigma);    
                                                         fy(i,j,k)=fy(i,j,k)-(psi_n(k2,l,l1)+0.0*phi_n(k2,l,l1))*(j-l)*potential(i,j,k)*dely/(sigma*sigma);
                                                         fz(i,j,k)=fz(i,j,k)-(psi_n(k2,l,l1)+0.0*phi_n(k2,l,l1))*(k-l1)*potential(i,j,k)*delz/(sigma*sigma);
                                                         l=l_copy;
                                                      
                                                  end
                                                  
                                                   if(k2>=1 && k2<=nx && l<1 && l1<1)            %1,2,2
                                                       l_copy=l;
                                                       l=ny+l;
                                                       l1_copy=l1;
                                                       l1=l1+nz;
                                                       gx(i,j,k)=gx(i,j,k)-(0.0*psi_n(k2,l,l1)+phi_n(k2,l,l1))*potential(i,j,k)*(i-k2)*delx/(sigma*sigma);   
                                                       gy(i,j,k)=gy(i,j,k)-(0.0*psi_n(k2,l,l1)+phi_n(k2,l,l1))*potential(i,j,k)*(j-l)*dely/(sigma*sigma);
                                                       gz(i,j,k)=gz(i,j,k)-(0.0*psi_n(k2,l,l1)+phi_n(k2,l,l1))*potential(i,j,k)*(k-l1)*delz/(sigma*sigma);
                                                       fx(i,j,k)=fx(i,j,k)-(psi_n(k2,l,l1)+0.0*phi_n(k2,l,l1))*(i-k2)*potential(i,j,k)*delx/(sigma*sigma);    
                                                       fy(i,j,k)=fy(i,j,k)-(psi_n(k2,l,l1)+0.0*phi_n(k2,l,l1))*(j-l)*potential(i,j,k)*dely/(sigma*sigma);
                                                       fz(i,j,k)=fz(i,j,k)-(psi_n(k2,l,l1)+0.0*phi_n(k2,l,l1))*(k-l1)*potential(i,j,k)*delz/(sigma*sigma);
                                                       l=l_copy;
                                                       l1=l1_copy;
                                                       
                                                   end
                                                   
                                                   if(k2>=1 && k2<=nx && l<1 && l1>nz)   %1,2,3
                                                       
                                                       l_copy=l;
                                                       l=ny+l;
                                                       l1_copy=l1;
                                                       l1=l1-nz;
                                                       gx(i,j,k)=gx(i,j,k)-(0.0*psi_n(k2,l,l1)+phi_n(k2,l,l1))*potential(i,j,k)*(i-k2)*delx/(sigma*sigma);   
                                                       gy(i,j,k)=gy(i,j,k)-(0.0*psi_n(k2,l,l1)+phi_n(k2,l,l1))*potential(i,j,k)*(j-l)*dely/(sigma*sigma);
                                                       gz(i,j,k)=gz(i,j,k)-(0.0*psi_n(k2,l,l1)+phi_n(k2,l,l1))*potential(i,j,k)*(k-l1)*delz/(sigma*sigma);
                                                       fx(i,j,k)=fx(i,j,k)-(psi_n(k2,l,l1)+0.0*phi_n(k2,l,l1))*(i-k2)*potential(i,j,k)*delx/(sigma*sigma);    
                                                       fy(i,j,k)=fy(i,j,k)-(psi_n(k2,l,l1)+0.0*phi_n(k2,l,l1))*(j-l)*potential(i,j,k)*dely/(sigma*sigma);
                                                       fz(i,j,k)=fz(i,j,k)-(psi_n(k2,l,l1)+0.0*phi_n(k2,l,l1))*(k-l1)*potential(i,j,k)*delz/(sigma*sigma);
                                                       l=l_copy;
                                                       l1=l1_copy;
                                                       
                                                       
                                                   end
                                                   
                                                    if(k2>=1 && k2<=nx && l>ny && l1>=1 && l1 <=nz)    %1,3,1
                                                       l_copy=l;
                                                       l=l-nz;
                                                       gx(i,j,k)=gx(i,j,k)-(0.0*psi_n(k2,l,l1)+phi_n(k2,l,l1))*potential(i,j,k)*(i-k2)*delx/(sigma*sigma);   
                                                       gy(i,j,k)=gy(i,j,k)-(0.0*psi_n(k2,l,l1)+phi_n(k2,l,l1))*potential(i,j,k)*(j-l)*dely/(sigma*sigma);
                                                       gz(i,j,k)=gz(i,j,k)-(0.0*psi_n(k2,l,l1)+phi_n(k2,l,l1))*potential(i,j,k)*(k-l1)*delz/(sigma*sigma);
                                                       fx(i,j,k)=fx(i,j,k)-(psi_n(k2,l,l1)+0.0*phi_n(k2,l,l1))*(i-k2)*potential(i,j,k)*delx/(sigma*sigma);    
                                                       fy(i,j,k)=fy(i,j,k)-(psi_n(k2,l,l1)+0.0*phi_n(k2,l,l1))*(j-l)*potential(i,j,k)*dely/(sigma*sigma);
                                                       fz(i,j,k)=fz(i,j,k)-(psi_n(k2,l,l1)+0.0*phi_n(k2,l,l1))*(k-l1)*potential(i,j,k)*delz/(sigma*sigma);
                                                       l=l_copy;
                                                        
                                                        
                                                    end
                                                    
                                                    if(k2>=1 && k2<=nx && l>ny && l1<1)     %1,3,2
                                                        l_copy=l;
                                                        l=l-nz;
                                                        l1_copy=l1;
                                                        l1=l1+nz;
                                                        gx(i,j,k)=gx(i,j,k)-(0.0*psi_n(k2,l,l1)+phi_n(k2,l,l1))*potential(i,j,k)*(i-k2)*delx/(sigma*sigma);   
                                                       gy(i,j,k)=gy(i,j,k)-(0.0*psi_n(k2,l,l1)+phi_n(k2,l,l1))*potential(i,j,k)*(j-l)*dely/(sigma*sigma);
                                                       gz(i,j,k)=gz(i,j,k)-(0.0*psi_n(k2,l,l1)+phi_n(k2,l,l1))*potential(i,j,k)*(k-l1)*delz/(sigma*sigma);
                                                       fx(i,j,k)=fx(i,j,k)-(psi_n(k2,l,l1)+0.0*phi_n(k2,l,l1))*(i-k2)*potential(i,j,k)*delx/(sigma*sigma);    
                                                       fy(i,j,k)=fy(i,j,k)-(psi_n(k2,l,l1)+0.0*phi_n(k2,l,l1))*(j-l)*potential(i,j,k)*dely/(sigma*sigma);
                                                       fz(i,j,k)=fz(i,j,k)-(psi_n(k2,l,l1)+0.0*phi_n(k2,l,l1))*(k-l1)*potential(i,j,k)*delz/(sigma*sigma);
                                                       l=l_copy;
                                                       l1=l1_copy;
                                                        
                                                        
                                                    end
                                                    
                                                     if(k2>=1 && k2<=nx && l>ny && l1>nz) %1,3,3
                                                         l_copy=l;
                                                         l=l-nz;
                                                         l1_copy=l1;
                                                         l1=l1-nz;
                                                          gx(i,j,k)=gx(i,j,k)-(0.0*psi_n(k2,l,l1)+phi_n(k2,l,l1))*potential(i,j,k)*(i-k2)*delx/(sigma*sigma);   
                                                         gy(i,j,k)=gy(i,j,k)-(0.0*psi_n(k2,l,l1)+phi_n(k2,l,l1))*potential(i,j,k)*(j-l)*dely/(sigma*sigma);
                                                         gz(i,j,k)=gz(i,j,k)-(0.0*psi_n(k2,l,l1)+phi_n(k2,l,l1))*potential(i,j,k)*(k-l1)*delz/(sigma*sigma);
                                                       fx(i,j,k)=fx(i,j,k)-(psi_n(k2,l,l1)+0.0*phi_n(k2,l,l1))*(i-k2)*potential(i,j,k)*delx/(sigma*sigma);    
                                                       fy(i,j,k)=fy(i,j,k)-(psi_n(k2,l,l1)+0.0*phi_n(k2,l,l1))*(j-l)*potential(i,j,k)*dely/(sigma*sigma);
                                                       fz(i,j,k)=fz(i,j,k)-(psi_n(k2,l,l1)+0.0*phi_n(k2,l,l1))*(k-l1)*potential(i,j,k)*delz/(sigma*sigma);
                                                       l=l_copy;
                                                       l1=l1_copy;
                                                         
                                                     end
                                                     
                                                     
                                                     
                                                      if(k2<1 && l>=1 && l<=ny && l1>=1 && l1 <=nz)  %2,1,1
                                                          
                                                           k2_copy=k2;
                                                           k2=k2+nz;
                                                           
                                                            gx(i,j,k)=gx(i,j,k)-(0.0*psi_n(k2,l,l1)+phi_n(k2,l,l1))*potential(i,j,k)*(i-k2)*delx/(sigma*sigma);   
                                                         gy(i,j,k)=gy(i,j,k)-(0.0*psi_n(k2,l,l1)+phi_n(k2,l,l1))*potential(i,j,k)*(j-l)*dely/(sigma*sigma);
                                                         gz(i,j,k)=gz(i,j,k)-(0.0*psi_n(k2,l,l1)+phi_n(k2,l,l1))*potential(i,j,k)*(k-l1)*delz/(sigma*sigma);
                                                       fx(i,j,k)=fx(i,j,k)-(psi_n(k2,l,l1)+0.0*phi_n(k2,l,l1))*(i-k2)*potential(i,j,k)*delx/(sigma*sigma);    
                                                       fy(i,j,k)=fy(i,j,k)-(psi_n(k2,l,l1)+0.0*phi_n(k2,l,l1))*(j-l)*potential(i,j,k)*dely/(sigma*sigma);
                                                       fz(i,j,k)=fz(i,j,k)-(psi_n(k2,l,l1)+0.0*phi_n(k2,l,l1))*(k-l1)*potential(i,j,k)*delz/(sigma*sigma);
                                                       k2=k2_copy;
                                                      
                                                          
                                                      end
                                                      
                                                       if(k2<1 && l>=1 && l<=ny && l1<1) %2,1,2
                                                            k2_copy=k2;
                                                           k2=k2+nz;
                                                            l1_copy=l1;
                                                           l1=l1+nz;
                                                            gx(i,j,k)=gx(i,j,k)-(0.0*psi_n(k2,l,l1)+phi_n(k2,l,l1))*potential(i,j,k)*(i-k2)*delx/(sigma*sigma);   
                                                         gy(i,j,k)=gy(i,j,k)-(0.0*psi_n(k2,l,l1)+phi_n(k2,l,l1))*potential(i,j,k)*(j-l)*dely/(sigma*sigma);
                                                         gz(i,j,k)=gz(i,j,k)-(0.0*psi_n(k2,l,l1)+phi_n(k2,l,l1))*potential(i,j,k)*(k-l1)*delz/(sigma*sigma);
                                                       fx(i,j,k)=fx(i,j,k)-(psi_n(k2,l,l1)+0.0*phi_n(k2,l,l1))*(i-k2)*potential(i,j,k)*delx/(sigma*sigma);    
                                                       fy(i,j,k)=fy(i,j,k)-(psi_n(k2,l,l1)+0.0*phi_n(k2,l,l1))*(j-l)*potential(i,j,k)*dely/(sigma*sigma);
                                                       fz(i,j,k)=fz(i,j,k)-(psi_n(k2,l,l1)+0.0*phi_n(k2,l,l1))*(k-l1)*potential(i,j,k)*delz/(sigma*sigma);
                                                       k2=k2_copy;
                                                       l1=l1_copy;
                                                           
                                                       end
                                                       
                                                       if(k2<1 && l>=1 && l<=ny && l1>nz)   %2,1,3
                                                            k2_copy=k2;
                                                           k2=k2+nz;
                                                           l1_copy=l1;
                                                           l1=l1-nz;
                                                           gx(i,j,k)=gx(i,j,k)-(0.0*psi_n(k2,l,l1)+phi_n(k2,l,l1))*potential(i,j,k)*(i-k2)*delx/(sigma*sigma);   
                                                         gy(i,j,k)=gy(i,j,k)-(0.0*psi_n(k2,l,l1)+phi_n(k2,l,l1))*potential(i,j,k)*(j-l)*dely/(sigma*sigma);
                                                         gz(i,j,k)=gz(i,j,k)-(0.0*psi_n(k2,l,l1)+phi_n(k2,l,l1))*potential(i,j,k)*(k-l1)*delz/(sigma*sigma);
                                                       fx(i,j,k)=fx(i,j,k)-(psi_n(k2,l,l1)+0.0*phi_n(k2,l,l1))*(i-k2)*potential(i,j,k)*delx/(sigma*sigma);    
                                                       fy(i,j,k)=fy(i,j,k)-(psi_n(k2,l,l1)+0.0*phi_n(k2,l,l1))*(j-l)*potential(i,j,k)*dely/(sigma*sigma);
                                                       fz(i,j,k)=fz(i,j,k)-(psi_n(k2,l,l1)+0.0*phi_n(k2,l,l1))*(k-l1)*potential(i,j,k)*delz/(sigma*sigma);
                                                       k2=k2_copy;
                                                       l1=l1_copy;
                                                           
                                                       end
                                                       
                                                       if(k2<1 && l<1 && l1>=1 && l1 <=nz) %2,2,1
                                                             k2_copy=k2;
                                                             k2=k2+nz;
                                                             l_copy=l;
                                                             l=l+nz;
                                                             gx(i,j,k)=gx(i,j,k)-(0.0*psi_n(k2,l,l1)+phi_n(k2,l,l1))*potential(i,j,k)*(i-k2)*delx/(sigma*sigma);   
                                                         gy(i,j,k)=gy(i,j,k)-(0.0*psi_n(k2,l,l1)+phi_n(k2,l,l1))*potential(i,j,k)*(j-l)*dely/(sigma*sigma);
                                                         gz(i,j,k)=gz(i,j,k)-(0.0*psi_n(k2,l,l1)+phi_n(k2,l,l1))*potential(i,j,k)*(k-l1)*delz/(sigma*sigma);
                                                       fx(i,j,k)=fx(i,j,k)-(psi_n(k2,l,l1)+0.0*phi_n(k2,l,l1))*(i-k2)*potential(i,j,k)*delx/(sigma*sigma);    
                                                       fy(i,j,k)=fy(i,j,k)-(psi_n(k2,l,l1)+0.0*phi_n(k2,l,l1))*(j-l)*potential(i,j,k)*dely/(sigma*sigma);
                                                       fz(i,j,k)=fz(i,j,k)-(psi_n(k2,l,l1)+0.0*phi_n(k2,l,l1))*(k-l1)*potential(i,j,k)*delz/(sigma*sigma);
                                                       k2=k2_copy;
                                                       l=l_copy;
                                                             
                                                           
                                                       end
                                                       
                                                       if(k2<1 && l<1 && l1<1)      %2,2,2
                                                           k2_copy=k2;
                                                             k2=k2+nz;
                                                             l_copy=l;
                                                             l=l+nz;
                                                              l1_copy=l1;
                                                             l1=l1+nz;
                                                             gx(i,j,k)=gx(i,j,k)-(0.0*psi_n(k2,l,l1)+phi_n(k2,l,l1))*potential(i,j,k)*(i-k2)*delx/(sigma*sigma);   
                                                         gy(i,j,k)=gy(i,j,k)-(0.0*psi_n(k2,l,l1)+phi_n(k2,l,l1))*potential(i,j,k)*(j-l)*dely/(sigma*sigma);
                                                         gz(i,j,k)=gz(i,j,k)-(0.0*psi_n(k2,l,l1)+phi_n(k2,l,l1))*potential(i,j,k)*(k-l1)*delz/(sigma*sigma);
                                                       fx(i,j,k)=fx(i,j,k)-(psi_n(k2,l,l1)+0.0*phi_n(k2,l,l1))*(i-k2)*potential(i,j,k)*delx/(sigma*sigma);    
                                                       fy(i,j,k)=fy(i,j,k)-(psi_n(k2,l,l1)+0.0*phi_n(k2,l,l1))*(j-l)*potential(i,j,k)*dely/(sigma*sigma);
                                                       fz(i,j,k)=fz(i,j,k)-(psi_n(k2,l,l1)+0.0*phi_n(k2,l,l1))*(k-l1)*potential(i,j,k)*delz/(sigma*sigma);
                                                       k2=k2_copy;
                                                       l=l_copy;
                                                       l1=l1_copy;
                                                           
                                                           
                                                       end
                                                       
                                                       if(k2<1 && l<1 && l1>nz)      %2,2,3
                                                           
                                                            k2_copy=k2;
                                                             k2=k2+nz;
                                                             l_copy=l;
                                                             l=l+nz;
                                                              l1_copy=l1;
                                                             l1=l1-nz;
                                                             gx(i,j,k)=gx(i,j,k)-(0.0*psi_n(k2,l,l1)+phi_n(k2,l,l1))*potential(i,j,k)*(i-k2)*delx/(sigma*sigma);   
                                                         gy(i,j,k)=gy(i,j,k)-(0.0*psi_n(k2,l,l1)+phi_n(k2,l,l1))*potential(i,j,k)*(j-l)*dely/(sigma*sigma);
                                                         gz(i,j,k)=gz(i,j,k)-(0.0*psi_n(k2,l,l1)+phi_n(k2,l,l1))*potential(i,j,k)*(k-l1)*delz/(sigma*sigma);
                                                       fx(i,j,k)=fx(i,j,k)-(psi_n(k2,l,l1)+0.0*phi_n(k2,l,l1))*(i-k2)*potential(i,j,k)*delx/(sigma*sigma);    
                                                       fy(i,j,k)=fy(i,j,k)-(psi_n(k2,l,l1)+0.0*phi_n(k2,l,l1))*(j-l)*potential(i,j,k)*dely/(sigma*sigma);
                                                       fz(i,j,k)=fz(i,j,k)-(psi_n(k2,l,l1)+0.0*phi_n(k2,l,l1))*(k-l1)*potential(i,j,k)*delz/(sigma*sigma);
                                                       k2=k2_copy;
                                                       l=l_copy;
                                                       l1=l1_copy;
                                                           
                                                       end
                                                       
                                                       if(k2<1 && l>ny && l1>=1 && l1 <=nz) %2,3,1
                                                           k2_copy=k2;
                                                             k2=k2+nz;
                                                             l_copy=l;
                                                             l=l-nz;
                                                              gx(i,j,k)=gx(i,j,k)-(0.0*psi_n(k2,l,l1)+phi_n(k2,l,l1))*potential(i,j,k)*(i-k2)*delx/(sigma*sigma);   
                                                         gy(i,j,k)=gy(i,j,k)-(0.0*psi_n(k2,l,l1)+phi_n(k2,l,l1))*potential(i,j,k)*(j-l)*dely/(sigma*sigma);
                                                         gz(i,j,k)=gz(i,j,k)-(0.0*psi_n(k2,l,l1)+phi_n(k2,l,l1))*potential(i,j,k)*(k-l1)*delz/(sigma*sigma);
                                                       fx(i,j,k)=fx(i,j,k)-(psi_n(k2,l,l1)+0.0*phi_n(k2,l,l1))*(i-k2)*potential(i,j,k)*delx/(sigma*sigma);    
                                                       fy(i,j,k)=fy(i,j,k)-(psi_n(k2,l,l1)+0.0*phi_n(k2,l,l1))*(j-l)*potential(i,j,k)*dely/(sigma*sigma);
                                                       fz(i,j,k)=fz(i,j,k)-(psi_n(k2,l,l1)+0.0*phi_n(k2,l,l1))*(k-l1)*potential(i,j,k)*delz/(sigma*sigma);
                                                       k2=k2_copy;
                                                       l=l_copy;
                                                           
                                                       end
                                                       
                                                       if(k2<1 && l>ny && l1<1)     %2,3,2
                                                            k2_copy=k2;
                                                             k2=k2+nz;
                                                             l_copy=l;
                                                             l=l-nz;
                                                             l1_copy=l1;
                                                             l1=l1+nz;
                                                             gx(i,j,k)=gx(i,j,k)-(0.0*psi_n(k2,l,l1)+phi_n(k2,l,l1))*potential(i,j,k)*(i-k2)*delx/(sigma*sigma);   
                                                         gy(i,j,k)=gy(i,j,k)-(0.0*psi_n(k2,l,l1)+phi_n(k2,l,l1))*potential(i,j,k)*(j-l)*dely/(sigma*sigma);
                                                         gz(i,j,k)=gz(i,j,k)-(0.0*psi_n(k2,l,l1)+phi_n(k2,l,l1))*potential(i,j,k)*(k-l1)*delz/(sigma*sigma);
                                                       fx(i,j,k)=fx(i,j,k)-(psi_n(k2,l,l1)+0.0*phi_n(k2,l,l1))*(i-k2)*potential(i,j,k)*delx/(sigma*sigma);    
                                                       fy(i,j,k)=fy(i,j,k)-(psi_n(k2,l,l1)+0.0*phi_n(k2,l,l1))*(j-l)*potential(i,j,k)*dely/(sigma*sigma);
                                                       fz(i,j,k)=fz(i,j,k)-(psi_n(k2,l,l1)+0.0*phi_n(k2,l,l1))*(k-l1)*potential(i,j,k)*delz/(sigma*sigma);
                                                       k2=k2_copy;
                                                       l=l_copy;
                                                       l1=l1_copy;
                                                       end
                                                       
                                                       if(k2<1 && l>ny && l1>nz)   %2,3,3
                                                           k2_copy=k2;
                                                             k2=k2+nz;
                                                             l_copy=l;
                                                             l=l-nz;
                                                             l1_copy=l1;
                                                             l1=l1-nz;
                                                             gx(i,j,k)=gx(i,j,k)-(0.0*psi_n(k2,l,l1)+phi_n(k2,l,l1))*potential(i,j,k)*(i-k2)*delx/(sigma*sigma);   
                                                         gy(i,j,k)=gy(i,j,k)-(0.0*psi_n(k2,l,l1)+phi_n(k2,l,l1))*potential(i,j,k)*(j-l)*dely/(sigma*sigma);
                                                         gz(i,j,k)=gz(i,j,k)-(0.0*psi_n(k2,l,l1)+phi_n(k2,l,l1))*potential(i,j,k)*(k-l1)*delz/(sigma*sigma);
                                                       fx(i,j,k)=fx(i,j,k)-(psi_n(k2,l,l1)+0.0*phi_n(k2,l,l1))*(i-k2)*potential(i,j,k)*delx/(sigma*sigma);    
                                                       fy(i,j,k)=fy(i,j,k)-(psi_n(k2,l,l1)+0.0*phi_n(k2,l,l1))*(j-l)*potential(i,j,k)*dely/(sigma*sigma);
                                                       fz(i,j,k)=fz(i,j,k)-(psi_n(k2,l,l1)+0.0*phi_n(k2,l,l1))*(k-l1)*potential(i,j,k)*delz/(sigma*sigma);
                                                       k2=k2_copy;
                                                       l=l_copy;
                                                       l1=l1_copy;
                                                           
                                                       end
                                                       
                                                       if(k2>nx && l>=1 && l<=ny && l1>=1 && l1 <=nz)   % 3,1,1
                                                           k2_copy=k2;
                                                           k2=k2-nz;
                                                           gx(i,j,k)=gx(i,j,k)-(0.0*psi_n(k2,l,l1)+phi_n(k2,l,l1))*potential(i,j,k)*(i-k2)*delx/(sigma*sigma);   
                                                         gy(i,j,k)=gy(i,j,k)-(0.0*psi_n(k2,l,l1)+phi_n(k2,l,l1))*potential(i,j,k)*(j-l)*dely/(sigma*sigma);
                                                         gz(i,j,k)=gz(i,j,k)-(0.0*psi_n(k2,l,l1)+phi_n(k2,l,l1))*potential(i,j,k)*(k-l1)*delz/(sigma*sigma);
                                                       fx(i,j,k)=fx(i,j,k)-(psi_n(k2,l,l1)+0.0*phi_n(k2,l,l1))*(i-k2)*potential(i,j,k)*delx/(sigma*sigma);    
                                                       fy(i,j,k)=fy(i,j,k)-(psi_n(k2,l,l1)+0.0*phi_n(k2,l,l1))*(j-l)*potential(i,j,k)*dely/(sigma*sigma);
                                                       fz(i,j,k)=fz(i,j,k)-(psi_n(k2,l,l1)+0.0*phi_n(k2,l,l1))*(k-l1)*potential(i,j,k)*delz/(sigma*sigma);
                                                       k2=k2_copy;
                                                           
                                                       end
                                                       
                                                       if(k2>nx && l>=1 && l<=ny && l1<1)   % 3,1,2
                                                           k2_copy=k2;
                                                           k2=k2-nz;
                                                           l1_copy=l1;
                                                             l1=l1+nz;
                                                             gx(i,j,k)=gx(i,j,k)-(0.0*psi_n(k2,l,l1)+phi_n(k2,l,l1))*potential(i,j,k)*(i-k2)*delx/(sigma*sigma);   
                                                         gy(i,j,k)=gy(i,j,k)-(0.0*psi_n(k2,l,l1)+phi_n(k2,l,l1))*potential(i,j,k)*(j-l)*dely/(sigma*sigma);
                                                         gz(i,j,k)=gz(i,j,k)-(0.0*psi_n(k2,l,l1)+phi_n(k2,l,l1))*potential(i,j,k)*(k-l1)*delz/(sigma*sigma);
                                                       fx(i,j,k)=fx(i,j,k)-(psi_n(k2,l,l1)+0.0*phi_n(k2,l,l1))*(i-k2)*potential(i,j,k)*delx/(sigma*sigma);    
                                                       fy(i,j,k)=fy(i,j,k)-(psi_n(k2,l,l1)+0.0*phi_n(k2,l,l1))*(j-l)*potential(i,j,k)*dely/(sigma*sigma);
                                                       fz(i,j,k)=fz(i,j,k)-(psi_n(k2,l,l1)+0.0*phi_n(k2,l,l1))*(k-l1)*potential(i,j,k)*delz/(sigma*sigma);
                                                       k2=k2_copy;
                                                       l1=l1_copy;
                                                       end 
                                                       
                                                       if(k2>nx && l>=1 && l<=ny && l1>nz)   % 3,1,3
                                                           
                                                           k2_copy=k2;
                                                           k2=k2-nz;
                                                           l1_copy=l1;
                                                             l1=l1-nz;
                                                             gx(i,j,k)=gx(i,j,k)-(0.0*psi_n(k2,l,l1)+phi_n(k2,l,l1))*potential(i,j,k)*(i-k2)*delx/(sigma*sigma);   
                                                         gy(i,j,k)=gy(i,j,k)-(0.0*psi_n(k2,l,l1)+phi_n(k2,l,l1))*potential(i,j,k)*(j-l)*dely/(sigma*sigma);
                                                         gz(i,j,k)=gz(i,j,k)-(0.0*psi_n(k2,l,l1)+phi_n(k2,l,l1))*potential(i,j,k)*(k-l1)*delz/(sigma*sigma);
                                                       fx(i,j,k)=fx(i,j,k)-(psi_n(k2,l,l1)+0.0*phi_n(k2,l,l1))*(i-k2)*potential(i,j,k)*delx/(sigma*sigma);    
                                                       fy(i,j,k)=fy(i,j,k)-(psi_n(k2,l,l1)+0.0*phi_n(k2,l,l1))*(j-l)*potential(i,j,k)*dely/(sigma*sigma);
                                                       fz(i,j,k)=fz(i,j,k)-(psi_n(k2,l,l1)+0.0*phi_n(k2,l,l1))*(k-l1)*potential(i,j,k)*delz/(sigma*sigma);
                                                       k2=k2_copy;
                                                       l1=l1_copy;
                                                           
                                                       end
                                                       
                                                       if(k2>nx && l<1 && l1>=1 && l1 <=nz)   % 3,2,1
                                                           
                                                           k2_copy=k2;
                                                           k2=k2-nz;
                                                           l_copy=l;
                                                             l=l+nz;
                                                             gx(i,j,k)=gx(i,j,k)-(0.0*psi_n(k2,l,l1)+phi_n(k2,l,l1))*potential(i,j,k)*(i-k2)*delx/(sigma*sigma);   
                                                         gy(i,j,k)=gy(i,j,k)-(0.0*psi_n(k2,l,l1)+phi_n(k2,l,l1))*potential(i,j,k)*(j-l)*dely/(sigma*sigma);
                                                         gz(i,j,k)=gz(i,j,k)-(0.0*psi_n(k2,l,l1)+phi_n(k2,l,l1))*potential(i,j,k)*(k-l1)*delz/(sigma*sigma);
                                                       fx(i,j,k)=fx(i,j,k)-(psi_n(k2,l,l1)+0.0*phi_n(k2,l,l1))*(i-k2)*potential(i,j,k)*delx/(sigma*sigma);    
                                                       fy(i,j,k)=fy(i,j,k)-(psi_n(k2,l,l1)+0.0*phi_n(k2,l,l1))*(j-l)*potential(i,j,k)*dely/(sigma*sigma);
                                                       fz(i,j,k)=fz(i,j,k)-(psi_n(k2,l,l1)+0.0*phi_n(k2,l,l1))*(k-l1)*potential(i,j,k)*delz/(sigma*sigma);
                                                       k2=k2_copy;
                                                       l=l_copy;
                                                       end
                                                       
                                                       if(k2>nx && l<1 && l1<1)   % 3,2,2
                                                           k2_copy=k2;
                                                             k2=k2-nz;
                                                             l_copy=l;
                                                             l=l+nz;
                                                             l1_copy=l1;
                                                             l1=l1+nz;
                                                             gx(i,j,k)=gx(i,j,k)-(0.0*psi_n(k2,l,l1)+phi_n(k2,l,l1))*potential(i,j,k)*(i-k2)*delx/(sigma*sigma);   
                                                         gy(i,j,k)=gy(i,j,k)-(0.0*psi_n(k2,l,l1)+phi_n(k2,l,l1))*potential(i,j,k)*(j-l)*dely/(sigma*sigma);
                                                         gz(i,j,k)=gz(i,j,k)-(0.0*psi_n(k2,l,l1)+phi_n(k2,l,l1))*potential(i,j,k)*(k-l1)*delz/(sigma*sigma);
                                                       fx(i,j,k)=fx(i,j,k)-(psi_n(k2,l,l1)+0.0*phi_n(k2,l,l1))*(i-k2)*potential(i,j,k)*delx/(sigma*sigma);    
                                                       fy(i,j,k)=fy(i,j,k)-(psi_n(k2,l,l1)+0.0*phi_n(k2,l,l1))*(j-l)*potential(i,j,k)*dely/(sigma*sigma);
                                                       fz(i,j,k)=fz(i,j,k)-(psi_n(k2,l,l1)+0.0*phi_n(k2,l,l1))*(k-l1)*potential(i,j,k)*delz/(sigma*sigma);
                                                       k2=k2_copy;
                                                       l=l_copy;
                                                       l1=l1_copy;
                                                           
                                                       end
                                                       
                                                       if(k2>nx && l<1 && l1>nz)   % 3,2,3
                                                           
                                                           k2_copy=k2;
                                                             k2=k2-nz;
                                                             l_copy=l;
                                                             l=l+nz;
                                                             l1_copy=l1;
                                                             l1=l1-nz;
                                                             gx(i,j,k)=gx(i,j,k)-(0.0*psi_n(k2,l,l1)+phi_n(k2,l,l1))*potential(i,j,k)*(i-k2)*delx/(sigma*sigma);   
                                                         gy(i,j,k)=gy(i,j,k)-(0.0*psi_n(k2,l,l1)+phi_n(k2,l,l1))*potential(i,j,k)*(j-l)*dely/(sigma*sigma);
                                                         gz(i,j,k)=gz(i,j,k)-(0.0*psi_n(k2,l,l1)+phi_n(k2,l,l1))*potential(i,j,k)*(k-l1)*delz/(sigma*sigma);
                                                       fx(i,j,k)=fx(i,j,k)-(psi_n(k2,l,l1)+0.0*phi_n(k2,l,l1))*(i-k2)*potential(i,j,k)*delx/(sigma*sigma);    
                                                       fy(i,j,k)=fy(i,j,k)-(psi_n(k2,l,l1)+0.0*phi_n(k2,l,l1))*(j-l)*potential(i,j,k)*dely/(sigma*sigma);
                                                       fz(i,j,k)=fz(i,j,k)-(psi_n(k2,l,l1)+0.0*phi_n(k2,l,l1))*(k-l1)*potential(i,j,k)*delz/(sigma*sigma);
                                                       k2=k2_copy;
                                                       l=l_copy;
                                                       l1=l1_copy;
                                                           
                                                       end
                                                       
                                                       if(k2>nx && l>ny && l1>=1 && l1 <=nz)   % 3,3,1
                                                           k2_copy=k2;
                                                             k2=k2-nz;
                                                             l_copy=l;
                                                             l=l-nz;
                                                            
                                                             gx(i,j,k)=gx(i,j,k)-(0.0*psi_n(k2,l,l1)+phi_n(k2,l,l1))*potential(i,j,k)*(i-k2)*delx/(sigma*sigma);   
                                                         gy(i,j,k)=gy(i,j,k)-(0.0*psi_n(k2,l,l1)+phi_n(k2,l,l1))*potential(i,j,k)*(j-l)*dely/(sigma*sigma);
                                                         gz(i,j,k)=gz(i,j,k)-(0.0*psi_n(k2,l,l1)+phi_n(k2,l,l1))*potential(i,j,k)*(k-l1)*delz/(sigma*sigma);
                                                       fx(i,j,k)=fx(i,j,k)-(psi_n(k2,l,l1)+0.0*phi_n(k2,l,l1))*(i-k2)*potential(i,j,k)*delx/(sigma*sigma);    
                                                       fy(i,j,k)=fy(i,j,k)-(psi_n(k2,l,l1)+0.0*phi_n(k2,l,l1))*(j-l)*potential(i,j,k)*dely/(sigma*sigma);
                                                       fz(i,j,k)=fz(i,j,k)-(psi_n(k2,l,l1)+0.0*phi_n(k2,l,l1))*(k-l1)*potential(i,j,k)*delz/(sigma*sigma);
                                                       k2=k2_copy;
                                                       l=l_copy;
                                                      
                                                       end
                                                       
                                                       if(k2>nx && l>ny && l1<1)   % 3,3,2
                                                           
                                                           k2_copy=k2;
                                                             k2=k2-nz;
                                                             l_copy=l;
                                                             l=l-nz;
                                                             l1_copy=l1;
                                                             l1=l1+nz;
                                                             gx(i,j,k)=gx(i,j,k)-(0.0*psi_n(k2,l,l1)+phi_n(k2,l,l1))*potential(i,j,k)*(i-k2)*delx/(sigma*sigma);   
                                                         gy(i,j,k)=gy(i,j,k)-(0.0*psi_n(k2,l,l1)+phi_n(k2,l,l1))*potential(i,j,k)*(j-l)*dely/(sigma*sigma);
                                                         gz(i,j,k)=gz(i,j,k)-(0.0*psi_n(k2,l,l1)+phi_n(k2,l,l1))*potential(i,j,k)*(k-l1)*delz/(sigma*sigma);
                                                       fx(i,j,k)=fx(i,j,k)-(psi_n(k2,l,l1)+0.0*phi_n(k2,l,l1))*(i-k2)*potential(i,j,k)*delx/(sigma*sigma);    
                                                       fy(i,j,k)=fy(i,j,k)-(psi_n(k2,l,l1)+0.0*phi_n(k2,l,l1))*(j-l)*potential(i,j,k)*dely/(sigma*sigma);
                                                       fz(i,j,k)=fz(i,j,k)-(psi_n(k2,l,l1)+0.0*phi_n(k2,l,l1))*(k-l1)*potential(i,j,k)*delz/(sigma*sigma);
                                                       k2=k2_copy;
                                                       l=l_copy;
                                                       l1=l1_copy;
                                                       end
                                                       
                                                       if(k2>nx && l>ny && l1>nz)   % 3,3,3
                                                           
                                                           k2_copy=k2;
                                                             k2=k2-nz;
                                                             l_copy=l;
                                                             l=l-nz;
                                                             l1_copy=l1;
                                                             l1=l1-nz;
                                                             gx(i,j,k)=gx(i,j,k)-(0.0*psi_n(k2,l,l1)+phi_n(k2,l,l1))*potential(i,j,k)*(i-k2)*delx/(sigma*sigma);   
                                                         gy(i,j,k)=gy(i,j,k)-(0.0*psi_n(k2,l,l1)+phi_n(k2,l,l1))*potential(i,j,k)*(j-l)*dely/(sigma*sigma);
                                                         gz(i,j,k)=gz(i,j,k)-(0.0*psi_n(k2,l,l1)+phi_n(k2,l,l1))*potential(i,j,k)*(k-l1)*delz/(sigma*sigma);
                                                       fx(i,j,k)=fx(i,j,k)-(psi_n(k2,l,l1)+0.0*phi_n(k2,l,l1))*(i-k2)*potential(i,j,k)*delx/(sigma*sigma);    
                                                       fy(i,j,k)=fy(i,j,k)-(psi_n(k2,l,l1)+0.0*phi_n(k2,l,l1))*(j-l)*potential(i,j,k)*dely/(sigma*sigma);
                                                       fz(i,j,k)=fz(i,j,k)-(psi_n(k2,l,l1)+0.0*phi_n(k2,l,l1))*(k-l1)*potential(i,j,k)*delz/(sigma*sigma);
                                                       k2=k2_copy;
                                                       l=l_copy;
                                                       l1=l1_copy;
                                                       
                                                       end
                                                       
                                                      
                                               
                                            
                                            
                                        end
                                        
                                 end
                            end
                            
                            
                            if(i>1 && i<nx && j>1 && j<ny && k>1 && k<nz)  %1,1,1
                                    x1=i;
                                    x2=i+1;
                                    x3=i-1;
                                    
                                    y1=j;
                                    y2=j+1;
                                    y3=j-1;
                                    
                                    z1=k;
                                    z2=k+1;
                                    z3=k-1;
                                    
                                    phi(x1,y1,z1) = phi_n(x1,y1,z1) + ff * ((phi_n(x2,y1,z1)+phi_n(x3,y1,z1)+phi_n(x1,y2,z1)+phi_n(x1,y3,z1)+phi_n(x1,y1,z2)+phi_n(x1,y1,z3))-6*phi_n(x1,y1,z1))+delt*k_a*phi_n(x1,y1,z1)*(k_b/k_a-phi_n(x1,y1,z1))+b_noise(x1,y1,z1)+(1/(delx))*(eta_phi(x2,y1,z1)-eta_phi(x1,y1,z1)+eta_phi(x1,y2,z1)-eta_phi(x1,y1,z1)+eta_phi(x1,y1,z2)-eta_phi(x1,y1,z1)) +(delt/delx)*(phi_n(x2,y1,z1)-phi_n(x1,y1,z1))*gx(x1,y1,z1)+(delt/dely)*(phi_n(x1,y2,z1)-phi_n(x1,y1,z1))*gy(x1,y1,z1)+(delt/delz)*(phi_n(x1,y1,z2)-phi_n(x1,y1,z1))*gz(x1,y1,z1)+(delt/delx)*phi_n(x1,y1,z1)*(gx(x2,y1,z1)-gx(x1,y1,z1))+(delt/delz)*phi_n(x1,y1,z1)*(gz(x1,y1,z2)-gz(x1,y1,z1)) + (delt/dely)*phi_n(x1,y1,z1)*(gy(x1,y2,z1)-gy(x1,y1,z1));
                                    psi(x1,y1,z1) = psi_n(x1,y1,z1) + ff * ((psi_n(x2,y1,z1)+psi_n(x3,y1,z1)+psi_n(x1,y2,z1)+psi_n(x1,y3,z1)+phi_n(x1,y1,z2)+phi_n(x1,y1,z3))-6*psi_n(x1,y1,z1))+(1/delx)*(eta(x2,y1,z1)-eta(x1,y1,z1)+eta(x1,y2,z1)-eta(x1,y1,z1)+eta(x1,y1,z2)-eta(x1,y1,z1)) + (delt/delx)*(psi_n(x2,y1,z1)-psi_n(x1,y1,z1))*fx(x1,y1,z1)+(delt/dely)*(psi_n(x1,y2,z1)-psi_n(x1,y1,z1))*fy(x1,y1,z1)+(delt/dely)*(psi_n(x1,y1,z2)-psi_n(x1,y1,z1))*fz(x1,y1,z1)+(delt/delx)*psi_n(x1,y1,z1)*(fx(x2,y1,z1)-fx(x1,y1,z1))+(delt/delz)*psi_n(x1,y1,z1)*(fz(x1,y1,z2)-fz(x1,y1,z1)) + (delt/dely)*psi_n(x1,y1,z1)*(fy(x1,y2,z1)-fy(x1,y1,z1));
                                    
                                    
                                    if(phi(x1,y1,z1)<0)
                                        phi(x1,y1,z1)= phi_n(x1,y1,z1);
                                        %psi(x1,y1,z1) = psi_n(x1,y1,z1);
                                    end
                                    if(psi(x1,y1,z1)<0)
                                        psi(x1,y1,z1)= psi_n(x1,y1,z1);
                                        %psi(x1,y1,z1) = psi_n(x1,y1,z1);
                                    end
                                    
                                 %   phi(i,j,k) = phi_n(i,j,k) + ff * ((phi_n(i+1,j,k)+phi_n(i-1,j,k)+phi_n(i,j+1,k)+phi_n(i,j-1,k)+phi_n(i,j,k+1)+phi_n(i,j,k-1))-6*phi_n(i,j,k))+delt* k_a*phi_n(i,j,k)*(k_b/k_a-phi_n(i,j,k))+delt*b_noise*randn + (delt/delx)*(phi_n(i+1,j,k)-phi_n(i,j,k))*gx(i,j,k)+(delt/dely)*(phi_n(i,j+1,k)-phi_n(i,j,k))*gy(i,j,k)+(delt/delz)*(phi_n(i,j,k+1)-phi_n(i,j,k))*gz(i,j,k)+(delt/delx)*phi_n(i,j,k)*(gx(i+1,j,k)-gx(i,j,k))+(delt/delz)*phi_n(i,j,k)*(gz(i,j,k+1)-gz(i,j,k)) + (delt/dely)*phi_n(i,j,k)*(gy(i,j+1,k)-gy(i,j,k));
                                  %  psi(i,j,k) = psi_n(i,j,k) + ff * ((psi_n(i+1,j,k)+psi_n(i-1,j,k)+psi_n(i,j+1,k)+psi_n(i,j-1,k)+phi_n(i,j,k+1)+phi_n(i,j,k-1))-6*psi_n(i,j,k))+(delt/delx)*(eta(i+1,j,k)-eta(i,j,k)+eta(i,j+1,k)-eta(i,j,k)+eta(i,j,k+1)-eta(i,j,k)) + (delt/delx)*(psi_n(i+1,j,k)-psi_n(i,j,k))*fx(i,j,k)+(delt/dely)*(psi_n(i,j+1,k)-psi_n(i,j,k))*fy(i,j,k)+(delt/dely)*(psi_n(i,j,k+1)-psi_n(i,j,k))*fz(i,j,k)+(delt/delx)*psi_n(i,j,k)*(fx(i+1,j,k)-fx(i,j,k))+(delt/delz)*psi_n(i,j,k)*(fz(i,j,k+1)-fz(i,j,k)) + (delt/dely)*psi_n(i,j,k)*(fy(i,j+1,k)-fy(i,j,k));
                            end

                            
                           if(i==1 && j>1 && j<ny && k>1 && k<nz) %2,1,1
                               
                                    x1=i;
                                    x2=i+1;
                                    x3=nx;
                                    
                                    y1=j;
                                    y2=j+1;
                                    y3=j-1;
                                    
                                    z1=k;
                                    z2=k+1;
                                    z3=k-1;
                                    
                                    phi(x1,y1,z1) = phi_n(x1,y1,z1) + ff * ((phi_n(x2,y1,z1)+phi_n(x3,y1,z1)+phi_n(x1,y2,z1)+phi_n(x1,y3,z1)+phi_n(x1,y1,z2)+phi_n(x1,y1,z3))-6*phi_n(x1,y1,z1))+delt*k_a*phi_n(x1,y1,z1)*(k_b/k_a-phi_n(x1,y1,z1))+b_noise(x1,y1,z1)+(1/(delx))*(eta_phi(x2,y1,z1)-eta_phi(x1,y1,z1)+eta_phi(x1,y2,z1)-eta_phi(x1,y1,z1)+eta_phi(x1,y1,z2)-eta_phi(x1,y1,z1)) + (delt/delx)*(phi_n(x2,y1,z1)-phi_n(x1,y1,z1))*gx(x1,y1,z1)+(delt/dely)*(phi_n(x1,y2,z1)-phi_n(x1,y1,z1))*gy(x1,y1,z1)+(delt/delz)*(phi_n(x1,y1,z2)-phi_n(x1,y1,z1))*gz(x1,y1,z1)+(delt/delx)*phi_n(x1,y1,z1)*(gx(x2,y1,z1)-gx(x1,y1,z1))+(delt/delz)*phi_n(x1,y1,z1)*(gz(x1,y1,z2)-gz(x1,y1,z1)) + (delt/dely)*phi_n(x1,y1,z1)*(gy(x1,y2,z1)-gy(x1,y1,z1));
                                    psi(x1,y1,z1) = psi_n(x1,y1,z1) + ff * ((psi_n(x2,y1,z1)+psi_n(x3,y1,z1)+psi_n(x1,y2,z1)+psi_n(x1,y3,z1)+phi_n(x1,y1,z2)+phi_n(x1,y1,z3))-6*psi_n(x1,y1,z1))+(1/delx)*(eta(x2,y1,z1)-eta(x1,y1,z1)+eta(x1,y2,z1)-eta(x1,y1,z1)+eta(x1,y1,z2)-eta(x1,y1,z1)) + (delt/delx)*(psi_n(x2,y1,z1)-psi_n(x1,y1,z1))*fx(x1,y1,z1)+(delt/dely)*(psi_n(x1,y2,z1)-psi_n(x1,y1,z1))*fy(x1,y1,z1)+(delt/dely)*(psi_n(x1,y1,z2)-psi_n(x1,y1,z1))*fz(x1,y1,z1)+(delt/delx)*psi_n(x1,y1,z1)*(fx(x2,y1,z1)-fx(x1,y1,z1))+(delt/delz)*psi_n(x1,y1,z1)*(fz(x1,y1,z2)-fz(x1,y1,z1)) + (delt/dely)*psi_n(x1,y1,z1)*(fy(x1,y2,z1)-fy(x1,y1,z1));
                                   
                               if(phi(x1,y1,z1)<0)
                                        phi(x1,y1,z1)= phi_n(x1,y1,z1);
                                        %psi(x1,y1,z1) = psi_n(x1,y1,z1);
                               end
                                    if(psi(x1,y1,z1)<0)
                                        psi(x1,y1,z1)= psi_n(x1,y1,z1);
                                        %psi(x1,y1,z1) = psi_n(x1,y1,z1);
                                    end
                               
                           end
                           
                           if(i==nx && j>1 && j<ny && k>1 && k<nz) %3,1,1
                                    x1=i;
                                    x2=1;
                                    x3=nx;
                                    
                                    y1=j;
                                    y2=j+1;
                                    y3=j-1;
                                    
                                    z1=k;
                                    z2=k+1;
                                    z3=k-1;
                                    
                                    phi(x1,y1,z1) = phi_n(x1,y1,z1) + ff * ((phi_n(x2,y1,z1)+phi_n(x3,y1,z1)+phi_n(x1,y2,z1)+phi_n(x1,y3,z1)+phi_n(x1,y1,z2)+phi_n(x1,y1,z3))-6*phi_n(x1,y1,z1))+delt*k_a*phi_n(x1,y1,z1)*(k_b/k_a-phi_n(x1,y1,z1))+b_noise(x1,y1,z1)+(1/(delx))*(eta_phi(x2,y1,z1)-eta_phi(x1,y1,z1)+eta_phi(x1,y2,z1)-eta_phi(x1,y1,z1)+eta_phi(x1,y1,z2)-eta_phi(x1,y1,z1)) + (delt/delx)*(phi_n(x2,y1,z1)-phi_n(x1,y1,z1))*gx(x1,y1,z1)+(delt/dely)*(phi_n(x1,y2,z1)-phi_n(x1,y1,z1))*gy(x1,y1,z1)+(delt/delz)*(phi_n(x1,y1,z2)-phi_n(x1,y1,z1))*gz(x1,y1,z1)+(delt/delx)*phi_n(x1,y1,z1)*(gx(x2,y1,z1)-gx(x1,y1,z1))+(delt/delz)*phi_n(x1,y1,z1)*(gz(x1,y1,z2)-gz(x1,y1,z1)) + (delt/dely)*phi_n(x1,y1,z1)*(gy(x1,y2,z1)-gy(x1,y1,z1));
                                    psi(x1,y1,z1) = psi_n(x1,y1,z1) + ff * ((psi_n(x2,y1,z1)+psi_n(x3,y1,z1)+psi_n(x1,y2,z1)+psi_n(x1,y3,z1)+phi_n(x1,y1,z2)+phi_n(x1,y1,z3))-6*psi_n(x1,y1,z1))+(1/delx)*(eta(x2,y1,z1)-eta(x1,y1,z1)+eta(x1,y2,z1)-eta(x1,y1,z1)+eta(x1,y1,z2)-eta(x1,y1,z1)) + (delt/delx)*(psi_n(x2,y1,z1)-psi_n(x1,y1,z1))*fx(x1,y1,z1)+(delt/dely)*(psi_n(x1,y2,z1)-psi_n(x1,y1,z1))*fy(x1,y1,z1)+(delt/dely)*(psi_n(x1,y1,z2)-psi_n(x1,y1,z1))*fz(x1,y1,z1)+(delt/delx)*psi_n(x1,y1,z1)*(fx(x2,y1,z1)-fx(x1,y1,z1))+(delt/delz)*psi_n(x1,y1,z1)*(fz(x1,y1,z2)-fz(x1,y1,z1)) + (delt/dely)*psi_n(x1,y1,z1)*(fy(x1,y2,z1)-fy(x1,y1,z1));
                                   
                               if(phi(x1,y1,z1)<0)
                                        phi(x1,y1,z1)= phi_n(x1,y1,z1);
                                        %psi(x1,y1,z1) = psi_n(x1,y1,z1);
                               end
                               if(psi(x1,y1,z1)<0)
                                        psi(x1,y1,z1)= psi_n(x1,y1,z1);
                                        %psi(x1,y1,z1) = psi_n(x1,y1,z1);
                               end
                           end
                           if(i>1 && i<nx && j==1  && k>1 && k<nz) % 1,2,1
                                    x1=i;
                                    x2=i+1;
                                    x3=i-1;
                                    
                                    y1=j;
                                    y2=j+1;
                                    y3=ny;
                                    
                                    z1=k;
                                    z2=k+1;
                                    z3=k-1;
                                    
                                    phi(x1,y1,z1) = phi_n(x1,y1,z1) + ff * ((phi_n(x2,y1,z1)+phi_n(x3,y1,z1)+phi_n(x1,y2,z1)+phi_n(x1,y3,z1)+phi_n(x1,y1,z2)+phi_n(x1,y1,z3))-6*phi_n(x1,y1,z1))+delt*k_a*phi_n(x1,y1,z1)*(k_b/k_a-phi_n(x1,y1,z1))+b_noise(x1,y1,z1)+(1/(delx))*(eta_phi(x2,y1,z1)-eta_phi(x1,y1,z1)+eta_phi(x1,y2,z1)-eta_phi(x1,y1,z1)+eta_phi(x1,y1,z2)-eta_phi(x1,y1,z1)) + (delt/delx)*(phi_n(x2,y1,z1)-phi_n(x1,y1,z1))*gx(x1,y1,z1)+(delt/dely)*(phi_n(x1,y2,z1)-phi_n(x1,y1,z1))*gy(x1,y1,z1)+(delt/delz)*(phi_n(x1,y1,z2)-phi_n(x1,y1,z1))*gz(x1,y1,z1)+(delt/delx)*phi_n(x1,y1,z1)*(gx(x2,y1,z1)-gx(x1,y1,z1))+(delt/delz)*phi_n(x1,y1,z1)*(gz(x1,y1,z2)-gz(x1,y1,z1)) + (delt/dely)*phi_n(x1,y1,z1)*(gy(x1,y2,z1)-gy(x1,y1,z1));
                                    psi(x1,y1,z1) = psi_n(x1,y1,z1) + ff * ((psi_n(x2,y1,z1)+psi_n(x3,y1,z1)+psi_n(x1,y2,z1)+psi_n(x1,y3,z1)+phi_n(x1,y1,z2)+phi_n(x1,y1,z3))-6*psi_n(x1,y1,z1))+(1/delx)*(eta(x2,y1,z1)-eta(x1,y1,z1)+eta(x1,y2,z1)-eta(x1,y1,z1)+eta(x1,y1,z2)-eta(x1,y1,z1)) + (delt/delx)*(psi_n(x2,y1,z1)-psi_n(x1,y1,z1))*fx(x1,y1,z1)+(delt/dely)*(psi_n(x1,y2,z1)-psi_n(x1,y1,z1))*fy(x1,y1,z1)+(delt/dely)*(psi_n(x1,y1,z2)-psi_n(x1,y1,z1))*fz(x1,y1,z1)+(delt/delx)*psi_n(x1,y1,z1)*(fx(x2,y1,z1)-fx(x1,y1,z1))+(delt/delz)*psi_n(x1,y1,z1)*(fz(x1,y1,z2)-fz(x1,y1,z1)) + (delt/dely)*psi_n(x1,y1,z1)*(fy(x1,y2,z1)-fy(x1,y1,z1));
                                    
                                    if(phi(x1,y1,z1)<0)
                                        phi(x1,y1,z1)= phi_n(x1,y1,z1);
                                        %psi(x1,y1,z1) = psi_n(x1,y1,z1);
                                    end
                                    if(psi(x1,y1,z1)<0)
                                        psi(x1,y1,z1)= psi_n(x1,y1,z1);
                                        %psi(x1,y1,z1) = psi_n(x1,y1,z1);
                                    end
                                    
                                    
                           end
                               
                               
                           if(i>1 && k==1 && i<nx && j>1 && j<ny) % 1,1,2 
                                    x1=i;
                                    x2=i+1;
                                    x3=i-1;
                                    
                                    y1=j;
                                    y2=j+1;
                                    y3=j-1;
                                    
                                    z1=k;
                                    z2=k+1;
                                    z3=nz;
                                    
                                    phi(x1,y1,z1) = phi_n(x1,y1,z1) + ff * ((phi_n(x2,y1,z1)+phi_n(x3,y1,z1)+phi_n(x1,y2,z1)+phi_n(x1,y3,z1)+phi_n(x1,y1,z2)+phi_n(x1,y1,z3))-6*phi_n(x1,y1,z1))+delt*k_a*phi_n(x1,y1,z1)*(k_b/k_a-phi_n(x1,y1,z1))+b_noise(x1,y1,z1)+(1/(delx))*(eta_phi(x2,y1,z1)-eta_phi(x1,y1,z1)+eta_phi(x1,y2,z1)-eta_phi(x1,y1,z1)+eta_phi(x1,y1,z2)-eta_phi(x1,y1,z1)) + (delt/delx)*(phi_n(x2,y1,z1)-phi_n(x1,y1,z1))*gx(x1,y1,z1)+(delt/dely)*(phi_n(x1,y2,z1)-phi_n(x1,y1,z1))*gy(x1,y1,z1)+(delt/delz)*(phi_n(x1,y1,z2)-phi_n(x1,y1,z1))*gz(x1,y1,z1)+(delt/delx)*phi_n(x1,y1,z1)*(gx(x2,y1,z1)-gx(x1,y1,z1))+(delt/delz)*phi_n(x1,y1,z1)*(gz(x1,y1,z2)-gz(x1,y1,z1)) + (delt/dely)*phi_n(x1,y1,z1)*(gy(x1,y2,z1)-gy(x1,y1,z1));
                                    psi(x1,y1,z1) = psi_n(x1,y1,z1) + ff * ((psi_n(x2,y1,z1)+psi_n(x3,y1,z1)+psi_n(x1,y2,z1)+psi_n(x1,y3,z1)+phi_n(x1,y1,z2)+phi_n(x1,y1,z3))-6*psi_n(x1,y1,z1))+(1/delx)*(eta(x2,y1,z1)-eta(x1,y1,z1)+eta(x1,y2,z1)-eta(x1,y1,z1)+eta(x1,y1,z2)-eta(x1,y1,z1)) + (delt/delx)*(psi_n(x2,y1,z1)-psi_n(x1,y1,z1))*fx(x1,y1,z1)+(delt/dely)*(psi_n(x1,y2,z1)-psi_n(x1,y1,z1))*fy(x1,y1,z1)+(delt/dely)*(psi_n(x1,y1,z2)-psi_n(x1,y1,z1))*fz(x1,y1,z1)+(delt/delx)*psi_n(x1,y1,z1)*(fx(x2,y1,z1)-fx(x1,y1,z1))+(delt/delz)*psi_n(x1,y1,z1)*(fz(x1,y1,z2)-fz(x1,y1,z1)) + (delt/dely)*psi_n(x1,y1,z1)*(fy(x1,y2,z1)-fy(x1,y1,z1));
                                    
                                    
                                    if(phi(x1,y1,z1)<0)
                                        phi(x1,y1,z1)= phi_n(x1,y1,z1);
                                        %psi(x1,y1,z1) = psi_n(x1,y1,z1);
                                    end
                                    if(psi(x1,y1,z1)<0)
                                        psi(x1,y1,z1)= psi_n(x1,y1,z1);
                                        %psi(x1,y1,z1) = psi_n(x1,y1,z1);
                                    end
                                    
                                    
                           end
                               
                               
                               
                           if(i>1 && k==nz && i<nx && j>1 && j<ny) % 1, 1,3
                                    x1=i;
                                    x2=i+1;
                                    x3=i-1;
                                    
                                    y1=j;
                                    y2=j+1;
                                    y3=j-1;
                                    
                                    z1=k;
                                    z2=1;
                                    z3=k-1;
                                    
                                    phi(x1,y1,z1) = phi_n(x1,y1,z1) + ff * ((phi_n(x2,y1,z1)+phi_n(x3,y1,z1)+phi_n(x1,y2,z1)+phi_n(x1,y3,z1)+phi_n(x1,y1,z2)+phi_n(x1,y1,z3))-6*phi_n(x1,y1,z1))+delt*k_a*phi_n(x1,y1,z1)*(k_b/k_a-phi_n(x1,y1,z1))+b_noise(x1,y1,z1)+(1/(delx))*(eta_phi(x2,y1,z1)-eta_phi(x1,y1,z1)+eta_phi(x1,y2,z1)-eta_phi(x1,y1,z1)+eta_phi(x1,y1,z2)-eta_phi(x1,y1,z1)) + (delt/delx)*(phi_n(x2,y1,z1)-phi_n(x1,y1,z1))*gx(x1,y1,z1)+(delt/dely)*(phi_n(x1,y2,z1)-phi_n(x1,y1,z1))*gy(x1,y1,z1)+(delt/delz)*(phi_n(x1,y1,z2)-phi_n(x1,y1,z1))*gz(x1,y1,z1)+(delt/delx)*phi_n(x1,y1,z1)*(gx(x2,y1,z1)-gx(x1,y1,z1))+(delt/delz)*phi_n(x1,y1,z1)*(gz(x1,y1,z2)-gz(x1,y1,z1)) + (delt/dely)*phi_n(x1,y1,z1)*(gy(x1,y2,z1)-gy(x1,y1,z1));
                                    psi(x1,y1,z1) = psi_n(x1,y1,z1) + ff * ((psi_n(x2,y1,z1)+psi_n(x3,y1,z1)+psi_n(x1,y2,z1)+psi_n(x1,y3,z1)+phi_n(x1,y1,z2)+phi_n(x1,y1,z3))-6*psi_n(x1,y1,z1))+(1/delx)*(eta(x2,y1,z1)-eta(x1,y1,z1)+eta(x1,y2,z1)-eta(x1,y1,z1)+eta(x1,y1,z2)-eta(x1,y1,z1)) + (delt/delx)*(psi_n(x2,y1,z1)-psi_n(x1,y1,z1))*fx(x1,y1,z1)+(delt/dely)*(psi_n(x1,y2,z1)-psi_n(x1,y1,z1))*fy(x1,y1,z1)+(delt/dely)*(psi_n(x1,y1,z2)-psi_n(x1,y1,z1))*fz(x1,y1,z1)+(delt/delx)*psi_n(x1,y1,z1)*(fx(x2,y1,z1)-fx(x1,y1,z1))+(delt/delz)*psi_n(x1,y1,z1)*(fz(x1,y1,z2)-fz(x1,y1,z1)) + (delt/dely)*psi_n(x1,y1,z1)*(fy(x1,y2,z1)-fy(x1,y1,z1));
                                    
                                    
                                    if(phi(x1,y1,z1)<0)
                                        phi(x1,y1,z1)= phi_n(x1,y1,z1);
                                        %psi(x1,y1,z1) = psi_n(x1,y1,z1);
                                    end
                                    if(psi(x1,y1,z1)<0)
                                        psi(x1,y1,z1)= psi_n(x1,y1,z1);
                                        %psi(x1,y1,z1) = psi_n(x1,y1,z1);
                                    end
                             
                           end
                           
                           if(i>1 && j==ny && i<nx && k>1 && k<nz) % 1, 3, 1
                                    x1=i;
                                    x2=i+1;
                                    x3=i-1;
                                    
                                    y1=j;
                                    y2=1;
                                    y3=j-1;
                                    
                                    z1=k;
                                    z2=k+1;
                                    z3=k-1;
                                    
                                    phi(x1,y1,z1) = phi_n(x1,y1,z1) + ff * ((phi_n(x2,y1,z1)+phi_n(x3,y1,z1)+phi_n(x1,y2,z1)+phi_n(x1,y3,z1)+phi_n(x1,y1,z2)+phi_n(x1,y1,z3))-6*phi_n(x1,y1,z1))+delt*k_a*phi_n(x1,y1,z1)*(k_b/k_a-phi_n(x1,y1,z1))+b_noise(x1,y1,z1)+(1/(delx))*(eta_phi(x2,y1,z1)-eta_phi(x1,y1,z1)+eta_phi(x1,y2,z1)-eta_phi(x1,y1,z1)+eta_phi(x1,y1,z2)-eta_phi(x1,y1,z1)) + (delt/delx)*(phi_n(x2,y1,z1)-phi_n(x1,y1,z1))*gx(x1,y1,z1)+(delt/dely)*(phi_n(x1,y2,z1)-phi_n(x1,y1,z1))*gy(x1,y1,z1)+(delt/delz)*(phi_n(x1,y1,z2)-phi_n(x1,y1,z1))*gz(x1,y1,z1)+(delt/delx)*phi_n(x1,y1,z1)*(gx(x2,y1,z1)-gx(x1,y1,z1))+(delt/delz)*phi_n(x1,y1,z1)*(gz(x1,y1,z2)-gz(x1,y1,z1)) + (delt/dely)*phi_n(x1,y1,z1)*(gy(x1,y2,z1)-gy(x1,y1,z1));
                                    psi(x1,y1,z1) = psi_n(x1,y1,z1) + ff * ((psi_n(x2,y1,z1)+psi_n(x3,y1,z1)+psi_n(x1,y2,z1)+psi_n(x1,y3,z1)+phi_n(x1,y1,z2)+phi_n(x1,y1,z3))-6*psi_n(x1,y1,z1))+(1/delx)*(eta(x2,y1,z1)-eta(x1,y1,z1)+eta(x1,y2,z1)-eta(x1,y1,z1)+eta(x1,y1,z2)-eta(x1,y1,z1)) + (delt/delx)*(psi_n(x2,y1,z1)-psi_n(x1,y1,z1))*fx(x1,y1,z1)+(delt/dely)*(psi_n(x1,y2,z1)-psi_n(x1,y1,z1))*fy(x1,y1,z1)+(delt/dely)*(psi_n(x1,y1,z2)-psi_n(x1,y1,z1))*fz(x1,y1,z1)+(delt/delx)*psi_n(x1,y1,z1)*(fx(x2,y1,z1)-fx(x1,y1,z1))+(delt/delz)*psi_n(x1,y1,z1)*(fz(x1,y1,z2)-fz(x1,y1,z1)) + (delt/dely)*psi_n(x1,y1,z1)*(fy(x1,y2,z1)-fy(x1,y1,z1));
                                    
                                    
                                    if(phi(x1,y1,z1)<0)
                                        phi(x1,y1,z1)= phi_n(x1,y1,z1);
                                        %psi(x1,y1,z1) = psi_n(x1,y1,z1);
                                    end
                                    if(psi(x1,y1,z1)<0)
                                        psi(x1,y1,z1)= psi_n(x1,y1,z1);
                                        %psi(x1,y1,z1) = psi_n(x1,y1,z1);
                                    end
                                    
                                    
                           end
                           
                           if(i==1 && j==1 && k>1 && k<nz) % 2,2,1
                                    x1=i;
                                    x2=i+1;
                                    x3=nx;
                                    
                                    y1=j;
                                    y2=j+1;
                                    y3=ny;
                                    
                                    z1=k;
                                    z2=k+1;
                                    z3=k-1;
                                    
                                    phi(x1,y1,z1) = phi_n(x1,y1,z1) + ff * ((phi_n(x2,y1,z1)+phi_n(x3,y1,z1)+phi_n(x1,y2,z1)+phi_n(x1,y3,z1)+phi_n(x1,y1,z2)+phi_n(x1,y1,z3))-6*phi_n(x1,y1,z1))+delt*k_a*phi_n(x1,y1,z1)*(k_b/k_a-phi_n(x1,y1,z1))+b_noise(x1,y1,z1)+(1/(delx))*(eta_phi(x2,y1,z1)-eta_phi(x1,y1,z1)+eta_phi(x1,y2,z1)-eta_phi(x1,y1,z1)+eta_phi(x1,y1,z2)-eta_phi(x1,y1,z1)) + (delt/delx)*(phi_n(x2,y1,z1)-phi_n(x1,y1,z1))*gx(x1,y1,z1)+(delt/dely)*(phi_n(x1,y2,z1)-phi_n(x1,y1,z1))*gy(x1,y1,z1)+(delt/delz)*(phi_n(x1,y1,z2)-phi_n(x1,y1,z1))*gz(x1,y1,z1)+(delt/delx)*phi_n(x1,y1,z1)*(gx(x2,y1,z1)-gx(x1,y1,z1))+(delt/delz)*phi_n(x1,y1,z1)*(gz(x1,y1,z2)-gz(x1,y1,z1)) + (delt/dely)*phi_n(x1,y1,z1)*(gy(x1,y2,z1)-gy(x1,y1,z1));
                                    psi(x1,y1,z1) = psi_n(x1,y1,z1) + ff * ((psi_n(x2,y1,z1)+psi_n(x3,y1,z1)+psi_n(x1,y2,z1)+psi_n(x1,y3,z1)+phi_n(x1,y1,z2)+phi_n(x1,y1,z3))-6*psi_n(x1,y1,z1))+(1/delx)*(eta(x2,y1,z1)-eta(x1,y1,z1)+eta(x1,y2,z1)-eta(x1,y1,z1)+eta(x1,y1,z2)-eta(x1,y1,z1)) + (delt/delx)*(psi_n(x2,y1,z1)-psi_n(x1,y1,z1))*fx(x1,y1,z1)+(delt/dely)*(psi_n(x1,y2,z1)-psi_n(x1,y1,z1))*fy(x1,y1,z1)+(delt/dely)*(psi_n(x1,y1,z2)-psi_n(x1,y1,z1))*fz(x1,y1,z1)+(delt/delx)*psi_n(x1,y1,z1)*(fx(x2,y1,z1)-fx(x1,y1,z1))+(delt/delz)*psi_n(x1,y1,z1)*(fz(x1,y1,z2)-fz(x1,y1,z1)) + (delt/dely)*psi_n(x1,y1,z1)*(fy(x1,y2,z1)-fy(x1,y1,z1));
                                    
                                    if(phi(x1,y1,z1)<0)
                                        phi(x1,y1,z1)= phi_n(x1,y1,z1);
                                        %psi(x1,y1,z1) = psi_n(x1,y1,z1);
                                    end
                                    if(psi(x1,y1,z1)<0)
                                        psi(x1,y1,z1)= psi_n(x1,y1,z1);
                                        %psi(x1,y1,z1) = psi_n(x1,y1,z1);
                                    end
                                    
                           end
                           
                           if(i==nx && j==1 && k>1 && k<nz) % 3,2,1
                                    x1=i;
                                    x2=1;
                                    x3=i-1;
                                    
                                    y1=j;
                                    y2=j+1;
                                    y3=ny;
                                    
                                    z1=k;
                                    z2=k+1;
                                    z3=k-1;
                                    
                                    phi(x1,y1,z1) = phi_n(x1,y1,z1) + ff * ((phi_n(x2,y1,z1)+phi_n(x3,y1,z1)+phi_n(x1,y2,z1)+phi_n(x1,y3,z1)+phi_n(x1,y1,z2)+phi_n(x1,y1,z3))-6*phi_n(x1,y1,z1))+delt*k_a*phi_n(x1,y1,z1)*(k_b/k_a-phi_n(x1,y1,z1))+b_noise(x1,y1,z1)+(1/(delx))*(eta_phi(x2,y1,z1)-eta_phi(x1,y1,z1)+eta_phi(x1,y2,z1)-eta_phi(x1,y1,z1)+eta_phi(x1,y1,z2)-eta_phi(x1,y1,z1)) + (delt/delx)*(phi_n(x2,y1,z1)-phi_n(x1,y1,z1))*gx(x1,y1,z1)+(delt/dely)*(phi_n(x1,y2,z1)-phi_n(x1,y1,z1))*gy(x1,y1,z1)+(delt/delz)*(phi_n(x1,y1,z2)-phi_n(x1,y1,z1))*gz(x1,y1,z1)+(delt/delx)*phi_n(x1,y1,z1)*(gx(x2,y1,z1)-gx(x1,y1,z1))+(delt/delz)*phi_n(x1,y1,z1)*(gz(x1,y1,z2)-gz(x1,y1,z1)) + (delt/dely)*phi_n(x1,y1,z1)*(gy(x1,y2,z1)-gy(x1,y1,z1));
                                    psi(x1,y1,z1) = psi_n(x1,y1,z1) + ff * ((psi_n(x2,y1,z1)+psi_n(x3,y1,z1)+psi_n(x1,y2,z1)+psi_n(x1,y3,z1)+phi_n(x1,y1,z2)+phi_n(x1,y1,z3))-6*psi_n(x1,y1,z1))+(1/delx)*(eta(x2,y1,z1)-eta(x1,y1,z1)+eta(x1,y2,z1)-eta(x1,y1,z1)+eta(x1,y1,z2)-eta(x1,y1,z1)) + (delt/delx)*(psi_n(x2,y1,z1)-psi_n(x1,y1,z1))*fx(x1,y1,z1)+(delt/dely)*(psi_n(x1,y2,z1)-psi_n(x1,y1,z1))*fy(x1,y1,z1)+(delt/dely)*(psi_n(x1,y1,z2)-psi_n(x1,y1,z1))*fz(x1,y1,z1)+(delt/delx)*psi_n(x1,y1,z1)*(fx(x2,y1,z1)-fx(x1,y1,z1))+(delt/delz)*psi_n(x1,y1,z1)*(fz(x1,y1,z2)-fz(x1,y1,z1)) + (delt/dely)*psi_n(x1,y1,z1)*(fy(x1,y2,z1)-fy(x1,y1,z1));
                                    
                                    if(phi(x1,y1,z1)<0)
                                        phi(x1,y1,z1)= phi_n(x1,y1,z1);
                                        %psi(x1,y1,z1) = psi_n(x1,y1,z1);
                                    end
                                    if(psi(x1,y1,z1)<0)
                                        psi(x1,y1,z1)= psi_n(x1,y1,z1);
                                        %psi(x1,y1,z1) = psi_n(x1,y1,z1);
                                    end
                                    
                                    
                           end
                           
                           if(i==1 && j==ny && k>1 && k<nz) %2,3,1
                                    x1=i;
                                    x2=i+1;
                                    x3=nx;
                                    
                                    y1=j;
                                    y2=1;
                                    y3=j-1;
                                    
                                    z1=k;
                                    z2=k+1;
                                    z3=k-1;
                                    
                                    phi(x1,y1,z1) = phi_n(x1,y1,z1) + ff * ((phi_n(x2,y1,z1)+phi_n(x3,y1,z1)+phi_n(x1,y2,z1)+phi_n(x1,y3,z1)+phi_n(x1,y1,z2)+phi_n(x1,y1,z3))-6*phi_n(x1,y1,z1))+delt*k_a*phi_n(x1,y1,z1)*(k_b/k_a-phi_n(x1,y1,z1))+b_noise(x1,y1,z1)+(1/(delx))*(eta_phi(x2,y1,z1)-eta_phi(x1,y1,z1)+eta_phi(x1,y2,z1)-eta_phi(x1,y1,z1)+eta_phi(x1,y1,z2)-eta_phi(x1,y1,z1)) + (delt/delx)*(phi_n(x2,y1,z1)-phi_n(x1,y1,z1))*gx(x1,y1,z1)+(delt/dely)*(phi_n(x1,y2,z1)-phi_n(x1,y1,z1))*gy(x1,y1,z1)+(delt/delz)*(phi_n(x1,y1,z2)-phi_n(x1,y1,z1))*gz(x1,y1,z1)+(delt/delx)*phi_n(x1,y1,z1)*(gx(x2,y1,z1)-gx(x1,y1,z1))+(delt/delz)*phi_n(x1,y1,z1)*(gz(x1,y1,z2)-gz(x1,y1,z1)) + (delt/dely)*phi_n(x1,y1,z1)*(gy(x1,y2,z1)-gy(x1,y1,z1));
                                    psi(x1,y1,z1) = psi_n(x1,y1,z1) + ff * ((psi_n(x2,y1,z1)+psi_n(x3,y1,z1)+psi_n(x1,y2,z1)+psi_n(x1,y3,z1)+phi_n(x1,y1,z2)+phi_n(x1,y1,z3))-6*psi_n(x1,y1,z1))+(1/delx)*(eta(x2,y1,z1)-eta(x1,y1,z1)+eta(x1,y2,z1)-eta(x1,y1,z1)+eta(x1,y1,z2)-eta(x1,y1,z1)) + (delt/delx)*(psi_n(x2,y1,z1)-psi_n(x1,y1,z1))*fx(x1,y1,z1)+(delt/dely)*(psi_n(x1,y2,z1)-psi_n(x1,y1,z1))*fy(x1,y1,z1)+(delt/dely)*(psi_n(x1,y1,z2)-psi_n(x1,y1,z1))*fz(x1,y1,z1)+(delt/delx)*psi_n(x1,y1,z1)*(fx(x2,y1,z1)-fx(x1,y1,z1))+(delt/delz)*psi_n(x1,y1,z1)*(fz(x1,y1,z2)-fz(x1,y1,z1)) + (delt/dely)*psi_n(x1,y1,z1)*(fy(x1,y2,z1)-fy(x1,y1,z1));
                                    
                                    if(phi(x1,y1,z1)<0)
                                        phi(x1,y1,z1)= phi_n(x1,y1,z1);
                                        %psi(x1,y1,z1) = psi_n(x1,y1,z1);
                                    end
                                    if(psi(x1,y1,z1)<0)
                                        psi(x1,y1,z1)= psi_n(x1,y1,z1);
                                        %psi(x1,y1,z1) = psi_n(x1,y1,z1);
                                    end
                                    
                           end
                           
                           if(i==nx && j==ny && k>1 && k<nz) %3,3,1
                                    x1=i;
                                    x2=1;
                                    x3=i-1;
                                    
                                    y1=j;
                                    y2=1;
                                    y3=j-1;
                                    
                                    z1=k;
                                    z2=k+1;
                                    z3=k-1;
                                    
                                    phi(x1,y1,z1) = phi_n(x1,y1,z1) + ff * ((phi_n(x2,y1,z1)+phi_n(x3,y1,z1)+phi_n(x1,y2,z1)+phi_n(x1,y3,z1)+phi_n(x1,y1,z2)+phi_n(x1,y1,z3))-6*phi_n(x1,y1,z1))+delt*k_a*phi_n(x1,y1,z1)*(k_b/k_a-phi_n(x1,y1,z1))+b_noise(x1,y1,z1)+(1/(delx))*(eta_phi(x2,y1,z1)-eta_phi(x1,y1,z1)+eta_phi(x1,y2,z1)-eta_phi(x1,y1,z1)+eta_phi(x1,y1,z2)-eta_phi(x1,y1,z1)) + (delt/delx)*(phi_n(x2,y1,z1)-phi_n(x1,y1,z1))*gx(x1,y1,z1)+(delt/dely)*(phi_n(x1,y2,z1)-phi_n(x1,y1,z1))*gy(x1,y1,z1)+(delt/delz)*(phi_n(x1,y1,z2)-phi_n(x1,y1,z1))*gz(x1,y1,z1)+(delt/delx)*phi_n(x1,y1,z1)*(gx(x2,y1,z1)-gx(x1,y1,z1))+(delt/delz)*phi_n(x1,y1,z1)*(gz(x1,y1,z2)-gz(x1,y1,z1)) + (delt/dely)*phi_n(x1,y1,z1)*(gy(x1,y2,z1)-gy(x1,y1,z1));
                                    psi(x1,y1,z1) = psi_n(x1,y1,z1) + ff * ((psi_n(x2,y1,z1)+psi_n(x3,y1,z1)+psi_n(x1,y2,z1)+psi_n(x1,y3,z1)+phi_n(x1,y1,z2)+phi_n(x1,y1,z3))-6*psi_n(x1,y1,z1))+(1/delx)*(eta(x2,y1,z1)-eta(x1,y1,z1)+eta(x1,y2,z1)-eta(x1,y1,z1)+eta(x1,y1,z2)-eta(x1,y1,z1)) + (delt/delx)*(psi_n(x2,y1,z1)-psi_n(x1,y1,z1))*fx(x1,y1,z1)+(delt/dely)*(psi_n(x1,y2,z1)-psi_n(x1,y1,z1))*fy(x1,y1,z1)+(delt/dely)*(psi_n(x1,y1,z2)-psi_n(x1,y1,z1))*fz(x1,y1,z1)+(delt/delx)*psi_n(x1,y1,z1)*(fx(x2,y1,z1)-fx(x1,y1,z1))+(delt/delz)*psi_n(x1,y1,z1)*(fz(x1,y1,z2)-fz(x1,y1,z1)) + (delt/dely)*psi_n(x1,y1,z1)*(fy(x1,y2,z1)-fy(x1,y1,z1));
                                    
                                    
                                    if(phi(x1,y1,z1)<0)
                                        phi(x1,y1,z1)= phi_n(x1,y1,z1);
                                        %psi(x1,y1,z1) = psi_n(x1,y1,z1);
                                    end
                                    if(psi(x1,y1,z1)<0)
                                        psi(x1,y1,z1)= psi_n(x1,y1,z1);
                                        %psi(x1,y1,z1) = psi_n(x1,y1,z1);
                                    end
                           end
                           
                           if(i==1 && k==1 && j>1 && j<ny) % 2,1,2
                                    x1=i;
                                    x2=i+1;
                                    x3=nx;
                                    
                                    y1=j;
                                    y2=j+1;
                                    y3=j-1;
                                    
                                    z1=k;
                                    z2=k+1;
                                    z3=nz;
                                    
                                    phi(x1,y1,z1) = phi_n(x1,y1,z1) + ff * ((phi_n(x2,y1,z1)+phi_n(x3,y1,z1)+phi_n(x1,y2,z1)+phi_n(x1,y3,z1)+phi_n(x1,y1,z2)+phi_n(x1,y1,z3))-6*phi_n(x1,y1,z1))+delt*k_a*phi_n(x1,y1,z1)*(k_b/k_a-phi_n(x1,y1,z1))+b_noise(x1,y1,z1)+(1/(delx))*(eta_phi(x2,y1,z1)-eta_phi(x1,y1,z1)+eta_phi(x1,y2,z1)-eta_phi(x1,y1,z1)+eta_phi(x1,y1,z2)-eta_phi(x1,y1,z1)) + (delt/delx)*(phi_n(x2,y1,z1)-phi_n(x1,y1,z1))*gx(x1,y1,z1)+(delt/dely)*(phi_n(x1,y2,z1)-phi_n(x1,y1,z1))*gy(x1,y1,z1)+(delt/delz)*(phi_n(x1,y1,z2)-phi_n(x1,y1,z1))*gz(x1,y1,z1)+(delt/delx)*phi_n(x1,y1,z1)*(gx(x2,y1,z1)-gx(x1,y1,z1))+(delt/delz)*phi_n(x1,y1,z1)*(gz(x1,y1,z2)-gz(x1,y1,z1)) + (delt/dely)*phi_n(x1,y1,z1)*(gy(x1,y2,z1)-gy(x1,y1,z1));
                                    psi(x1,y1,z1) = psi_n(x1,y1,z1) + ff * ((psi_n(x2,y1,z1)+psi_n(x3,y1,z1)+psi_n(x1,y2,z1)+psi_n(x1,y3,z1)+phi_n(x1,y1,z2)+phi_n(x1,y1,z3))-6*psi_n(x1,y1,z1))+(1/delx)*(eta(x2,y1,z1)-eta(x1,y1,z1)+eta(x1,y2,z1)-eta(x1,y1,z1)+eta(x1,y1,z2)-eta(x1,y1,z1)) + (delt/delx)*(psi_n(x2,y1,z1)-psi_n(x1,y1,z1))*fx(x1,y1,z1)+(delt/dely)*(psi_n(x1,y2,z1)-psi_n(x1,y1,z1))*fy(x1,y1,z1)+(delt/dely)*(psi_n(x1,y1,z2)-psi_n(x1,y1,z1))*fz(x1,y1,z1)+(delt/delx)*psi_n(x1,y1,z1)*(fx(x2,y1,z1)-fx(x1,y1,z1))+(delt/delz)*psi_n(x1,y1,z1)*(fz(x1,y1,z2)-fz(x1,y1,z1)) + (delt/dely)*psi_n(x1,y1,z1)*(fy(x1,y2,z1)-fy(x1,y1,z1));
                                    
                                    
                                    if(phi(x1,y1,z1)<0)
                                        phi(x1,y1,z1)= phi_n(x1,y1,z1);
                                        %psi(x1,y1,z1) = psi_n(x1,y1,z1);
                                    end
                                    if(psi(x1,y1,z1)<0)
                                        psi(x1,y1,z1)= psi_n(x1,y1,z1);
                                        %psi(x1,y1,z1) = psi_n(x1,y1,z1);
                                    end
                           end
                           
                           if(i==nx && k==1 && j>1 && j<ny) % 3,1,2
                                    x1=i;
                                    x2=1;
                                    x3=i-1;
                                    
                                    y1=j;
                                    y2=j+1;
                                    y3=j-1;
                                    
                                    z1=k;
                                    z2=k+1;
                                    z3=nz;
                                    
                                    phi(x1,y1,z1) = phi_n(x1,y1,z1) + ff * ((phi_n(x2,y1,z1)+phi_n(x3,y1,z1)+phi_n(x1,y2,z1)+phi_n(x1,y3,z1)+phi_n(x1,y1,z2)+phi_n(x1,y1,z3))-6*phi_n(x1,y1,z1))+delt*k_a*phi_n(x1,y1,z1)*(k_b/k_a-phi_n(x1,y1,z1))+b_noise(x1,y1,z1)+(1/(delx))*(eta_phi(x2,y1,z1)-eta_phi(x1,y1,z1)+eta_phi(x1,y2,z1)-eta_phi(x1,y1,z1)+eta_phi(x1,y1,z2)-eta_phi(x1,y1,z1)) + (delt/delx)*(phi_n(x2,y1,z1)-phi_n(x1,y1,z1))*gx(x1,y1,z1)+(delt/dely)*(phi_n(x1,y2,z1)-phi_n(x1,y1,z1))*gy(x1,y1,z1)+(delt/delz)*(phi_n(x1,y1,z2)-phi_n(x1,y1,z1))*gz(x1,y1,z1)+(delt/delx)*phi_n(x1,y1,z1)*(gx(x2,y1,z1)-gx(x1,y1,z1))+(delt/delz)*phi_n(x1,y1,z1)*(gz(x1,y1,z2)-gz(x1,y1,z1)) + (delt/dely)*phi_n(x1,y1,z1)*(gy(x1,y2,z1)-gy(x1,y1,z1));
                                    psi(x1,y1,z1) = psi_n(x1,y1,z1) + ff * ((psi_n(x2,y1,z1)+psi_n(x3,y1,z1)+psi_n(x1,y2,z1)+psi_n(x1,y3,z1)+phi_n(x1,y1,z2)+phi_n(x1,y1,z3))-6*psi_n(x1,y1,z1))+(1/delx)*(eta(x2,y1,z1)-eta(x1,y1,z1)+eta(x1,y2,z1)-eta(x1,y1,z1)+eta(x1,y1,z2)-eta(x1,y1,z1)) + (delt/delx)*(psi_n(x2,y1,z1)-psi_n(x1,y1,z1))*fx(x1,y1,z1)+(delt/dely)*(psi_n(x1,y2,z1)-psi_n(x1,y1,z1))*fy(x1,y1,z1)+(delt/dely)*(psi_n(x1,y1,z2)-psi_n(x1,y1,z1))*fz(x1,y1,z1)+(delt/delx)*psi_n(x1,y1,z1)*(fx(x2,y1,z1)-fx(x1,y1,z1))+(delt/delz)*psi_n(x1,y1,z1)*(fz(x1,y1,z2)-fz(x1,y1,z1)) + (delt/dely)*psi_n(x1,y1,z1)*(fy(x1,y2,z1)-fy(x1,y1,z1));
                                    
                                    if(phi(x1,y1,z1)<0)
                                        phi(x1,y1,z1)= phi_n(x1,y1,z1);
                                        %psi(x1,y1,z1) = psi_n(x1,y1,z1);
                                    end
                                    if(psi(x1,y1,z1)<0)
                                        psi(x1,y1,z1)= psi_n(x1,y1,z1);
                                        %psi(x1,y1,z1) = psi_n(x1,y1,z1);
                                    end
                                    
                                    
                           end
                           
                           if(i==1 && k==nz && j>1 && j<ny) % 2,1,3
                                    x1=i;
                                    x2=i+1;
                                    x3=nx;
                                    
                                    y1=j;
                                    y2=j+1;
                                    y3=j-1;
                                    
                                    z1=k;
                                    z2=1;
                                    z3=k-1;
                                    
                                    phi(x1,y1,z1) = phi_n(x1,y1,z1) + ff * ((phi_n(x2,y1,z1)+phi_n(x3,y1,z1)+phi_n(x1,y2,z1)+phi_n(x1,y3,z1)+phi_n(x1,y1,z2)+phi_n(x1,y1,z3))-6*phi_n(x1,y1,z1))+delt*k_a*phi_n(x1,y1,z1)*(k_b/k_a-phi_n(x1,y1,z1))+b_noise(x1,y1,z1)+(1/(delx))*(eta_phi(x2,y1,z1)-eta_phi(x1,y1,z1)+eta_phi(x1,y2,z1)-eta_phi(x1,y1,z1)+eta_phi(x1,y1,z2)-eta_phi(x1,y1,z1)) + (delt/delx)*(phi_n(x2,y1,z1)-phi_n(x1,y1,z1))*gx(x1,y1,z1)+(delt/dely)*(phi_n(x1,y2,z1)-phi_n(x1,y1,z1))*gy(x1,y1,z1)+(delt/delz)*(phi_n(x1,y1,z2)-phi_n(x1,y1,z1))*gz(x1,y1,z1)+(delt/delx)*phi_n(x1,y1,z1)*(gx(x2,y1,z1)-gx(x1,y1,z1))+(delt/delz)*phi_n(x1,y1,z1)*(gz(x1,y1,z2)-gz(x1,y1,z1)) + (delt/dely)*phi_n(x1,y1,z1)*(gy(x1,y2,z1)-gy(x1,y1,z1));
                                    psi(x1,y1,z1) = psi_n(x1,y1,z1) + ff * ((psi_n(x2,y1,z1)+psi_n(x3,y1,z1)+psi_n(x1,y2,z1)+psi_n(x1,y3,z1)+phi_n(x1,y1,z2)+phi_n(x1,y1,z3))-6*psi_n(x1,y1,z1))+(1/delx)*(eta(x2,y1,z1)-eta(x1,y1,z1)+eta(x1,y2,z1)-eta(x1,y1,z1)+eta(x1,y1,z2)-eta(x1,y1,z1)) + (delt/delx)*(psi_n(x2,y1,z1)-psi_n(x1,y1,z1))*fx(x1,y1,z1)+(delt/dely)*(psi_n(x1,y2,z1)-psi_n(x1,y1,z1))*fy(x1,y1,z1)+(delt/dely)*(psi_n(x1,y1,z2)-psi_n(x1,y1,z1))*fz(x1,y1,z1)+(delt/delx)*psi_n(x1,y1,z1)*(fx(x2,y1,z1)-fx(x1,y1,z1))+(delt/delz)*psi_n(x1,y1,z1)*(fz(x1,y1,z2)-fz(x1,y1,z1)) + (delt/dely)*psi_n(x1,y1,z1)*(fy(x1,y2,z1)-fy(x1,y1,z1));
                                    
                                    
                                    if(phi(x1,y1,z1)<0)
                                        phi(x1,y1,z1)= phi_n(x1,y1,z1);
                                        %psi(x1,y1,z1) = psi_n(x1,y1,z1);
                                    end
                                    if(psi(x1,y1,z1)<0)
                                        psi(x1,y1,z1)= psi_n(x1,y1,z1);
                                        %psi(x1,y1,z1) = psi_n(x1,y1,z1);
                                    end
                                    
                           end
                           
                           if(i==nx && k==nz && j>1 && j<ny) %3,1,3
                                    x1=i;
                                    x2=1;
                                    x3=i-1;
                                    
                                    y1=j;
                                    y2=j+1;
                                    y3=j-1;
                                    
                                    z1=k;
                                    z2=1;
                                    z3=k-1;
                                    
                                    phi(x1,y1,z1) = phi_n(x1,y1,z1) + ff * ((phi_n(x2,y1,z1)+phi_n(x3,y1,z1)+phi_n(x1,y2,z1)+phi_n(x1,y3,z1)+phi_n(x1,y1,z2)+phi_n(x1,y1,z3))-6*phi_n(x1,y1,z1))+delt*k_a*phi_n(x1,y1,z1)*(k_b/k_a-phi_n(x1,y1,z1))+b_noise(x1,y1,z1)+(1/(delx))*(eta_phi(x2,y1,z1)-eta_phi(x1,y1,z1)+eta_phi(x1,y2,z1)-eta_phi(x1,y1,z1)+eta_phi(x1,y1,z2)-eta_phi(x1,y1,z1)) + (delt/delx)*(phi_n(x2,y1,z1)-phi_n(x1,y1,z1))*gx(x1,y1,z1)+(delt/dely)*(phi_n(x1,y2,z1)-phi_n(x1,y1,z1))*gy(x1,y1,z1)+(delt/delz)*(phi_n(x1,y1,z2)-phi_n(x1,y1,z1))*gz(x1,y1,z1)+(delt/delx)*phi_n(x1,y1,z1)*(gx(x2,y1,z1)-gx(x1,y1,z1))+(delt/delz)*phi_n(x1,y1,z1)*(gz(x1,y1,z2)-gz(x1,y1,z1)) + (delt/dely)*phi_n(x1,y1,z1)*(gy(x1,y2,z1)-gy(x1,y1,z1));
                                    psi(x1,y1,z1) = psi_n(x1,y1,z1) + ff * ((psi_n(x2,y1,z1)+psi_n(x3,y1,z1)+psi_n(x1,y2,z1)+psi_n(x1,y3,z1)+phi_n(x1,y1,z2)+phi_n(x1,y1,z3))-6*psi_n(x1,y1,z1))+(1/delx)*(eta(x2,y1,z1)-eta(x1,y1,z1)+eta(x1,y2,z1)-eta(x1,y1,z1)+eta(x1,y1,z2)-eta(x1,y1,z1)) + (delt/delx)*(psi_n(x2,y1,z1)-psi_n(x1,y1,z1))*fx(x1,y1,z1)+(delt/dely)*(psi_n(x1,y2,z1)-psi_n(x1,y1,z1))*fy(x1,y1,z1)+(delt/dely)*(psi_n(x1,y1,z2)-psi_n(x1,y1,z1))*fz(x1,y1,z1)+(delt/delx)*psi_n(x1,y1,z1)*(fx(x2,y1,z1)-fx(x1,y1,z1))+(delt/delz)*psi_n(x1,y1,z1)*(fz(x1,y1,z2)-fz(x1,y1,z1)) + (delt/dely)*psi_n(x1,y1,z1)*(fy(x1,y2,z1)-fy(x1,y1,z1));
                                    
                                    
                                    if(phi(x1,y1,z1)<0)
                                        phi(x1,y1,z1)= phi_n(x1,y1,z1);
                                        %psi(x1,y1,z1) = psi_n(x1,y1,z1);
                                    end
                                    if(psi(x1,y1,z1)<0)
                                        psi(x1,y1,z1)= psi_n(x1,y1,z1);
                                        %psi(x1,y1,z1) = psi_n(x1,y1,z1);
                                    end
                                    
                           end
                           
                           if(j==1 && k==1 && i>1 && i<nx) %1,2,2
                                    x1=i;
                                    x2=i+1;
                                    x3=i-1;
                                    
                                    y1=j;
                                    y2=j+1;
                                    y3=ny;
                                    
                                    z1=k;
                                    z2=k+1;
                                    z3=nz;
                                    
                                    phi(x1,y1,z1) = phi_n(x1,y1,z1) + ff * ((phi_n(x2,y1,z1)+phi_n(x3,y1,z1)+phi_n(x1,y2,z1)+phi_n(x1,y3,z1)+phi_n(x1,y1,z2)+phi_n(x1,y1,z3))-6*phi_n(x1,y1,z1))+delt*k_a*phi_n(x1,y1,z1)*(k_b/k_a-phi_n(x1,y1,z1))+b_noise(x1,y1,z1)+(1/(delx))*(eta_phi(x2,y1,z1)-eta_phi(x1,y1,z1)+eta_phi(x1,y2,z1)-eta_phi(x1,y1,z1)+eta_phi(x1,y1,z2)-eta_phi(x1,y1,z1)) + (delt/delx)*(phi_n(x2,y1,z1)-phi_n(x1,y1,z1))*gx(x1,y1,z1)+(delt/dely)*(phi_n(x1,y2,z1)-phi_n(x1,y1,z1))*gy(x1,y1,z1)+(delt/delz)*(phi_n(x1,y1,z2)-phi_n(x1,y1,z1))*gz(x1,y1,z1)+(delt/delx)*phi_n(x1,y1,z1)*(gx(x2,y1,z1)-gx(x1,y1,z1))+(delt/delz)*phi_n(x1,y1,z1)*(gz(x1,y1,z2)-gz(x1,y1,z1)) + (delt/dely)*phi_n(x1,y1,z1)*(gy(x1,y2,z1)-gy(x1,y1,z1));
                                    psi(x1,y1,z1) = psi_n(x1,y1,z1) + ff * ((psi_n(x2,y1,z1)+psi_n(x3,y1,z1)+psi_n(x1,y2,z1)+psi_n(x1,y3,z1)+phi_n(x1,y1,z2)+phi_n(x1,y1,z3))-6*psi_n(x1,y1,z1))+(1/delx)*(eta(x2,y1,z1)-eta(x1,y1,z1)+eta(x1,y2,z1)-eta(x1,y1,z1)+eta(x1,y1,z2)-eta(x1,y1,z1)) + (delt/delx)*(psi_n(x2,y1,z1)-psi_n(x1,y1,z1))*fx(x1,y1,z1)+(delt/dely)*(psi_n(x1,y2,z1)-psi_n(x1,y1,z1))*fy(x1,y1,z1)+(delt/dely)*(psi_n(x1,y1,z2)-psi_n(x1,y1,z1))*fz(x1,y1,z1)+(delt/delx)*psi_n(x1,y1,z1)*(fx(x2,y1,z1)-fx(x1,y1,z1))+(delt/delz)*psi_n(x1,y1,z1)*(fz(x1,y1,z2)-fz(x1,y1,z1)) + (delt/dely)*psi_n(x1,y1,z1)*(fy(x1,y2,z1)-fy(x1,y1,z1));
                                    
                                    
                                    if(phi(x1,y1,z1)<0)
                                        phi(x1,y1,z1)= phi_n(x1,y1,z1);
                                        %psi(x1,y1,z1) = psi_n(x1,y1,z1);
                                    end
                                    if(psi(x1,y1,z1)<0)
                                        psi(x1,y1,z1)= psi_n(x1,y1,z1);
                                        %psi(x1,y1,z1) = psi_n(x1,y1,z1);
                                    end
                                    
                           end
                           
                           if(j==ny && k==1 && i>1 && i<nx) %1,3,2
                                    x1=i;
                                    x2=i+1;
                                    x3=i-1;
                                    
                                    y1=j;
                                    y2=1;
                                    y3=j-1;
                                    
                                    z1=k;
                                    z2=k+1;
                                    z3=nz;
                                    
                                    phi(x1,y1,z1) = phi_n(x1,y1,z1) + ff * ((phi_n(x2,y1,z1)+phi_n(x3,y1,z1)+phi_n(x1,y2,z1)+phi_n(x1,y3,z1)+phi_n(x1,y1,z2)+phi_n(x1,y1,z3))-6*phi_n(x1,y1,z1))+delt*k_a*phi_n(x1,y1,z1)*(k_b/k_a-phi_n(x1,y1,z1))+b_noise(x1,y1,z1)+(1/(delx))*(eta_phi(x2,y1,z1)-eta_phi(x1,y1,z1)+eta_phi(x1,y2,z1)-eta_phi(x1,y1,z1)+eta_phi(x1,y1,z2)-eta_phi(x1,y1,z1)) + (delt/delx)*(phi_n(x2,y1,z1)-phi_n(x1,y1,z1))*gx(x1,y1,z1)+(delt/dely)*(phi_n(x1,y2,z1)-phi_n(x1,y1,z1))*gy(x1,y1,z1)+(delt/delz)*(phi_n(x1,y1,z2)-phi_n(x1,y1,z1))*gz(x1,y1,z1)+(delt/delx)*phi_n(x1,y1,z1)*(gx(x2,y1,z1)-gx(x1,y1,z1))+(delt/delz)*phi_n(x1,y1,z1)*(gz(x1,y1,z2)-gz(x1,y1,z1)) + (delt/dely)*phi_n(x1,y1,z1)*(gy(x1,y2,z1)-gy(x1,y1,z1));
                                    psi(x1,y1,z1) = psi_n(x1,y1,z1) + ff * ((psi_n(x2,y1,z1)+psi_n(x3,y1,z1)+psi_n(x1,y2,z1)+psi_n(x1,y3,z1)+phi_n(x1,y1,z2)+phi_n(x1,y1,z3))-6*psi_n(x1,y1,z1))+(1/delx)*(eta(x2,y1,z1)-eta(x1,y1,z1)+eta(x1,y2,z1)-eta(x1,y1,z1)+eta(x1,y1,z2)-eta(x1,y1,z1)) + (delt/delx)*(psi_n(x2,y1,z1)-psi_n(x1,y1,z1))*fx(x1,y1,z1)+(delt/dely)*(psi_n(x1,y2,z1)-psi_n(x1,y1,z1))*fy(x1,y1,z1)+(delt/dely)*(psi_n(x1,y1,z2)-psi_n(x1,y1,z1))*fz(x1,y1,z1)+(delt/delx)*psi_n(x1,y1,z1)*(fx(x2,y1,z1)-fx(x1,y1,z1))+(delt/delz)*psi_n(x1,y1,z1)*(fz(x1,y1,z2)-fz(x1,y1,z1)) + (delt/dely)*psi_n(x1,y1,z1)*(fy(x1,y2,z1)-fy(x1,y1,z1));
                                    
                                    
                                    if(phi(x1,y1,z1)<0)
                                        phi(x1,y1,z1)= phi_n(x1,y1,z1);
                                        %psi(x1,y1,z1) = psi_n(x1,y1,z1);
                                    end
                                    if(psi(x1,y1,z1)<0)
                                        psi(x1,y1,z1)= psi_n(x1,y1,z1);
                                        %psi(x1,y1,z1) = psi_n(x1,y1,z1);
                                    end
                           end
                           
                           if(j==1 && k==nz && i>1 && i<nx) %1,2,3
                                    x1=i;
                                    x2=i+1;
                                    x3=i-1;
                                    
                                    y1=j;
                                    y2=j+1;
                                    y3=ny;
                                    
                                    z1=k;
                                    z2=1;
                                    z3=k-1;
                                    
                                    phi(x1,y1,z1) = phi_n(x1,y1,z1) + ff * ((phi_n(x2,y1,z1)+phi_n(x3,y1,z1)+phi_n(x1,y2,z1)+phi_n(x1,y3,z1)+phi_n(x1,y1,z2)+phi_n(x1,y1,z3))-6*phi_n(x1,y1,z1))+delt*k_a*phi_n(x1,y1,z1)*(k_b/k_a-phi_n(x1,y1,z1))+b_noise(x1,y1,z1)+(1/(delx))*(eta_phi(x2,y1,z1)-eta_phi(x1,y1,z1)+eta_phi(x1,y2,z1)-eta_phi(x1,y1,z1)+eta_phi(x1,y1,z2)-eta_phi(x1,y1,z1)) + (delt/delx)*(phi_n(x2,y1,z1)-phi_n(x1,y1,z1))*gx(x1,y1,z1)+(delt/dely)*(phi_n(x1,y2,z1)-phi_n(x1,y1,z1))*gy(x1,y1,z1)+(delt/delz)*(phi_n(x1,y1,z2)-phi_n(x1,y1,z1))*gz(x1,y1,z1)+(delt/delx)*phi_n(x1,y1,z1)*(gx(x2,y1,z1)-gx(x1,y1,z1))+(delt/delz)*phi_n(x1,y1,z1)*(gz(x1,y1,z2)-gz(x1,y1,z1)) + (delt/dely)*phi_n(x1,y1,z1)*(gy(x1,y2,z1)-gy(x1,y1,z1));
                                    psi(x1,y1,z1) = psi_n(x1,y1,z1) + ff * ((psi_n(x2,y1,z1)+psi_n(x3,y1,z1)+psi_n(x1,y2,z1)+psi_n(x1,y3,z1)+phi_n(x1,y1,z2)+phi_n(x1,y1,z3))-6*psi_n(x1,y1,z1))+(1/delx)*(eta(x2,y1,z1)-eta(x1,y1,z1)+eta(x1,y2,z1)-eta(x1,y1,z1)+eta(x1,y1,z2)-eta(x1,y1,z1)) + (delt/delx)*(psi_n(x2,y1,z1)-psi_n(x1,y1,z1))*fx(x1,y1,z1)+(delt/dely)*(psi_n(x1,y2,z1)-psi_n(x1,y1,z1))*fy(x1,y1,z1)+(delt/dely)*(psi_n(x1,y1,z2)-psi_n(x1,y1,z1))*fz(x1,y1,z1)+(delt/delx)*psi_n(x1,y1,z1)*(fx(x2,y1,z1)-fx(x1,y1,z1))+(delt/delz)*psi_n(x1,y1,z1)*(fz(x1,y1,z2)-fz(x1,y1,z1)) + (delt/dely)*psi_n(x1,y1,z1)*(fy(x1,y2,z1)-fy(x1,y1,z1));
                                    
                                    
                                    if(phi(x1,y1,z1)<0)
                                        phi(x1,y1,z1)= phi_n(x1,y1,z1);
                                        %psi(x1,y1,z1) = psi_n(x1,y1,z1);
                                    end
                                    if(psi(x1,y1,z1)<0)
                                        psi(x1,y1,z1)= psi_n(x1,y1,z1);
                                        %psi(x1,y1,z1) = psi_n(x1,y1,z1);
                                    end
                                    
                           end
                           
                           if(j==ny && k==nz && i>1 && i<nx) %1,3,3
                                    x1=i;
                                    x2=i+1;
                                    x3=i-1;
                                    
                                    y1=j;
                                    y2=1;
                                    y3=j-1;
                                    
                                    z1=k;
                                    z2=1;
                                    z3=k-1;
                                    
                                    phi(x1,y1,z1) = phi_n(x1,y1,z1) + ff * ((phi_n(x2,y1,z1)+phi_n(x3,y1,z1)+phi_n(x1,y2,z1)+phi_n(x1,y3,z1)+phi_n(x1,y1,z2)+phi_n(x1,y1,z3))-6*phi_n(x1,y1,z1))+delt*k_a*phi_n(x1,y1,z1)*(k_b/k_a-phi_n(x1,y1,z1))+b_noise(x1,y1,z1)+(1/(delx))*(eta_phi(x2,y1,z1)-eta_phi(x1,y1,z1)+eta_phi(x1,y2,z1)-eta_phi(x1,y1,z1)+eta_phi(x1,y1,z2)-eta_phi(x1,y1,z1))+ (delt/delx)*(phi_n(x2,y1,z1)-phi_n(x1,y1,z1))*gx(x1,y1,z1)+(delt/dely)*(phi_n(x1,y2,z1)-phi_n(x1,y1,z1))*gy(x1,y1,z1)+(delt/delz)*(phi_n(x1,y1,z2)-phi_n(x1,y1,z1))*gz(x1,y1,z1)+(delt/delx)*phi_n(x1,y1,z1)*(gx(x2,y1,z1)-gx(x1,y1,z1))+(delt/delz)*phi_n(x1,y1,z1)*(gz(x1,y1,z2)-gz(x1,y1,z1)) + (delt/dely)*phi_n(x1,y1,z1)*(gy(x1,y2,z1)-gy(x1,y1,z1));
                                    psi(x1,y1,z1) = psi_n(x1,y1,z1) + ff * ((psi_n(x2,y1,z1)+psi_n(x3,y1,z1)+psi_n(x1,y2,z1)+psi_n(x1,y3,z1)+phi_n(x1,y1,z2)+phi_n(x1,y1,z3))-6*psi_n(x1,y1,z1))+(1/delx)*(eta(x2,y1,z1)-eta(x1,y1,z1)+eta(x1,y2,z1)-eta(x1,y1,z1)+eta(x1,y1,z2)-eta(x1,y1,z1)) + (delt/delx)*(psi_n(x2,y1,z1)-psi_n(x1,y1,z1))*fx(x1,y1,z1)+(delt/dely)*(psi_n(x1,y2,z1)-psi_n(x1,y1,z1))*fy(x1,y1,z1)+(delt/dely)*(psi_n(x1,y1,z2)-psi_n(x1,y1,z1))*fz(x1,y1,z1)+(delt/delx)*psi_n(x1,y1,z1)*(fx(x2,y1,z1)-fx(x1,y1,z1))+(delt/delz)*psi_n(x1,y1,z1)*(fz(x1,y1,z2)-fz(x1,y1,z1)) + (delt/dely)*psi_n(x1,y1,z1)*(fy(x1,y2,z1)-fy(x1,y1,z1));
                                    
                                    
                                    if(phi(x1,y1,z1)<0)
                                        phi(x1,y1,z1)= phi_n(x1,y1,z1);
                                        %psi(x1,y1,z1) = psi_n(x1,y1,z1);
                                    end
                                    if(psi(x1,y1,z1)<0)
                                        psi(x1,y1,z1)= psi_n(x1,y1,z1);
                                        %psi(x1,y1,z1) = psi_n(x1,y1,z1);
                                    end
                                    
                                    
                           end
                           
                           if(j==1 && k==1 && i==1)  % 2,2,2
                                    x1=i;
                                    x2=i+1;
                                    x3=nx;
                                    
                                    y1=j;
                                    y2=j+1;
                                    y3=ny;
                                    
                                    z1=k;
                                    z2=k+1;
                                    z3=nz;
                                    
                                    phi(x1,y1,z1) = phi_n(x1,y1,z1) + ff * ((phi_n(x2,y1,z1)+phi_n(x3,y1,z1)+phi_n(x1,y2,z1)+phi_n(x1,y3,z1)+phi_n(x1,y1,z2)+phi_n(x1,y1,z3))-6*phi_n(x1,y1,z1))+delt*k_a*phi_n(x1,y1,z1)*(k_b/k_a-phi_n(x1,y1,z1))+b_noise(x1,y1,z1)+(1/(delx))*(eta_phi(x2,y1,z1)-eta_phi(x1,y1,z1)+eta_phi(x1,y2,z1)-eta_phi(x1,y1,z1)+eta_phi(x1,y1,z2)-eta_phi(x1,y1,z1)) + (delt/delx)*(phi_n(x2,y1,z1)-phi_n(x1,y1,z1))*gx(x1,y1,z1)+(delt/dely)*(phi_n(x1,y2,z1)-phi_n(x1,y1,z1))*gy(x1,y1,z1)+(delt/delz)*(phi_n(x1,y1,z2)-phi_n(x1,y1,z1))*gz(x1,y1,z1)+(delt/delx)*phi_n(x1,y1,z1)*(gx(x2,y1,z1)-gx(x1,y1,z1))+(delt/delz)*phi_n(x1,y1,z1)*(gz(x1,y1,z2)-gz(x1,y1,z1)) + (delt/dely)*phi_n(x1,y1,z1)*(gy(x1,y2,z1)-gy(x1,y1,z1));
                                    psi(x1,y1,z1) = psi_n(x1,y1,z1) + ff * ((psi_n(x2,y1,z1)+psi_n(x3,y1,z1)+psi_n(x1,y2,z1)+psi_n(x1,y3,z1)+phi_n(x1,y1,z2)+phi_n(x1,y1,z3))-6*psi_n(x1,y1,z1))+(1/delx)*(eta(x2,y1,z1)-eta(x1,y1,z1)+eta(x1,y2,z1)-eta(x1,y1,z1)+eta(x1,y1,z2)-eta(x1,y1,z1)) + (delt/delx)*(psi_n(x2,y1,z1)-psi_n(x1,y1,z1))*fx(x1,y1,z1)+(delt/dely)*(psi_n(x1,y2,z1)-psi_n(x1,y1,z1))*fy(x1,y1,z1)+(delt/dely)*(psi_n(x1,y1,z2)-psi_n(x1,y1,z1))*fz(x1,y1,z1)+(delt/delx)*psi_n(x1,y1,z1)*(fx(x2,y1,z1)-fx(x1,y1,z1))+(delt/delz)*psi_n(x1,y1,z1)*(fz(x1,y1,z2)-fz(x1,y1,z1)) + (delt/dely)*psi_n(x1,y1,z1)*(fy(x1,y2,z1)-fy(x1,y1,z1));
                                    
                                    
                                    if(phi(x1,y1,z1)<0)
                                        phi(x1,y1,z1)= phi_n(x1,y1,z1);
                                        %psi(x1,y1,z1) = psi_n(x1,y1,z1);
                                    end
                                    if(psi(x1,y1,z1)<0)
                                        psi(x1,y1,z1)= psi_n(x1,y1,z1);
                                        %psi(x1,y1,z1) = psi_n(x1,y1,z1);
                                    end
                                    
                           end
                           
                           if(j==ny && k==1 && i==1) %2,3,2
                                    x1=i;
                                    x2=i+1;
                                    x3=nx;
                                    
                                    y1=j;
                                    y2=1;
                                    y3=j-1;
                                    
                                    z1=k;
                                    z2=k+1;
                                    z3=nz;
                                    
                                    phi(x1,y1,z1) = phi_n(x1,y1,z1) + ff * ((phi_n(x2,y1,z1)+phi_n(x3,y1,z1)+phi_n(x1,y2,z1)+phi_n(x1,y3,z1)+phi_n(x1,y1,z2)+phi_n(x1,y1,z3))-6*phi_n(x1,y1,z1))+delt*k_a*phi_n(x1,y1,z1)*(k_b/k_a-phi_n(x1,y1,z1))+b_noise(x1,y1,z1)+(1/(delx))*(eta_phi(x2,y1,z1)-eta_phi(x1,y1,z1)+eta_phi(x1,y2,z1)-eta_phi(x1,y1,z1)+eta_phi(x1,y1,z2)-eta_phi(x1,y1,z1)) + (delt/delx)*(phi_n(x2,y1,z1)-phi_n(x1,y1,z1))*gx(x1,y1,z1)+(delt/dely)*(phi_n(x1,y2,z1)-phi_n(x1,y1,z1))*gy(x1,y1,z1)+(delt/delz)*(phi_n(x1,y1,z2)-phi_n(x1,y1,z1))*gz(x1,y1,z1)+(delt/delx)*phi_n(x1,y1,z1)*(gx(x2,y1,z1)-gx(x1,y1,z1))+(delt/delz)*phi_n(x1,y1,z1)*(gz(x1,y1,z2)-gz(x1,y1,z1)) + (delt/dely)*phi_n(x1,y1,z1)*(gy(x1,y2,z1)-gy(x1,y1,z1));
                                    psi(x1,y1,z1) = psi_n(x1,y1,z1) + ff * ((psi_n(x2,y1,z1)+psi_n(x3,y1,z1)+psi_n(x1,y2,z1)+psi_n(x1,y3,z1)+phi_n(x1,y1,z2)+phi_n(x1,y1,z3))-6*psi_n(x1,y1,z1))+(1/delx)*(eta(x2,y1,z1)-eta(x1,y1,z1)+eta(x1,y2,z1)-eta(x1,y1,z1)+eta(x1,y1,z2)-eta(x1,y1,z1)) + (delt/delx)*(psi_n(x2,y1,z1)-psi_n(x1,y1,z1))*fx(x1,y1,z1)+(delt/dely)*(psi_n(x1,y2,z1)-psi_n(x1,y1,z1))*fy(x1,y1,z1)+(delt/dely)*(psi_n(x1,y1,z2)-psi_n(x1,y1,z1))*fz(x1,y1,z1)+(delt/delx)*psi_n(x1,y1,z1)*(fx(x2,y1,z1)-fx(x1,y1,z1))+(delt/delz)*psi_n(x1,y1,z1)*(fz(x1,y1,z2)-fz(x1,y1,z1)) + (delt/dely)*psi_n(x1,y1,z1)*(fy(x1,y2,z1)-fy(x1,y1,z1));
                                    
                                    
                                    if(phi(x1,y1,z1)<0)
                                        phi(x1,y1,z1)= phi_n(x1,y1,z1);
                                        %psi(x1,y1,z1) = psi_n(x1,y1,z1);
                                    end
                                    if(psi(x1,y1,z1)<0)
                                        psi(x1,y1,z1)= psi_n(x1,y1,z1);
                                        %psi(x1,y1,z1) = psi_n(x1,y1,z1);
                                    end
                                    
                                    
                           end
                           
                           if(j==1 && k==nz && i==1) %2,2,3
                                    x1=i;
                                    x2=i+1;
                                    x3=nx;
                                    
                                    y1=j;
                                    y2=j+1;
                                    y3=ny;
                                    
                                    z1=k;
                                    z2=1;
                                    z3=k-1;
                                    
                                    phi(x1,y1,z1) = phi_n(x1,y1,z1) + ff * ((phi_n(x2,y1,z1)+phi_n(x3,y1,z1)+phi_n(x1,y2,z1)+phi_n(x1,y3,z1)+phi_n(x1,y1,z2)+phi_n(x1,y1,z3))-6*phi_n(x1,y1,z1))+delt*k_a*phi_n(x1,y1,z1)*(k_b/k_a-phi_n(x1,y1,z1))+b_noise(x1,y1,z1)+(1/(delx))*(eta_phi(x2,y1,z1)-eta_phi(x1,y1,z1)+eta_phi(x1,y2,z1)-eta_phi(x1,y1,z1)+eta_phi(x1,y1,z2)-eta_phi(x1,y1,z1)) + (delt/delx)*(phi_n(x2,y1,z1)-phi_n(x1,y1,z1))*gx(x1,y1,z1)+(delt/dely)*(phi_n(x1,y2,z1)-phi_n(x1,y1,z1))*gy(x1,y1,z1)+(delt/delz)*(phi_n(x1,y1,z2)-phi_n(x1,y1,z1))*gz(x1,y1,z1)+(delt/delx)*phi_n(x1,y1,z1)*(gx(x2,y1,z1)-gx(x1,y1,z1))+(delt/delz)*phi_n(x1,y1,z1)*(gz(x1,y1,z2)-gz(x1,y1,z1)) + (delt/dely)*phi_n(x1,y1,z1)*(gy(x1,y2,z1)-gy(x1,y1,z1));
                                    psi(x1,y1,z1) = psi_n(x1,y1,z1) + ff * ((psi_n(x2,y1,z1)+psi_n(x3,y1,z1)+psi_n(x1,y2,z1)+psi_n(x1,y3,z1)+phi_n(x1,y1,z2)+phi_n(x1,y1,z3))-6*psi_n(x1,y1,z1))+(1/delx)*(eta(x2,y1,z1)-eta(x1,y1,z1)+eta(x1,y2,z1)-eta(x1,y1,z1)+eta(x1,y1,z2)-eta(x1,y1,z1)) + (delt/delx)*(psi_n(x2,y1,z1)-psi_n(x1,y1,z1))*fx(x1,y1,z1)+(delt/dely)*(psi_n(x1,y2,z1)-psi_n(x1,y1,z1))*fy(x1,y1,z1)+(delt/dely)*(psi_n(x1,y1,z2)-psi_n(x1,y1,z1))*fz(x1,y1,z1)+(delt/delx)*psi_n(x1,y1,z1)*(fx(x2,y1,z1)-fx(x1,y1,z1))+(delt/delz)*psi_n(x1,y1,z1)*(fz(x1,y1,z2)-fz(x1,y1,z1)) + (delt/dely)*psi_n(x1,y1,z1)*(fy(x1,y2,z1)-fy(x1,y1,z1));
                                    
                                    
                                    if(phi(x1,y1,z1)<0)
                                        phi(x1,y1,z1)= phi_n(x1,y1,z1);
                                        %psi(x1,y1,z1) = psi_n(x1,y1,z1);
                                    end
                                    if(psi(x1,y1,z1)<0)
                                        psi(x1,y1,z1)= psi_n(x1,y1,z1);
                                        %psi(x1,y1,z1) = psi_n(x1,y1,z1);
                                    end
                                    
                           end
                           
                           if(j==ny && k==nz && i==1)  % 2,3,3
                                    x1=i;
                                    x2=i+1;
                                    x3=nx;
                                    
                                    y1=j;
                                    y2=1;
                                    y3=j-1;
                                    
                                    z1=k;
                                    z2=1;
                                    z3=k-1;
                                    
                                    phi(x1,y1,z1) = phi_n(x1,y1,z1) + ff * ((phi_n(x2,y1,z1)+phi_n(x3,y1,z1)+phi_n(x1,y2,z1)+phi_n(x1,y3,z1)+phi_n(x1,y1,z2)+phi_n(x1,y1,z3))-6*phi_n(x1,y1,z1))+delt*k_a*phi_n(x1,y1,z1)*(k_b/k_a-phi_n(x1,y1,z1))+b_noise(x1,y1,z1)+(1/(delx))*(eta_phi(x2,y1,z1)-eta_phi(x1,y1,z1)+eta_phi(x1,y2,z1)-eta_phi(x1,y1,z1)+eta_phi(x1,y1,z2)-eta_phi(x1,y1,z1)) + (delt/delx)*(phi_n(x2,y1,z1)-phi_n(x1,y1,z1))*gx(x1,y1,z1)+(delt/dely)*(phi_n(x1,y2,z1)-phi_n(x1,y1,z1))*gy(x1,y1,z1)+(delt/delz)*(phi_n(x1,y1,z2)-phi_n(x1,y1,z1))*gz(x1,y1,z1)+(delt/delx)*phi_n(x1,y1,z1)*(gx(x2,y1,z1)-gx(x1,y1,z1))+(delt/delz)*phi_n(x1,y1,z1)*(gz(x1,y1,z2)-gz(x1,y1,z1)) + (delt/dely)*phi_n(x1,y1,z1)*(gy(x1,y2,z1)-gy(x1,y1,z1));
                                    psi(x1,y1,z1) = psi_n(x1,y1,z1) + ff * ((psi_n(x2,y1,z1)+psi_n(x3,y1,z1)+psi_n(x1,y2,z1)+psi_n(x1,y3,z1)+phi_n(x1,y1,z2)+phi_n(x1,y1,z3))-6*psi_n(x1,y1,z1))+(1/delx)*(eta(x2,y1,z1)-eta(x1,y1,z1)+eta(x1,y2,z1)-eta(x1,y1,z1)+eta(x1,y1,z2)-eta(x1,y1,z1)) + (delt/delx)*(psi_n(x2,y1,z1)-psi_n(x1,y1,z1))*fx(x1,y1,z1)+(delt/dely)*(psi_n(x1,y2,z1)-psi_n(x1,y1,z1))*fy(x1,y1,z1)+(delt/dely)*(psi_n(x1,y1,z2)-psi_n(x1,y1,z1))*fz(x1,y1,z1)+(delt/delx)*psi_n(x1,y1,z1)*(fx(x2,y1,z1)-fx(x1,y1,z1))+(delt/delz)*psi_n(x1,y1,z1)*(fz(x1,y1,z2)-fz(x1,y1,z1)) + (delt/dely)*psi_n(x1,y1,z1)*(fy(x1,y2,z1)-fy(x1,y1,z1));
                                    
                                    
                                    if(phi(x1,y1,z1)<0)
                                        phi(x1,y1,z1)= phi_n(x1,y1,z1);
                                        %psi(x1,y1,z1) = psi_n(x1,y1,z1);
                                    end
                                    if(psi(x1,y1,z1)<0)
                                        psi(x1,y1,z1)= psi_n(x1,y1,z1);
                                        %psi(x1,y1,z1) = psi_n(x1,y1,z1);
                                    end
                                    
                           end
                           
                           if(i==nx && k==1 && j==1) %3,2,2
                                    x1=i;
                                    x2=1;
                                    x3=i-1;
                                    
                                    y1=j;
                                    y2=j+1;
                                    y3=ny;
                                    
                                    z1=k;
                                    z2=k+1;
                                    z3=nz;
                                    
                                    phi(x1,y1,z1) = phi_n(x1,y1,z1) + ff * ((phi_n(x2,y1,z1)+phi_n(x3,y1,z1)+phi_n(x1,y2,z1)+phi_n(x1,y3,z1)+phi_n(x1,y1,z2)+phi_n(x1,y1,z3))-6*phi_n(x1,y1,z1))+delt*k_a*phi_n(x1,y1,z1)*(k_b/k_a-phi_n(x1,y1,z1))+b_noise(x1,y1,z1)+(1/(delx))*(eta_phi(x2,y1,z1)-eta_phi(x1,y1,z1)+eta_phi(x1,y2,z1)-eta_phi(x1,y1,z1)+eta_phi(x1,y1,z2)-eta_phi(x1,y1,z1)) + (delt/delx)*(phi_n(x2,y1,z1)-phi_n(x1,y1,z1))*gx(x1,y1,z1)+(delt/dely)*(phi_n(x1,y2,z1)-phi_n(x1,y1,z1))*gy(x1,y1,z1)+(delt/delz)*(phi_n(x1,y1,z2)-phi_n(x1,y1,z1))*gz(x1,y1,z1)+(delt/delx)*phi_n(x1,y1,z1)*(gx(x2,y1,z1)-gx(x1,y1,z1))+(delt/delz)*phi_n(x1,y1,z1)*(gz(x1,y1,z2)-gz(x1,y1,z1)) + (delt/dely)*phi_n(x1,y1,z1)*(gy(x1,y2,z1)-gy(x1,y1,z1));
                                    psi(x1,y1,z1) = psi_n(x1,y1,z1) + ff * ((psi_n(x2,y1,z1)+psi_n(x3,y1,z1)+psi_n(x1,y2,z1)+psi_n(x1,y3,z1)+phi_n(x1,y1,z2)+phi_n(x1,y1,z3))-6*psi_n(x1,y1,z1))+(1/delx)*(eta(x2,y1,z1)-eta(x1,y1,z1)+eta(x1,y2,z1)-eta(x1,y1,z1)+eta(x1,y1,z2)-eta(x1,y1,z1)) + (delt/delx)*(psi_n(x2,y1,z1)-psi_n(x1,y1,z1))*fx(x1,y1,z1)+(delt/dely)*(psi_n(x1,y2,z1)-psi_n(x1,y1,z1))*fy(x1,y1,z1)+(delt/dely)*(psi_n(x1,y1,z2)-psi_n(x1,y1,z1))*fz(x1,y1,z1)+(delt/delx)*psi_n(x1,y1,z1)*(fx(x2,y1,z1)-fx(x1,y1,z1))+(delt/delz)*psi_n(x1,y1,z1)*(fz(x1,y1,z2)-fz(x1,y1,z1)) + (delt/dely)*psi_n(x1,y1,z1)*(fy(x1,y2,z1)-fy(x1,y1,z1));
                                    
                                    
                                    if(phi(x1,y1,z1)<0)
                                        phi(x1,y1,z1)= phi_n(x1,y1,z1);
                                        %psi(x1,y1,z1) = psi_n(x1,y1,z1);
                                    end
                                    if(psi(x1,y1,z1)<0)
                                        psi(x1,y1,z1)= psi_n(x1,y1,z1);
                                        %psi(x1,y1,z1) = psi_n(x1,y1,z1);
                                    end
                                    
                           end
                           
                           if(i==nx && k==nz && j==1) %3,2,3
                                    x1=i;
                                    x2=1;
                                    x3=i-1;
                                    
                                    y1=j;
                                    y2=j+1;
                                    y3=ny;
                                    
                                    z1=k;
                                    z2=1;
                                    z3=k-1;
                                    
                                    phi(x1,y1,z1) = phi_n(x1,y1,z1) + ff * ((phi_n(x2,y1,z1)+phi_n(x3,y1,z1)+phi_n(x1,y2,z1)+phi_n(x1,y3,z1)+phi_n(x1,y1,z2)+phi_n(x1,y1,z3))-6*phi_n(x1,y1,z1))+delt*k_a*phi_n(x1,y1,z1)*(k_b/k_a-phi_n(x1,y1,z1))+b_noise(x1,y1,z1)+(1/(delx))*(eta_phi(x2,y1,z1)-eta_phi(x1,y1,z1)+eta_phi(x1,y2,z1)-eta_phi(x1,y1,z1)+eta_phi(x1,y1,z2)-eta_phi(x1,y1,z1)) + (delt/delx)*(phi_n(x2,y1,z1)-phi_n(x1,y1,z1))*gx(x1,y1,z1)+(delt/dely)*(phi_n(x1,y2,z1)-phi_n(x1,y1,z1))*gy(x1,y1,z1)+(delt/delz)*(phi_n(x1,y1,z2)-phi_n(x1,y1,z1))*gz(x1,y1,z1)+(delt/delx)*phi_n(x1,y1,z1)*(gx(x2,y1,z1)-gx(x1,y1,z1))+(delt/delz)*phi_n(x1,y1,z1)*(gz(x1,y1,z2)-gz(x1,y1,z1)) + (delt/dely)*phi_n(x1,y1,z1)*(gy(x1,y2,z1)-gy(x1,y1,z1));
                                    psi(x1,y1,z1) = psi_n(x1,y1,z1) + ff * ((psi_n(x2,y1,z1)+psi_n(x3,y1,z1)+psi_n(x1,y2,z1)+psi_n(x1,y3,z1)+phi_n(x1,y1,z2)+phi_n(x1,y1,z3))-6*psi_n(x1,y1,z1))+(1/delx)*(eta(x2,y1,z1)-eta(x1,y1,z1)+eta(x1,y2,z1)-eta(x1,y1,z1)+eta(x1,y1,z2)-eta(x1,y1,z1)) + (delt/delx)*(psi_n(x2,y1,z1)-psi_n(x1,y1,z1))*fx(x1,y1,z1)+(delt/dely)*(psi_n(x1,y2,z1)-psi_n(x1,y1,z1))*fy(x1,y1,z1)+(delt/dely)*(psi_n(x1,y1,z2)-psi_n(x1,y1,z1))*fz(x1,y1,z1)+(delt/delx)*psi_n(x1,y1,z1)*(fx(x2,y1,z1)-fx(x1,y1,z1))+(delt/delz)*psi_n(x1,y1,z1)*(fz(x1,y1,z2)-fz(x1,y1,z1)) + (delt/dely)*psi_n(x1,y1,z1)*(fy(x1,y2,z1)-fy(x1,y1,z1));
                                    
                                    if(phi(x1,y1,z1)<0)
                                        phi(x1,y1,z1)= phi_n(x1,y1,z1);
                                        %psi(x1,y1,z1) = psi_n(x1,y1,z1);
                                    end
                                    if(psi(x1,y1,z1)<0)
                                        psi(x1,y1,z1)= psi_n(x1,y1,z1);
                                        %psi(x1,y1,z1) = psi_n(x1,y1,z1);
                                    end
                                    
                                    
                           end
                           
                           if(j==ny && i==nx && k==1) %3,3,2
                                     x1=i;
                                    x2=1;
                                    x3=i-1;
                                    
                                    y1=j;
                                    y2=1;
                                    y3=j-1;
                                    
                                    z1=k;
                                    z2=k+1;
                                    z3=nz;
                                    
                                    phi(x1,y1,z1) = phi_n(x1,y1,z1) + ff * ((phi_n(x2,y1,z1)+phi_n(x3,y1,z1)+phi_n(x1,y2,z1)+phi_n(x1,y3,z1)+phi_n(x1,y1,z2)+phi_n(x1,y1,z3))-6*phi_n(x1,y1,z1))+delt*k_a*phi_n(x1,y1,z1)*(k_b/k_a-phi_n(x1,y1,z1))+b_noise(x1,y1,z1)+(1/(delx))*(eta_phi(x2,y1,z1)-eta_phi(x1,y1,z1)+eta_phi(x1,y2,z1)-eta_phi(x1,y1,z1)+eta_phi(x1,y1,z2)-eta_phi(x1,y1,z1)) + (delt/delx)*(phi_n(x2,y1,z1)-phi_n(x1,y1,z1))*gx(x1,y1,z1)+(delt/dely)*(phi_n(x1,y2,z1)-phi_n(x1,y1,z1))*gy(x1,y1,z1)+(delt/delz)*(phi_n(x1,y1,z2)-phi_n(x1,y1,z1))*gz(x1,y1,z1)+(delt/delx)*phi_n(x1,y1,z1)*(gx(x2,y1,z1)-gx(x1,y1,z1))+(delt/delz)*phi_n(x1,y1,z1)*(gz(x1,y1,z2)-gz(x1,y1,z1)) + (delt/dely)*phi_n(x1,y1,z1)*(gy(x1,y2,z1)-gy(x1,y1,z1));
                                    psi(x1,y1,z1) = psi_n(x1,y1,z1) + ff * ((psi_n(x2,y1,z1)+psi_n(x3,y1,z1)+psi_n(x1,y2,z1)+psi_n(x1,y3,z1)+phi_n(x1,y1,z2)+phi_n(x1,y1,z3))-6*psi_n(x1,y1,z1))+(1/delx)*(eta(x2,y1,z1)-eta(x1,y1,z1)+eta(x1,y2,z1)-eta(x1,y1,z1)+eta(x1,y1,z2)-eta(x1,y1,z1)) + (delt/delx)*(psi_n(x2,y1,z1)-psi_n(x1,y1,z1))*fx(x1,y1,z1)+(delt/dely)*(psi_n(x1,y2,z1)-psi_n(x1,y1,z1))*fy(x1,y1,z1)+(delt/dely)*(psi_n(x1,y1,z2)-psi_n(x1,y1,z1))*fz(x1,y1,z1)+(delt/delx)*psi_n(x1,y1,z1)*(fx(x2,y1,z1)-fx(x1,y1,z1))+(delt/delz)*psi_n(x1,y1,z1)*(fz(x1,y1,z2)-fz(x1,y1,z1)) + (delt/dely)*psi_n(x1,y1,z1)*(fy(x1,y2,z1)-fy(x1,y1,z1));
                                    
                                    
                                    if(phi(x1,y1,z1)<0)
                                        phi(x1,y1,z1)= phi_n(x1,y1,z1);
                                        %psi(x1,y1,z1) = psi_n(x1,y1,z1);
                                    end
                                    if(psi(x1,y1,z1)<0)
                                        psi(x1,y1,z1)= psi_n(x1,y1,z1);
                                        %psi(x1,y1,z1) = psi_n(x1,y1,z1);
                                    end
                               
                           end
                           
                            if(j==ny && i==nx && k==nz) %3,3,3
                                
                                   x1=i;
                                    x2=1;
                                    x3=i-1;
                                    
                                    y1=j;
                                    y2=1;
                                    y3=j-1;
                                    
                                    z1=k;
                                    z2=1;
                                    z3=k-1;
                                    
                                    phi(x1,y1,z1) = phi_n(x1,y1,z1) + ff * ((phi_n(x2,y1,z1)+phi_n(x3,y1,z1)+phi_n(x1,y2,z1)+phi_n(x1,y3,z1)+phi_n(x1,y1,z2)+phi_n(x1,y1,z3))-6*phi_n(x1,y1,z1))+delt*k_a*phi_n(x1,y1,z1)*(k_b/k_a-phi_n(x1,y1,z1))+b_noise(x1,y1,z1)+(1/(delx))*(eta_phi(x2,y1,z1)-eta_phi(x1,y1,z1)+eta_phi(x1,y2,z1)-eta_phi(x1,y1,z1)+eta_phi(x1,y1,z2)-eta_phi(x1,y1,z1)) + (delt/delx)*(phi_n(x2,y1,z1)-phi_n(x1,y1,z1))*gx(x1,y1,z1)+(delt/dely)*(phi_n(x1,y2,z1)-phi_n(x1,y1,z1))*gy(x1,y1,z1)+(delt/delz)*(phi_n(x1,y1,z2)-phi_n(x1,y1,z1))*gz(x1,y1,z1)+(delt/delx)*phi_n(x1,y1,z1)*(gx(x2,y1,z1)-gx(x1,y1,z1))+(delt/delz)*phi_n(x1,y1,z1)*(gz(x1,y1,z2)-gz(x1,y1,z1)) + (delt/dely)*phi_n(x1,y1,z1)*(gy(x1,y2,z1)-gy(x1,y1,z1));
                                    psi(x1,y1,z1) = psi_n(x1,y1,z1) + ff * ((psi_n(x2,y1,z1)+psi_n(x3,y1,z1)+psi_n(x1,y2,z1)+psi_n(x1,y3,z1)+phi_n(x1,y1,z2)+phi_n(x1,y1,z3))-6*psi_n(x1,y1,z1))+(1/delx)*(eta(x2,y1,z1)-eta(x1,y1,z1)+eta(x1,y2,z1)-eta(x1,y1,z1)+eta(x1,y1,z2)-eta(x1,y1,z1)) + (delt/delx)*(psi_n(x2,y1,z1)-psi_n(x1,y1,z1))*fx(x1,y1,z1)+(delt/dely)*(psi_n(x1,y2,z1)-psi_n(x1,y1,z1))*fy(x1,y1,z1)+(delt/dely)*(psi_n(x1,y1,z2)-psi_n(x1,y1,z1))*fz(x1,y1,z1)+(delt/delx)*psi_n(x1,y1,z1)*(fx(x2,y1,z1)-fx(x1,y1,z1))+(delt/delz)*psi_n(x1,y1,z1)*(fz(x1,y1,z2)-fz(x1,y1,z1)) + (delt/dely)*psi_n(x1,y1,z1)*(fy(x1,y2,z1)-fy(x1,y1,z1));
                                    
                                    
                                    if(phi(x1,y1,z1)<0)
                                        phi(x1,y1,z1)= phi_n(x1,y1,z1);
                                        %psi(x1,y1,z1) = psi_n(x1,y1,z1);
                                    end
                                    if(psi(x1,y1,z1)<0)
                                        psi(x1,y1,z1)= psi_n(x1,y1,z1);
                                        %psi(x1,y1,z1) = psi_n(x1,y1,z1);
                                    end
                                    
                            end
                            
                            
                
            
           
            end
            
        end
        
     end
     
 norm_t_phi = sum(sum(sum(phi)));
 phi = phi/norm_t_phi;
 
 norm_t_psi = sum(sum(sum(psi)));
 psi = psi/norm_t_psi;
 
     
    total_phi(:,:,:,k1)=phi;
    total_psi(:,:,:,k1)=psi;
  
  
   phi_n=phi;
   psi_n=psi;
  
end
save('total_phi3d_tpcc20_4000_nt_test_cc_300_e.mat')
save('total_psi3d_tpcc20_4000_nt_test_tp_300_e.mat')
