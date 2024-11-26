
clc
clear
%% Array factor caculation  
 M=30;
 G=6;
 N=M;
 L=16;
j=1;  
c=3e8;
f=10e9;
f0=0.5e6;
landa=c./f;
dx=landa/3;
k=2*pi*f/c;
dy=dx;
A=M*dx*N*dy;
eta=120*pi;
r=200;
phi=linspace(0,2*pi,r);
theta=linspace(0,pi/2,r);
dtheta=theta(2)-theta(1);
dphi=phi(2)-phi(1);
% v1=exp(1i*pi);
% v2=exp(1i*0);
% s1=[v1 v2 v2 v2 v2 v2 v2 v2];
% s2=[v2 v1 v2 v2 v2 v2 v2 v2];
% s3=[v2 v2 v1 v2 v2 v2 v2 v2];
% s4=[v2 v2 v2 v1 v2 v2 v2 v2];
% s5=[v2 v2 v2 v2 v1 v2 v2 v2];
% s6=[v2 v2 v2 v2 v2 v1 v2 v2];
% s7=[v2 v2 v2 v2 v2 v2 v1 v2];
% s8=[v2 v2 v2 v2 v2 v2 v2 v1];
% phase1=v2*ones(1,N*L);
% for i=1:N
%     phase1(i+(i-1)*N)=v1;
% end
% % phase1=[s1 s2 s3 s4 s5 s6 s7 s8];
% %% cat
% Ga=phase1;
% for h=1:N-1
%     Ga=cat(1,Ga,phase1);
% end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% generation of phase/amplitude matrix
tetades1=20;
tetades2=50;
tetades3=25;
a=1;
b=0.5;
cc=0;
kx1=k*sind(tetades1);
kx2=k*sind(tetades2);
kx3=k*sind(tetades3);
x=dx*(1:N);
phase1=exp(1i*kx1*x);
phase11=phase1;
for h=1:N-1
    phase11=cat(1,phase11,phase1);
end
phase2=exp(1i*kx2*x);

phase22=phase2;
for h=1:N-1
    phase22=cat(1,phase22,phase2);
end
phase3=exp(1i*kx3*x);
phase33=phase3;
for h=1:N-1
    phase33=cat(1,phase33,phase3);
end
g1=rot90(phase11);
g11=rot90(rot90(rot90(phase11)));
g2=rot90(phase22);
g22=rot90(rot90(rot90(phase22)));
g3=rot90(phase33);
g33=rot90(rot90(rot90(phase33)));
A=phase1;
NN=M;
for ke=1:NN
    B11(1,ke)=A(1,ke);
end
for i=1:NN-1
    for j=1:NN
        if j==1
            B11(i+1,j+NN-1)=B11(i,j);
        else
        B11(i+1,j-1)=B11(i,j);
        end
    end
end
 
%%nimsaz phase2%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A=phase2;
NN=M;
for ke=1:NN
    B22(1,ke)=A(1,ke);
end
for i=1:NN-1
    for j=1:NN
        if j==1
            B22(i+1,j+NN-1)=B22(i,j);
        else
        B22(i+1,j-1)=B22(i,j);
        end
    end
end
%% nimsaz phase3%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A=phase3;
NN=M;
for ke=1:NN
    B33(1,ke)=A(1,ke);
end
for i=1:NN-1
    for j=1:NN
        if j==1
            B33(i+1,j+NN-1)=B33(i,j);
        else
        B33(i+1,j-1)=B33(i,j);
        end
    end
end
C1=rot90(B11);
V1=rot90(C1);
F1=rot90(V1);
C2=rot90(B22);
V2=rot90(C2);
F2=rot90(V2);
C3=rot90(B33);
V3=rot90(C3);
F3=rot90(V3);

totalph1=(a*phase11)+(b*g22)+(cc*phase33);
MM=max(abs(totalph1(:)));
totalph1=totalph1/max(abs(totalph1(:)));
 AC11=abs(totalph1);
 phimn11=angle(totalph1);
 %%%%%%%%%%%%%%%%%%%%%%3bit
%%%%%%%%%%%%%%%%%%%%%%discrete%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 for i=1:N
     for j=1:N
          if phimn11(i,j)>=-pi && phimn11(i,j)<-3*pi/4
           phimn11(i,j)=-7*pi/8;     
         elseif phimn11(i,j)>=-3*pi/4 && phimn11(i,j)<-pi/2
         phimn11(i,j)=-5*pi/8;
            
          elseif phimn11(i,j)>=-pi/2 && phimn11(i,j)<-pi/4
           phimn11(i,j)=-3*pi/8;
            elseif phimn11(i,j)>=-pi/4 && phimn11(i,j)<0
           phimn11(i,j)=-pi/8;
            elseif phimn11(i,j)>=0 && phimn11(i,j)<pi/4
           phimn11(i,j)=pi/8;
            elseif phimn11(i,j)>=pi/4 && phimn11(i,j)<pi/2
          phimn11(i,j)=3*pi/8;
            elseif phimn11(i,j)>=pi/2 && phimn11(i,j)<3*pi/4
           phimn11(i,j)=5*pi/8;
          else 
           phimn11(i,j)=7*pi/8;
         
          end
     end
 end
  for i=1:N
      for j=1:N
  
          if AC11(i,j)>=0 && AC11(i,j)<0.125
           AC11(i,j)=0.0625;     
         elseif AC11(i,j)>=0.125 && AC11(i,j)<0.25
          AC11(i,j)=0.1875;
          elseif AC11(i,j)>=0.25 && AC11(i,j)<0.375
          AC11(i,j)=0.3125;
           elseif AC11(i,j)>=0.375 && AC11(i,j)<0.5
          AC11(i,j)=0.4375;
           elseif AC11(i,j)>=0.5 && AC11(i,j)<0.625
         AC11(i,j)=0.5625;
           elseif AC11(i,j)>=0.625 && AC11(i,j)<0.75
          AC11(i,j)=0.6875;
            
          elseif AC11(i,j)>=0.75 && AC11(i,j)<0.875
           AC11(i,j)=0.8125;
         
           else
           AC11(i,j)=0.9375;
          end
      end
  end
  AC11=AC11+(1/16);
%   AC11=ones(30);
  phimn11=180/pi.*phimn11+157.5;
  % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  64 %%%%%%%%%%%%%%%%%%%%%
  g=0;
for i=1:N
    for j=1:N
        if AC11(i,j)==1/8
            if phimn11(i,j)==0
                B=[exp(1i*0) exp(1i*0) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2)];
              for n=1:L
            ap=(B(n)/L);
            B1(n)=ap;
              end
        C=sum(B1);
        apq(i,j)=C;
            elseif phimn11(i,j)==45
                 B=[exp(1i*pi/4) exp(1i*pi/4) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2)];
               for n=1:L
            ap=(B(n)/L);
            B1(n)=ap;
              end
        C=sum(B1);
        apq(i,j)=C;
            elseif phimn11(i,j)==90
                 B=[exp(1i*2*pi/4) exp(1i*2*pi/4) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2)];
                    for n=1:L
            ap=(B(n)/L);
            B1(n)=ap;
              end
        C=sum(B1);
        apq(i,j)=C;
            elseif phimn11(i,j)==135
                 B=[exp(1i*3*pi/4) exp(1i*3*pi/4) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2)];
              for n=1:L
            ap=(B(n)/L);
            B1(n)=ap;
              end
        C=sum(B1);
        apq(i,j)=C; 
            elseif phimn11(i,j)==180
                 B=[exp(1i*4*pi/4) exp(1i*4*pi/4) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2)];
                    for n=1:L
            ap=(B(n)/L);
            B1(n)=ap;
              end
        C=sum(B1);
        apq(i,j)=C;
            elseif phimn11(i,j)==225
                 B=[exp(1i*5*pi/4) exp(1i*5*pi/4) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2)];
                   for n=1:L
            ap=(B(n)/L);
            B1(n)=ap;
              end
        C=sum(B1);
        apq(i,j)=C;
            elseif phimn11(i,j)==270
                 B=[exp(1i*6*pi/4) exp(1i*6*pi/4) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2)];
                  for n=1:L
            ap=(B(n)/L);
            B1(n)=ap;
              end
        C=sum(B1);
        apq(i,j)=C;
            else
                 B=[exp(1i*7*pi/4) exp(1i*7*pi/4) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2)];
               for n=1:L
            ap=(B(n)/L);
            B1(n)=ap;
              end
        C=sum(B1);
        apq(i,j)=C;
            end
        end
%%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 if AC11(i,j)==2/8
            if phimn11(i,j)==0
                B=[exp(1i*0*pi/4) exp(1i*0*pi/4) exp(1i*0*pi/4) exp(1i*0*pi/4) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2)];
               for n=1:L
            ap=(B(n)/L);
            B1(n)=ap;
              end
        C=sum(B1);
        apq(i,j)=C;
            elseif phimn11(i,j)==45
                 B=[exp(1i*1*pi/4) exp(1i*1*pi/4) exp(1i*1*pi/4) exp(1i*1*pi/4) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2)];
                for n=1:L
            ap=(B(n)/L);
            B1(n)=ap;
              end
        C=sum(B1);
        apq(i,j)=C;
            elseif phimn11(i,j)==90
                 B=[exp(1i*2*pi/4) exp(1i*2*pi/4) exp(1i*2*pi/4) exp(1i*2*pi/4) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2)];
                   for n=1:L
            ap=(B(n)/L);
            B1(n)=ap;
              end
        C=sum(B1);
        apq(i,j)=C; 
            elseif phimn11(i,j)==135
                 B=[exp(1i*3*pi/4) exp(1i*3*pi/4) exp(1i*3*pi/4) exp(1i*3*pi/4) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2)];
                   for n=1:L
            ap=(B(n)/L);
            B1(n)=ap;
              end
        C=sum(B1);
        apq(i,j)=C;
            elseif phimn11(i,j)==180
                 B=[exp(1i*4*pi/4) exp(1i*4*pi/4) exp(1i*4*pi/4) exp(1i*4*pi/4) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2)];
                    for n=1:L
            ap=(B(n)/L);
            B1(n)=ap;
              end
        C=sum(B1);
        apq(i,j)=C;
            elseif phimn11(i,j)==225
                 B=[exp(1i*5*pi/4) exp(1i*5*pi/4) exp(1i*5*pi/4) exp(1i*5*pi/4) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2)];
                   for n=1:L
            ap=(B(n)/L);
            B1(n)=ap;
              end
        C=sum(B1);
        apq(i,j)=C;
            elseif phimn11(i,j)==270
                 B=[exp(1i*6*pi/4) exp(1i*6*pi/4) exp(1i*6*pi/4) exp(1i*6*pi/4) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2)];
                  for n=1:L
            ap=(B(n)/L);
            B1(n)=ap;
              end
        C=sum(B1);
        apq(i,j)=C;
            else
                 B=[exp(1i*7*pi/4) exp(1i*7*pi/4) exp(1i*7*pi/4) exp(1i*7*pi/4) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2)];
               for n=1:L
            ap=(B(n)/L);
            B1(n)=ap;
              end
        C=sum(B1);
        apq(i,j)=C;
            end
 end
        %%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
         if AC11(i,j)==3/8
            if phimn11(i,j)==0
                B=[exp(1i*0*pi/4) exp(1i*0*pi/4) exp(1i*0*pi/4) exp(1i*0*pi/4) exp(1i*0*pi/4) exp(1i*0*pi/4) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2)];
               for n=1:L
            ap=(B(n)/L);
            B1(n)=ap;
              end
        C=sum(B1);
        apq(i,j)=C;
            elseif phimn11(i,j)==45
                 B=[exp(1i*1*pi/4) exp(1i*1*pi/4) exp(1i*1*pi/4) exp(1i*1*pi/4) exp(1i*1*pi/4) exp(1i*1*pi/4) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2)];
                for n=1:L
            ap=(B(n)/L);
            B1(n)=ap;
              end
        C=sum(B1);
        apq(i,j)=C;
            elseif phimn11(i,j)==90
                 B=[exp(1i*2*pi/4) exp(1i*2*pi/4) exp(1i*2*pi/4) exp(1i*2*pi/4) exp(1i*2*pi/4) exp(1i*2*pi/4) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2)];
                   for n=1:L
            ap=(B(n)/L);
            B1(n)=ap;
              end
        C=sum(B1);
        apq(i,j)=C;
            elseif phimn11(i,j)==135
                 B=[exp(1i*3*pi/4) exp(1i*3*pi/4) exp(1i*3*pi/4) exp(1i*3*pi/4) exp(1i*3*pi/4) exp(1i*3*pi/4) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2)];
                   for n=1:L
            ap=(B(n)/L);
            B1(n)=ap;
              end
        C=sum(B1);
        apq(i,j)=C;
            elseif phimn11(i,j)==180
                 B=[exp(1i*4*pi/4) exp(1i*4*pi/4) exp(1i*4*pi/4) exp(1i*4*pi/4) exp(1i*4*pi/4) exp(1i*4*pi/4) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2)];
                for n=1:L
            ap=(B(n)/L);
            B1(n)=ap;
              end
        C=sum(B1);
        apq(i,j)=C;   
            elseif phimn11(i,j)==225
                 B=[exp(1i*5*pi/4) exp(1i*5*pi/4) exp(1i*5*pi/4) exp(1i*5*pi/4) exp(1i*5*pi/4) exp(1i*5*pi/4) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2)];
                for n=1:L
            ap=(B(n)/L);
            B1(n)=ap;
              end
        C=sum(B1);
        apq(i,j)=C;   
            elseif phimn11(i,j)==270
                 B=[exp(1i*6*pi/4) exp(1i*6*pi/4) exp(1i*6*pi/4) exp(1i*6*pi/4) exp(1i*6*pi/4) exp(1i*6*pi/4) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2)];
                  for n=1:L
            ap=(B(n)/L);
            B1(n)=ap;
              end
        C=sum(B1);
        apq(i,j)=C;
            else
                 B=[exp(1i*7*pi/4) exp(1i*7*pi/4) exp(1i*7*pi/4) exp(1i*7*pi/4) exp(1i*7*pi/4) exp(1i*7*pi/4) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2)];
               for n=1:L
            ap=(B(n)/L);
            B1(n)=ap;
              end
        C=sum(B1);
        apq(i,j)=C;
            end
        end


            %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
                
    if AC11(i,j)==4/8
            if phimn11(i,j)==0
                B=[exp(1i*0*pi/4) exp(1i*0*pi/4) exp(1i*0*pi/4) exp(1i*0*pi/4) exp(1i*0*pi/4) exp(1i*0*pi/4) exp(1i*0*pi/4) exp(1i*0*pi/4) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2)];
               for n=1:L
            ap=(B(n)/L);
            B1(n)=ap;
              end
        C=sum(B1);
        apq(i,j)=C;
            elseif phimn11(i,j)==45
                 B=[exp(1i*1*pi/4) exp(1i*1*pi/4) exp(1i*1*pi/4) exp(1i*1*pi/4) exp(1i*1*pi/4) exp(1i*1*pi/4) exp(1i*1*pi/4) exp(1i*1*pi/4) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2)];
                for n=1:L
            ap=(B(n)/L);
            B1(n)=ap;
              end
        C=sum(B1);
        apq(i,j)=C;
            elseif phimn11(i,j)==90
                 B=[exp(1i*2*pi/4) exp(1i*2*pi/4) exp(1i*2*pi/4) exp(1i*2*pi/4) exp(1i*2*pi/4) exp(1i*2*pi/4) exp(1i*2*pi/4) exp(1i*2*pi/4) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2)];
                for n=1:L
            ap=(B(n)/L);
            B1(n)=ap;
              end
        C=sum(B1);
        apq(i,j)=C;   
            elseif phimn11(i,j)==135
                 B=[exp(1i*3*pi/4) exp(1i*3*pi/4) exp(1i*3*pi/4) exp(1i*3*pi/4) exp(1i*3*pi/4) exp(1i*3*pi/4) exp(1i*3*pi/4) exp(1i*3*pi/4) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2)];
                for n=1:L
            ap=(B(n)/L);
            B1(n)=ap;
              end
        C=sum(B1);
        apq(i,j)=C;    
            elseif phimn11(i,j)==180
                 B=[exp(1i*4*pi/4) exp(1i*4*pi/4) exp(1i*4*pi/4) exp(1i*4*pi/4) exp(1i*4*pi/4) exp(1i*4*pi/4) exp(1i*4*pi/4) exp(1i*4*pi/4) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2)];
                for n=1:L
            ap=(B(n)/L);
            B1(n)=ap;
              end
        C=sum(B1);
        apq(i,j)=C;   
            elseif phimn11(i,j)==225
                 B=[exp(1i*5*pi/4) exp(1i*5*pi/4) exp(1i*5*pi/4) exp(1i*5*pi/4) exp(1i*5*pi/4) exp(1i*5*pi/4) exp(1i*5*pi/4) exp(1i*5*pi/4) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2)];
                for n=1:L
            ap=(B(n)/L);
            B1(n)=ap;
              end
        C=sum(B1);
        apq(i,j)=C;   
            elseif phimn11(i,j)==270
                 B=[exp(1i*6*pi/4) exp(1i*6*pi/4) exp(1i*6*pi/4) exp(1i*6*pi/4) exp(1i*6*pi/4) exp(1i*6*pi/4) exp(1i*6*pi/4) exp(1i*6*pi/4) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2)];
                  for n=1:L
            ap=(B(n)/L);
            B1(n)=ap;
              end
        C=sum(B1);
        apq(i,j)=C;
            else
                 B=[exp(1i*7*pi/4) exp(1i*7*pi/4) exp(1i*7*pi/4) exp(1i*7*pi/4) exp(1i*7*pi/4) exp(1i*7*pi/4) exp(1i*7*pi/4) exp(1i*7*pi/4) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2)];
              for n=1:L
            ap=(B(n)/L);
            B1(n)=ap;
              end
        C=sum(B1);
        apq(i,j)=C;
            end
    end     
        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         if AC11(i,j)==5/8
            if phimn11(i,j)==0
                B=[exp(1i*0*pi/4) exp(1i*0*pi/4) exp(1i*0*pi/4) exp(1i*0*pi/4) exp(1i*0*pi/4) exp(1i*0*pi/4) exp(1i*0*pi/4) exp(1i*0*pi/4) exp(1i*0*pi/4) exp(1i*0*pi/4) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2)];
               for n=1:L
            ap=(B(n)/L);
            B1(n)=ap;
              end
        C=sum(B1);
        apq(i,j)=C;
            elseif phimn11(i,j)==45
                 B=[exp(1i*1*pi/4) exp(1i*1*pi/4) exp(1i*1*pi/4) exp(1i*1*pi/4) exp(1i*1*pi/4) exp(1i*1*pi/4) exp(1i*1*pi/4) exp(1i*1*pi/4) exp(1i*1*pi/4) exp(1i*1*pi/4) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2)];
              for n=1:L
            ap=(B(n)/L);
            B1(n)=ap;
              end
        C=sum(B1);
        apq(i,j)=C;
            elseif phimn11(i,j)==90
                 B=[exp(1i*2*pi/4) exp(1i*2*pi/4) exp(1i*2*pi/4) exp(1i*2*pi/4) exp(1i*2*pi/4) exp(1i*2*pi/4) exp(1i*2*pi/4) exp(1i*2*pi/4) exp(1i*2*pi/4) exp(1i*2*pi/4) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2)];
                for n=1:L
            ap=(B(n)/L);
            B1(n)=ap;
              end
        C=sum(B1);
        apq(i,j)=C;  
            elseif phimn11(i,j)==135
                 B=[exp(1i*3*pi/4) exp(1i*3*pi/4) exp(1i*3*pi/4) exp(1i*3*pi/4) exp(1i*3*pi/4) exp(1i*3*pi/4) exp(1i*3*pi/4) exp(1i*3*pi/4) exp(1i*3*pi/4) exp(1i*3*pi/4) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2)];
                for n=1:L
            ap=(B(n)/L);
            B1(n)=ap;
              end
        C=sum(B1);
        apq(i,j)=C;  
            elseif phimn11(i,j)==180
                 B=[exp(1i*4*pi/4) exp(1i*4*pi/4) exp(1i*4*pi/4) exp(1i*4*pi/4) exp(1i*4*pi/4) exp(1i*4*pi/4) exp(1i*4*pi/4) exp(1i*4*pi/4) exp(1i*4*pi/4) exp(1i*4*pi/4) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2)];
                for n=1:L
            ap=(B(n)/L);
            B1(n)=ap;
              end
        C=sum(B1);
        apq(i,j)=C;  
            elseif phimn11(i,j)==225
                 B=[exp(1i*5*pi/4) exp(1i*5*pi/4) exp(1i*5*pi/4) exp(1i*5*pi/4) exp(1i*5*pi/4) exp(1i*5*pi/4) exp(1i*5*pi/4) exp(1i*5*pi/4) exp(1i*5*pi/4) exp(1i*5*pi/4) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2)];
                for n=1:L
            ap=(B(n)/L);
            B1(n)=ap;
              end
        C=sum(B1);
        apq(i,j)=C;
            elseif phimn11(i,j)==270
                 B=[exp(1i*6*pi/4) exp(1i*6*pi/4) exp(1i*6*pi/4) exp(1i*6*pi/4) exp(1i*6*pi/4) exp(1i*6*pi/4) exp(1i*6*pi/4) exp(1i*6*pi/4) exp(1i*6*pi/4) exp(1i*6*pi/4) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2)];
                 for n=1:L
            ap=(B(n)/L);
            B1(n)=ap;
              end
        C=sum(B1);
        apq(i,j)=C; 
            else
                 B=[exp(1i*7*pi/4) exp(1i*7*pi/4) exp(1i*7*pi/4) exp(1i*7*pi/4) exp(1i*7*pi/4) exp(1i*7*pi/4) exp(1i*7*pi/4) exp(1i*7*pi/4) exp(1i*7*pi/4) exp(1i*7*pi/4) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2)];
              for n=1:L
            ap=(B(n)/L);
            B1(n)=ap;
              end
        C=sum(B1);
        apq(i,j)=C;
            end
         end
        %%  5%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
         if AC11(i,j)==6/8
            if phimn11(i,j)==0
                B=[exp(1i*0*pi/4) exp(1i*0*pi/4) exp(1i*0*pi/4) exp(1i*0*pi/4) exp(1i*0*pi/4) exp(1i*0*pi/4) exp(1i*0*pi/4) exp(1i*0*pi/4) exp(1i*0*pi/4) exp(1i*0*pi/4) exp(1i*0*pi/4) exp(1i*0*pi/4) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2)];
              for n=1:L
            ap=(B(n)/L);
            B1(n)=ap;
              end
        C=sum(B1);
        apq(i,j)=C;
            elseif phimn11(i,j)==45
                 B=[exp(1i*1*pi/4) exp(1i*1*pi/4) exp(1i*1*pi/4) exp(1i*1*pi/4) exp(1i*1*pi/4) exp(1i*1*pi/4) exp(1i*1*pi/4) exp(1i*1*pi/4) exp(1i*1*pi/4) exp(1i*1*pi/4) exp(1i*1*pi/4) exp(1i*1*pi/4) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2)];
               for n=1:L
            ap=(B(n)/L);
            B1(n)=ap;
              end
        C=sum(B1);
        apq(i,j)=C;
            elseif phimn11(i,j)==90
                 B=[exp(1i*2*pi/4) exp(1i*2*pi/4) exp(1i*2*pi/4) exp(1i*2*pi/4) exp(1i*2*pi/4) exp(1i*2*pi/4) exp(1i*2*pi/4) exp(1i*2*pi/4) exp(1i*2*pi/4) exp(1i*2*pi/4) exp(1i*2*pi/4) exp(1i*2*pi/4) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2)];
                for n=1:L
            ap=(B(n)/L);
            B1(n)=ap;
              end
        C=sum(B1);
        apq(i,j)=C;  
            elseif phimn11(i,j)==135
                 B=[exp(1i*3*pi/4) exp(1i*3*pi/4) exp(1i*3*pi/4) exp(1i*3*pi/4) exp(1i*3*pi/4) exp(1i*3*pi/4) exp(1i*3*pi/4) exp(1i*3*pi/4) exp(1i*3*pi/4) exp(1i*3*pi/4) exp(1i*3*pi/4) exp(1i*3*pi/4) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2)];
                for n=1:L
            ap=(B(n)/L);
            B1(n)=ap;
              end
        C=sum(B1);
        apq(i,j)=C;  
            elseif phimn11(i,j)==180
                 B=[exp(1i*4*pi/4) exp(1i*4*pi/4) exp(1i*4*pi/4) exp(1i*4*pi/4) exp(1i*4*pi/4) exp(1i*4*pi/4) exp(1i*4*pi/4) exp(1i*4*pi/4) exp(1i*4*pi/4) exp(1i*4*pi/4) exp(1i*4*pi/4) exp(1i*4*pi/4) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2)];
                for n=1:L
            ap=(B(n)/L);
            B1(n)=ap;
              end
        C=sum(B1);
        apq(i,j)=C;   
            elseif phimn11(i,j)==225
                 B=[exp(1i*5*pi/4) exp(1i*5*pi/4) exp(1i*5*pi/4) exp(1i*5*pi/4) exp(1i*5*pi/4) exp(1i*5*pi/4) exp(1i*5*pi/4) exp(1i*5*pi/4) exp(1i*5*pi/4) exp(1i*5*pi/4) exp(1i*5*pi/4) exp(1i*5*pi/4) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2)];
                for n=1:L
            ap=(B(n)/L);
            B1(n)=ap;
              end
        C=sum(B1);
        apq(i,j)=C; 
            elseif phimn11(i,j)==270
                 B=[exp(1i*6*pi/4) exp(1i*6*pi/4) exp(1i*6*pi/4) exp(1i*6*pi/4) exp(1i*6*pi/4) exp(1i*6*pi/4) exp(1i*6*pi/4) exp(1i*6*pi/4) exp(1i*6*pi/4) exp(1i*6*pi/4) exp(1i*6*pi/4) exp(1i*6*pi/4) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2)];
                  for n=1:L
            ap=(B(n)/L);
            B1(n)=ap;
              end
        C=sum(B1);
        apq(i,j)=C;
            else
                 B=[exp(1i*7*pi/4) exp(1i*7*pi/4) exp(1i*7*pi/4) exp(1i*7*pi/4) exp(1i*7*pi/4) exp(1i*7*pi/4) exp(1i*7*pi/4) exp(1i*7*pi/4) exp(1i*7*pi/4) exp(1i*7*pi/4) exp(1i*7*pi/4) exp(1i*7*pi/4) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2)];
             for n=1:L
            ap=(B(n)/L);
            B1(n)=ap;
              end
        C=sum(B1);
        apq(i,j)=C;
            end
         end
        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         if AC11(i,j)==7/8
            if phimn11(i,j)==0
                B=[exp(1i*0*pi/4) exp(1i*0*pi/4) exp(1i*0*pi/4) exp(1i*0*pi/4) exp(1i*0*pi/4) exp(1i*0*pi/4) exp(1i*0*pi/4) exp(1i*0*pi/4) exp(1i*0*pi/4) exp(1i*0*pi/4) exp(1i*0*pi/4) exp(1i*0*pi/4) exp(1i*0*pi/4) exp(1i*0*pi/4) exp(1i*1*pi/2) exp(1i*3*pi/2)];
              for n=1:L
            ap=(B(n)/L);
            B1(n)=ap;
              end
        C=sum(B1);
        apq(i,j)=C;
            elseif phimn11(i,j)==45
                 B=[exp(1i*1*pi/4) exp(1i*1*pi/4) exp(1i*1*pi/4) exp(1i*1*pi/4) exp(1i*1*pi/4) exp(1i*1*pi/4) exp(1i*1*pi/4) exp(1i*1*pi/4) exp(1i*1*pi/4) exp(1i*1*pi/4) exp(1i*1*pi/4) exp(1i*1*pi/4) exp(1i*1*pi/4) exp(1i*1*pi/4) exp(1i*1*pi/2) exp(1i*3*pi/2)];
               for n=1:L
            ap=(B(n)/L);
            B1(n)=ap;
              end
        C=sum(B1);
        apq(i,j)=C;
            elseif phimn11(i,j)==90
                 B=[exp(1i*2*pi/4) exp(1i*2*pi/4) exp(1i*2*pi/4) exp(1i*2*pi/4) exp(1i*2*pi/4) exp(1i*2*pi/4) exp(1i*2*pi/4) exp(1i*2*pi/4) exp(1i*2*pi/4) exp(1i*2*pi/4) exp(1i*2*pi/4) exp(1i*2*pi/4) exp(1i*2*pi/4) exp(1i*2*pi/4) exp(1i*1*pi/2) exp(1i*3*pi/2)];
                for n=1:L
            ap=(B(n)/L);
            B1(n)=ap;
              end
        C=sum(B1);
        apq(i,j)=C; 
            elseif phimn11(i,j)==135
                 B=[exp(1i*3*pi/4) exp(1i*3*pi/4) exp(1i*3*pi/4) exp(1i*3*pi/4) exp(1i*3*pi/4) exp(1i*3*pi/4) exp(1i*3*pi/4) exp(1i*3*pi/4) exp(1i*3*pi/4) exp(1i*3*pi/4) exp(1i*3*pi/4) exp(1i*3*pi/4) exp(1i*3*pi/4) exp(1i*3*pi/4) exp(1i*1*pi/2) exp(1i*3*pi/2)];
                for n=1:L
            ap=(B(n)/L);
            B1(n)=ap;
              end
        C=sum(B1);
        apq(i,j)=C;  
            elseif phimn11(i,j)==180
                 B=[exp(1i*4*pi/4) exp(1i*4*pi/4) exp(1i*4*pi/4) exp(1i*4*pi/4) exp(1i*4*pi/4) exp(1i*4*pi/4) exp(1i*4*pi/4) exp(1i*4*pi/4) exp(1i*4*pi/4) exp(1i*4*pi/4) exp(1i*4*pi/4) exp(1i*0*pi/4) exp(1i*4*pi/4) exp(1i*4*pi/4) exp(1i*1*pi/2) exp(1i*3*pi/2)];
                for n=1:L
            ap=(B(n)/L);
            B1(n)=ap;
              end
        C=sum(B1);
        apq(i,j)=C;  
            elseif phimn11(i,j)==225
                 B=[exp(1i*5*pi/4) exp(1i*5*pi/4) exp(1i*5*pi/4) exp(1i*5*pi/4) exp(1i*5*pi/4) exp(1i*5*pi/4) exp(1i*5*pi/4) exp(1i*5*pi/4) exp(1i*5*pi/4) exp(1i*5*pi/4) exp(1i*5*pi/4) exp(1i*5*pi/4) exp(1i*5*pi/4) exp(1i*5*pi/4) exp(1i*1*pi/2) exp(1i*3*pi/2)];
                for n=1:L
            ap=(B(n)/L);
            B1(n)=ap;
              end
        C=sum(B1);
        apq(i,j)=C;   
            elseif phimn11(i,j)==270
                 B=[exp(1i*6*pi/4) exp(1i*6*pi/4) exp(1i*6*pi/4) exp(1i*6*pi/4) exp(1i*6*pi/4) exp(1i*6*pi/4) exp(1i*6*pi/4) exp(1i*6*pi/4) exp(1i*6*pi/4) exp(1i*6*pi/4) exp(1i*6*pi/4) exp(1i*6*pi/4) exp(1i*6*pi/4) exp(1i*6*pi/4) exp(1i*1*pi/2) exp(1i*3*pi/2)];
                  for n=1:L
            ap=(B(n)/L);
            B1(n)=ap;
              end
        C=sum(B1);
        apq(i,j)=C;
            else
                 B=[exp(1i*7*pi/4) exp(1i*7*pi/4) exp(1i*7*pi/4) exp(1i*7*pi/4) exp(1i*7*pi/4) exp(1i*7*pi/4) exp(1i*7*pi/4) exp(1i*7*pi/4) exp(1i*7*pi/4) exp(1i*7*pi/4) exp(1i*7*pi/4) exp(1i*7*pi/4) exp(1i*7*pi/4) exp(1i*7*pi/4) exp(1i*1*pi/2) exp(1i*3*pi/2)];
               for n=1:L
            ap=(B(n)/L);
            B1(n)=ap;
              end
        C=sum(B1);
        apq(i,j)=C;
            end
         end
        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         if AC11(i,j)==8/8
            if phimn11(i,j)==0
                B=[exp(1i*0*pi/4) exp(1i*0*pi/4) exp(1i*0*pi/4) exp(1i*0*pi/4) exp(1i*0*pi/4) exp(1i*0*pi/4) exp(1i*0*pi/4) exp(1i*0*pi/4) exp(1i*0*pi/4) exp(1i*0*pi/4) exp(1i*0*pi/4) exp(1i*0*pi/4) exp(1i*0*pi/4) exp(1i*0*pi/4) exp(1i*0*pi/4) exp(1i*0*pi/4)];
              for n=1:L
            ap=(B(n)/L);
            B1(n)=ap;
              end
        C=sum(B1);
        apq(i,j)=C;
            elseif phimn11(i,j)==45
                 B=[exp(1i*1*pi/4) exp(1i*1*pi/4) exp(1i*1*pi/4) exp(1i*1*pi/4) exp(1i*1*pi/4) exp(1i*1*pi/4) exp(1i*1*pi/4) exp(1i*1*pi/4) exp(1i*1*pi/4) exp(1i*1*pi/4) exp(1i*1*pi/4) exp(1i*1*pi/4) exp(1i*1*pi/4) exp(1i*1*pi/4) exp(1i*1*pi/4) exp(1i*1*pi/4)];
               for n=1:L
            ap=(B(n)/L);
            B1(n)=ap;
              end
        C=sum(B1);
        apq(i,j)=C;
            elseif phimn11(i,j)==90
                 B=[exp(1i*2*pi/4) exp(1i*2*pi/4) exp(1i*2*pi/4) exp(1i*2*pi/4) exp(1i*2*pi/4) exp(1i*2*pi/4) exp(1i*2*pi/4) exp(1i*2*pi/4) exp(1i*2*pi/4) exp(1i*2*pi/4) exp(1i*2*pi/4) exp(1i*2*pi/4) exp(1i*2*pi/4) exp(1i*2*pi/4) exp(1i*2*pi/4) exp(1i*2*pi/4)];
                for n=1:L
            ap=(B(n)/L);
            B1(n)=ap;
              end
        C=sum(B1);
        apq(i,j)=C; 
            elseif phimn11(i,j)==135
                 B=[exp(1i*3*pi/4) exp(1i*3*pi/4) exp(1i*3*pi/4) exp(1i*3*pi/4) exp(1i*3*pi/4) exp(1i*3*pi/4) exp(1i*3*pi/4) exp(1i*3*pi/4) exp(1i*3*pi/4) exp(1i*3*pi/4) exp(1i*3*pi/4) exp(1i*3*pi/4) exp(1i*3*pi/4) exp(1i*3*pi/4) exp(1i*3*pi/4) exp(1i*3*pi/4)];
                for n=1:L
            ap=(B(n)/L);
            B1(n)=ap;
              end
        C=sum(B1);
        apq(i,j)=C;  
            elseif phimn11(i,j)==180
                 B=[exp(1i*4*pi/4) exp(1i*4*pi/4) exp(1i*4*pi/4) exp(1i*4*pi/4) exp(1i*4*pi/4) exp(1i*4*pi/4) exp(1i*4*pi/4) exp(1i*4*pi/4) exp(1i*4*pi/4) exp(1i*4*pi/4) exp(1i*4*pi/4) exp(1i*4*pi/4) exp(1i*4*pi/4) exp(1i*4*pi/4) exp(1i*4*pi/4) exp(1i*4*pi/4)];
                for n=1:L
            ap=(B(n)/L);
            B1(n)=ap;
              end
        C=sum(B1);
        apq(i,j)=C;  
            elseif phimn11(i,j)==225
                 B=[exp(1i*5*pi/4) exp(1i*5*pi/4) exp(1i*5*pi/4) exp(1i*5*pi/4) exp(1i*5*pi/4) exp(1i*5*pi/4) exp(1i*5*pi/4) exp(1i*5*pi/4) exp(1i*5*pi/4) exp(1i*5*pi/4) exp(1i*5*pi/4) exp(1i*5*pi/4) exp(1i*5*pi/4) exp(1i*5*pi/4) exp(1i*5*pi/4) exp(1i*5*pi/4)];
                for n=1:L
            ap=(B(n)/L);
            B1(n)=ap;
              end
        C=sum(B1);
        apq(i,j)=C;   
            elseif phimn11(i,j)==270
                 B=[exp(1i*6*pi/4) exp(1i*6*pi/4) exp(1i*6*pi/4) exp(1i*6*pi/4) exp(1i*6*pi/4) exp(1i*6*pi/4) exp(1i*6*pi/4) exp(1i*6*pi/4) exp(1i*6*pi/4) exp(1i*6*pi/4) exp(1i*6*pi/4) exp(1i*6*pi/4) exp(1i*6*pi/4) exp(1i*6*pi/4) exp(1i*6*pi/4) exp(1i*6*pi/4)];
                  for n=1:L
            ap=(B(n)/L);
            B1(n)=ap;
              end
        C=sum(B1);
        apq(i,j)=C;
            else
                 B=[exp(1i*7*pi/4) exp(1i*7*pi/4) exp(1i*7*pi/4) exp(1i*7*pi/4) exp(1i*7*pi/4) exp(1i*7*pi/4) exp(1i*7*pi/4) exp(1i*7*pi/4) exp(1i*7*pi/4) exp(1i*7*pi/4) exp(1i*7*pi/4) exp(1i*7*pi/4) exp(1i*7*pi/4) exp(1i*7*pi/4) exp(1i*7*pi/4) exp(1i*7*pi/4)];
               for n=1:L
            ap=(B(n)/L);
            B1(n)=ap;
              end
        C=sum(B1);
        apq(i,j)=C;
            end
         end
    end
end
AC1=abs(apq)*MM;
phimn1=angle(apq);
k=(2*pi*(f+g*f0))/c;
for u=1:length(phi)
    for v=1:length(theta)
        E=0;
        for m=1:M
            for n=1:N      
                E=E+AC1(m,n)*exp(+1i*k*dx*(n-1)*sin(theta(v))*cos(phi(u))+1i*k*dy*(m-1)*sin(theta(v))*sin(phi(u))+1i*phimn1(m,n));
            end
        end
%         AF(u,v)=E*cos(theta(v));
        AF1(u,v)=E;
    end
end
% %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Directivity
Prad=0;
for u=1:length(phi)
    for v=1:length(theta)
        Prad=Prad+(abs(AF1(u,v)))^2*sin(theta(v));
    end
end
Pradiation0=dphi*dtheta*Prad;
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%64 64 64 64  g=======,,,,,,,====1,2,3...
for g=1:G
for i=1:N
    for j=1:N
        if AC11(i,j)==1/8
            if phimn11(i,j)==0
                B=[exp(1i*0) exp(1i*0) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2)];
              for n=1:L
            ap=(B(n)/L)*sin(g*pi/L)/(g*pi/L)*exp(-1i*pi*g*((2*n-1))/L);
            B1(n)=ap;
              end
        C=sum(B1);
        apq(i,j)=C;
            elseif phimn11(i,j)==45
                 B=[exp(1i*pi/4) exp(1i*pi/4) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2)];
               for n=1:L
            ap=(B(n)/L)*sin(g*pi/L)/(g*pi/L)*exp(-1i*pi*g*((2*n-1))/L);
            B1(n)=ap;
              end
        C=sum(B1);
        apq(i,j)=C;
            elseif phimn11(i,j)==90
                 B=[exp(1i*2*pi/4) exp(1i*2*pi/4) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2)];
                    for n=1:L
            ap=(B(n)/L)*sin(g*pi/L)/(g*pi/L)*exp(-1i*pi*g*((2*n-1))/L);
            B1(n)=ap;
              end
        C=sum(B1);
        apq(i,j)=C;
            elseif phimn11(i,j)==135
                 B=[exp(1i*3*pi/4) exp(1i*3*pi/4) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2)];
              for n=1:L
           ap=(B(n)/L)*sin(g*pi/L)/(g*pi/L)*exp(-1i*pi*g*((2*n-1))/L);
            B1(n)=ap;
              end
        C=sum(B1);
        apq(i,j)=C; 
            elseif phimn11(i,j)==180
                 B=[exp(1i*4*pi/4) exp(1i*4*pi/4) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2)];
                    for n=1:L
            ap=(B(n)/L)*sin(g*pi/L)/(g*pi/L)*exp(-1i*pi*g*((2*n-1))/L);
            B1(n)=ap;
              end
        C=sum(B1);
        apq(i,j)=C;
            elseif phimn11(i,j)==225
                 B=[exp(1i*5*pi/4) exp(1i*5*pi/4) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2)];
                   for n=1:L
            ap=(B(n)/L)*sin(g*pi/L)/(g*pi/L)*exp(-1i*pi*g*((2*n-1))/L);
            B1(n)=ap;
              end
        C=sum(B1);
        apq(i,j)=C;
            elseif phimn11(i,j)==270
                 B=[exp(1i*6*pi/4) exp(1i*6*pi/4) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2)];
                  for n=1:L
            ap=(B(n)/L)*sin(g*pi/L)/(g*pi/L)*exp(-1i*pi*g*((2*n-1))/L);
            B1(n)=ap;
              end
        C=sum(B1);
        apq(i,j)=C;
            else
                 B=[exp(1i*7*pi/4) exp(1i*7*pi/4) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2)];
               for n=1:L
            ap=(B(n)/L)*sin(g*pi/L)/(g*pi/L)*exp(-1i*pi*g*((2*n-1))/L);
            B1(n)=ap;
              end
        C=sum(B1);
        apq(i,j)=C;
            end
        end
%%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 if AC11(i,j)==2/8
            if phimn11(i,j)==0
                B=[exp(1i*0*pi/4) exp(1i*0*pi/4) exp(1i*0*pi/4) exp(1i*0*pi/4) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2)];
               for n=1:L
            ap=(B(n)/L)*sin(g*pi/L)/(g*pi/L)*exp(-1i*pi*g*((2*n-1))/L);
            B1(n)=ap;
              end
        C=sum(B1);
        apq(i,j)=C;
            elseif phimn11(i,j)==45
                 B=[exp(1i*1*pi/4) exp(1i*1*pi/4) exp(1i*1*pi/4) exp(1i*1*pi/4) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2)];
                for n=1:L
            ap=(B(n)/L)*sin(g*pi/L)/(g*pi/L)*exp(-1i*pi*g*((2*n-1))/L);
            B1(n)=ap;
              end
        C=sum(B1);
        apq(i,j)=C;
            elseif phimn11(i,j)==90
                 B=[exp(1i*2*pi/4) exp(1i*2*pi/4) exp(1i*2*pi/4) exp(1i*2*pi/4) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2)];
                   for n=1:L
            ap=(B(n)/L)*sin(g*pi/L)/(g*pi/L)*exp(-1i*pi*g*((2*n-1))/L);
            B1(n)=ap;
              end
        C=sum(B1);
        apq(i,j)=C; 
            elseif phimn11(i,j)==135
                 B=[exp(1i*3*pi/4) exp(1i*3*pi/4) exp(1i*3*pi/4) exp(1i*3*pi/4) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2)];
                   for n=1:L
            ap=(B(n)/L)*sin(g*pi/L)/(g*pi/L)*exp(-1i*pi*g*((2*n-1))/L);
            B1(n)=ap;
              end
        C=sum(B1);
        apq(i,j)=C;
            elseif phimn11(i,j)==180
                 B=[exp(1i*4*pi/4) exp(1i*4*pi/4) exp(1i*4*pi/4) exp(1i*4*pi/4) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2)];
                    for n=1:L
            ap=(B(n)/L)*sin(g*pi/L)/(g*pi/L)*exp(-1i*pi*g*((2*n-1))/L);
            B1(n)=ap;
              end
        C=sum(B1);
        apq(i,j)=C;
            elseif phimn11(i,j)==225
                 B=[exp(1i*5*pi/4) exp(1i*5*pi/4) exp(1i*5*pi/4) exp(1i*5*pi/4) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2)];
                   for n=1:L
            ap=(B(n)/L)*sin(g*pi/L)/(g*pi/L)*exp(-1i*pi*g*((2*n-1))/L);
            B1(n)=ap;
              end
        C=sum(B1);
        apq(i,j)=C;
            elseif phimn11(i,j)==270
                 B=[exp(1i*6*pi/4) exp(1i*6*pi/4) exp(1i*6*pi/4) exp(1i*6*pi/4) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2)];
                  for n=1:L
            ap=(B(n)/L)*sin(g*pi/L)/(g*pi/L)*exp(-1i*pi*g*((2*n-1))/L);
            B1(n)=ap;
              end
        C=sum(B1);
        apq(i,j)=C;
            else
                 B=[exp(1i*7*pi/4) exp(1i*7*pi/4) exp(1i*7*pi/4) exp(1i*7*pi/4) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2)];
               for n=1:L
            ap=(B(n)/L)*sin(g*pi/L)/(g*pi/L)*exp(-1i*pi*g*((2*n-1))/L);
            B1(n)=ap;
              end
        C=sum(B1);
        apq(i,j)=C;
            end
 end
        %%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
         if AC11(i,j)==3/8
            if phimn11(i,j)==0
                B=[exp(1i*0*pi/4) exp(1i*0*pi/4) exp(1i*0*pi/4) exp(1i*0*pi/4) exp(1i*0*pi/4) exp(1i*0*pi/4) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2)];
               for n=1:L
            ap=(B(n)/L)*sin(g*pi/L)/(g*pi/L)*exp(-1i*pi*g*((2*n-1))/L);
            B1(n)=ap;
              end
        C=sum(B1);
        apq(i,j)=C;
            elseif phimn11(i,j)==45
                 B=[exp(1i*1*pi/4) exp(1i*1*pi/4) exp(1i*1*pi/4) exp(1i*1*pi/4) exp(1i*1*pi/4) exp(1i*1*pi/4) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2)];
                for n=1:L
            ap=(B(n)/L)*sin(g*pi/L)/(g*pi/L)*exp(-1i*pi*g*((2*n-1))/L);
            B1(n)=ap;
              end
        C=sum(B1);
        apq(i,j)=C;
            elseif phimn11(i,j)==90
                 B=[exp(1i*2*pi/4) exp(1i*2*pi/4) exp(1i*2*pi/4) exp(1i*2*pi/4) exp(1i*2*pi/4) exp(1i*2*pi/4) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2)];
                   for n=1:L
            ap=(B(n)/L)*sin(g*pi/L)/(g*pi/L)*exp(-1i*pi*g*((2*n-1))/L);
            B1(n)=ap;
              end
        C=sum(B1);
        apq(i,j)=C;
            elseif phimn11(i,j)==135
                 B=[exp(1i*3*pi/4) exp(1i*3*pi/4) exp(1i*3*pi/4) exp(1i*3*pi/4) exp(1i*3*pi/4) exp(1i*3*pi/4) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2)];
                   for n=1:L
            ap=(B(n)/L)*sin(g*pi/L)/(g*pi/L)*exp(-1i*pi*g*((2*n-1))/L);
            B1(n)=ap;
              end
        C=sum(B1);
        apq(i,j)=C;
            elseif phimn11(i,j)==180
                 B=[exp(1i*4*pi/4) exp(1i*4*pi/4) exp(1i*4*pi/4) exp(1i*4*pi/4) exp(1i*4*pi/4) exp(1i*4*pi/4) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2)];
                for n=1:L
            ap=(B(n)/L)*sin(g*pi/L)/(g*pi/L)*exp(-1i*pi*g*((2*n-1))/L);
            B1(n)=ap;
              end
        C=sum(B1);
        apq(i,j)=C;   
            elseif phimn11(i,j)==225
                 B=[exp(1i*5*pi/4) exp(1i*5*pi/4) exp(1i*5*pi/4) exp(1i*5*pi/4) exp(1i*5*pi/4) exp(1i*5*pi/4) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2)];
                for n=1:L
            ap=(B(n)/L)*sin(g*pi/L)/(g*pi/L)*exp(-1i*pi*g*((2*n-1))/L);
            B1(n)=ap;
              end
        C=sum(B1);
        apq(i,j)=C;   
            elseif phimn11(i,j)==270
                 B=[exp(1i*6*pi/4) exp(1i*6*pi/4) exp(1i*6*pi/4) exp(1i*6*pi/4) exp(1i*6*pi/4) exp(1i*6*pi/4) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2)];
                  for n=1:L
            ap=(B(n)/L)*sin(g*pi/L)/(g*pi/L)*exp(-1i*pi*g*((2*n-1))/L);
            B1(n)=ap;
              end
        C=sum(B1);
        apq(i,j)=C;
            else
                 B=[exp(1i*7*pi/4) exp(1i*7*pi/4) exp(1i*7*pi/4) exp(1i*7*pi/4) exp(1i*7*pi/4) exp(1i*7*pi/4) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2)];
               for n=1:L
            ap=(B(n)/L)*sin(g*pi/L)/(g*pi/L)*exp(-1i*pi*g*((2*n-1))/L);
            B1(n)=ap;
              end
        C=sum(B1);
        apq(i,j)=C;
            end
        end


            %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
                
    if AC11(i,j)==4/8
            if phimn11(i,j)==0
                B=[exp(1i*0*pi/4) exp(1i*0*pi/4) exp(1i*0*pi/4) exp(1i*0*pi/4) exp(1i*0*pi/4) exp(1i*0*pi/4) exp(1i*0*pi/4) exp(1i*0*pi/4) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2)];
               for n=1:L
            ap=(B(n)/L)*sin(g*pi/L)/(g*pi/L)*exp(-1i*pi*g*((2*n-1))/L);
            B1(n)=ap;
              end
        C=sum(B1);
        apq(i,j)=C;
            elseif phimn11(i,j)==45
                 B=[exp(1i*1*pi/4) exp(1i*1*pi/4) exp(1i*1*pi/4) exp(1i*1*pi/4) exp(1i*1*pi/4) exp(1i*1*pi/4) exp(1i*1*pi/4) exp(1i*1*pi/4) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2)];
                for n=1:L
            ap=(B(n)/L)*sin(g*pi/L)/(g*pi/L)*exp(-1i*pi*g*((2*n-1))/L);
            B1(n)=ap;
              end
        C=sum(B1);
        apq(i,j)=C;
            elseif phimn11(i,j)==90
                 B=[exp(1i*2*pi/4) exp(1i*2*pi/4) exp(1i*2*pi/4) exp(1i*2*pi/4) exp(1i*2*pi/4) exp(1i*2*pi/4) exp(1i*2*pi/4) exp(1i*2*pi/4) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2)];
                for n=1:L
            ap=(B(n)/L)*sin(g*pi/L)/(g*pi/L)*exp(-1i*pi*g*((2*n-1))/L);
            B1(n)=ap;
              end
        C=sum(B1);
        apq(i,j)=C;   
            elseif phimn11(i,j)==135
                 B=[exp(1i*3*pi/4) exp(1i*3*pi/4) exp(1i*3*pi/4) exp(1i*3*pi/4) exp(1i*3*pi/4) exp(1i*3*pi/4) exp(1i*3*pi/4) exp(1i*3*pi/4) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2)];
                for n=1:L
            ap=(B(n)/L)*sin(g*pi/L)/(g*pi/L)*exp(-1i*pi*g*((2*n-1))/L);
            B1(n)=ap;
              end
        C=sum(B1);
        apq(i,j)=C;    
            elseif phimn11(i,j)==180
                 B=[exp(1i*4*pi/4) exp(1i*4*pi/4) exp(1i*4*pi/4) exp(1i*4*pi/4) exp(1i*4*pi/4) exp(1i*4*pi/4) exp(1i*4*pi/4) exp(1i*4*pi/4) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2)];
                for n=1:L
            ap=(B(n)/L)*sin(g*pi/L)/(g*pi/L)*exp(-1i*pi*g*((2*n-1))/L);
            B1(n)=ap;
              end
        C=sum(B1);
        apq(i,j)=C;   
            elseif phimn11(i,j)==225
                 B=[exp(1i*5*pi/4) exp(1i*5*pi/4) exp(1i*5*pi/4) exp(1i*5*pi/4) exp(1i*5*pi/4) exp(1i*5*pi/4) exp(1i*5*pi/4) exp(1i*5*pi/4) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2)];
                for n=1:L
            ap=(B(n)/L)*sin(g*pi/L)/(g*pi/L)*exp(-1i*pi*g*((2*n-1))/L);
            B1(n)=ap;
              end
        C=sum(B1);
        apq(i,j)=C;   
            elseif phimn11(i,j)==270
                 B=[exp(1i*6*pi/4) exp(1i*6*pi/4) exp(1i*6*pi/4) exp(1i*6*pi/4) exp(1i*6*pi/4) exp(1i*6*pi/4) exp(1i*6*pi/4) exp(1i*6*pi/4) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2)];
                  for n=1:L
            ap=(B(n)/L)*sin(g*pi/L)/(g*pi/L)*exp(-1i*pi*g*((2*n-1))/L);
            B1(n)=ap;
              end
        C=sum(B1);
        apq(i,j)=C;
            else
                 B=[exp(1i*7*pi/4) exp(1i*7*pi/4) exp(1i*7*pi/4) exp(1i*7*pi/4) exp(1i*7*pi/4) exp(1i*7*pi/4) exp(1i*7*pi/4) exp(1i*7*pi/4) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2)];
              for n=1:L
           ap=(B(n)/L)*sin(g*pi/L)/(g*pi/L)*exp(-1i*pi*g*((2*n-1))/L);
            B1(n)=ap;
              end
        C=sum(B1);
        apq(i,j)=C;
            end
    end     
        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         if AC11(i,j)==5/8
            if phimn11(i,j)==0
                B=[exp(1i*0*pi/4) exp(1i*0*pi/4) exp(1i*0*pi/4) exp(1i*0*pi/4) exp(1i*0*pi/4) exp(1i*0*pi/4) exp(1i*0*pi/4) exp(1i*0*pi/4) exp(1i*0*pi/4) exp(1i*0*pi/4) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2)];
               for n=1:L
           ap=(B(n)/L)*sin(g*pi/L)/(g*pi/L)*exp(-1i*pi*g*((2*n-1))/L);
            B1(n)=ap;
              end
        C=sum(B1);
        apq(i,j)=C;
            elseif phimn11(i,j)==45
                 B=[exp(1i*1*pi/4) exp(1i*1*pi/4) exp(1i*1*pi/4) exp(1i*1*pi/4) exp(1i*1*pi/4) exp(1i*1*pi/4) exp(1i*1*pi/4) exp(1i*1*pi/4) exp(1i*1*pi/4) exp(1i*1*pi/4) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2)];
              for n=1:L
            ap=(B(n)/L)*sin(g*pi/L)/(g*pi/L)*exp(-1i*pi*g*((2*n-1))/L);
            B1(n)=ap;
              end
        C=sum(B1);
        apq(i,j)=C;
            elseif phimn11(i,j)==90
                 B=[exp(1i*2*pi/4) exp(1i*2*pi/4) exp(1i*2*pi/4) exp(1i*2*pi/4) exp(1i*2*pi/4) exp(1i*2*pi/4) exp(1i*2*pi/4) exp(1i*2*pi/4) exp(1i*2*pi/4) exp(1i*2*pi/4) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2)];
                for n=1:L
           ap=(B(n)/L)*sin(g*pi/L)/(g*pi/L)*exp(-1i*pi*g*((2*n-1))/L);
            B1(n)=ap;
              end
        C=sum(B1);
        apq(i,j)=C;  
            elseif phimn11(i,j)==135
                 B=[exp(1i*3*pi/4) exp(1i*3*pi/4) exp(1i*3*pi/4) exp(1i*3*pi/4) exp(1i*3*pi/4) exp(1i*3*pi/4) exp(1i*3*pi/4) exp(1i*3*pi/4) exp(1i*3*pi/4) exp(1i*3*pi/4) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2)];
                for n=1:L
            ap=(B(n)/L)*sin(g*pi/L)/(g*pi/L)*exp(-1i*pi*g*((2*n-1))/L);
            B1(n)=ap;
              end
        C=sum(B1);
        apq(i,j)=C;  
            elseif phimn11(i,j)==180
                 B=[exp(1i*4*pi/4) exp(1i*4*pi/4) exp(1i*4*pi/4) exp(1i*4*pi/4) exp(1i*4*pi/4) exp(1i*4*pi/4) exp(1i*4*pi/4) exp(1i*4*pi/4) exp(1i*4*pi/4) exp(1i*4*pi/4) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2)];
                for n=1:L
           ap=(B(n)/L)*sin(g*pi/L)/(g*pi/L)*exp(-1i*pi*g*((2*n-1))/L);
            B1(n)=ap;
              end
        C=sum(B1);
        apq(i,j)=C;  
            elseif phimn11(i,j)==225
                 B=[exp(1i*5*pi/4) exp(1i*5*pi/4) exp(1i*5*pi/4) exp(1i*5*pi/4) exp(1i*5*pi/4) exp(1i*5*pi/4) exp(1i*5*pi/4) exp(1i*5*pi/4) exp(1i*5*pi/4) exp(1i*5*pi/4) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2)];
                for n=1:L
           ap=(B(n)/L)*sin(g*pi/L)/(g*pi/L)*exp(-1i*pi*g*((2*n-1))/L);
            B1(n)=ap;
              end
        C=sum(B1);
        apq(i,j)=C;
            elseif phimn11(i,j)==270
                 B=[exp(1i*6*pi/4) exp(1i*6*pi/4) exp(1i*6*pi/4) exp(1i*6*pi/4) exp(1i*6*pi/4) exp(1i*6*pi/4) exp(1i*6*pi/4) exp(1i*6*pi/4) exp(1i*6*pi/4) exp(1i*6*pi/4) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2)];
                 for n=1:L
           ap=(B(n)/L)*sin(g*pi/L)/(g*pi/L)*exp(-1i*pi*g*((2*n-1))/L);
            B1(n)=ap;
              end
        C=sum(B1);
        apq(i,j)=C; 
            else
                 B=[exp(1i*7*pi/4) exp(1i*7*pi/4) exp(1i*7*pi/4) exp(1i*7*pi/4) exp(1i*7*pi/4) exp(1i*7*pi/4) exp(1i*7*pi/4) exp(1i*7*pi/4) exp(1i*7*pi/4) exp(1i*7*pi/4) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2)];
              for n=1:L
           ap=(B(n)/L)*sin(g*pi/L)/(g*pi/L)*exp(-1i*pi*g*((2*n-1))/L);
            B1(n)=ap;
              end
        C=sum(B1);
        apq(i,j)=C;
            end
         end
        %%  5%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
         if AC11(i,j)==6/8
            if phimn11(i,j)==0
                B=[exp(1i*0*pi/4) exp(1i*0*pi/4) exp(1i*0*pi/4) exp(1i*0*pi/4) exp(1i*0*pi/4) exp(1i*0*pi/4) exp(1i*0*pi/4) exp(1i*0*pi/4) exp(1i*0*pi/4) exp(1i*0*pi/4) exp(1i*0*pi/4) exp(1i*0*pi/4) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2)];
              for n=1:L
           ap=(B(n)/L)*sin(g*pi/L)/(g*pi/L)*exp(-1i*pi*g*((2*n-1))/L);
            B1(n)=ap;
              end
        C=sum(B1);
        apq(i,j)=C;
            elseif phimn11(i,j)==45
                 B=[exp(1i*1*pi/4) exp(1i*1*pi/4) exp(1i*1*pi/4) exp(1i*1*pi/4) exp(1i*1*pi/4) exp(1i*1*pi/4) exp(1i*1*pi/4) exp(1i*1*pi/4) exp(1i*1*pi/4) exp(1i*1*pi/4) exp(1i*1*pi/4) exp(1i*1*pi/4) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2)];
               for n=1:L
           ap=(B(n)/L)*sin(g*pi/L)/(g*pi/L)*exp(-1i*pi*g*((2*n-1))/L);
            B1(n)=ap;
              end
        C=sum(B1);
        apq(i,j)=C;
            elseif phimn11(i,j)==90
                 B=[exp(1i*2*pi/4) exp(1i*2*pi/4) exp(1i*2*pi/4) exp(1i*2*pi/4) exp(1i*2*pi/4) exp(1i*2*pi/4) exp(1i*2*pi/4) exp(1i*2*pi/4) exp(1i*2*pi/4) exp(1i*2*pi/4) exp(1i*2*pi/4) exp(1i*2*pi/4) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2)];
                for n=1:L
          ap=(B(n)/L)*sin(g*pi/L)/(g*pi/L)*exp(-1i*pi*g*((2*n-1))/L);
            B1(n)=ap;
              end
        C=sum(B1);
        apq(i,j)=C;  
            elseif phimn11(i,j)==135
                 B=[exp(1i*3*pi/4) exp(1i*3*pi/4) exp(1i*3*pi/4) exp(1i*3*pi/4) exp(1i*3*pi/4) exp(1i*3*pi/4) exp(1i*3*pi/4) exp(1i*3*pi/4) exp(1i*3*pi/4) exp(1i*3*pi/4) exp(1i*3*pi/4) exp(1i*3*pi/4) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2)];
                for n=1:L
           ap=(B(n)/L)*sin(g*pi/L)/(g*pi/L)*exp(-1i*pi*g*((2*n-1))/L);
            B1(n)=ap;
              end
        C=sum(B1);
        apq(i,j)=C;  
            elseif phimn11(i,j)==180
                 B=[exp(1i*4*pi/4) exp(1i*4*pi/4) exp(1i*4*pi/4) exp(1i*4*pi/4) exp(1i*4*pi/4) exp(1i*4*pi/4) exp(1i*4*pi/4) exp(1i*4*pi/4) exp(1i*4*pi/4) exp(1i*4*pi/4) exp(1i*4*pi/4) exp(1i*4*pi/4) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2)];
                for n=1:L
           ap=(B(n)/L)*sin(g*pi/L)/(g*pi/L)*exp(-1i*pi*g*((2*n-1))/L);
            B1(n)=ap;
              end
        C=sum(B1);
        apq(i,j)=C;   
            elseif phimn11(i,j)==225
                 B=[exp(1i*5*pi/4) exp(1i*5*pi/4) exp(1i*5*pi/4) exp(1i*5*pi/4) exp(1i*5*pi/4) exp(1i*5*pi/4) exp(1i*5*pi/4) exp(1i*5*pi/4) exp(1i*5*pi/4) exp(1i*5*pi/4) exp(1i*5*pi/4) exp(1i*5*pi/4) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2)];
                for n=1:L
           ap=(B(n)/L)*sin(g*pi/L)/(g*pi/L)*exp(-1i*pi*g*((2*n-1))/L);
            B1(n)=ap;
              end
        C=sum(B1);
        apq(i,j)=C; 
            elseif phimn11(i,j)==270
                 B=[exp(1i*6*pi/4) exp(1i*6*pi/4) exp(1i*6*pi/4) exp(1i*6*pi/4) exp(1i*6*pi/4) exp(1i*6*pi/4) exp(1i*6*pi/4) exp(1i*6*pi/4) exp(1i*6*pi/4) exp(1i*6*pi/4) exp(1i*6*pi/4) exp(1i*6*pi/4) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2)];
                  for n=1:L
            ap=(B(n)/L)*sin(g*pi/L)/(g*pi/L)*exp(-1i*pi*g*((2*n-1))/L);
            B1(n)=ap;
              end
        C=sum(B1);
        apq(i,j)=C;
            else
                 B=[exp(1i*7*pi/4) exp(1i*7*pi/4) exp(1i*7*pi/4) exp(1i*7*pi/4) exp(1i*7*pi/4) exp(1i*7*pi/4) exp(1i*7*pi/4) exp(1i*7*pi/4) exp(1i*7*pi/4) exp(1i*7*pi/4) exp(1i*7*pi/4) exp(1i*7*pi/4) exp(1i*1*pi/2) exp(1i*3*pi/2) exp(1i*1*pi/2) exp(1i*3*pi/2)];
             for n=1:L
           ap=(B(n)/L)*sin(g*pi/L)/(g*pi/L)*exp(-1i*pi*g*((2*n-1))/L);
            B1(n)=ap;
              end
        C=sum(B1);
        apq(i,j)=C;
            end
         end
        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         if AC11(i,j)==7/8
            if phimn11(i,j)==0
                B=[exp(1i*0*pi/4) exp(1i*0*pi/4) exp(1i*0*pi/4) exp(1i*0*pi/4) exp(1i*0*pi/4) exp(1i*0*pi/4) exp(1i*0*pi/4) exp(1i*0*pi/4) exp(1i*0*pi/4) exp(1i*0*pi/4) exp(1i*0*pi/4) exp(1i*0*pi/4) exp(1i*0*pi/4) exp(1i*0*pi/4) exp(1i*1*pi/2) exp(1i*3*pi/2)];
              for n=1:L
           ap=(B(n)/L)*sin(g*pi/L)/(g*pi/L)*exp(-1i*pi*g*((2*n-1))/L);
            B1(n)=ap;
              end
        C=sum(B1);
        apq(i,j)=C;
            elseif phimn11(i,j)==45
                 B=[exp(1i*1*pi/4) exp(1i*1*pi/4) exp(1i*1*pi/4) exp(1i*1*pi/4) exp(1i*1*pi/4) exp(1i*1*pi/4) exp(1i*1*pi/4) exp(1i*1*pi/4) exp(1i*1*pi/4) exp(1i*1*pi/4) exp(1i*1*pi/4) exp(1i*1*pi/4) exp(1i*1*pi/4) exp(1i*1*pi/4) exp(1i*1*pi/2) exp(1i*3*pi/2)];
               for n=1:L
           ap=(B(n)/L)*sin(g*pi/L)/(g*pi/L)*exp(-1i*pi*g*((2*n-1))/L);
            B1(n)=ap;
              end
        C=sum(B1);
        apq(i,j)=C;
            elseif phimn11(i,j)==90
                 B=[exp(1i*2*pi/4) exp(1i*2*pi/4) exp(1i*2*pi/4) exp(1i*2*pi/4) exp(1i*2*pi/4) exp(1i*2*pi/4) exp(1i*2*pi/4) exp(1i*2*pi/4) exp(1i*2*pi/4) exp(1i*2*pi/4) exp(1i*2*pi/4) exp(1i*2*pi/4) exp(1i*2*pi/4) exp(1i*2*pi/4) exp(1i*1*pi/2) exp(1i*3*pi/2)];
                for n=1:L
           ap=(B(n)/L)*sin(g*pi/L)/(g*pi/L)*exp(-1i*pi*g*((2*n-1))/L);
            B1(n)=ap;
              end
        C=sum(B1);
        apq(i,j)=C; 
            elseif phimn11(i,j)==135
                 B=[exp(1i*3*pi/4) exp(1i*3*pi/4) exp(1i*3*pi/4) exp(1i*3*pi/4) exp(1i*3*pi/4) exp(1i*3*pi/4) exp(1i*3*pi/4) exp(1i*3*pi/4) exp(1i*3*pi/4) exp(1i*3*pi/4) exp(1i*3*pi/4) exp(1i*3*pi/4) exp(1i*3*pi/4) exp(1i*3*pi/4) exp(1i*1*pi/2) exp(1i*3*pi/2)];
                for n=1:L
           ap=(B(n)/L)*sin(g*pi/L)/(g*pi/L)*exp(-1i*pi*g*((2*n-1))/L);
            B1(n)=ap;
              end
        C=sum(B1);
        apq(i,j)=C;  
            elseif phimn11(i,j)==180
                 B=[exp(1i*4*pi/4) exp(1i*4*pi/4) exp(1i*4*pi/4) exp(1i*4*pi/4) exp(1i*4*pi/4) exp(1i*4*pi/4) exp(1i*4*pi/4) exp(1i*4*pi/4) exp(1i*4*pi/4) exp(1i*4*pi/4) exp(1i*4*pi/4) exp(1i*0*pi/4) exp(1i*4*pi/4) exp(1i*4*pi/4) exp(1i*1*pi/2) exp(1i*3*pi/2)];
                for n=1:L
           ap=(B(n)/L)*sin(g*pi/L)/(g*pi/L)*exp(-1i*pi*g*((2*n-1))/L);
            B1(n)=ap;
              end
        C=sum(B1);
        apq(i,j)=C;  
            elseif phimn11(i,j)==225
                 B=[exp(1i*5*pi/4) exp(1i*5*pi/4) exp(1i*5*pi/4) exp(1i*5*pi/4) exp(1i*5*pi/4) exp(1i*5*pi/4) exp(1i*5*pi/4) exp(1i*5*pi/4) exp(1i*5*pi/4) exp(1i*5*pi/4) exp(1i*5*pi/4) exp(1i*5*pi/4) exp(1i*5*pi/4) exp(1i*5*pi/4) exp(1i*1*pi/2) exp(1i*3*pi/2)];
                for n=1:L
           ap=(B(n)/L)*sin(g*pi/L)/(g*pi/L)*exp(-1i*pi*g*((2*n-1))/L);
            B1(n)=ap;
              end
        C=sum(B1);
        apq(i,j)=C;   
            elseif phimn11(i,j)==270
                 B=[exp(1i*6*pi/4) exp(1i*6*pi/4) exp(1i*6*pi/4) exp(1i*6*pi/4) exp(1i*6*pi/4) exp(1i*6*pi/4) exp(1i*6*pi/4) exp(1i*6*pi/4) exp(1i*6*pi/4) exp(1i*6*pi/4) exp(1i*6*pi/4) exp(1i*6*pi/4) exp(1i*6*pi/4) exp(1i*6*pi/4) exp(1i*1*pi/2) exp(1i*3*pi/2)];
                  for n=1:L
           ap=(B(n)/L)*sin(g*pi/L)/(g*pi/L)*exp(-1i*pi*g*((2*n-1))/L);
            B1(n)=ap;
              end
        C=sum(B1);
        apq(i,j)=C;
            else
                 B=[exp(1i*7*pi/4) exp(1i*7*pi/4) exp(1i*7*pi/4) exp(1i*7*pi/4) exp(1i*7*pi/4) exp(1i*7*pi/4) exp(1i*7*pi/4) exp(1i*7*pi/4) exp(1i*7*pi/4) exp(1i*7*pi/4) exp(1i*7*pi/4) exp(1i*7*pi/4) exp(1i*7*pi/4) exp(1i*7*pi/4) exp(1i*1*pi/2) exp(1i*3*pi/2)];
               for n=1:L
           ap=(B(n)/L)*sin(g*pi/L)/(g*pi/L)*exp(-1i*pi*g*((2*n-1))/L);
            B1(n)=ap;
              end
        C=sum(B1);
        apq(i,j)=C;
            end
         end
        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         if AC11(i,j)==8/8
            if phimn11(i,j)==0
                B=[exp(1i*0*pi/4) exp(1i*0*pi/4) exp(1i*0*pi/4) exp(1i*0*pi/4) exp(1i*0*pi/4) exp(1i*0*pi/4) exp(1i*0*pi/4) exp(1i*0*pi/4) exp(1i*0*pi/4) exp(1i*0*pi/4) exp(1i*0*pi/4) exp(1i*0*pi/4) exp(1i*0*pi/4) exp(1i*0*pi/4) exp(1i*0*pi/4) exp(1i*0*pi/4)];
              for n=1:L
           ap=(B(n)/L)*sin(g*pi/L)/(g*pi/L)*exp(-1i*pi*g*((2*n-1))/L);
            B1(n)=ap;
              end
        C=sum(B1);
        apq(i,j)=C;
            elseif phimn11(i,j)==45
                 B=[exp(1i*1*pi/4) exp(1i*1*pi/4) exp(1i*1*pi/4) exp(1i*1*pi/4) exp(1i*1*pi/4) exp(1i*1*pi/4) exp(1i*1*pi/4) exp(1i*1*pi/4) exp(1i*1*pi/4) exp(1i*1*pi/4) exp(1i*1*pi/4) exp(1i*1*pi/4) exp(1i*1*pi/4) exp(1i*1*pi/4) exp(1i*1*pi/4) exp(1i*1*pi/4)];
               for n=1:L
           ap=(B(n)/L)*sin(g*pi/L)/(g*pi/L)*exp(-1i*pi*g*((2*n-1))/L);
            B1(n)=ap;
              end
        C=sum(B1);
        apq(i,j)=C;
            elseif phimn11(i,j)==90
                 B=[exp(1i*2*pi/4) exp(1i*2*pi/4) exp(1i*2*pi/4) exp(1i*2*pi/4) exp(1i*2*pi/4) exp(1i*2*pi/4) exp(1i*2*pi/4) exp(1i*2*pi/4) exp(1i*2*pi/4) exp(1i*2*pi/4) exp(1i*2*pi/4) exp(1i*2*pi/4) exp(1i*2*pi/4) exp(1i*2*pi/4) exp(1i*2*pi/4) exp(1i*2*pi/4)];
                for n=1:L
           ap=(B(n)/L)*sin(g*pi/L)/(g*pi/L)*exp(-1i*pi*g*((2*n-1))/L);
            B1(n)=ap;
              end
        C=sum(B1);
        apq(i,j)=C; 
            elseif phimn11(i,j)==135
                 B=[exp(1i*3*pi/4) exp(1i*3*pi/4) exp(1i*3*pi/4) exp(1i*3*pi/4) exp(1i*3*pi/4) exp(1i*3*pi/4) exp(1i*3*pi/4) exp(1i*3*pi/4) exp(1i*3*pi/4) exp(1i*3*pi/4) exp(1i*3*pi/4) exp(1i*3*pi/4) exp(1i*3*pi/4) exp(1i*3*pi/4) exp(1i*3*pi/4) exp(1i*3*pi/4)];
                for n=1:L
           ap=(B(n)/L)*sin(g*pi/L)/(g*pi/L)*exp(-1i*pi*g*((2*n-1))/L);
            B1(n)=ap;
              end
        C=sum(B1);
        apq(i,j)=C;  
            elseif phimn11(i,j)==180
                 B=[exp(1i*4*pi/4) exp(1i*4*pi/4) exp(1i*4*pi/4) exp(1i*4*pi/4) exp(1i*4*pi/4) exp(1i*4*pi/4) exp(1i*4*pi/4) exp(1i*4*pi/4) exp(1i*4*pi/4) exp(1i*4*pi/4) exp(1i*4*pi/4) exp(1i*4*pi/4) exp(1i*4*pi/4) exp(1i*4*pi/4) exp(1i*4*pi/4) exp(1i*4*pi/4)];
                for n=1:L
           ap=(B(n)/L)*sin(g*pi/L)/(g*pi/L)*exp(-1i*pi*g*((2*n-1))/L);
            B1(n)=ap;
              end
        C=sum(B1);
        apq(i,j)=C;  
            elseif phimn11(i,j)==225
                 B=[exp(1i*5*pi/4) exp(1i*5*pi/4) exp(1i*5*pi/4) exp(1i*5*pi/4) exp(1i*5*pi/4) exp(1i*5*pi/4) exp(1i*5*pi/4) exp(1i*5*pi/4) exp(1i*5*pi/4) exp(1i*5*pi/4) exp(1i*5*pi/4) exp(1i*5*pi/4) exp(1i*5*pi/4) exp(1i*5*pi/4) exp(1i*5*pi/4) exp(1i*5*pi/4)];
                for n=1:L
           ap=(B(n)/L)*sin(g*pi/L)/(g*pi/L)*exp(-1i*pi*g*((2*n-1))/L);
            B1(n)=ap;
              end
        C=sum(B1);
        apq(i,j)=C;   
            elseif phimn11(i,j)==270
                 B=[exp(1i*6*pi/4) exp(1i*6*pi/4) exp(1i*6*pi/4) exp(1i*6*pi/4) exp(1i*6*pi/4) exp(1i*6*pi/4) exp(1i*6*pi/4) exp(1i*6*pi/4) exp(1i*6*pi/4) exp(1i*6*pi/4) exp(1i*6*pi/4) exp(1i*6*pi/4) exp(1i*6*pi/4) exp(1i*6*pi/4) exp(1i*6*pi/4) exp(1i*6*pi/4)];
                  for n=1:L
           ap=(B(n)/L)*sin(g*pi/L)/(g*pi/L)*exp(-1i*pi*g*((2*n-1))/L);
            B1(n)=ap;
              end
        C=sum(B1);
        apq(i,j)=C;
            else
                 B=[exp(1i*7*pi/4) exp(1i*7*pi/4) exp(1i*7*pi/4) exp(1i*7*pi/4) exp(1i*7*pi/4) exp(1i*7*pi/4) exp(1i*7*pi/4) exp(1i*7*pi/4) exp(1i*7*pi/4) exp(1i*7*pi/4) exp(1i*7*pi/4) exp(1i*7*pi/4) exp(1i*7*pi/4) exp(1i*7*pi/4) exp(1i*7*pi/4) exp(1i*7*pi/4)];
               for n=1:L
           ap=(B(n)/L)*sin(g*pi/L)/(g*pi/L)*exp(-1i*pi*g*((2*n-1))/L);
            B1(n)=ap;
              end
        C=sum(B1);
        apq(i,j)=C;
            end
         end
    end
end
AC=abs(apq)*MM;
phimn=angle(apq);
%%  %% %%%%%%%%%%%%%%%%%%%?Array factor
k=(2*pi*(f+g*f0))/c;
for u=1:length(phi)
    for v=1:length(theta)
        E=0;
        for m=1:M
            for n=1:N      
                E=E+AC(m,n)*exp(+1i*k*dx*(n-1)*sin(theta(v))*cos(phi(u))+1i*k*dy*(m-1)*sin(theta(v))*sin(phi(u))+1i*phimn(m,n));
            end
        end
%         AF(u,v)=E*cos(theta(v));
        AF(u,v)=E;
    end
end
% %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Directivity
Prad=0;
for u=1:length(phi)
    for v=1:length(theta)
        Prad=Prad+(abs(AF(u,v)))^2*sin(theta(v));
    end
end
Pradiation=dphi*dtheta*Prad;
PP(g)=Pradiation;
end
PM=sum(PP);
Ptotal=Pradiation0+(2*PM);
for u=1:length(phi)
    for v=1:length(theta)
        dir(u,v)=4*pi*abs(AF1(u,v))^2/Ptotal;
    end
end
% %% %%%%%%%%%%%%%%%%%%%%%%%%plot
figure(111)
imagesc(theta*180/pi,phi*180/pi,abs(dir))
colormap(jet);
phimn=phimn.*180/pi;
figure 
imagesc(phimn1)
colorbar
title('Surface Phase Distribution')
figure 
imagesc(AC1)
colorbar
title('Surface amplitute Distribution')
figure
plot(theta*180/pi,abs(dir(180*r/360,:)),'linewidth',2)
hold on
plot(theta*180/pi,abs(dir(270*r/360,:)),'linewidth',2)
title('Directivity')
figure
plot(theta*180/pi,10*log10(abs(dir(180*r/360,:))),'linewidth',2)
hold on
plot(theta*180/pi,10*log10(abs(dir(270*r/360,:))),'linewidth',2)
title('Directivity')
set(gca,'YLim',[-10 30])
figure
plot(theta*180/pi,abs(AF1(180*r/360,:)),'linewidth',2)
hold on
plot(theta*180/pi,abs(AF1(270*r/360,:)),'linewidth',2)
title('Directivity')
% figure
% plot(theta*180/pi,abs(dir(360*r/360,:)),'linewidth',2)
% title('Directivity')
% set(gca,'YLim',[-30 0])