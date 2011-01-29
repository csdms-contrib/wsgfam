%WSGFAM v.1.3.1 July 29, 2009
%Code by Malcolm E. Scully, Yanxia Ma, and Carl T. Friedrichs,
%applied to the Waiapu Shelf as documented by Yanxia Ma, 
%Carl T. Friedrichs, Courtney K. Harris, & L. Donelson Wright (2009)
%Deposition on the Waiapu, New Zealand, continental shelf by 
%seasonal wave- and current-supported sediment gravity flows interacting
%with spatially varying bathymetry. To be submitted to Marine Geology.
%
%WSGFAM v.1.2 Code by Malcolm E. Scully and Carl T. Friedrichs,
%applied to Po Shelf as documented by Friedrichs & Scully (2007)
%Continental Shelf Research, 27: 322-337.
%
%WSGFAM v.1.1 Code by Malcolm E. Scully and Carl T. Friedrichs,
%applied to Eel Shelf as documented by Scully et al. (2003)
%Journal of Geophysical Research, 108(C4), doi:10.1029/2002JC001467.
%
%WSGFAM v.1.3.1 is for the pb210 constrained, "High Energy Period" run
clear
load wavepb210; % n*3 matrix:(time, r.m.s. wave height, period)
%
load sedpb210; %n*2 matris: (time, Q_s)
%sed(:,2)=sed(:,2)*2; %sensitive test, double sediment input 

load depthnew; %model domain bathymetry% 
load rotlon; %model domain longitude
load rotlat;  %model domain latitude
load outmodel; %model domain boundary 
[l,m,n]=size(outmodel);

load sed_linear+3  %along-shelf sediment distribution 
load critical35; %critical wave orbital velocity

load ucpb210; %current

%%

dx=193.9;%across-shelf length of grid cell , 
dy=178.1; %along-shelf length of grid cell

area=dx*dy;
g=9.8;
psed=2650;
s=1.65;
Ri=.25;
Cd=.003;
%limit steep slope to prevent beta >= 1
totalslope=sqrt(X_slope.^2+Y_slope.^2);
steep=find(totalslope >= 0.011);
totalslope(steep)=0.011;
beta=totalslope*Ri/Cd;

Load=sed(:,2)*0.23; %reduce the sediment input

totalsedDw=zeros(1,m,n);
maxconc=zeros(1,m,n);
Deposit=zeros(1,m,n);
eroded=zeros(1,m,n);
Dwconc=zeros(1,m,n);
volconc=zeros(1,m,n);
B=zeros(1,m,n);
uplumeX=zeros(1,m,n);
vplumeY=zeros(1,m,n);
offshelflux=zeros(1,m,n);
outfluxX=zeros(1,m,n);
outfluxY=zeros(1,m,n);
fluxinX=zeros(1,m,n+1);
Northfluxin=zeros(1,m+2,n);
Southfluxin=zeros(1,m+2,n);
Dwconcno0=zeros(1,m,n);
uplumeno0=zeros(1,m,n);
vplumeno0=zeros(1,m,n);
totaloutflux=zeros(1,m,n+1);
Xpercent=zeros(1,m,n);
Ypercent=zeros(1,m,n);
outmodelflux=zeros(1,m,n);
a=zeros(m,n);
u_c=zeros(1,m,n);


A1=0;B1=0;C1=0;D1=0;X1=0;
x=1:n;
y=1:m;
y2=2:m+1;

ifX= (outmodel(1,y,x)==1);
North= Y_slope(1,y,x)>=0;
South= Y_slope(1,y,x)<0;

totsed_input=0;
t1=0;


for t=1:10461  

   t1=t1+1;
   
   if (t/100)==round(t/100)
	   t10=t/100
   end   
   %calculate max wave orbital velocity at all points in grid
    a(:)=wave(t,2)/2;
    gamma=.4*depth;
    bad=find((2*a)>gamma);
    a(bad)=gamma(bad)/2;
%    sigma=(2*pi)./wave(t,3);
     sigma2=((2*pi)/wave(t,3))^2;   
     kh=1;
    %if depth=0 at some grid points,there would be Warning: Divide by zero.
    ind=find(depth==0);
    depth(ind)=0.5;
    for j=1:10
        kh=(depth*sigma2)./(9.8*tanh(kh));
    end

     u_wave(1,:,:)=(sqrt(sigma2)*a)./sinh(kh);
 
 
 % u_wave=u_wave*(1/sqrt(2)); %if a is significant wave height
   u_c(1,:,:)=uc(t);
   
    newdeposit(1,y,x)=((900*Load(t)*sed_norm(1,y,x)));
    totsed_input=totsed_input+sum(sum(newdeposit));
    Deposit(1,y,x)=Deposit(1,y,x)+newdeposit(1,y,x);

   totalsedDw(1,y,x)=totalsedDw(1,y,x)+fluxinX(1,y,x)+Northfluxin(1,y2,x)+Southfluxin(1,y2,x);
   outmodelflux(1,y,x)=outmodelflux(1,y,x)+(totalsedDw(1,y,x).*(outmodel));
   totalsedDw(1,y,x)=totalsedDw(1,y,x)-(totalsedDw(1,y,x).*(outmodel));
 
    %
	umaxconc=sqrt((u_wave.^2+u_c.^2)./(1-(beta.^2)));
    %
	maxconc(1,y,x)=(((umaxconc(1,y,x)).^2)*Ri*psed)/(g*s);
    neg=maxconc<0;
    maxconc(neg)=0;

	ifA= ifX==0 & totalsedDw(1,y,x)>=(maxconc(1,y,x)*area);
	ifB= ifX==0 & totalsedDw(1,y,x)<(maxconc(1,y,x)*area)&u_wave(1,y,x)>critical&((Deposit(:,y,x)))>((maxconc(1,y,x)*area)-totalsedDw(1,y,x));
	ifC= ifX==0 & totalsedDw(1,y,x)<(maxconc(1,y,x)*area)&u_wave(1,y,x)>critical&((Deposit(:,y,x)))<=((maxconc(1,y,x)*area)-totalsedDw(1,y,x));
	ifD= ifX==0 & totalsedDw(1,y,x)<(maxconc(1,y,x)*area)&u_wave(1,y,x)<=critical;
      

   Deposit(1,y,x)=Deposit(1,y,x)+((totalsedDw(1,y,x)-(maxconc(1,y,x)*area)).*ifA);

   eroded=(((maxconc(1,y,x)*area)-totalsedDw(1,y,x)).*ifB)+(Deposit.*ifC);
   Deposit=Deposit-(eroded);
   totalsedDw(1,y,x)=((maxconc(1,y,x)*area).*ifA)+(((totalsedDw(1,y,x)+eroded(1,y,x)).*(ifB+ifC)))+(totalsedDw(1,y,x).*ifD);
 
   
   Dwconc(1,y,x)=((totalsedDw(1,y,x))/area);
   volconc(1,y,x)=((Dwconc(1,y,x)/2650));
   B(1,y,x)=((g*s*volconc(1,y,x)));

   umax=sqrt(u_wave.^2+u_c.^2);
   for i=1:5 
      uplumeX(1,y,x)=(((X_slope(1,y,x).*B(1,y,x))./(umax(1,y,x).*Cd)).*(ifA+ifB+ifC+ifD));
      vplumeY(1,y,x)=(((abs(Y_slope(1,y,x)).*B(1,y,x))./(umax(1,y,x).*Cd)).*(ifA+ifB+ifC+ifD));
      umax=sqrt(uplumeX.^2+vplumeY.^2+u_wave(1,y,x).^2+u_c.^2);
  end

%  
   ifplume(1,y,x)=  uplumeX(1,y,x)>0.0001;
   ifnoplume(1,y,x)= ifplume(1,y,x)==0;
   uplumeno0(1,y,x)=uplumeX+ifnoplume;
   vplumeno0(1,y,x)=vplumeY+ifnoplume;
   Dwconcno0(1,y,x)=Dwconc+ifnoplume;
   Xpercent(1,y,x)=(uplumeno0(1,y,x).*Dwconcno0(1,y,x)*900*dy)./(((uplumeno0(1,y,x).*Dwconcno0(1,y,x)*900*dy))+((vplumeno0(1,y,x).*Dwconcno0(1,y,x)*900*dx)));
   Ypercent(1,y,x)=(vplumeno0(1,y,x).*Dwconcno0(1,y,x)*900*dx)./(((uplumeno0(1,y,x).*Dwconcno0(1,y,x)*900*dy))+((vplumeno0(1,y,x).*Dwconcno0(1,y,x)*900*dx)));
  
  
   ifA1= ifA==1 & (((uplumeX(1,y,x).*Dwconc(1,y,x)*900*dy))+((vplumeY(1,y,x).*Dwconc(1,y,x)*900*dx)))<(Dwconc(1,y,x)*area);
   ifA2= ifA==1 & (((uplumeX(1,y,x).*Dwconc(1,y,x)*900*dy))+((vplumeY(1,y,x).*Dwconc(1,y,x)*900*dx)))>=(Dwconc(1,y,x)*area);
   ifB1= ifB==1 & (((uplumeX(1,y,x).*Dwconc(1,y,x)*900*dy))+((vplumeY(1,y,x).*Dwconc(1,y,x)*900*dx)))<(Dwconc(1,y,x)*area);
   ifB2= ifB==1 & (((uplumeX(1,y,x).*Dwconc(1,y,x)*900*dy))+((vplumeY(1,y,x).*Dwconc(1,y,x)*900*dx)))>=(Dwconc(1,y,x)*area);
   ifC1= ifC==1 & (((uplumeX(1,y,x).*Dwconc(1,y,x)*900*dy))+((vplumeY(1,y,x).*Dwconc(1,y,x)*900*dx)))<(Dwconc(1,y,x)*area);
   ifC2= ifC==1 & (((uplumeX(1,y,x).*Dwconc(1,y,x)*900*dy))+((vplumeY(1,y,x).*Dwconc(1,y,x)*900*dx)))>=(Dwconc(1,y,x)*area);
   ifD1= ifD==1 & (((uplumeX(1,y,x).*Dwconc(1,y,x)*900*dy))+((vplumeY(1,y,x).*Dwconc(1,y,x)*900*dx)))<(Dwconc(1,y,x)*area);
   ifD2= ifD==1 &(((uplumeX(1,y,x).*Dwconc(1,y,x)*900*dy))+((vplumeY(1,y,x).*Dwconc(1,y,x)*900*dx)))>=(Dwconc(1,y,x)*area);
   
   outfluxX(1,y,x)=((totalsedDw(1,y,x)).*(Xpercent(1,y,x)).*(ifA2+ifB2+ifC2+ifD2))+((ifA1+ifB1+ifC1+ifD1).*(uplumeX(1,y,x).*Dwconc(1,y,x)*900*dy));
   outfluxY(1,y,x)=((totalsedDw(1,y,x)).*(Ypercent(1,y,x)).*(ifA2+ifB2+ifC2+ifD2))+((ifA1+ifB1+ifC1+ifD1).*(vplumeY(1,y,x).*Dwconc(1,y,x)*900*dx));

   totalsedDw(1,y,x)=(totalsedDw(1,y,x)-(outfluxX(1,y,x)+outfluxY(1,y,x)));


	X1=X1+sum(sum(ifX));
	A1=A1+sum(sum(ifA));
   B1=B1+sum(sum(ifB));
   C1=C1+sum(sum(ifC));
   D1=D1+sum(sum(ifD));
   
   fluxinX(1,y,x+1)=(outfluxX(1,y,x));
   Northfluxin(1,y2+1,x)=outfluxY(1,y,x).*(North);
   Southfluxin(1,y2-1,x)=outfluxY(1,y,x).*(South);
  
   dep60(t)=Deposit(1,103,57);
   dep40(t)=Deposit(1,101,47);

end

(sum(sum(sum(Deposit)))+sum(sum(outmodelflux))+sum(sum(totalsedDw))+sum(sum(fluxinX))+sum(sum(Northfluxin))+sum(sum(Southfluxin)))/totsed_input

%
Final_Deposit_uni=squeeze(Deposit);
save Final_Deposit_uni Final_Deposit_uni


