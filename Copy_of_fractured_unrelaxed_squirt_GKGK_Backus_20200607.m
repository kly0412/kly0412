clc;
clear;

fractured=load('fractured.txt');
[n1 n2]=size(fractured);
pi=3.1415926;
w0=500e3*2*pi;  
nSw=101;
nInc=3;

%%%%%%%%%% correct anisotropic coefficients from fractured data %%%%%%%%%%%%%%%%%%%%%%
mult0 = 1.0347;% average(fractured(:,8)/fractured(:,6))
mult45 = 1.0088;% average(fractured(:,8)/fractured(:,7))
mult0V = 1.0172;% average(fractured(:,5)/fractured(:,3))
mult45V = 1.0044;% average(fractured(:,5)/fractured(:,4))

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% anisotropy adjustment %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% shift all valuses in certain lines including Sw=100% %%%%%%%%%%%%%
% fractured_correct=fractured;
% fractured_correct(:,3)=fractured(:,3)*mult0V;
% fractured_correct(:,4)=fractured(:,4)*mult45V;
% fractured_correct(:,6)=fractured(:,6)*mult0;
% fractured_correct(:,7)=fractured(:,7)*mult45;
% 
%%%%%%%%% shift all valuses in certain lines except Sw=100% %%%%%%%%%%%%%
fractured_correct=fractured;
fractured_correct(1:n1-1,3)=fractured(1:n1-1,3)*mult0V;
fractured_correct(1:n1-1,4)=fractured(1:n1-1,4)*mult45V;
fractured_correct(1:n1-1,6)=fractured(1:n1-1,6)*mult0;
fractured_correct(1:n1-1,7)=fractured(1:n1-1,7)*mult45;

%%%%%%%%%%% plot data %%%%%%%%%%%%%%%%%%%%%%
figure() %%%%% original data  %%%%%%%%%%%%
subplot(1,2,1)
set(gcf,'outerposition',get(0,'screensize'));
h1=plot(fractured_correct(:,1)*0.01,fractured(:,8),'ws','markerfacecolor','r','markersize',20);
hold on;  
h2=plot(fractured_correct(:,1)*0.01,fractured(:,7),'wo','markerfacecolor','g','markersize',20);
hold on; 
h3=plot(fractured_correct(:,1)*0.01,fractured(:,6),'wd','markerfacecolor','b','markersize',20);
hold on;
set(gca,'FontSize',30);
set(gca,'xTick',[0 0.25 0.5 0.75 1]);
xlabel('S_w','fontsize',40);
ylabel('M_p(Pa)','fontsize',40);    
set(gca,'ylim',[14e9 25e9]);
lgd1=legend([h1,h2,h3],'Data-90бу','Data-45бу','Data-0бу','orientation','vertical','location','northwest');
set(lgd1,'FontName','Time New Roman','FontSize',30);

subplot(1,2,2)
set(gcf,'outerposition',get(0,'screensize'));
h1=plot(fractured_correct(:,1)*0.01,fractured(:,5),'ws','markerfacecolor','r','markersize',20);
hold on;  
h2=plot(fractured_correct(:,1)*0.01,fractured(:,4),'wo','markerfacecolor','g','markersize',20);
hold on; 
h3=plot(fractured_correct(:,1)*0.01,fractured(:,3),'wd','markerfacecolor','b','markersize',20);
hold on;
set(gca,'FontSize',30);
set(gca,'yTick',[2700 2900 3100 3300 3500]);
set(gca,'xTick',[0 0.25 0.5 0.75 1]);
xlabel('S_w','fontsize',40);
ylabel('V_p(m/s)','fontsize',40);    
set(gca,'ylim',[2700 3500]);

figure()%%%%% data after adjustment %%%%%%%%%%%%
subplot(1,2,1)
set(gcf,'outerposition',get(0,'screensize'));
h1=plot(fractured_correct(:,1)*0.01,fractured_correct(:,8),'ws','markerfacecolor','r','markersize',20);
hold on;  
h2=plot(fractured_correct(:,1)*0.01,fractured_correct(:,7),'wo','markerfacecolor','g','markersize',20);
hold on; 
h3=plot(fractured_correct(:,1)*0.01,fractured_correct(:,6),'wd','markerfacecolor','b','markersize',20);
hold on;
set(gca,'FontSize',30);
set(gca,'xTick',[0 0.25 0.5 0.75 1]);
xlabel('S_w','fontsize',40);
ylabel('M_p(Pa)','fontsize',40);    
set(gca,'ylim',[14e9 25e9]);
lgd1=legend([h1,h2,h3],'Data-90бу','Data-45бу','Data-0бу','orientation','vertical','location','northwest');
set(lgd1,'FontName','Time New Roman','FontSize',30);

subplot(1,2,2)
set(gcf,'outerposition',get(0,'screensize'));
h1=plot(fractured_correct(:,1)*0.01,fractured_correct(:,5),'ws','markerfacecolor','r','markersize',20);
hold on;  
h2=plot(fractured_correct(:,1)*0.01,fractured_correct(:,4),'wo','markerfacecolor','g','markersize',20);
hold on; 
h3=plot(fractured_correct(:,1)*0.01,fractured_correct(:,3),'wd','markerfacecolor','b','markersize',20);
hold on;
set(gca,'FontSize',30);
set(gca,'yTick',[2700 2900 3100 3300 3500]);
set(gca,'xTick',[0 0.25 0.5 0.75 1]);
xlabel('S_w','fontsize',40);
ylabel('V_p(m/s)','fontsize',40);    
set(gca,'ylim',[2700 3500]);

%% %%%%%%%%%%%%%  parameters input %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% BASIC PARAMETERS %%%%%%%%%%%%%%%%%%%%%%
Kg = 38e9;% table from Boris
Gg = 44e9;% table from Boris
Den_g = 2590;% table from Boris       

Kwater=2.25e9;% table from Boris
Dwater=1000;% table from Boris
Ywater=1e-3; % table 2 from Tillotson et al's paper (2014)

Kgas=0.0001e9;% table from Boris
Dgas=1.2;% table from Boris
Ygas=0.00182e-3; % table 2 from Tillotson et al's paper (2014)        

Por = 0.33; %%data from Kelvin et al's paper (2015)
Por_m = 0.3244; % data computed from Por_t(0.33)-Por_f(0.0056);
Por_p = 21e-15;% data from Kelvin et al's paper (2015)
phi_c=0.001;% microcrack porosity given by Boris

hc=0.0314;% data from Kelvin et al's paper (2015)
H0=6e-3; % fracture period from Tillotson et al, figure 2, page 1240
Por_f=0.0056; % fracture porosity from Tillotson et al, page 1241

Den_pgw=Den_g*(1-Por)+Por*Dwater;% fractured sample density for water saturation
Den_pgg=Den_g*(1-Por)+Por*Dgas;% fractured sample density for air saturation      

Por_t = Por_m + Por_f;% total porosity for fractured sample
Por_m_re = Por_m-phi_c;% relax porosity for fractured sample
Por_re = Por_t-phi_c;% relax porosity for fractured sample

Den_gw=Den_g*(1-Por_t)+Por_t*Dwater;% density for water saturated fractured sample
Den_gg=Den_g*(1-Por_t)+Por_t*Dgas;% density for air saturated fractured sample

V90s1=2226;% S1 wave velocity at 90 degrees for fractured sample with Sw=0
Gdry=Den_gg*V90s1*V90s1;% dry shear modulus for fractured sample 

%% %%%%%%%%%%%%%  computing codes for three Kdry methods %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for Kdry_method=1:3
    
    if Kdry_method==1 %%%%%%%%%%%%% computing Kdry from values at Sw=0
        
        V90pp=3242;% P-wave velocity at 90 degrees for fractured sample with Sw=0
        Mdry=Den_gg*V90pp*V90pp;% dry P wave modulus for fractured sample 
        Kdry=Mdry-4./3*Gdry;% dry bulk modulus for fractured sample    
        
    elseif Kdry_method==2  %%%%%%%%%%%%% computing Kdry from values at Sw=100
        
        Msat=fractured_correct(10,8);        
        Kdry =(3*Kg^2*Kwater - 3*Kg*Kwater*Msat + 4*Kg*Kwater*Gdry - 3*Kg^2*Msat*Por + 4*Kg^2*Por*Gdry + 3*Kg*Kwater*Msat*Por - 4*Kg*Kwater*Por*Gdry)...
                    /(4*Kwater*Gdry - 3*Kg^2*Por + 3*Kg*Kwater - 3*Kwater*Msat + 3*Kg*Kwater*Por);
                
    elseif Kdry_method==3 %%%%%%%%%%%%% computing Kdry from values at Sw=100 considering the unrelaxation at the beginning        
        %%%%%%%%%%%%%        (1)  Msat=Ksat+4./3*Usat;   %% Msat is known from data , fractured.txt     %%%%%%%%%%%
        %%%%%%%%%%%%%        (2) Ksat/(Kg-Ksat)=Kwater/(Por_t*(Kg-Kwater))+Kun/(Kg-Kun)  
        %%%%%%%%%%%%             Usat=Uun                                               %% Gassmann equation %%%%%%%%%%%
        %%%%%%%%%%%% (3) 1./Uun = 1./U0 - 4/15 * (1./K0 - 1/Kun)  %% equation 2 in Gurevich et al's paper(2009) %%%%%%%%%%
        %%%%%%%%%%%% Combining (1)(2)and (3), I have the expressions of Kun and Gun as follows %%%
        
        U0=Gdry;% 2226 is the velocity from SV data for fractured sample with Sw=0
        M0=fractured_correct(1,8);% M value from P-wave data for fractured sample with Sw=0 at Inc=90
        K0=M0-4./3*U0;              
        Msat=fractured_correct(10,8); %% Msat is known from data when Sw=100% at 90 incidence angle, fractured.txt     %%%%%%%%%%%
        Sw=fractured_correct(10,1)*0.01;
        Sg=1-Sw;
        Kf=1./(Sw/Kwater+Sg/Kgas);

        Kun1=((225*K0^2*Kg^4*Kf^2 - 450*K0^2*Kg^4*Kf*Msat*Por_t + 480*K0^2*Kg^4*Kf*Por_t*U0 + 225*K0^2*Kg^4*Msat^2*Por_t^2 - 480*K0^2*Kg^4*Msat*Por_t^2*U0 ...
        + 576*K0^2*Kg^4*Por_t^2*U0^2 + 450*K0^2*Kg^3*Kf^2*Msat*Por_t - 450*K0^2*Kg^3*Kf^2*Msat - 480*K0^2*Kg^3*Kf^2*Por_t*U0 + 720*K0^2*Kg^3*Kf^2*U0 ...
        - 450*K0^2*Kg^3*Kf*Msat^2*Por_t^2 + 450*K0^2*Kg^3*Kf*Msat^2*Por_t + 960*K0^2*Kg^3*Kf*Msat*Por_t^2*U0 - 1200*K0^2*Kg^3*Kf*Msat*Por_t*U0 ...
        - 1152*K0^2*Kg^3*Kf*Por_t^2*U0^2 + 768*K0^2*Kg^3*Kf*Por_t*U0^2 + 225*K0^2*Kg^2*Kf^2*Msat^2*Por_t^2 - 450*K0^2*Kg^2*Kf^2*Msat^2*Por_t ...
        + 225*K0^2*Kg^2*Kf^2*Msat^2 - 480*K0^2*Kg^2*Kf^2*Msat*Por_t^2*U0 + 1200*K0^2*Kg^2*Kf^2*Msat*Por_t*U0 - 840*K0^2*Kg^2*Kf^2*Msat*U0 ...
        + 576*K0^2*Kg^2*Kf^2*Por_t^2*U0^2 - 768*K0^2*Kg^2*Kf^2*Por_t*U0^2 + 576*K0^2*Kg^2*Kf^2*U0^2 + 120*K0^2*Kg^2*Kf*Msat^2*Por_t*U0 ...
        - 128*K0^2*Kg^2*Kf*Msat*Por_t*U0^2 - 120*K0^2*Kg*Kf^2*Msat^2*Por_t*U0 + 120*K0^2*Kg*Kf^2*Msat^2*U0 + 128*K0^2*Kg*Kf^2*Msat*Por_t*U0^2 ...
        - 192*K0^2*Kg*Kf^2*Msat*U0^2 + 16*K0^2*Kf^2*Msat^2*U0^2 - 120*K0*Kg^4*Kf^2*U0 + 240*K0*Kg^4*Kf*Msat*Por_t*U0 - 128*K0*Kg^4*Kf*Por_t*U0^2 ...
        - 120*K0*Kg^4*Msat^2*Por_t^2*U0 + 128*K0*Kg^4*Msat*Por_t^2*U0^2 - 240*K0*Kg^3*Kf^2*Msat*Por_t*U0 + 240*K0*Kg^3*Kf^2*Msat*U0 + 128*K0*Kg^3*Kf^2*Por_t*U0^2 ...
        - 192*K0*Kg^3*Kf^2*U0^2 + 240*K0*Kg^3*Kf*Msat^2*Por_t^2*U0 - 240*K0*Kg^3*Kf*Msat^2*Por_t*U0 - 256*K0*Kg^3*Kf*Msat*Por_t^2*U0^2 ...
        + 320*K0*Kg^3*Kf*Msat*Por_t*U0^2 - 120*K0*Kg^2*Kf^2*Msat^2*Por_t^2*U0 + 240*K0*Kg^2*Kf^2*Msat^2*Por_t*U0 - 120*K0*Kg^2*Kf^2*Msat^2*U0 ...
        + 128*K0*Kg^2*Kf^2*Msat*Por_t^2*U0^2 - 320*K0*Kg^2*Kf^2*Msat*Por_t*U0^2 + 224*K0*Kg^2*Kf^2*Msat*U0^2 - 32*K0*Kg^2*Kf*Msat^2*Por_t*U0^2 ...
        + 32*K0*Kg*Kf^2*Msat^2*Por_t*U0^2 - 32*K0*Kg*Kf^2*Msat^2*U0^2 + 16*Kg^4*Kf^2*U0^2 - 32*Kg^4*Kf*Msat*Por_t*U0^2 + 16*Kg^4*Msat^2*Por_t^2*U0^2 ...
        + 32*Kg^3*Kf^2*Msat*Por_t*U0^2 - 32*Kg^3*Kf^2*Msat*U0^2 - 32*Kg^3*Kf*Msat^2*Por_t^2*U0^2 + 32*Kg^3*Kf*Msat^2*Por_t*U0^2 + 16*Kg^2*Kf^2*Msat^2*Por_t^2*U0^2 ...
        - 32*Kg^2*Kf^2*Msat^2*Por_t*U0^2 + 16*Kg^2*Kf^2*Msat^2*U0^2)^(1/2) + 15*K0*Kg^2*Kf - 4*Kg^2*Kf*U0 - 15*K0*Kg*Kf*Msat + 16*K0*Kg*Kf*U0 + 4*K0*Kf*Msat*U0 ...
        + 4*Kg*Kf*Msat*U0 - 15*K0*Kg^2*Msat*Por_t + 24*K0*Kg^2*Por_t*U0 + 4*Kg^2*Msat*Por_t*U0 + 15*K0*Kg*Kf*Msat*Por_t - 24*K0*Kg*Kf*Por_t*U0 - 4*Kg*Kf*Msat*Por_t*U0)...
        /(2*(15*K0*Kg*Kf - 15*K0*Kf*Msat + 20*K0*Kf*U0 - 4*Kg*Kf*U0 + 4*Kf*Msat*U0 - 15*K0*Kg^2*Por_t + 4*Kg^2*Por_t*U0 + 15*K0*Kg*Kf*Por_t - 4*Kg*Kf*Por_t*U0));

        Kun2=(15*K0*Kg^2*Kf - (225*K0^2*Kg^4*Kf^2 - 450*K0^2*Kg^4*Kf*Msat*Por_t + 480*K0^2*Kg^4*Kf*Por_t*U0 + 225*K0^2*Kg^4*Msat^2*Por_t^2 - 480*K0^2*Kg^4*Msat*Por_t^2*U0 ...
        + 576*K0^2*Kg^4*Por_t^2*U0^2 + 450*K0^2*Kg^3*Kf^2*Msat*Por_t - 450*K0^2*Kg^3*Kf^2*Msat - 480*K0^2*Kg^3*Kf^2*Por_t*U0 + 720*K0^2*Kg^3*Kf^2*U0 ...
        - 450*K0^2*Kg^3*Kf*Msat^2*Por_t^2 + 450*K0^2*Kg^3*Kf*Msat^2*Por_t + 960*K0^2*Kg^3*Kf*Msat*Por_t^2*U0 - 1200*K0^2*Kg^3*Kf*Msat*Por_t*U0 ...
        - 1152*K0^2*Kg^3*Kf*Por_t^2*U0^2 + 768*K0^2*Kg^3*Kf*Por_t*U0^2 + 225*K0^2*Kg^2*Kf^2*Msat^2*Por_t^2 - 450*K0^2*Kg^2*Kf^2*Msat^2*Por_t ...
        + 225*K0^2*Kg^2*Kf^2*Msat^2 - 480*K0^2*Kg^2*Kf^2*Msat*Por_t^2*U0 + 1200*K0^2*Kg^2*Kf^2*Msat*Por_t*U0 - 840*K0^2*Kg^2*Kf^2*Msat*U0 ...
        + 576*K0^2*Kg^2*Kf^2*Por_t^2*U0^2 - 768*K0^2*Kg^2*Kf^2*Por_t*U0^2 + 576*K0^2*Kg^2*Kf^2*U0^2 + 120*K0^2*Kg^2*Kf*Msat^2*Por_t*U0 ...
        - 128*K0^2*Kg^2*Kf*Msat*Por_t*U0^2 - 120*K0^2*Kg*Kf^2*Msat^2*Por_t*U0 + 120*K0^2*Kg*Kf^2*Msat^2*U0 + 128*K0^2*Kg*Kf^2*Msat*Por_t*U0^2 ...
        - 192*K0^2*Kg*Kf^2*Msat*U0^2 + 16*K0^2*Kf^2*Msat^2*U0^2 - 120*K0*Kg^4*Kf^2*U0 + 240*K0*Kg^4*Kf*Msat*Por_t*U0 - 128*K0*Kg^4*Kf*Por_t*U0^2 ...
        - 120*K0*Kg^4*Msat^2*Por_t^2*U0 + 128*K0*Kg^4*Msat*Por_t^2*U0^2 - 240*K0*Kg^3*Kf^2*Msat*Por_t*U0 + 240*K0*Kg^3*Kf^2*Msat*U0 + 128*K0*Kg^3*Kf^2*Por_t*U0^2 ...
        - 192*K0*Kg^3*Kf^2*U0^2 + 240*K0*Kg^3*Kf*Msat^2*Por_t^2*U0 - 240*K0*Kg^3*Kf*Msat^2*Por_t*U0 - 256*K0*Kg^3*Kf*Msat*Por_t^2*U0^2 + 320*K0*Kg^3*Kf*Msat*Por_t*U0^2 ...
        - 120*K0*Kg^2*Kf^2*Msat^2*Por_t^2*U0 + 240*K0*Kg^2*Kf^2*Msat^2*Por_t*U0 - 120*K0*Kg^2*Kf^2*Msat^2*U0 + 128*K0*Kg^2*Kf^2*Msat*Por_t^2*U0^2 ...
        - 320*K0*Kg^2*Kf^2*Msat*Por_t*U0^2 + 224*K0*Kg^2*Kf^2*Msat*U0^2 - 32*K0*Kg^2*Kf*Msat^2*Por_t*U0^2 + 32*K0*Kg*Kf^2*Msat^2*Por_t*U0^2 ...
        - 32*K0*Kg*Kf^2*Msat^2*U0^2 + 16*Kg^4*Kf^2*U0^2 - 32*Kg^4*Kf*Msat*Por_t*U0^2 + 16*Kg^4*Msat^2*Por_t^2*U0^2 + 32*Kg^3*Kf^2*Msat*Por_t*U0^2 ...
        - 32*Kg^3*Kf^2*Msat*U0^2 - 32*Kg^3*Kf*Msat^2*Por_t^2*U0^2 + 32*Kg^3*Kf*Msat^2*Por_t*U0^2 + 16*Kg^2*Kf^2*Msat^2*Por_t^2*U0^2 - 32*Kg^2*Kf^2*Msat^2*Por_t*U0^2 ...
        + 16*Kg^2*Kf^2*Msat^2*U0^2)^(1/2) - 4*Kg^2*Kf*U0 - 15*K0*Kg*Kf*Msat + 16*K0*Kg*Kf*U0 + 4*K0*Kf*Msat*U0 + 4*Kg*Kf*Msat*U0 - 15*K0*Kg^2*Msat*Por_t ...
        + 24*K0*Kg^2*Por_t*U0 + 4*Kg^2*Msat*Por_t*U0 + 15*K0*Kg*Kf*Msat*Por_t - 24*K0*Kg*Kf*Por_t*U0 - 4*Kg*Kf*Msat*Por_t*U0)/(2*(15*K0*Kg*Kf - 15*K0*Kf*Msat ...
        + 20*K0*Kf*U0 - 4*Kg*Kf*U0 + 4*Kf*Msat*U0 - 15*K0*Kg^2*Por_t + 4*Kg^2*Por_t*U0 + 15*K0*Kg*Kf*Por_t - 4*Kg*Kf*Por_t*U0));


        Uun1= (3*((225*K0^2*Kg^4*Kf^2 - 450*K0^2*Kg^4*Kf*Msat*Por_t + 480*K0^2*Kg^4*Kf*Por_t*U0 + 225*K0^2*Kg^4*Msat^2*Por_t^2 - 480*K0^2*Kg^4*Msat*Por_t^2*U0 ...
        + 576*K0^2*Kg^4*Por_t^2*U0^2 + 450*K0^2*Kg^3*Kf^2*Msat*Por_t - 450*K0^2*Kg^3*Kf^2*Msat - 480*K0^2*Kg^3*Kf^2*Por_t*U0 + 720*K0^2*Kg^3*Kf^2*U0 ...
        - 450*K0^2*Kg^3*Kf*Msat^2*Por_t^2 + 450*K0^2*Kg^3*Kf*Msat^2*Por_t + 960*K0^2*Kg^3*Kf*Msat*Por_t^2*U0 - 1200*K0^2*Kg^3*Kf*Msat*Por_t*U0 ...
        - 1152*K0^2*Kg^3*Kf*Por_t^2*U0^2 + 768*K0^2*Kg^3*Kf*Por_t*U0^2 + 225*K0^2*Kg^2*Kf^2*Msat^2*Por_t^2 - 450*K0^2*Kg^2*Kf^2*Msat^2*Por_t ...
        + 225*K0^2*Kg^2*Kf^2*Msat^2 - 480*K0^2*Kg^2*Kf^2*Msat*Por_t^2*U0 + 1200*K0^2*Kg^2*Kf^2*Msat*Por_t*U0 - 840*K0^2*Kg^2*Kf^2*Msat*U0 ...
        + 576*K0^2*Kg^2*Kf^2*Por_t^2*U0^2 - 768*K0^2*Kg^2*Kf^2*Por_t*U0^2 + 576*K0^2*Kg^2*Kf^2*U0^2 + 120*K0^2*Kg^2*Kf*Msat^2*Por_t*U0 ...
        - 128*K0^2*Kg^2*Kf*Msat*Por_t*U0^2 - 120*K0^2*Kg*Kf^2*Msat^2*Por_t*U0 + 120*K0^2*Kg*Kf^2*Msat^2*U0 + 128*K0^2*Kg*Kf^2*Msat*Por_t*U0^2 ...
        - 192*K0^2*Kg*Kf^2*Msat*U0^2 + 16*K0^2*Kf^2*Msat^2*U0^2 - 120*K0*Kg^4*Kf^2*U0 + 240*K0*Kg^4*Kf*Msat*Por_t*U0 - 128*K0*Kg^4*Kf*Por_t*U0^2 ...
        - 120*K0*Kg^4*Msat^2*Por_t^2*U0 + 128*K0*Kg^4*Msat*Por_t^2*U0^2 - 240*K0*Kg^3*Kf^2*Msat*Por_t*U0 + 240*K0*Kg^3*Kf^2*Msat*U0 + 128*K0*Kg^3*Kf^2*Por_t*U0^2 ...
        - 192*K0*Kg^3*Kf^2*U0^2 + 240*K0*Kg^3*Kf*Msat^2*Por_t^2*U0 - 240*K0*Kg^3*Kf*Msat^2*Por_t*U0 - 256*K0*Kg^3*Kf*Msat*Por_t^2*U0^2 + 320*K0*Kg^3*Kf*Msat*Por_t*U0^2 ...
        - 120*K0*Kg^2*Kf^2*Msat^2*Por_t^2*U0 + 240*K0*Kg^2*Kf^2*Msat^2*Por_t*U0 - 120*K0*Kg^2*Kf^2*Msat^2*U0 + 128*K0*Kg^2*Kf^2*Msat*Por_t^2*U0^2 ...
        - 320*K0*Kg^2*Kf^2*Msat*Por_t*U0^2 + 224*K0*Kg^2*Kf^2*Msat*U0^2 - 32*K0*Kg^2*Kf*Msat^2*Por_t*U0^2 + 32*K0*Kg*Kf^2*Msat^2*Por_t*U0^2 - 32*K0*Kg*Kf^2*Msat^2*U0^2 ...
        + 16*Kg^4*Kf^2*U0^2 - 32*Kg^4*Kf*Msat*Por_t*U0^2 + 16*Kg^4*Msat^2*Por_t^2*U0^2 + 32*Kg^3*Kf^2*Msat*Por_t*U0^2 - 32*Kg^3*Kf^2*Msat*U0^2 - 32*Kg^3*Kf*Msat^2*Por_t^2*U0^2 ...
        + 32*Kg^3*Kf*Msat^2*Por_t*U0^2 + 16*Kg^2*Kf^2*Msat^2*Por_t^2*U0^2 - 32*Kg^2*Kf^2*Msat^2*Por_t*U0^2 + 16*Kg^2*Kf^2*Msat^2*U0^2)^(1/2) - 15*K0*Kg^2*Kf + 4*Kg^2*Kf*U0 ...
        + 15*K0*Kg*Kf*Msat + 16*K0*Kg*Kf*U0 + 4*K0*Kf*Msat*U0 - 4*Kg*Kf*Msat*U0 + 15*K0*Kg^2*Msat*Por_t + 24*K0*Kg^2*Por_t*U0 - 4*Kg^2*Msat*Por_t*U0 - 15*K0*Kg*Kf*Msat*Por_t ...
        - 24*K0*Kg*Kf*Por_t*U0 + 4*Kg*Kf*Msat*Por_t*U0))/(8*(15*K0*Kg*Kf + 4*K0*Kf*U0 - 4*Kg*Kf*U0 + 15*K0*Kg^2*Por_t - 4*Kg^2*Por_t*U0 - 15*K0*Kg*Kf*Por_t + 4*Kg*Kf*Por_t*U0));

        Uun2=(3*(4*Kg^2*Kf*U0 - 15*K0*Kg^2*Kf - (225*K0^2*Kg^4*Kf^2 - 450*K0^2*Kg^4*Kf*Msat*Por_t + 480*K0^2*Kg^4*Kf*Por_t*U0 + 225*K0^2*Kg^4*Msat^2*Por_t^2 - 480*K0^2*Kg^4*Msat*Por_t^2*U0 ...
         + 576*K0^2*Kg^4*Por_t^2*U0^2 + 450*K0^2*Kg^3*Kf^2*Msat*Por_t - 450*K0^2*Kg^3*Kf^2*Msat - 480*K0^2*Kg^3*Kf^2*Por_t*U0 + 720*K0^2*Kg^3*Kf^2*U0 - 450*K0^2*Kg^3*Kf*Msat^2*Por_t^2 ...
         + 450*K0^2*Kg^3*Kf*Msat^2*Por_t + 960*K0^2*Kg^3*Kf*Msat*Por_t^2*U0 - 1200*K0^2*Kg^3*Kf*Msat*Por_t*U0 - 1152*K0^2*Kg^3*Kf*Por_t^2*U0^2 + 768*K0^2*Kg^3*Kf*Por_t*U0^2 ...
         + 225*K0^2*Kg^2*Kf^2*Msat^2*Por_t^2 - 450*K0^2*Kg^2*Kf^2*Msat^2*Por_t + 225*K0^2*Kg^2*Kf^2*Msat^2 - 480*K0^2*Kg^2*Kf^2*Msat*Por_t^2*U0 + 1200*K0^2*Kg^2*Kf^2*Msat*Por_t*U0 ...
         - 840*K0^2*Kg^2*Kf^2*Msat*U0 + 576*K0^2*Kg^2*Kf^2*Por_t^2*U0^2 - 768*K0^2*Kg^2*Kf^2*Por_t*U0^2 + 576*K0^2*Kg^2*Kf^2*U0^2 + 120*K0^2*Kg^2*Kf*Msat^2*Por_t*U0 ...
         - 128*K0^2*Kg^2*Kf*Msat*Por_t*U0^2 - 120*K0^2*Kg*Kf^2*Msat^2*Por_t*U0 + 120*K0^2*Kg*Kf^2*Msat^2*U0 + 128*K0^2*Kg*Kf^2*Msat*Por_t*U0^2 - 192*K0^2*Kg*Kf^2*Msat*U0^2 ...
         + 16*K0^2*Kf^2*Msat^2*U0^2 - 120*K0*Kg^4*Kf^2*U0 + 240*K0*Kg^4*Kf*Msat*Por_t*U0 - 128*K0*Kg^4*Kf*Por_t*U0^2 - 120*K0*Kg^4*Msat^2*Por_t^2*U0 + 128*K0*Kg^4*Msat*Por_t^2*U0^2 ...
         - 240*K0*Kg^3*Kf^2*Msat*Por_t*U0 + 240*K0*Kg^3*Kf^2*Msat*U0 + 128*K0*Kg^3*Kf^2*Por_t*U0^2 - 192*K0*Kg^3*Kf^2*U0^2 + 240*K0*Kg^3*Kf*Msat^2*Por_t^2*U0 - 240*K0*Kg^3*Kf*Msat^2*Por_t*U0 ...
         - 256*K0*Kg^3*Kf*Msat*Por_t^2*U0^2 + 320*K0*Kg^3*Kf*Msat*Por_t*U0^2 - 120*K0*Kg^2*Kf^2*Msat^2*Por_t^2*U0 + 240*K0*Kg^2*Kf^2*Msat^2*Por_t*U0 - 120*K0*Kg^2*Kf^2*Msat^2*U0 ...
         + 128*K0*Kg^2*Kf^2*Msat*Por_t^2*U0^2 - 320*K0*Kg^2*Kf^2*Msat*Por_t*U0^2 + 224*K0*Kg^2*Kf^2*Msat*U0^2 - 32*K0*Kg^2*Kf*Msat^2*Por_t*U0^2 + 32*K0*Kg*Kf^2*Msat^2*Por_t*U0^2 ...
         - 32*K0*Kg*Kf^2*Msat^2*U0^2 + 16*Kg^4*Kf^2*U0^2 - 32*Kg^4*Kf*Msat*Por_t*U0^2 + 16*Kg^4*Msat^2*Por_t^2*U0^2 + 32*Kg^3*Kf^2*Msat*Por_t*U0^2 - 32*Kg^3*Kf^2*Msat*U0^2 ...
         - 32*Kg^3*Kf*Msat^2*Por_t^2*U0^2 + 32*Kg^3*Kf*Msat^2*Por_t*U0^2 + 16*Kg^2*Kf^2*Msat^2*Por_t^2*U0^2 - 32*Kg^2*Kf^2*Msat^2*Por_t*U0^2 + 16*Kg^2*Kf^2*Msat^2*U0^2)^(1/2) ...
         + 15*K0*Kg*Kf*Msat + 16*K0*Kg*Kf*U0 + 4*K0*Kf*Msat*U0 - 4*Kg*Kf*Msat*U0 + 15*K0*Kg^2*Msat*Por_t + 24*K0*Kg^2*Por_t*U0 - 4*Kg^2*Msat*Por_t*U0 - 15*K0*Kg*Kf*Msat*Por_t ...
         - 24*K0*Kg*Kf*Por_t*U0 + 4*Kg*Kf*Msat*Por_t*U0))/(8*(15*K0*Kg*Kf + 4*K0*Kf*U0 - 4*Kg*Kf*U0 + 15*K0*Kg^2*Por_t - 4*Kg^2*Por_t*U0 - 15*K0*Kg*Kf*Por_t + 4*Kg*Kf*Por_t*U0));

        if(Kun1>0)
         Kdry0=Kun1;
        else
         Kdry0=Kun2;
        end  

        if(Uun1*Uun2>0)
         Gdry0=min(Uun1,Uun2);
        else
         Gdry0=max(Uun1,Uun2);
        end 

        Kh=Kdry0;%%% approximation according to equation 24 in Gurevich et al's paper(2009) when Sw>0.1  
    end



%%%%%%%%%%%% COMPUTE C TENSORS AS FUNCTIONS OF SATURATION at INCIDENT ANGLE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    SSw=linspace(0,1.00,nSw);
    Inc=linspace(0,90,nInc);

    %%% modeling results of patchy model
    vp=zeros(nSw,nInc);
    mp=zeros(nSw,nInc);
    vp_uniform=zeros(nSw,nInc);
    mp_uniform=zeros(nSw,nInc);

    %%% modeling results of high-frequency limit
    vpr=zeros(nSw,nInc);
    mpr=zeros(nSw,nInc);

    %%% modeling results of low-frequency limit
    vpv=zeros(nSw,nInc);
    mpv=zeros(nSw,nInc);

    %%% modeling results of combining models
    vp_final=zeros(nSw,nInc);
    mp_final=zeros(nSw,nInc);

    C_low11=zeros(1,nSw);
    C_low12=zeros(1,nSw);
    C_low13=zeros(1,nSw);
    C_low22=zeros(1,nSw);
    C_low23=zeros(1,nSw);
    C_low33=zeros(1,nSw);
    C_low44=zeros(1,nSw);
    C_low55=zeros(1,nSw);
    C_low66=zeros(1,nSw);

    for iSw=1:nSw

        Sw=SSw(iSw);
        H=H0;%/Sw;
        Sg=1-Sw;            
        Den=Sw*Den_gw+Sg*Den_gg;
        Den_gwg=Sw*Den_gw+Sg*Den_gg;    
        
        if Kdry_method==3  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% squirt model %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
            Kfl=Sw*Kwater+(1-Sw)*Kgas;
            Kuf=1./(1./Kh+1./(1./(1./K0-1./Kh)+1./((1./Kfl-1./Kg)*phi_c)));  
            Kdry=Kuf;            
            Gdry=1./(1./U0-4./15*(1./K0-1./Kdry)); % equation 2 in Gurevich et al's paper(2009)    
            Mdry=Kdry+4./3*Gdry;% dry bulk modulus for fractured sample      
        end        
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% dry model %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
         %%% Hudson or linear slip model
        dN =(4*hc*(Kdry+4/3*Gdry)*(Kdry+4/3*Gdry))/(3*Gdry*(Kdry+1/3*Gdry));
        dT =(16*hc*(Kdry+4/3*Gdry))/(9*(Kdry+2/3*Gdry));

        Cp11 = Kdry + 4./3 * Gdry;
        Cp12 = Kdry - 2./3 * Gdry;
        Cp13 = Kdry - 2./3 * Gdry;
        Cp22 = Kdry + 4./3 * Gdry;
        Cp23 = Kdry - 2./3 * Gdry;
        Cp33 = Kdry + 4./3 * Gdry;
        Cp44 = Gdry;
        Cp55 = Gdry;
        Cp66 = Gdry;

        rb = Cp13 / Cp11;
        Cdry11 = Cp11 * ( 1 - rb * rb * dN );
        Cdry12 = Cp12 * ( 1 - rb * dN );
        Cdry13 = Cp13 * ( 1 - dN);
        Cdry22 = Cp22 * ( 1 - rb * rb * dN );
        Cdry23 = Cp23 * ( 1 - dN);
        Cdry33 = Cp33 * ( 1 - dN);
        Cdry44 = Cp44 * ( 1 - dT );
        Cdry55 = Cp55 * ( 1 - dT );
        Cdry66 = Cp66;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% uniform model %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%         alpha1=1-(Cdry11+Cdry12+Cdry13)/3/Kg;
%         alpha2=1-(Cdry12+Cdry22+Cdry23)/3/Kg;
%         alpha3=1-(Cdry13+Cdry23+Cdry33)/3/Kg;
%         Kstar=(Cdry11+Cdry12+Cdry13+Cdry12+Cdry22+Cdry23+Cdry13+Cdry23+Cdry33)/9;
        Kf=1./(Sw/Kwater+Sg/Kgas);
%         Mf=Kg/((1-Kstar/Kg)-Por*(1-Kg/Kf));
% 
%         C_low11(iSw)=Cdry11+alpha1*alpha1*Mf;
%         C_low12(iSw)=Cdry12+alpha1*alpha2*Mf;
%         C_low13(iSw)=Cdry13+alpha1*alpha3*Mf;
%         C_low22(iSw)=Cdry22+alpha2*alpha2*Mf;
%         C_low23(iSw)=Cdry23+alpha2*alpha3*Mf;
%         C_low33(iSw)=Cdry33+alpha3*alpha3*Mf;
%         C_low44(iSw)=Cdry44;
%         C_low55(iSw)=Cdry55;
%         C_low66(iSw)=Cdry66;  

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%% patchy model %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
        %%%% KG_model %%%%%%%%%%%%%%  
        [C_low11(iSw) C_low12(iSw) C_low13(iSw) C_low22(iSw) C_low23(iSw) C_low33(iSw) C_low44(iSw) C_low55(iSw) Cw_low66(iSw)]=KG_model(Por_f,H,Kg,Gg,Kdry,Gdry,Kf,Ywater,Por_p,Por_m_re,Kf,dN,dT,w0);
        [Cwater11 Cwater12 Cwater13 Cwater22 Cwater23 Cwater33 Cwater44 Cwater55 Cwater66]=KG_model(Por_f,H,Kg,Gg,Kdry,Gdry,Kwater,Ywater,Por_p,Por_m_re,Kwater,dN,dT,w0);
        [Cgas11 Cgas12 Cgas13 Cgas22 Cgas23 Cgas33 Cgas44 Cgas55 Cgas66]=KG_model(Por_f,H,Kg,Gg,Kdry,Gdry,Kgas,Ygas,Por_p,Por_m_re,Kgas,dN,dT,w0);
        
%         alpha1=1-(Cdry11+Cdry12+Cdry13)/3/Kg;
%         alpha2=1-(Cdry12+Cdry22+Cdry23)/3/Kg;
%         alpha3=1-(Cdry13+Cdry23+Cdry33)/3/Kg;
%         Kstar=(Cdry11+Cdry12+Cdry13+Cdry12+Cdry22+Cdry23+Cdry13+Cdry23+Cdry33)/9;
%         Kf=Kwater;
%         Mf=Kg/((1-Kstar/Kg)-Por*(1-Kg/Kf));
% 
%         Cwater11=Cdry11+alpha1*alpha1*Mf;
%         Cwater12=Cdry12+alpha1*alpha2*Mf;
%         Cwater13=Cdry13+alpha1*alpha3*Mf;
%         Cwater22=Cdry22+alpha2*alpha2*Mf;
%         Cwater23=Cdry23+alpha2*alpha3*Mf;
%         Cwater33=Cdry33+alpha3*alpha3*Mf;
%         Cwater44=Cdry44;
%         Cwater55=Cdry55;
%         Cwater66=Cdry66;  
%         
%         Kf=Kgas;
%         Mf=Kg/((1-Kstar/Kg)-Por*(1-Kg/Kf));
% 
%         Cgas11=Cdry11+alpha1*alpha1*Mf;
%         Cgas12=Cdry12+alpha1*alpha2*Mf;
%         Cgas13=Cdry13+alpha1*alpha3*Mf;
%         Cgas22=Cdry22+alpha2*alpha2*Mf;
%         Cgas23=Cdry23+alpha2*alpha3*Mf;
%         Cgas33=Cdry33+alpha3*alpha3*Mf;
%         Cgas44=Cdry44;
%         Cgas55=Cdry55;
%         Cgas66=Cdry66;  
        
        %%%%%%%%%%%%%%  Backus Average %%%%%%%%%%%%%%%%
        Cwg11=Sw*(Cwater11-Cwater13*Cwater13/Cwater33)+Sg*(Cgas11-Cgas13*Cgas13/Cgas33)+(Sw*(Cwater13/Cwater33)+Sg*(Cgas13/Cgas33))*(Sw*(Cwater13/Cwater33)+Sg*(Cgas13/Cgas33))*1./(Sw/Cwater33+Sg/Cgas33);
        Cwg12=Sw*(Cwater12-Cwater13*Cwater13/Cwater33)+Sg*(Cgas12-Cgas13*Cgas13/Cgas33)+(Sw*(Cwater13/Cwater33)+Sg*(Cgas13/Cgas33))*(Sw*(Cwater13/Cwater33)+Sg*(Cgas13/Cgas33))*1./(Sw/Cwater33+Sg/Cgas33);
        Cwg13=(Sw*(Cwater13/Cwater33)+Sg*(Cgas13/Cgas33))*1./(Sw/Cwater33+Sg/Cgas33);

        Cwg22=Cwg11;
        Cwg23=Cwg13;
        Cwg33=1./(Sw/Cwater33+Sg/Cgas33);

        Cwg44=1./(Sw/Cwater44+Sg/Cgas44);
        Cwg55=Cwg44;
        Cwg66=Sw*Cwater66+Sg*Cgas66;   

        %%%%%%%%%%%%%%%%% Ruess bound for mixture of KG_water and KG_gas %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
        Cwgr11=Sw*Cwater11+Sg*Cgas11;
        Cwgr12=Sw*Cwater12+Sg*Cgas12;
        Cwgr13=Sw*Cwater13+Sg*Cgas13;
        Cwgr22=Sw*Cwater22+Sg*Cgas22;
        Cwgr23=Sw*Cwater23+Sg*Cgas23;
        Cwgr33=Sw*Cwater33+Sg*Cgas33;
        Cwgr44=Sw*Cwater44+Sg*Cgas44;
        Cwgr55=Sw*Cwater55+Sg*Cgas55;
        Cwgr66=Sw*Cwater66+Sg*Cgas66;

        %%%%%%%%%%%%%%%%% Voigt bound for mixture of KG_water and KG_gas %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
        M=Cwater11*Cwater23*Cwater23+Cwater22*Cwater13*Cwater13+Cwater33*Cwater12*Cwater12-2*Cwater12*Cwater13*Cwater23-Cwater11*Cwater22*Cwater33;
        Swater11=(Cwater23*Cwater23-Cwater22*Cwater33)/M;
        Swater12=(Cwater12*Cwater33-Cwater13*Cwater23)/M;
        Swater13=(Cwater13*Cwater22-Cwater12*Cwater23)/M;
        Swater22=(Cwater13*Cwater13-Cwater11*Cwater33)/M;
        Swater23=(Cwater23*Cwater11-Cwater12*Cwater13)/M;
        Swater33=(Cwater12*Cwater12-Cwater11*Cwater22)/M;
        Swater44=1./Cwater44;
        Swater55=1./Cwater55;
        Swater66=1./Cwater66;

        M=Cgas11*Cgas23*Cgas23+Cgas22*Cgas13*Cgas13+Cgas33*Cgas12*Cgas12-2*Cgas12*Cgas13*Cgas23-Cgas11*Cgas22*Cgas33;
        Sgas11=(Cgas23*Cgas23-Cgas22*Cgas33)/M;
        Sgas12=(Cgas12*Cgas33-Cgas13*Cgas23)/M;
        Sgas13=(Cgas13*Cgas22-Cgas12*Cgas23)/M;
        Sgas22=(Cgas13*Cgas13-Cgas11*Cgas33)/M;
        Sgas23=(Cgas23*Cgas11-Cgas12*Cgas13)/M;
        Sgas33=(Cgas12*Cgas12-Cgas11*Cgas22)/M;
        Sgas44=1./Cgas44;
        Sgas55=1./Cgas55;
        Sgas66=1./Cgas66;

        Swgv11=Sw*Swater11+Sg*Sgas11;
        Swgv12=Sw*Swater12+Sg*Sgas12;
        Swgv13=Sw*Swater13+Sg*Sgas13;
        Swgv22=Sw*Swater22+Sg*Sgas22;
        Swgv23=Sw*Swater23+Sg*Sgas23;
        Swgv33=Sw*Swater33+Sg*Sgas33;
        Swgv44=Sw*Swater44+Sg*Sgas44;
        Swgv55=Sw*Swater55+Sg*Sgas55;
        Swgv66=Sw*Swater66+Sg*Sgas66;

        M=Swgv11*Swgv23*Swgv23+Swgv22*Swgv13*Swgv13+Swgv33*Swgv12*Swgv12-2*Swgv12*Swgv13*Swgv23-Swgv11*Swgv22*Swgv33;
        Cwgv11=(Swgv23*Swgv23-Swgv22*Swgv33)/M;
        Cwgv12=(Swgv12*Swgv33-Swgv13*Swgv23)/M;
        Cwgv13=(Swgv13*Swgv22-Swgv12*Swgv23)/M;
        Cwgv22=(Swgv13*Swgv13-Swgv11*Swgv33)/M;
        Cwgv23=(Swgv23*Swgv11-Swgv12*Swgv13)/M;
        Cwgv33=(Swgv12*Swgv12-Swgv11*Swgv22)/M;
        Cwgv44=1./Swgv44;
        Cwgv55=1./Swgv55;
        Cwgv66=1./Swgv66;


 %%%%%%%%%%%% COMPUTE VELOCITIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
        for iInc=1:nInc

            theta=Inc(iInc)*pi/180;                
            phi=0;
            Dip=0;

            sinT=sin(theta);
            cosT=cos(theta);
            sinA=sin(phi);
            cosA=cos(phi);
            sinD=sin(Dip);
            cosD=cos(Dip);

            E=(-sinT*cosA*sinD+cosT*cosD)*(-sinT*cosA*sinD+cosT*cosD);
            G=(sinT*cosA*cosD+cosT*sinD)*(sinT*cosA*cosD+cosT*sinD)+sinT*sinT*sinA*sinA;

            D=((Cwg11-Cwg44)*G-(Cwg33-Cwg44)*E)*((Cwg11-Cwg44)*G-(Cwg33-Cwg44)*E)+4*(Cwg13+Cwg44)*(Cwg13+Cwg44)*G*E;
            vv=sqrt((Cwg44+Cwg11*G+Cwg33*E+sqrt(D))/(2*Den_gwg));
            vp(iSw,iInc)=1./(real(1./vv));
            mp(iSw,iInc)=Den_gwg*vp(iSw,iInc)*vp(iSw,iInc);
              
            Dr=((Cwgr11-Cwgr44)*G-(Cwgr33-Cwgr44)*E)*((Cwgr11-Cwgr44)*G-(Cwgr33-Cwgr44)*E)+4*(Cwgr13+Cwgr44)*(Cwgr13+Cwgr44)*G*E;
            vv=sqrt((Cwgr44+Cwgr11*G+Cwgr33*E+sqrt(Dr))/(2*Den_gwg));
            vpr(iSw,iInc)=1./(real(1./vv));
            mpr(iSw,iInc)=Den_gwg*vpr(iSw,iInc)*vpr(iSw,iInc);

            Dv=((Cwgv11-Cwgv44)*G-(Cwgv33-Cwgv44)*E)*((Cwgv11-Cwgv44)*G-(Cwgv33-Cwgv44)*E)+4*(Cwgv13+Cwgv44)*(Cwgv13+Cwgv44)*G*E;
            vv=sqrt((Cwgv44+Cwgv11*G+Cwgv33*E+sqrt(Dv))/(2*Den_gwg));
            vpv(iSw,iInc)=1./(real(1./vv));
            mpv(iSw,iInc)=Den_gwg*vpv(iSw,iInc)*vpv(iSw,iInc);  
            
            D_unique=((C_low11(iSw)-C_low44(iSw))*G-(C_low33(iSw)-C_low44(iSw))*E)*((C_low11(iSw)-C_low44(iSw))*G-(C_low33(iSw)-C_low44(iSw))*E)+4*(C_low13(iSw)+C_low44(iSw))*(C_low13(iSw)+C_low44(iSw))*G*E;
            vv=sqrt((C_low44(iSw)+C_low11(iSw)*G+C_low33(iSw)*E+sqrt(D_unique))/(2*Den_gwg));
            vp_uniform(iSw,iInc)=1./(real(1./vv));
            mp_uniform(iSw,iInc)=Den_gwg*vp_uniform(iSw,iInc)*vp_uniform(iSw,iInc);
        end 
        
        %%%%%%%%% final modeling results of combining models %%%%%%%%%%%%%
        if Sw<0.5
            vp_final(iSw,:)=vp_uniform(iSw,:);
            mp_final(iSw,:)=mp_uniform(iSw,:);
        else
            vp_final(iSw,:)=vp(iSw,:);
            mp_final(iSw,:)=mp(iSw,:);
        end
    end


%%%%%%%%%%%% PLOT RESUTLS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   

    if Kdry_method==1 % Kdry from data at Sw=0
%         figure()
%         subplot(1,2,1)
%         set(gcf,'outerposition',get(0,'screensize'));
%         h7=plot(SSw,mpr(:,3),'r:.','LineWidth',4);
%         hold on; 
%         h8=plot(SSw,mpr(:,2),'g:.','LineWidth',4);
%         hold on;  
%         h9=plot(SSw,mpr(:,1),'b:.','LineWidth',4);
%         hold on; 
%         h10=plot(SSw,mpv(:,3),'r:.','LineWidth',4);
%         hold on;
%         h11=plot(SSw,mpv(:,2),'g:.','LineWidth',4);
%         hold on;   
%         h12=plot(SSw,mpv(:,1),'b:.','LineWidth',4);
%         hold on;  
%         h1=plot(fractured_correct(:,1)*0.01,fractured_correct(:,8),'ws','markerfacecolor','r','markersize',20);
%         hold on;  
%         h2=plot(SSw,mp(:,3),'r','LineWidth',4);
%         hold on; 
%         h3=plot(fractured_correct(:,1)*0.01,fractured_correct(:,7),'wo','markerfacecolor','g','markersize',20);
%         hold on; 
%         h4=plot(SSw,mp(:,2),'g','LineWidth',4);
%         hold on; 
%         h5=plot(fractured_correct(:,1)*0.01,fractured_correct(:,6),'wd','markerfacecolor','b','markersize',20);
%         hold on;
%         h6=plot(SSw,mp(:,1),'b','LineWidth',4);
%         hold on;
%         set(gca,'FontSize',30);
%         set(gca,'xTick',[0 0.25 0.5 0.75 1]);
%         xlabel('Sw','fontsize',40);
%         ylabel('Mp(Pa)','fontsize',40);    
%         set(gca,'ylim',[14e9 25e9]);
% 
%         subplot(1,2,2)
%         set(gcf,'outerposition',get(0,'screensize'));
%         h7=plot(SSw,vpr(:,3),'r:.','LineWidth',4);
%         hold on; 
%         h8=plot(SSw,vpr(:,2),'g:.','LineWidth',4);
%         hold on;  
%         h9=plot(SSw,vpr(:,1),'b:.','LineWidth',4);
%         hold on; 
%         h10=plot(SSw,vpv(:,3),'r:.','LineWidth',4);
%         hold on;
%         h11=plot(SSw,vpv(:,2),'g:.','LineWidth',4);
%         hold on;   
%         h12=plot(SSw,vpv(:,1),'b:.','LineWidth',4);
%         hold on;  
%         h1=plot(fractured_correct(:,1)*0.01,fractured_correct(:,5),'ws','markerfacecolor','r','markersize',20);
%         hold on;  
%         h2=plot(SSw,vp(:,3),'r','LineWidth',4);
%         hold on; 
%         h3=plot(fractured_correct(:,1)*0.01,fractured_correct(:,4),'wo','markerfacecolor','g','markersize',20);
%         hold on; 
%         h4=plot(SSw,vp(:,2),'g','LineWidth',4);
%         hold on; 
%         h5=plot(fractured_correct(:,1)*0.01,fractured_correct(:,3),'wd','markerfacecolor','b','markersize',20);
%         hold on;
%         h6=plot(SSw,vp(:,1),'b','LineWidth',4);
%         hold on;
%         set(gca,'FontSize',30);
%         set(gca,'yTick',[2700 2900 3100 3300 3500]);
%         set(gca,'xTick',[0 0.25 0.5 0.75 1]);
%         xlabel('Sw','fontsize',40);
%         ylabel('Vp(m/s)','fontsize',40);    
%         set(gca,'ylim',[2700 3500]);

    elseif Kdry_method==2  % Kdry from data at Sw=100

%         figure()
%         subplot(1,2,1)
%         set(gcf,'outerposition',get(0,'screensize'));
%         h7=plot(SSw,mpr(:,3),'r:.','LineWidth',4);
%         hold on; 
%         h8=plot(SSw,mpr(:,2),'g:.','LineWidth',4);
%         hold on;  
%         h9=plot(SSw,mpr(:,1),'b:.','LineWidth',4);
%         hold on; 
%         h10=plot(SSw,mpv(:,3),'r:.','LineWidth',4);
%         hold on;
%         h11=plot(SSw,mpv(:,2),'g:.','LineWidth',4);
%         hold on;   
%         h12=plot(SSw,mpv(:,1),'b:.','LineWidth',4);
%         hold on;  
%         h1=plot(fractured_correct(:,1)*0.01,fractured_correct(:,8),'ws','markerfacecolor','r','markersize',20);
%         hold on;  
%         h2=plot(SSw,mp(:,3),'r','LineWidth',4);
%         hold on; 
%         h3=plot(fractured_correct(:,1)*0.01,fractured_correct(:,7),'wo','markerfacecolor','g','markersize',20);
%         hold on; 
%         h4=plot(SSw,mp(:,2),'g','LineWidth',4);
%         hold on; 
%         h5=plot(fractured_correct(:,1)*0.01,fractured_correct(:,6),'wd','markerfacecolor','b','markersize',20);
%         hold on;
%         h6=plot(SSw,mp(:,1),'b','LineWidth',4);
%         hold on;
%         set(gca,'FontSize',30);
%         set(gca,'xTick',[0 0.25 0.5 0.75 1]);
%         xlabel('Sw','fontsize',40);
%         ylabel('Mp(Pa)','fontsize',40);    
%         set(gca,'ylim',[14e9 25e9]);
%         subplot(1,2,2)
%         set(gcf,'outerposition',get(0,'screensize'));
% 
%         h7=plot(SSw,vpr(:,3),'r:.','LineWidth',4);
%         hold on; 
%         h8=plot(SSw,vpr(:,2),'g:.','LineWidth',4);
%         hold on;  
%         h9=plot(SSw,vpr(:,1),'b:.','LineWidth',4);
%         hold on; 
%         h10=plot(SSw,vpv(:,3),'r:.','LineWidth',4);
%         hold on;
%         h11=plot(SSw,vpv(:,2),'g:.','LineWidth',4);
%         hold on;   
%         h12=plot(SSw,vpv(:,1),'b:.','LineWidth',4);
%         hold on;  
%         h1=plot(fractured_correct(:,1)*0.01,fractured_correct(:,5),'ws','markerfacecolor','r','markersize',20);
%         hold on;  
%         h2=plot(SSw,vp(:,3),'r','LineWidth',4);
%         hold on; 
%         h3=plot(fractured_correct(:,1)*0.01,fractured_correct(:,4),'wo','markerfacecolor','g','markersize',20);
%         hold on; 
%         h4=plot(SSw,vp(:,2),'g','LineWidth',4);
%         hold on; 
%         h5=plot(fractured_correct(:,1)*0.01,fractured_correct(:,3),'wd','markerfacecolor','b','markersize',20);
%         hold on;
%         h6=plot(SSw,vp(:,1),'b','LineWidth',4);
%         hold on;
%         set(gca,'FontSize',30);
%         set(gca,'yTick',[2700 2900 3100 3300 3500]);
%         set(gca,'xTick',[0 0.25 0.5 0.75 1]);
%         xlabel('Sw','fontsize',40);
%         ylabel('Vp(m/s)','fontsize',40);    
%         set(gca,'ylim',[2700 3500]);
% 
%         figure()
%         subplot(1,2,1)
%         set(gcf,'outerposition',get(0,'screensize'));
%         plot(SSw,mp_unique(:,3),'r','LineWidth',4);
%         hold on;  
%         plot(SSw,mp_unique(:,2),'g','LineWidth',4);
%         hold on;  
%         plot(SSw,mp_unique(:,1),'b','LineWidth',4);
%         hold on;  
%         h7=plot(SSw,mpr(:,3),'r:.','LineWidth',4);
%         hold on; 
%         h8=plot(SSw,mpr(:,2),'g:.','LineWidth',4);
%         hold on;  
%         h9=plot(SSw,mpr(:,1),'b:.','LineWidth',4);
%         hold on; 
%         h10=plot(SSw,mpv(:,3),'r:.','LineWidth',4);
%         hold on;
%         h11=plot(SSw,mpv(:,2),'g:.','LineWidth',4);
%         hold on;   
%         h12=plot(SSw,mpv(:,1),'b:.','LineWidth',4);
%         hold on;  
%         h1=plot(fractured_correct(:,1)*0.01,fractured_correct(:,8),'ws','markerfacecolor','r','markersize',20);
%         hold on;  
%         h2=plot(SSw,mp(:,3),'r','LineWidth',4);
%         hold on; 
%         h3=plot(fractured_correct(:,1)*0.01,fractured_correct(:,7),'wo','markerfacecolor','g','markersize',20);
%         hold on; 
%         h4=plot(SSw,mp(:,2),'g','LineWidth',4);
%         hold on; 
%         h5=plot(fractured_correct(:,1)*0.01,fractured_correct(:,6),'wd','markerfacecolor','b','markersize',20);
%         hold on;
%         h6=plot(SSw,mp(:,1),'b','LineWidth',4);
%         hold on;
%         set(gca,'FontSize',30);
%         set(gca,'xTick',[0 0.25 0.5 0.75 1]);
%         xlabel('Sw','fontsize',40);
%         ylabel('Mp(Pa)','fontsize',40);    
%         set(gca,'ylim',[14e9 25e9]);
% 
%         subplot(1,2,2)
%         set(gcf,'outerposition',get(0,'screensize'));
%         plot(SSw,vp_unique(:,3),'r','LineWidth',4);
%         hold on;  
%         plot(SSw,vp_unique(:,2),'g','LineWidth',4);
%         hold on;  
%         plot(SSw,vp_unique(:,1),'b','LineWidth',4);
%         hold on;  
%         h7=plot(SSw,vpr(:,3),'r:.','LineWidth',4);
%         hold on; 
%         h8=plot(SSw,vpr(:,2),'g:.','LineWidth',4);
%         hold on;  
%         h9=plot(SSw,vpr(:,1),'b:.','LineWidth',4);
%         hold on; 
%         h10=plot(SSw,vpv(:,3),'r:.','LineWidth',4);
%         hold on;
%         h11=plot(SSw,vpv(:,2),'g:.','LineWidth',4);
%         hold on;   
%         h12=plot(SSw,vpv(:,1),'b:.','LineWidth',4);
%         hold on;  
%         h1=plot(fractured_correct(:,1)*0.01,fractured_correct(:,5),'ws','markerfacecolor','r','markersize',20);
%         hold on;  
%         h2=plot(SSw,vp(:,3),'r','LineWidth',4);
%         hold on; 
%         h3=plot(fractured_correct(:,1)*0.01,fractured_correct(:,4),'wo','markerfacecolor','g','markersize',20);
%         hold on; 
%         h4=plot(SSw,vp(:,2),'g','LineWidth',4);
%         hold on; 
%         h5=plot(fractured_correct(:,1)*0.01,fractured_correct(:,3),'wd','markerfacecolor','b','markersize',20);
%         hold on;
%         h6=plot(SSw,vp(:,1),'b','LineWidth',4);
%         hold on;
%         set(gca,'FontSize',30);
%         set(gca,'yTick',[2700 2900 3100 3300 3500]);
%         set(gca,'xTick',[0 0.25 0.5 0.75 1]);
%         xlabel('Sw','fontsize',40);
%         ylabel('Vp(m/s)','fontsize',40);    
%         set(gca,'ylim',[2700 3500]);
        
    elseif Kdry_method==3  % Kdry from data at Sw=100 considering squirt process
%        
%         figure()
%         subplot(1,2,1)
%         set(gcf,'outerposition',get(0,'screensize'));
%           
% %         h7=plot(SSw,mpr(:,3),'r:.','LineWidth',4);
% %         hold on; 
% %         h8=plot(SSw,mpr(:,2),'g:.','LineWidth',4);
% %         hold on;  
% %         h9=plot(SSw,mpr(:,1),'b:.','LineWidth',4);
% %         hold on; 
% %         h10=plot(SSw,mpv(:,3),'r:.','LineWidth',4);
% %         hold on;
% %         h11=plot(SSw,mpv(:,2),'g:.','LineWidth',4);
% %         hold on;   
% %         h12=plot(SSw,mpv(:,1),'b:.','LineWidth',4);
% %         hold on;  
%         h1=plot(fractured_correct(:,1)*0.01,fractured_correct(:,8),'ws','markerfacecolor','r','markersize',20);
%         hold on;  
% %         h2=plot(SSw,mp(:,3),'r','LineWidth',4);
% %         hold on; 
%         h3=plot(fractured_correct(:,1)*0.01,fractured_correct(:,7),'wo','markerfacecolor','g','markersize',20);
%         hold on; 
% %         h4=plot(SSw,mp(:,2),'g','LineWidth',4);
% %         hold on; 
%         h5=plot(fractured_correct(:,1)*0.01,fractured_correct(:,6),'wd','markerfacecolor','b','markersize',20);
%         hold on;
%         
%         plot(SSw(nSw:nSw),mp_final(nSw:nSw,3),'or','markersize',30,'LineWidth',5);
%         hold on;  
%         plot(SSw(nSw:nSw),mp_final(nSw:nSw,2),'og','markersize',30,'LineWidth',5);
%         hold on;  
%         plot(SSw(nSw:nSw),mp_final(nSw:nSw,1),'ob','markersize',30,'LineWidth',5);
%         hold on;
% %         h6=plot(SSw,mp(:,1),'b','LineWidth',4);
% %         hold on;        
%         set(gca,'FontSize',30);
%         set(gca,'xTick',[0 0.25 0.5 0.75 1]);
%         xlabel('Sw','fontsize',40);
%         ylabel('Mp(Pa)','fontsize',40);    
%         set(gca,'ylim',[14e9 25e9]);
% 
%         
%         subplot(1,2,2)
%         set(gcf,'outerposition',get(0,'screensize'));
%        
% %         h7=plot(SSw,vpr(:,3),'r:.','LineWidth',4);
% %         hold on; 
% %         h8=plot(SSw,vpr(:,2),'g:.','LineWidth',4);
% %         hold on;  
% %         h9=plot(SSw,vpr(:,1),'b:.','LineWidth',4);
% %         hold on; 
% %         h10=plot(SSw,vpv(:,3),'r:.','LineWidth',4);
% %         hold on;
% %         h11=plot(SSw,vpv(:,2),'g:.','LineWidth',4);
% %         hold on;   
% %         h12=plot(SSw,vpv(:,1),'b:.','LineWidth',4);
% %         hold on;  
%         h1=plot(fractured_correct(:,1)*0.01,fractured_correct(:,5),'ws','markerfacecolor','r','markersize',20);
%         hold on;  
% %         h2=plot(SSw,vp(:,3),'r','LineWidth',4);
% %         hold on; 
%         h3=plot(fractured_correct(:,1)*0.01,fractured_correct(:,4),'wo','markerfacecolor','g','markersize',20);
%         hold on; 
% %         h4=plot(SSw,vp(:,2),'g','LineWidth',4);
% %         hold on; 
%         h5=plot(fractured_correct(:,1)*0.01,fractured_correct(:,3),'wd','markerfacecolor','b','markersize',20);
%         hold on;
% %         h6=plot(SSw,vp(:,1),'b','LineWidth',4);
% %         hold on;
%         plot(SSw(nSw:nSw),vp_final(nSw:nSw,3),'or','markersize',30,'LineWidth',5);
%         hold on;  
%         plot(SSw(nSw:nSw),vp_final(nSw:nSw,2),'og','markersize',30,'LineWidth',5);
%         hold on;  
%         plot(SSw(nSw:nSw),vp_final(nSw:nSw,1),'ob','markersize',30,'LineWidth',5);
%         hold on;
%         set(gca,'FontSize',30);
%         set(gca,'yTick',[2700 2900 3100 3300 3500]);
%         set(gca,'xTick',[0 0.25 0.5 0.75 1]);
%         xlabel('Sw','fontsize',40);
%         ylabel('Vp(m/s)','fontsize',40);    
%         set(gca,'ylim',[2700 3500]);
%         
%         figure()
%         subplot(1,2,1)
%         set(gcf,'outerposition',get(0,'screensize'));
%         plot(SSw(1:9),mp_final(1:9,3),'r','LineWidth',4);
%         hold on;  
%         plot(SSw(1:9),mp_final(1:9,2),'g','LineWidth',4);
%         hold on;  
%         plot(SSw(1:9),mp_final(1:9,1),'b','LineWidth',4);
%         hold on;  
% %         h7=plot(SSw,mpr(:,3),'r:.','LineWidth',4);
% %         hold on; 
% %         h8=plot(SSw,mpr(:,2),'g:.','LineWidth',4);
% %         hold on;  
% %         h9=plot(SSw,mpr(:,1),'b:.','LineWidth',4);
% %         hold on; 
% %         h10=plot(SSw,mpv(:,3),'r:.','LineWidth',4);
% %         hold on;
% %         h11=plot(SSw,mpv(:,2),'g:.','LineWidth',4);
% %         hold on;   
% %         h12=plot(SSw,mpv(:,1),'b:.','LineWidth',4);
% %         hold on;  
%         h1=plot(fractured_correct(:,1)*0.01,fractured_correct(:,8),'ws','markerfacecolor','r','markersize',20);
%         hold on;  
% %         h2=plot(SSw,mp(:,3),'r','LineWidth',4);
% %         hold on; 
%         h3=plot(fractured_correct(:,1)*0.01,fractured_correct(:,7),'wo','markerfacecolor','g','markersize',20);
%         hold on; 
% %         h4=plot(SSw,mp(:,2),'g','LineWidth',4);
% %         hold on; 
%         h5=plot(fractured_correct(:,1)*0.01,fractured_correct(:,6),'wd','markerfacecolor','b','markersize',20);
%         hold on;
% %         h6=plot(SSw,mp(:,1),'b','LineWidth',4);
% %         hold on;        
%         set(gca,'FontSize',30);
%         set(gca,'xTick',[0 0.25 0.5 0.75 1]);
%         xlabel('Sw','fontsize',40);
%         ylabel('Mp(Pa)','fontsize',40);    
%         set(gca,'ylim',[14e9 25e9]);
% 
%         
%         subplot(1,2,2)
%         set(gcf,'outerposition',get(0,'screensize'));
%         plot(SSw(1:9),vp_final(1:9,3),'r','LineWidth',4);
%         hold on;  
%         plot(SSw(1:9),vp_final(1:9,2),'g','LineWidth',4);
%         hold on;  
%         plot(SSw(1:9),vp_final(1:9,1),'b','LineWidth',4);
%         hold on;  
% %         h7=plot(SSw,vpr(:,3),'r:.','LineWidth',4);
% %         hold on; 
% %         h8=plot(SSw,vpr(:,2),'g:.','LineWidth',4);
% %         hold on;  
% %         h9=plot(SSw,vpr(:,1),'b:.','LineWidth',4);
% %         hold on; 
% %         h10=plot(SSw,vpv(:,3),'r:.','LineWidth',4);
% %         hold on;
% %         h11=plot(SSw,vpv(:,2),'g:.','LineWidth',4);
% %         hold on;   
% %         h12=plot(SSw,vpv(:,1),'b:.','LineWidth',4);
% %         hold on;  
%         h1=plot(fractured_correct(:,1)*0.01,fractured_correct(:,5),'ws','markerfacecolor','r','markersize',20);
%         hold on;  
% %         h2=plot(SSw,vp(:,3),'r','LineWidth',4);
% %         hold on; 
%         h3=plot(fractured_correct(:,1)*0.01,fractured_correct(:,4),'wo','markerfacecolor','g','markersize',20);
%         hold on; 
% %         h4=plot(SSw,vp(:,2),'g','LineWidth',4);
% %         hold on; 
%         h5=plot(fractured_correct(:,1)*0.01,fractured_correct(:,3),'wd','markerfacecolor','b','markersize',20);
%         hold on;
% %         h6=plot(SSw,vp(:,1),'b','LineWidth',4);
% %         hold on;
%         set(gca,'FontSize',30);
%         set(gca,'yTick',[2700 2900 3100 3300 3500]);
%         set(gca,'xTick',[0 0.25 0.5 0.75 1]);
%         xlabel('Sw','fontsize',40);
%         ylabel('Vp(m/s)','fontsize',40);    
%         set(gca,'ylim',[2700 3500]);
%         
%         figure()
%         subplot(1,2,1)
%         set(gcf,'outerposition',get(0,'screensize'));
%         plot(SSw(9:50),mp_final(9:50,3),'r','LineWidth',4);
%         hold on;  
%         plot(SSw(9:50),mp_final(9:50,2),'g','LineWidth',4);
%         hold on;  
%         plot(SSw(9:50),mp_final(9:50,1),'b','LineWidth',4);
%         hold on;  
% %         h7=plot(SSw,mpr(:,3),'r:.','LineWidth',4);
% %         hold on; 
% %         h8=plot(SSw,mpr(:,2),'g:.','LineWidth',4);
% %         hold on;  
% %         h9=plot(SSw,mpr(:,1),'b:.','LineWidth',4);
% %         hold on; 
% %         h10=plot(SSw,mpv(:,3),'r:.','LineWidth',4);
% %         hold on;
% %         h11=plot(SSw,mpv(:,2),'g:.','LineWidth',4);
% %         hold on;   
% %         h12=plot(SSw,mpv(:,1),'b:.','LineWidth',4);
% %         hold on;  
%         h1=plot(fractured_correct(:,1)*0.01,fractured_correct(:,8),'ws','markerfacecolor','r','markersize',20);
%         hold on;  
% %         h2=plot(SSw,mp(:,3),'r','LineWidth',4);
% %         hold on; 
%         h3=plot(fractured_correct(:,1)*0.01,fractured_correct(:,7),'wo','markerfacecolor','g','markersize',20);
%         hold on; 
% %         h4=plot(SSw,mp(:,2),'g','LineWidth',4);
% %         hold on; 
%         h5=plot(fractured_correct(:,1)*0.01,fractured_correct(:,6),'wd','markerfacecolor','b','markersize',20);
%         hold on;
% %         h6=plot(SSw,mp(:,1),'b','LineWidth',4);
% %         hold on;        
%         set(gca,'FontSize',30);
%         set(gca,'xTick',[0 0.25 0.5 0.75 1]);
%         xlabel('Sw','fontsize',40);
%         ylabel('Mp(Pa)','fontsize',40);    
%         set(gca,'ylim',[14e9 25e9]);
% 
%         
%         subplot(1,2,2)
%         set(gcf,'outerposition',get(0,'screensize'));
%         plot(SSw(9:50),vp_final(9:50,3),'r','LineWidth',4);
%         hold on;  
%         plot(SSw(9:50),vp_final(9:50,2),'g','LineWidth',4);
%         hold on;  
%         plot(SSw(9:50),vp_final(9:50,1),'b','LineWidth',4);
%         hold on;  
% %         h7=plot(SSw,vpr(:,3),'r:.','LineWidth',4);
% %         hold on; 
% %         h8=plot(SSw,vpr(:,2),'g:.','LineWidth',4);
% %         hold on;  
% %         h9=plot(SSw,vpr(:,1),'b:.','LineWidth',4);
% %         hold on; 
% %         h10=plot(SSw,vpv(:,3),'r:.','LineWidth',4);
% %         hold on;
% %         h11=plot(SSw,vpv(:,2),'g:.','LineWidth',4);
% %         hold on;   
% %         h12=plot(SSw,vpv(:,1),'b:.','LineWidth',4);
% %         hold on;  
%         h1=plot(fractured_correct(:,1)*0.01,fractured_correct(:,5),'ws','markerfacecolor','r','markersize',20);
%         hold on;  
% %         h2=plot(SSw,vp(:,3),'r','LineWidth',4);
% %         hold on; 
%         h3=plot(fractured_correct(:,1)*0.01,fractured_correct(:,4),'wo','markerfacecolor','g','markersize',20);
%         hold on; 
% %         h4=plot(SSw,vp(:,2),'g','LineWidth',4);
% %         hold on; 
%         h5=plot(fractured_correct(:,1)*0.01,fractured_correct(:,3),'wd','markerfacecolor','b','markersize',20);
%         hold on;
% %         h6=plot(SSw,vp(:,1),'b','LineWidth',4);
% %         hold on;
%         set(gca,'FontSize',30);
%         set(gca,'yTick',[2700 2900 3100 3300 3500]);
%         set(gca,'xTick',[0 0.25 0.5 0.75 1]);
%         xlabel('Sw','fontsize',40);
%         ylabel('Vp(m/s)','fontsize',40);    
%         set(gca,'ylim',[2700 3500]);
%         
%         figure()
%         subplot(1,2,1)
%         set(gcf,'outerposition',get(0,'screensize'));
%         plot(SSw(51:nSw-1),mp_final(51:nSw-1,3),'r','LineWidth',4);
%         hold on;  
%         plot(SSw(51:nSw-1),mp_final(51:nSw-1,2),'g','LineWidth',4);
%         hold on;  
%         plot(SSw(51:nSw-1),mp_final(51:nSw-1,1),'b','LineWidth',4);
%         hold on;  
% %         h7=plot(SSw,mpr(:,3),'r:.','LineWidth',4);
% %         hold on; 
% %         h8=plot(SSw,mpr(:,2),'g:.','LineWidth',4);
% %         hold on;  
% %         h9=plot(SSw,mpr(:,1),'b:.','LineWidth',4);
% %         hold on; 
% %         h10=plot(SSw,mpv(:,3),'r:.','LineWidth',4);
% %         hold on;
% %         h11=plot(SSw,mpv(:,2),'g:.','LineWidth',4);
% %         hold on;   
% %         h12=plot(SSw,mpv(:,1),'b:.','LineWidth',4);
% %         hold on;  
%         h1=plot(fractured_correct(:,1)*0.01,fractured_correct(:,8),'ws','markerfacecolor','r','markersize',20);
%         hold on;  
% %         h2=plot(SSw,mp(:,3),'r','LineWidth',4);
% %         hold on; 
%         h3=plot(fractured_correct(:,1)*0.01,fractured_correct(:,7),'wo','markerfacecolor','g','markersize',20);
%         hold on; 
% %         h4=plot(SSw,mp(:,2),'g','LineWidth',4);
% %         hold on; 
%         h5=plot(fractured_correct(:,1)*0.01,fractured_correct(:,6),'wd','markerfacecolor','b','markersize',20);
%         hold on;
% %         h6=plot(SSw,mp(:,1),'b','LineWidth',4);
% %         hold on;        
%         set(gca,'FontSize',30);
%         set(gca,'xTick',[0 0.25 0.5 0.75 1]);
%         xlabel('Sw','fontsize',40);
%         ylabel('Mp(Pa)','fontsize',40);    
%         set(gca,'ylim',[14e9 25e9]);
% 
%         
%         subplot(1,2,2)
%         set(gcf,'outerposition',get(0,'screensize'));
%         plot(SSw(51:nSw-1),vp_final(51:nSw-1,3),'r','LineWidth',4);
%         hold on;  
%         plot(SSw(51:nSw-1),vp_final(51:nSw-1,2),'g','LineWidth',4);
%         hold on;  
%         plot(SSw(51:nSw-1),vp_final(51:nSw-1,1),'b','LineWidth',4);
%         hold on;  
% %         h7=plot(SSw,vpr(:,3),'r:.','LineWidth',4);
% %         hold on; 
% %         h8=plot(SSw,vpr(:,2),'g:.','LineWidth',4);
% %         hold on;  
% %         h9=plot(SSw,vpr(:,1),'b:.','LineWidth',4);
% %         hold on; 
% %         h10=plot(SSw,vpv(:,3),'r:.','LineWidth',4);
% %         hold on;
% %         h11=plot(SSw,vpv(:,2),'g:.','LineWidth',4);
% %         hold on;   
% %         h12=plot(SSw,vpv(:,1),'b:.','LineWidth',4);
% %         hold on;  
%         h1=plot(fractured_correct(:,1)*0.01,fractured_correct(:,5),'ws','markerfacecolor','r','markersize',20);
%         hold on;  
% %         h2=plot(SSw,vp(:,3),'r','LineWidth',4);
% %         hold on; 
%         h3=plot(fractured_correct(:,1)*0.01,fractured_correct(:,4),'wo','markerfacecolor','g','markersize',20);
%         hold on; 
% %         h4=plot(SSw,vp(:,2),'g','LineWidth',4);
% %         hold on; 
%         h5=plot(fractured_correct(:,1)*0.01,fractured_correct(:,3),'wd','markerfacecolor','b','markersize',20);
%         hold on;
% %         h6=plot(SSw,vp(:,1),'b','LineWidth',4);
% %         hold on;
%         set(gca,'FontSize',30);
%         set(gca,'yTick',[2700 2900 3100 3300 3500]);
%         set(gca,'xTick',[0 0.25 0.5 0.75 1]);
%         xlabel('Sw','fontsize',40);
%         ylabel('Vp(m/s)','fontsize',40);    
%         set(gca,'ylim',[2700 3500]);
        
        figure()
%         subplot(2,1,1)
        set(gcf,'outerposition',get(0,'screensize'));
        h1=plot(fractured_correct(:,1)*0.01,fractured_correct(:,8),'ws','markerfacecolor','r','markersize',20);
        hold on;          
        h2=plot(fractured_correct(:,1)*0.01,fractured_correct(:,7),'wo','markerfacecolor','g','markersize',20);
        hold on;         
        h3=plot(fractured_correct(:,1)*0.01,fractured_correct(:,6),'wd','markerfacecolor','b','markersize',20);
        hold on;
        h4=plot(SSw(1:20),mp_final(1:20,3),'r','LineWidth',4);
        hold on;
        h5=plot(SSw(1:20),mp_final(1:20,2),'g','LineWidth',4);
        hold on; 
        h6=plot(SSw(1:20),mp_final(1:20,1),'b','LineWidth',4);
        hold on;
        h7=plot(SSw(20:50),mp_final(20:50,3),'r--','LineWidth',4);
        hold on;
        h8=plot(SSw(20:50),mp_final(20:50,2),'g--','LineWidth',4);
        hold on; 
        h9=plot(SSw(20:50),mp_final(20:50,1),'b--','LineWidth',4);
        hold on;
        h10=plot(SSw(51:nSw-1),mp_final(51:nSw-1,3),'r.','markersize',25,'LineWidth',3);
        hold on;
        h11=plot(SSw(51:nSw-1),mp_final(51:nSw-1,2),'g.','markersize',25,'LineWidth',3);
        hold on; 
        h12=plot(SSw(51:nSw-1),mp_final(51:nSw-1,1),'b.','markersize',25,'LineWidth',3);
        hold on;
        h13=plot(SSw(nSw:nSw),mp_final(nSw:nSw,3),'or','markersize',25,'LineWidth',3);
        hold on;
        h14=plot(SSw(nSw:nSw),mp_final(nSw:nSw,2),'og','markersize',25,'LineWidth',3);
        hold on; 
        h15=plot(SSw(nSw:nSw),mp_final(nSw:nSw,1),'ob','markersize',25,'LineWidth',3);   
        hold on; 
        h16=plot(SSw,mp_uniform(:,3),'r');  
        hold on; 
        h17=plot(SSw,mp_uniform(:,2),'g');  
        hold on; 
        h18=plot(SSw,mp_uniform(:,1),'b');  
        
        set(gca,'FontSize',30);
        set(gca,'FontSize',30);
        set(gca,'xTick',[0 0.25 0.5 0.75 1]);
        xlabel('S_w','fontsize',40);
        ylabel('M_p(Pa)','fontsize',40);    
        set(gca,'ylim',[14e9 25e9]);
        lgd1=legend([h1,h2,h3],'Data-90бу','Data-45бу','Data-0бу','orientation','vertical','location','south');
        set(lgd1,'FontName','Time New Roman','FontSize',30);
        legend boxoff;
        ah1=axes('position',get(gca,'position'),'visible','off');
        lgd2=legend(ah1,[h4,h5,h6],'Squirt-90бу','Squirt-45бу','Squirt-0бу','orientation','vertical','location','south');
        set(lgd2,'FontName','Time New Roman','FontSize',30);
        legend boxoff;
        ah2=axes('position',get(gca,'position'),'visible','off');
        lgd3=legend(ah2,[h7,h8,h9],'Uniform-90бу','Uniform-45бу','Uniform-0бу','orientation','vertical','location','south');
        set(lgd3,'FontName','Time New Roman','FontSize',30);
        legend boxoff;
        ah3=axes('position',get(gca,'position'),'visible','off');
        lgd4=legend(ah3,[h10,h11,h12],'Patchy-90бу','Patchy-45бу','Patchy-0бу','orientation','vertical','location','south');
        set(lgd4,'FontName','Time New Roman','FontSize',30);
        legend boxoff;
        ah4=axes('position',get(gca,'position'),'visible','off');
        lgd5=legend(ah4,[h13,h14,h15],'Sat-90бу','Sat-45бу','Sat-0бу','orientation','vertical','location','south');
        set(lgd5,'FontName','Time New Roman','FontSize',30);
        legend boxoff;
        
%         legend('Data-90бу','Data-45бу','Data-0бу','Imbibition-90бу','Imbibition-45бу','Imbibition-0бу','Uniform-90бу','Uniform-45бу','Uniform-0бу','Patchy-90бу','Patchy-45бу','Patchy-0бу','Saturated-90бу','Saturated-45бу','Saturated-0бу');

%         subplot(2,1,2)
 figure()
        set(gcf,'outerposition',get(0,'screensize'));
        h1=plot(fractured_correct(:,1)*0.01,fractured_correct(:,5),'ws','markerfacecolor','r','markersize',20);
        hold on;         
        h2=plot(fractured_correct(:,1)*0.01,fractured_correct(:,4),'wo','markerfacecolor','g','markersize',20);
        hold on;       
        h3=plot(fractured_correct(:,1)*0.01,fractured_correct(:,3),'wd','markerfacecolor','b','markersize',20);
        hold on;
        
        h4=plot(SSw(1:20),vp_final(1:20,3),'r','LineWidth',4);
        hold on;
        h5=plot(SSw(1:20),vp_final(1:20,2),'g','LineWidth',4);
        hold on; 
        h6=plot(SSw(1:20),vp_final(1:20,1),'b','LineWidth',4);
        hold on;
        h7=plot(SSw(20:50),vp_final(20:50,3),'r--','LineWidth',4);
        hold on;
        h8=plot(SSw(20:50),vp_final(20:50,2),'g--','LineWidth',4);
        hold on; 
        h9=plot(SSw(20:50),vp_final(20:50,1),'b--','LineWidth',4);
        hold on;
        h10=plot(SSw(51:nSw-1),vp_final(51:nSw-1,3),'r.','markersize',25,'LineWidth',3);
        hold on;
        h11=plot(SSw(51:nSw-1),vp_final(51:nSw-1,2),'g.','markersize',25,'LineWidth',3);
        hold on; 
        h12=plot(SSw(51:nSw-1),vp_final(51:nSw-1,1),'b.','markersize',25,'LineWidth',3);
        hold on;
        h13=plot(SSw(nSw:nSw),vp_final(nSw:nSw,3),'or','markersize',25,'LineWidth',3);
        hold on;
        h14=plot(SSw(nSw:nSw),vp_final(nSw:nSw,2),'og','markersize',25,'LineWidth',3);
        hold on; 
        h15=plot(SSw(nSw:nSw),vp_final(nSw:nSw,1),'ob','markersize',25,'LineWidth',3);  
        hold on;
        hold on; 
        h16=plot(SSw,vp_uniform(:,3),'r');  
        hold on; 
        h17=plot(SSw,vp_uniform(:,2),'g');  
        hold on; 
        h18=plot(SSw,vp_uniform(:,1),'b');  
        set(gca,'FontSize',30);
        set(gca,'yTick',[2700 2900 3100 3300 3500]);
        set(gca,'xTick',[0 0.25 0.5 0.75 1]);
        xlabel('S_w','fontsize',40);
        ylabel('V_p(m/s)','fontsize',40);    
        set(gca,'ylim',[2700 3500]);
    end    
end

