clc;
clear;

blank=load('blank.txt');

pi=3.1415926;
flag_method=2;
nSw=101;
nInc=3;

%%%%%%%%%% correct anisotropic coefficients from blank data %%%%%%%%%%%%%%%%%%%%%%
mult0 = 1.0347;% average(blank(:,8)/blank(:,6))
mult45 = 1.0088;% average(blank(:,8)/blank(:,7))
mult0V = 1.0172;% average(blank(:,5)/blank(:,3))
mult45V = 1.0044;% average(blank(:,5)/blank(:,4))

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% anisotropy adjustment %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% shift all valuses in certain lines including Sw=100% %%%%%%%%%%%%%
blank_correct=blank;
blank_correct(:,3)=blank(:,3)*mult0V;
blank_correct(:,4)=blank(:,4)*mult45V;
blank_correct(:,6)=blank(:,6)*mult0;
blank_correct(:,7)=blank(:,7)*mult45;

%%%%%%%%%%% plot data %%%%%%%%%%%%%%%%%%%%%%
figure() %%%% original data  %%%%%%%%%%%%
subplot(1,2,1)
set(gcf,'outerposition',get(0,'screensize'));
h1=plot(blank_correct(:,1)*0.01,blank(:,8),'ws','markerfacecolor','r','markersize',20);
hold on;  
h2=plot(blank_correct(:,1)*0.01,blank(:,7),'wo','markerfacecolor','g','markersize',20);
hold on; 
h3=plot(blank_correct(:,1)*0.01,blank(:,6),'wd','markerfacecolor','b','markersize',20);
hold on;
set(gca,'FontSize',30);
set(gca,'xTick',[0 0.25 0.5 0.75 1]);
xlabel('S_w','fontsize',40);
ylabel('M_p(Pa)','fontsize',40);    
set(gca,'ylim',[14e9 25e9]);
lgd1=legend([h1,h2,h3],'Data-90бу','Data-45бу','Data-0бу','orientation','vertical','location','northwest');
set(lgd1,'FontName','Time New Roman','FontSize',30);
% legend boxoff;


subplot(1,2,2)
set(gcf,'outerposition',get(0,'screensize'));
h1=plot(blank_correct(:,1)*0.01,blank(:,5),'ws','markerfacecolor','r','markersize',20);
hold on;  
h2=plot(blank_correct(:,1)*0.01,blank(:,4),'wo','markerfacecolor','g','markersize',20);
hold on; 
h3=plot(blank_correct(:,1)*0.01,blank(:,3),'wd','markerfacecolor','b','markersize',20);
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
h1=plot(blank_correct(:,1)*0.01,blank_correct(:,8),'ws','markerfacecolor','r','markersize',20);
hold on;  
h2=plot(blank_correct(:,1)*0.01,blank_correct(:,7),'wo','markerfacecolor','g','markersize',20);
hold on; 
h3=plot(blank_correct(:,1)*0.01,blank_correct(:,6),'wd','markerfacecolor','b','markersize',20);
hold on;
set(gca,'FontSize',30);
set(gca,'xTick',[0 0.25 0.5 0.75 1]);
xlabel('S_w','fontsize',40);
ylabel('M_p(Pa)','fontsize',40);    
set(gca,'ylim',[14e9 25e9]);
lgd1=legend([h1,h2,h3],'Data-90бу','Data-45бу','Data-0бу','orientation','vertical','location','northwest');
set(lgd1,'FontName','Time New Roman','FontSize',30);
% legend boxoff;

subplot(1,2,2)
set(gcf,'outerposition',get(0,'screensize'));
h1=plot(blank_correct(:,1)*0.01,blank_correct(:,5),'ws','markerfacecolor','r','markersize',20);
hold on;  
h2=plot(blank_correct(:,1)*0.01,blank_correct(:,4),'wo','markerfacecolor','g','markersize',20);
hold on; 
h3=plot(blank_correct(:,1)*0.01,blank_correct(:,3),'wd','markerfacecolor','b','markersize',20);
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

Por=0.3;
Por_t = 0.3;% data from Kelvin et al's paper (2015)   
Por_p = 21e-15;% data from Kelvin et al's paper (2015)        

H=6e-3; % fracture period from Tillotson et al, figure 2, page 1240
phi_c=0.001;% microcrack porosity given by Boris
Por_re = Por_t-phi_c;% relax porosity for fractured sample

Den_gw=Den_g*(1-Por_t)+Por_t*Dwater;% density for water saturated sample
Den_gg=Den_g*(1-Por_t)+Por_t*Dgas;% density for air saturated sample 

V90s1=2226;% S1 wave velocity at 90 degrees for fractured sample with Sw=0
Udry=Den_gg*V90s1*V90s1;% dry shear modulus for blank sample 

%% %%%%%%%%%%%%%  computing codes for three Kdry methods %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for Kdry_method=1:3
    
    if Kdry_method==1 %%%%%%%%%%%%% computing Kdry from values at Sw=0
        
        V90pp=3242;% P-wave velocity at 90 degrees for blank sample with Sw=0
        Mdry=Den_gg*V90pp*V90pp;% dry P wave modulus for blank sample 
        Kdry=Mdry-4./3*Udry% dry bulk modulus for blank sample    
        
    elseif Kdry_method==2  %%%%%%%%%%%%% computing Kdry from values at Sw=100
        
        Msat=blank(10,8);        
        Kdry =(3*Kg^2*Kwater - 3*Kg*Kwater*Msat + 4*Kg*Kwater*Udry - 3*Kg^2*Msat*Por + 4*Kg^2*Por*Udry + 3*Kg*Kwater*Msat*Por - 4*Kg*Kwater*Por*Udry)...
                    /(4*Kwater*Udry - 3*Kg^2*Por + 3*Kg*Kwater - 3*Kwater*Msat + 3*Kg*Kwater*Por);
                
    elseif Kdry_method==3 %%%%%%%%%%%%% computing Kdry from values at Sw=100 considering the unrelaxation at the beginning        
        %%%%%%%%%%%%%        (1)  Msat=Ksat+4./3*Usat;   %% Msat is known from data , fractured.txt     %%%%%%%%%%%
        %%%%%%%%%%%%%        (2) Ksat/(Kg-Ksat)=Kwater/(Por_t*(Kg-Kwater))+Kun/(Kg-Kun)  
        %%%%%%%%%%%%             Usat=Uun                                               %% Gassmann equation %%%%%%%%%%%
        %%%%%%%%%%%% (3) 1./Uun = 1./U0 - 4/15 * (1./K0 - 1/Kun)  %% equation 2 in Gurevich et al's paper(2009) %%%%%%%%%%
        %%%%%%%%%%%% Combining (1)(2)and (3), I have the expressions of Kun and Gun as follows %%%

        U0=7.6e9;% table from Boris for matrix when Sw=0
        K0=8.9e9;% table from Boris for matrix when Sw=0
        Msat=blank(10,8)   %% Msat is known from data when Sw=100% at 90 incidence angle, blank.txt     %%%%%%%%%%%

        Kun1=((225*K0^2*Kg^4*Kwater^2 - 450*K0^2*Kg^4*Kwater*Msat*Por_t + 480*K0^2*Kg^4*Kwater*Por_t*U0 + 225*K0^2*Kg^4*Msat^2*Por_t^2 - 480*K0^2*Kg^4*Msat*Por_t^2*U0 ...
            + 576*K0^2*Kg^4*Por_t^2*U0^2 + 450*K0^2*Kg^3*Kwater^2*Msat*Por_t - 450*K0^2*Kg^3*Kwater^2*Msat - 480*K0^2*Kg^3*Kwater^2*Por_t*U0 + 720*K0^2*Kg^3*Kwater^2*U0 ...
            - 450*K0^2*Kg^3*Kwater*Msat^2*Por_t^2 + 450*K0^2*Kg^3*Kwater*Msat^2*Por_t + 960*K0^2*Kg^3*Kwater*Msat*Por_t^2*U0 - 1200*K0^2*Kg^3*Kwater*Msat*Por_t*U0 ...
            - 1152*K0^2*Kg^3*Kwater*Por_t^2*U0^2 + 768*K0^2*Kg^3*Kwater*Por_t*U0^2 + 225*K0^2*Kg^2*Kwater^2*Msat^2*Por_t^2 - 450*K0^2*Kg^2*Kwater^2*Msat^2*Por_t ...
            + 225*K0^2*Kg^2*Kwater^2*Msat^2 - 480*K0^2*Kg^2*Kwater^2*Msat*Por_t^2*U0 + 1200*K0^2*Kg^2*Kwater^2*Msat*Por_t*U0 - 840*K0^2*Kg^2*Kwater^2*Msat*U0 ...
            + 576*K0^2*Kg^2*Kwater^2*Por_t^2*U0^2 - 768*K0^2*Kg^2*Kwater^2*Por_t*U0^2 + 576*K0^2*Kg^2*Kwater^2*U0^2 + 120*K0^2*Kg^2*Kwater*Msat^2*Por_t*U0 ...
            - 128*K0^2*Kg^2*Kwater*Msat*Por_t*U0^2 - 120*K0^2*Kg*Kwater^2*Msat^2*Por_t*U0 + 120*K0^2*Kg*Kwater^2*Msat^2*U0 + 128*K0^2*Kg*Kwater^2*Msat*Por_t*U0^2 ...
            - 192*K0^2*Kg*Kwater^2*Msat*U0^2 + 16*K0^2*Kwater^2*Msat^2*U0^2 - 120*K0*Kg^4*Kwater^2*U0 + 240*K0*Kg^4*Kwater*Msat*Por_t*U0 - 128*K0*Kg^4*Kwater*Por_t*U0^2 ...
            - 120*K0*Kg^4*Msat^2*Por_t^2*U0 + 128*K0*Kg^4*Msat*Por_t^2*U0^2 - 240*K0*Kg^3*Kwater^2*Msat*Por_t*U0 + 240*K0*Kg^3*Kwater^2*Msat*U0 + 128*K0*Kg^3*Kwater^2*Por_t*U0^2 ...
            - 192*K0*Kg^3*Kwater^2*U0^2 + 240*K0*Kg^3*Kwater*Msat^2*Por_t^2*U0 - 240*K0*Kg^3*Kwater*Msat^2*Por_t*U0 - 256*K0*Kg^3*Kwater*Msat*Por_t^2*U0^2 ...
            + 320*K0*Kg^3*Kwater*Msat*Por_t*U0^2 - 120*K0*Kg^2*Kwater^2*Msat^2*Por_t^2*U0 + 240*K0*Kg^2*Kwater^2*Msat^2*Por_t*U0 - 120*K0*Kg^2*Kwater^2*Msat^2*U0 ...
            + 128*K0*Kg^2*Kwater^2*Msat*Por_t^2*U0^2 - 320*K0*Kg^2*Kwater^2*Msat*Por_t*U0^2 + 224*K0*Kg^2*Kwater^2*Msat*U0^2 - 32*K0*Kg^2*Kwater*Msat^2*Por_t*U0^2 ...
            + 32*K0*Kg*Kwater^2*Msat^2*Por_t*U0^2 - 32*K0*Kg*Kwater^2*Msat^2*U0^2 + 16*Kg^4*Kwater^2*U0^2 - 32*Kg^4*Kwater*Msat*Por_t*U0^2 + 16*Kg^4*Msat^2*Por_t^2*U0^2 ...
            + 32*Kg^3*Kwater^2*Msat*Por_t*U0^2 - 32*Kg^3*Kwater^2*Msat*U0^2 - 32*Kg^3*Kwater*Msat^2*Por_t^2*U0^2 + 32*Kg^3*Kwater*Msat^2*Por_t*U0^2 + 16*Kg^2*Kwater^2*Msat^2*Por_t^2*U0^2 ...
            - 32*Kg^2*Kwater^2*Msat^2*Por_t*U0^2 + 16*Kg^2*Kwater^2*Msat^2*U0^2)^(1/2) + 15*K0*Kg^2*Kwater - 4*Kg^2*Kwater*U0 - 15*K0*Kg*Kwater*Msat + 16*K0*Kg*Kwater*U0 + 4*K0*Kwater*Msat*U0 ...
            + 4*Kg*Kwater*Msat*U0 - 15*K0*Kg^2*Msat*Por_t + 24*K0*Kg^2*Por_t*U0 + 4*Kg^2*Msat*Por_t*U0 + 15*K0*Kg*Kwater*Msat*Por_t - 24*K0*Kg*Kwater*Por_t*U0 - 4*Kg*Kwater*Msat*Por_t*U0)...
            /(2*(15*K0*Kg*Kwater - 15*K0*Kwater*Msat + 20*K0*Kwater*U0 - 4*Kg*Kwater*U0 + 4*Kwater*Msat*U0 - 15*K0*Kg^2*Por_t + 4*Kg^2*Por_t*U0 + 15*K0*Kg*Kwater*Por_t - 4*Kg*Kwater*Por_t*U0));

        Kun2=(15*K0*Kg^2*Kwater - (225*K0^2*Kg^4*Kwater^2 - 450*K0^2*Kg^4*Kwater*Msat*Por_t + 480*K0^2*Kg^4*Kwater*Por_t*U0 + 225*K0^2*Kg^4*Msat^2*Por_t^2 - 480*K0^2*Kg^4*Msat*Por_t^2*U0 ...
            + 576*K0^2*Kg^4*Por_t^2*U0^2 + 450*K0^2*Kg^3*Kwater^2*Msat*Por_t - 450*K0^2*Kg^3*Kwater^2*Msat - 480*K0^2*Kg^3*Kwater^2*Por_t*U0 + 720*K0^2*Kg^3*Kwater^2*U0 ...
            - 450*K0^2*Kg^3*Kwater*Msat^2*Por_t^2 + 450*K0^2*Kg^3*Kwater*Msat^2*Por_t + 960*K0^2*Kg^3*Kwater*Msat*Por_t^2*U0 - 1200*K0^2*Kg^3*Kwater*Msat*Por_t*U0 ...
            - 1152*K0^2*Kg^3*Kwater*Por_t^2*U0^2 + 768*K0^2*Kg^3*Kwater*Por_t*U0^2 + 225*K0^2*Kg^2*Kwater^2*Msat^2*Por_t^2 - 450*K0^2*Kg^2*Kwater^2*Msat^2*Por_t ...
            + 225*K0^2*Kg^2*Kwater^2*Msat^2 - 480*K0^2*Kg^2*Kwater^2*Msat*Por_t^2*U0 + 1200*K0^2*Kg^2*Kwater^2*Msat*Por_t*U0 - 840*K0^2*Kg^2*Kwater^2*Msat*U0 ...
            + 576*K0^2*Kg^2*Kwater^2*Por_t^2*U0^2 - 768*K0^2*Kg^2*Kwater^2*Por_t*U0^2 + 576*K0^2*Kg^2*Kwater^2*U0^2 + 120*K0^2*Kg^2*Kwater*Msat^2*Por_t*U0 ...
            - 128*K0^2*Kg^2*Kwater*Msat*Por_t*U0^2 - 120*K0^2*Kg*Kwater^2*Msat^2*Por_t*U0 + 120*K0^2*Kg*Kwater^2*Msat^2*U0 + 128*K0^2*Kg*Kwater^2*Msat*Por_t*U0^2 ...
            - 192*K0^2*Kg*Kwater^2*Msat*U0^2 + 16*K0^2*Kwater^2*Msat^2*U0^2 - 120*K0*Kg^4*Kwater^2*U0 + 240*K0*Kg^4*Kwater*Msat*Por_t*U0 - 128*K0*Kg^4*Kwater*Por_t*U0^2 ...
            - 120*K0*Kg^4*Msat^2*Por_t^2*U0 + 128*K0*Kg^4*Msat*Por_t^2*U0^2 - 240*K0*Kg^3*Kwater^2*Msat*Por_t*U0 + 240*K0*Kg^3*Kwater^2*Msat*U0 + 128*K0*Kg^3*Kwater^2*Por_t*U0^2 ...
            - 192*K0*Kg^3*Kwater^2*U0^2 + 240*K0*Kg^3*Kwater*Msat^2*Por_t^2*U0 - 240*K0*Kg^3*Kwater*Msat^2*Por_t*U0 - 256*K0*Kg^3*Kwater*Msat*Por_t^2*U0^2 + 320*K0*Kg^3*Kwater*Msat*Por_t*U0^2 ...
            - 120*K0*Kg^2*Kwater^2*Msat^2*Por_t^2*U0 + 240*K0*Kg^2*Kwater^2*Msat^2*Por_t*U0 - 120*K0*Kg^2*Kwater^2*Msat^2*U0 + 128*K0*Kg^2*Kwater^2*Msat*Por_t^2*U0^2 ...
            - 320*K0*Kg^2*Kwater^2*Msat*Por_t*U0^2 + 224*K0*Kg^2*Kwater^2*Msat*U0^2 - 32*K0*Kg^2*Kwater*Msat^2*Por_t*U0^2 + 32*K0*Kg*Kwater^2*Msat^2*Por_t*U0^2 ...
            - 32*K0*Kg*Kwater^2*Msat^2*U0^2 + 16*Kg^4*Kwater^2*U0^2 - 32*Kg^4*Kwater*Msat*Por_t*U0^2 + 16*Kg^4*Msat^2*Por_t^2*U0^2 + 32*Kg^3*Kwater^2*Msat*Por_t*U0^2 ...
            - 32*Kg^3*Kwater^2*Msat*U0^2 - 32*Kg^3*Kwater*Msat^2*Por_t^2*U0^2 + 32*Kg^3*Kwater*Msat^2*Por_t*U0^2 + 16*Kg^2*Kwater^2*Msat^2*Por_t^2*U0^2 - 32*Kg^2*Kwater^2*Msat^2*Por_t*U0^2 ...
            + 16*Kg^2*Kwater^2*Msat^2*U0^2)^(1/2) - 4*Kg^2*Kwater*U0 - 15*K0*Kg*Kwater*Msat + 16*K0*Kg*Kwater*U0 + 4*K0*Kwater*Msat*U0 + 4*Kg*Kwater*Msat*U0 - 15*K0*Kg^2*Msat*Por_t ...
            + 24*K0*Kg^2*Por_t*U0 + 4*Kg^2*Msat*Por_t*U0 + 15*K0*Kg*Kwater*Msat*Por_t - 24*K0*Kg*Kwater*Por_t*U0 - 4*Kg*Kwater*Msat*Por_t*U0)/(2*(15*K0*Kg*Kwater - 15*K0*Kwater*Msat ...
            + 20*K0*Kwater*U0 - 4*Kg*Kwater*U0 + 4*Kwater*Msat*U0 - 15*K0*Kg^2*Por_t + 4*Kg^2*Por_t*U0 + 15*K0*Kg*Kwater*Por_t - 4*Kg*Kwater*Por_t*U0));


        Uun1= (3*((225*K0^2*Kg^4*Kwater^2 - 450*K0^2*Kg^4*Kwater*Msat*Por_t + 480*K0^2*Kg^4*Kwater*Por_t*U0 + 225*K0^2*Kg^4*Msat^2*Por_t^2 - 480*K0^2*Kg^4*Msat*Por_t^2*U0 ...
            + 576*K0^2*Kg^4*Por_t^2*U0^2 + 450*K0^2*Kg^3*Kwater^2*Msat*Por_t - 450*K0^2*Kg^3*Kwater^2*Msat - 480*K0^2*Kg^3*Kwater^2*Por_t*U0 + 720*K0^2*Kg^3*Kwater^2*U0 ...
            - 450*K0^2*Kg^3*Kwater*Msat^2*Por_t^2 + 450*K0^2*Kg^3*Kwater*Msat^2*Por_t + 960*K0^2*Kg^3*Kwater*Msat*Por_t^2*U0 - 1200*K0^2*Kg^3*Kwater*Msat*Por_t*U0 ...
            - 1152*K0^2*Kg^3*Kwater*Por_t^2*U0^2 + 768*K0^2*Kg^3*Kwater*Por_t*U0^2 + 225*K0^2*Kg^2*Kwater^2*Msat^2*Por_t^2 - 450*K0^2*Kg^2*Kwater^2*Msat^2*Por_t ...
            + 225*K0^2*Kg^2*Kwater^2*Msat^2 - 480*K0^2*Kg^2*Kwater^2*Msat*Por_t^2*U0 + 1200*K0^2*Kg^2*Kwater^2*Msat*Por_t*U0 - 840*K0^2*Kg^2*Kwater^2*Msat*U0 ...
            + 576*K0^2*Kg^2*Kwater^2*Por_t^2*U0^2 - 768*K0^2*Kg^2*Kwater^2*Por_t*U0^2 + 576*K0^2*Kg^2*Kwater^2*U0^2 + 120*K0^2*Kg^2*Kwater*Msat^2*Por_t*U0 ...
            - 128*K0^2*Kg^2*Kwater*Msat*Por_t*U0^2 - 120*K0^2*Kg*Kwater^2*Msat^2*Por_t*U0 + 120*K0^2*Kg*Kwater^2*Msat^2*U0 + 128*K0^2*Kg*Kwater^2*Msat*Por_t*U0^2 ...
            - 192*K0^2*Kg*Kwater^2*Msat*U0^2 + 16*K0^2*Kwater^2*Msat^2*U0^2 - 120*K0*Kg^4*Kwater^2*U0 + 240*K0*Kg^4*Kwater*Msat*Por_t*U0 - 128*K0*Kg^4*Kwater*Por_t*U0^2 ...
            - 120*K0*Kg^4*Msat^2*Por_t^2*U0 + 128*K0*Kg^4*Msat*Por_t^2*U0^2 - 240*K0*Kg^3*Kwater^2*Msat*Por_t*U0 + 240*K0*Kg^3*Kwater^2*Msat*U0 + 128*K0*Kg^3*Kwater^2*Por_t*U0^2 ...
            - 192*K0*Kg^3*Kwater^2*U0^2 + 240*K0*Kg^3*Kwater*Msat^2*Por_t^2*U0 - 240*K0*Kg^3*Kwater*Msat^2*Por_t*U0 - 256*K0*Kg^3*Kwater*Msat*Por_t^2*U0^2 + 320*K0*Kg^3*Kwater*Msat*Por_t*U0^2 ...
            - 120*K0*Kg^2*Kwater^2*Msat^2*Por_t^2*U0 + 240*K0*Kg^2*Kwater^2*Msat^2*Por_t*U0 - 120*K0*Kg^2*Kwater^2*Msat^2*U0 + 128*K0*Kg^2*Kwater^2*Msat*Por_t^2*U0^2 ...
            - 320*K0*Kg^2*Kwater^2*Msat*Por_t*U0^2 + 224*K0*Kg^2*Kwater^2*Msat*U0^2 - 32*K0*Kg^2*Kwater*Msat^2*Por_t*U0^2 + 32*K0*Kg*Kwater^2*Msat^2*Por_t*U0^2 - 32*K0*Kg*Kwater^2*Msat^2*U0^2 ...
            + 16*Kg^4*Kwater^2*U0^2 - 32*Kg^4*Kwater*Msat*Por_t*U0^2 + 16*Kg^4*Msat^2*Por_t^2*U0^2 + 32*Kg^3*Kwater^2*Msat*Por_t*U0^2 - 32*Kg^3*Kwater^2*Msat*U0^2 - 32*Kg^3*Kwater*Msat^2*Por_t^2*U0^2 ...
            + 32*Kg^3*Kwater*Msat^2*Por_t*U0^2 + 16*Kg^2*Kwater^2*Msat^2*Por_t^2*U0^2 - 32*Kg^2*Kwater^2*Msat^2*Por_t*U0^2 + 16*Kg^2*Kwater^2*Msat^2*U0^2)^(1/2) - 15*K0*Kg^2*Kwater + 4*Kg^2*Kwater*U0 ...
            + 15*K0*Kg*Kwater*Msat + 16*K0*Kg*Kwater*U0 + 4*K0*Kwater*Msat*U0 - 4*Kg*Kwater*Msat*U0 + 15*K0*Kg^2*Msat*Por_t + 24*K0*Kg^2*Por_t*U0 - 4*Kg^2*Msat*Por_t*U0 - 15*K0*Kg*Kwater*Msat*Por_t ...
            - 24*K0*Kg*Kwater*Por_t*U0 + 4*Kg*Kwater*Msat*Por_t*U0))/(8*(15*K0*Kg*Kwater + 4*K0*Kwater*U0 - 4*Kg*Kwater*U0 + 15*K0*Kg^2*Por_t - 4*Kg^2*Por_t*U0 - 15*K0*Kg*Kwater*Por_t + 4*Kg*Kwater*Por_t*U0));

         Uun2=(3*(4*Kg^2*Kwater*U0 - 15*K0*Kg^2*Kwater - (225*K0^2*Kg^4*Kwater^2 - 450*K0^2*Kg^4*Kwater*Msat*Por_t + 480*K0^2*Kg^4*Kwater*Por_t*U0 + 225*K0^2*Kg^4*Msat^2*Por_t^2 - 480*K0^2*Kg^4*Msat*Por_t^2*U0 ...
             + 576*K0^2*Kg^4*Por_t^2*U0^2 + 450*K0^2*Kg^3*Kwater^2*Msat*Por_t - 450*K0^2*Kg^3*Kwater^2*Msat - 480*K0^2*Kg^3*Kwater^2*Por_t*U0 + 720*K0^2*Kg^3*Kwater^2*U0 - 450*K0^2*Kg^3*Kwater*Msat^2*Por_t^2 ...
             + 450*K0^2*Kg^3*Kwater*Msat^2*Por_t + 960*K0^2*Kg^3*Kwater*Msat*Por_t^2*U0 - 1200*K0^2*Kg^3*Kwater*Msat*Por_t*U0 - 1152*K0^2*Kg^3*Kwater*Por_t^2*U0^2 + 768*K0^2*Kg^3*Kwater*Por_t*U0^2 ...
             + 225*K0^2*Kg^2*Kwater^2*Msat^2*Por_t^2 - 450*K0^2*Kg^2*Kwater^2*Msat^2*Por_t + 225*K0^2*Kg^2*Kwater^2*Msat^2 - 480*K0^2*Kg^2*Kwater^2*Msat*Por_t^2*U0 + 1200*K0^2*Kg^2*Kwater^2*Msat*Por_t*U0 ...
             - 840*K0^2*Kg^2*Kwater^2*Msat*U0 + 576*K0^2*Kg^2*Kwater^2*Por_t^2*U0^2 - 768*K0^2*Kg^2*Kwater^2*Por_t*U0^2 + 576*K0^2*Kg^2*Kwater^2*U0^2 + 120*K0^2*Kg^2*Kwater*Msat^2*Por_t*U0 ...
             - 128*K0^2*Kg^2*Kwater*Msat*Por_t*U0^2 - 120*K0^2*Kg*Kwater^2*Msat^2*Por_t*U0 + 120*K0^2*Kg*Kwater^2*Msat^2*U0 + 128*K0^2*Kg*Kwater^2*Msat*Por_t*U0^2 - 192*K0^2*Kg*Kwater^2*Msat*U0^2 ...
             + 16*K0^2*Kwater^2*Msat^2*U0^2 - 120*K0*Kg^4*Kwater^2*U0 + 240*K0*Kg^4*Kwater*Msat*Por_t*U0 - 128*K0*Kg^4*Kwater*Por_t*U0^2 - 120*K0*Kg^4*Msat^2*Por_t^2*U0 + 128*K0*Kg^4*Msat*Por_t^2*U0^2 ...
             - 240*K0*Kg^3*Kwater^2*Msat*Por_t*U0 + 240*K0*Kg^3*Kwater^2*Msat*U0 + 128*K0*Kg^3*Kwater^2*Por_t*U0^2 - 192*K0*Kg^3*Kwater^2*U0^2 + 240*K0*Kg^3*Kwater*Msat^2*Por_t^2*U0 - 240*K0*Kg^3*Kwater*Msat^2*Por_t*U0 ...
             - 256*K0*Kg^3*Kwater*Msat*Por_t^2*U0^2 + 320*K0*Kg^3*Kwater*Msat*Por_t*U0^2 - 120*K0*Kg^2*Kwater^2*Msat^2*Por_t^2*U0 + 240*K0*Kg^2*Kwater^2*Msat^2*Por_t*U0 - 120*K0*Kg^2*Kwater^2*Msat^2*U0 ...
             + 128*K0*Kg^2*Kwater^2*Msat*Por_t^2*U0^2 - 320*K0*Kg^2*Kwater^2*Msat*Por_t*U0^2 + 224*K0*Kg^2*Kwater^2*Msat*U0^2 - 32*K0*Kg^2*Kwater*Msat^2*Por_t*U0^2 + 32*K0*Kg*Kwater^2*Msat^2*Por_t*U0^2 ...
             - 32*K0*Kg*Kwater^2*Msat^2*U0^2 + 16*Kg^4*Kwater^2*U0^2 - 32*Kg^4*Kwater*Msat*Por_t*U0^2 + 16*Kg^4*Msat^2*Por_t^2*U0^2 + 32*Kg^3*Kwater^2*Msat*Por_t*U0^2 - 32*Kg^3*Kwater^2*Msat*U0^2 ...
             - 32*Kg^3*Kwater*Msat^2*Por_t^2*U0^2 + 32*Kg^3*Kwater*Msat^2*Por_t*U0^2 + 16*Kg^2*Kwater^2*Msat^2*Por_t^2*U0^2 - 32*Kg^2*Kwater^2*Msat^2*Por_t*U0^2 + 16*Kg^2*Kwater^2*Msat^2*U0^2)^(1/2) ...
             + 15*K0*Kg*Kwater*Msat + 16*K0*Kg*Kwater*U0 + 4*K0*Kwater*Msat*U0 - 4*Kg*Kwater*Msat*U0 + 15*K0*Kg^2*Msat*Por_t + 24*K0*Kg^2*Por_t*U0 - 4*Kg^2*Msat*Por_t*U0 - 15*K0*Kg*Kwater*Msat*Por_t ...
             - 24*K0*Kg*Kwater*Por_t*U0 + 4*Kg*Kwater*Msat*Por_t*U0))/(8*(15*K0*Kg*Kwater + 4*K0*Kwater*U0 - 4*Kg*Kwater*U0 + 15*K0*Kg^2*Por_t - 4*Kg^2*Por_t*U0 - 15*K0*Kg*Kwater*Por_t + 4*Kg*Kwater*Por_t*U0));


         if(Kun1>0)
             Kun=Kun1;
         else
             Kun=Kun2;
         end  

         if(Uun1*Uun2>0)
             Uun=min(Uun1,Uun2);
         else
             Uun=max(Uun1,Uun2);
         end 
         Kh= Kun; 
    end



%%%%%%%%%%%% COMPUTE C TENSORS AS FUNCTIONS OF SATURATION SNA INCIDENT ANGLE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    SSw=linspace(0,1,nSw);
    Inc=linspace(0,90,nInc);

    vp=zeros(nSw,nInc);
    mp=zeros(nSw,nInc);
    vp_r=zeros(nSw,nInc);
    mp_r=zeros(nSw,nInc);
    vp_v=zeros(nSw,nInc);
    mp_v=zeros(nSw,nInc);
    Den=zeros(1,nSw);


    vp_low=zeros(nSw,nInc);
    mp_low=zeros(nSw,nInc);
    vp_r_low=zeros(nSw,nInc);
    mp_r_low=zeros(nSw,nInc);
    vp_v_low=zeros(nSw,nInc);
    mp_v_low=zeros(nSw,nInc);  
    
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
        Sg=1-Sw;            
        Den(iSw)=Sw*Den_gw+Sg*Den_gg;
        Den_gwg=Sw*Den_gw+Sg*Den_gg;


        if Kdry_method==3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% squirt model %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

            Kfl=Sw*Kwater+(1-Sw)*Kgas;
            Kuf=1/(1/Kh+1/(1/(1/K0-1/Kh)+1/((1/Kfl-1/Kg)*phi_c)));% equation 24 in Gurevich et al's paper(2009)      
            Kdry=Kuf;
            Udry=1./(1./U0-4./15*(1./K0-1./Kdry)); % equation 2 in Gurevich et al's paper(2009)     
        end

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% dry model %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
         %%% Hudson or linear slip model
        Cdry11 = Kdry + 4./3 * Udry;
        Cdry12 = Kdry - 2./3 * Udry;
        Cdry13 = Kdry - 2./3 * Udry;
        Cdry22 = Kdry + 4./3 * Udry;
        Cdry23 = Kdry - 2./3 * Udry;
        Cdry33 = Kdry + 4./3 * Udry;
        Cdry44 = Udry;
        Cdry55 = Udry;
        Cdry66 = Udry;

        Kf=1./(Sw/Kwater+Sg/Kgas);

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%% uniform model %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
        alpha1=1-(Cdry11+Cdry12+Cdry13)/3/Kg;
        alpha2=1-(Cdry12+Cdry22+Cdry23)/3/Kg;
        alpha3=1-(Cdry13+Cdry23+Cdry33)/3/Kg;

        Kstar=(Cdry11+Cdry12+Cdry13+Cdry12+Cdry22+Cdry23+Cdry13+Cdry23+Cdry33)/9;
        Mf=Kg/((1-Kstar/Kg)-Por*(1-Kg/Kf));

        C_low11(iSw)=Cdry11+alpha1*alpha1*Mf;
        C_low12(iSw)=Cdry12+alpha1*alpha2*Mf;
        C_low13(iSw)=Cdry13+alpha1*alpha3*Mf;
        C_low22(iSw)=Cdry22+alpha2*alpha2*Mf;
        C_low23(iSw)=Cdry23+alpha2*alpha3*Mf;
        C_low33(iSw)=Cdry33+alpha3*alpha3*Mf;
        C_low44(iSw)=Cdry44;
        C_low55(iSw)=Cdry55;
        C_low66(iSw)=Cdry66;  
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% patchy model %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

        %%%% Gassmann_model %%%%%%%%%%%%%% 
        alpha1=1-(Cdry11+Cdry12+Cdry13)/3/Kg;
        alpha2=1-(Cdry12+Cdry22+Cdry23)/3/Kg;
        alpha3=1-(Cdry13+Cdry23+Cdry33)/3/Kg;

        Kstar=(Cdry11+Cdry12+Cdry13+Cdry12+Cdry22+Cdry23+Cdry13+Cdry23+Cdry33)/9;
        Mwater=Kg/((1-Kstar/Kg)-Por*(1-Kg/Kwater));
        Mgas=Kg/((1-Kstar/Kg)-Por*(1-Kg/Kgas));

        Cwater11=Cdry11+alpha1*alpha1*Mwater;
        Cwater12=Cdry12+alpha1*alpha2*Mwater;
        Cwater13=Cdry13+alpha1*alpha3*Mwater;
        Cwater22=Cdry22+alpha2*alpha2*Mwater;
        Cwater23=Cdry23+alpha2*alpha3*Mwater;
        Cwater33=Cdry33+alpha3*alpha3*Mwater;
        Cwater44=Cdry44;
        Cwater55=Cdry55;
        Cwater66=Cdry66;

        Cgas11=Cdry11+alpha1*alpha1*Mgas;
        Cgas12=Cdry12+alpha1*alpha2*Mgas;
        Cgas13=Cdry13+alpha1*alpha3*Mgas;
        Cgas22=Cdry22+alpha2*alpha2*Mgas;
        Cgas23=Cdry23+alpha2*alpha3*Mgas;
        Cgas33=Cdry33+alpha3*alpha3*Mgas;
        Cgas44=Cdry44;
        Cgas55=Cdry55;
        Cgas66=Cdry66;   
        
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

        %%%%%%%%%%%%%%%%% Ruess bound %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
        Cwgr11=Sw*Cwater11+Sg*Cgas11;
        Cwgr12=Sw*Cwater12+Sg*Cgas12;
        Cwgr13=Sw*Cwater13+Sg*Cgas13;
        Cwgr22=Sw*Cwater22+Sg*Cgas22;
        Cwgr23=Sw*Cwater23+Sg*Cgas23;
        Cwgr33=Sw*Cwater33+Sg*Cgas33;
        Cwgr44=Sw*Cwater44+Sg*Cgas44;
        Cwgr55=Sw*Cwater55+Sg*Cgas55;
        Cwgr66=Sw*Cwater66+Sg*Cgas66;

        %%%%%%%%%%%%%%%%% Voigt bound %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
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
            vp(iSw,iInc)=sqrt((Cwg44+Cwg11*G+Cwg33*E+sqrt(D))/(2*Den_gwg));
            mp(iSw,iInc)=Den_gwg*vp(iSw,iInc)*vp(iSw,iInc);

            Dr=((Cwgr11-Cwgr44)*G-(Cwgr33-Cwgr44)*E)*((Cwgr11-Cwgr44)*G-(Cwgr33-Cwgr44)*E)+4*(Cwgr13+Cwgr44)*(Cwgr13+Cwgr44)*G*E;
            vpr(iSw,iInc)=sqrt((Cwgr44+Cwgr11*G+Cwgr33*E+sqrt(D))/(2*Den_gwg));
            mpr(iSw,iInc)=Den_gwg*vpr(iSw,iInc)*vpr(iSw,iInc);

            Dv=((Cwgv11-Cwgv44)*G-(Cwgv33-Cwgv44)*E)*((Cwgv11-Cwgv44)*G-(Cwgv33-Cwgv44)*E)+4*(Cwgv13+Cwgv44)*(Cwgv13+Cwgv44)*G*E;
            vpv(iSw,iInc)=sqrt((Cwgv44+Cwgv11*G+Cwgv33*E+sqrt(Dv))/(2*Den_gwg));
            mpv(iSw,iInc)=Den_gwg*vpv(iSw,iInc)*vpv(iSw,iInc);

            D=((C_low11(iSw)-C_low44(iSw))*G-(C_low33(iSw)-C_low44(iSw))*E)*((C_low11(iSw)-C_low44(iSw))*G-(C_low33(iSw)-C_low44(iSw))*E)+4*(C_low13(iSw)+C_low44(iSw))*(C_low13(iSw)+C_low44(iSw))*G*E;
            vv=sqrt((C_low44(iSw)+C_low11(iSw)*G+C_low33(iSw)*E+sqrt(D))/(2*Den_gwg));
            vp_low(iSw,iInc)=1./(real(1./vv));
            mp_low(iSw,iInc)=Den_gwg*vp_low(iSw,iInc)*vp_low(iSw,iInc);
        end
        
        %%%%%%%%% final modeling results of combining models %%%%%%%%%%%%%
        if Sw<0.5
            vp_final(iSw,:)=vp_low(iSw,:);
            mp_final(iSw,:)=mp_low(iSw,:);
        else
            vp_final(iSw,:)=vp(iSw,:);
            mp_final(iSw,:)=mp(iSw,:);
        end       
    end

%%%%%%%%%%%% PLOT RESUTLS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   

    if Kdry_method==1  % Kdry from data at Sw=0
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
%         h1=plot(blank_correct(:,1)*0.01,blank_correct(:,8),'ws','markerfacecolor','r','markersize',20);
%         hold on;  
%         h2=plot(SSw,mp(:,3),'r','LineWidth',4);
%         hold on; 
%         h3=plot(blank_correct(:,1)*0.01,blank_correct(:,7),'wo','markerfacecolor','g','markersize',20);
%         hold on; 
%         h4=plot(SSw,mp(:,2),'g','LineWidth',4);
%         hold on; 
%         h5=plot(blank_correct(:,1)*0.01,blank_correct(:,6),'wd','markerfacecolor','b','markersize',20);
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
%         h1=plot(blank_correct(:,1)*0.01,blank_correct(:,5),'ws','markerfacecolor','r','markersize',20);
%         hold on;  
%         h2=plot(SSw,vp(:,3),'r','LineWidth',4);
%         hold on; 
%         h3=plot(blank_correct(:,1)*0.01,blank_correct(:,4),'wo','markerfacecolor','g','markersize',20);
%         hold on; 
%         h4=plot(SSw,vp(:,2),'g','LineWidth',4);
%         hold on; 
%         h5=plot(blank_correct(:,1)*0.01,blank_correct(:,3),'wd','markerfacecolor','b','markersize',20);
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
% %         h7=plot(SSw(11:nSw),mpr(11:nSw,3),'r:.','LineWidth',4);
% %         hold on; 
% %         h8=plot(SSw(11:nSw),mpr(11:nSw,2),'g:.','LineWidth',4);
% %         hold on;  
% %         h9=plot(SSw(11:nSw),mpr(11:nSw,1),'b:.','LineWidth',4);
% %         hold on; 
% %         h10=plot(SSw(11:nSw),mpv(11:nSw,3),'r:.','LineWidth',4);
% %         hold on;
% %         h11=plot(SSw(11:nSw),mpv(11:nSw,2),'g:.','LineWidth',4);
% %         hold on;   
% %         h12=plot(SSw(11:nSw),mpv(11:nSw,1),'b:.','LineWidth',4);
% %         hold on;  
%         h1=plot(blank_correct(:,1)*0.01,blank_correct(:,8),'ws','markerfacecolor','r','markersize',20);
%         hold on;  
%         h2=plot(SSw(11:nSw),mp(11:nSw,3),'r','LineWidth',4);
%         hold on; 
%         h3=plot(blank_correct(:,1)*0.01,blank_correct(:,7),'wo','markerfacecolor','g','markersize',20);
%         hold on; 
%         h4=plot(SSw(11:nSw),mp(11:nSw,2),'g','LineWidth',4);
%         hold on; 
%         h5=plot(blank_correct(:,1)*0.01,blank_correct(:,6),'wd','markerfacecolor','b','markersize',20);
%         hold on;
%         h6=plot(SSw(11:nSw),mp(11:nSw,1),'b','LineWidth',4);
%         hold on;
%         set(gca,'FontSize',30);
%         set(gca,'xTick',[0 0.25 0.5 0.75 1]);
%         xlabel('Sw','fontsize',40);
%         ylabel('Mp(Pa)','fontsize',40);    
%         set(gca,'ylim',[14e9 25e9]);
%         
%         subplot(1,2,2)
%         set(gcf,'outerposition',get(0,'screensize'));         
% %         h7=plot(SSw(11:nSw),vpr(11:nSw,3),'r:.','LineWidth',4);
% %         hold on; 
% %         h8=plot(SSw(11:nSw),vpr(11:nSw,2),'g:.','LineWidth',4);
% %         hold on;  
% %         h9=plot(SSw(11:nSw),vpr(11:nSw,1),'b:.','LineWidth',4);
% %         hold on; 
% %         h10=plot(SSw(11:nSw),vpv(11:nSw,3),'r:.','LineWidth',4);
% %         hold on;
% %         h11=plot(SSw(11:nSw),vpv(11:nSw,2),'g:.','LineWidth',4);
% %         hold on;   
% %         h12=plot(SSw(11:nSw),vpv(11:nSw,1),'b:.','LineWidth',4);
% %         hold on;  
%         h1=plot(blank_correct(:,1)*0.01,blank_correct(:,5),'ws','markerfacecolor','r','markersize',20);
%         hold on;  
%         h2=plot(SSw(11:nSw),vp(11:nSw,3),'r','LineWidth',4);
%         hold on; 
%         h3=plot(blank_correct(:,1)*0.01,blank_correct(:,4),'wo','markerfacecolor','g','markersize',20);
%         hold on; 
%         h4=plot(SSw(11:nSw),vp(11:nSw,2),'g','LineWidth',4);
%         hold on; 
%         h5=plot(blank_correct(:,1)*0.01,blank_correct(:,3),'wd','markerfacecolor','b','markersize',20);
%         hold on;
%         h6=plot(SSw(11:nSw),vp(11:nSw,1),'b','LineWidth',4);
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
%         plot(SSw(2:11),mp_low(2:11,3),'r','LineWidth',4);
%         hold on;  
%         plot(SSw(2:11),mp_low(2:11,2),'g','LineWidth',4);
%         hold on;  
%         plot(SSw(2:11),mp_low(2:11,1),'b','LineWidth',4);
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
%         h1=plot(blank_correct(:,1)*0.01,blank_correct(:,8),'ws','markerfacecolor','r','markersize',20);
%         hold on;  
% %         h2=plot(SSw,mp(:,3),'r','LineWidth',4);
% %         hold on; 
%         h3=plot(blank_correct(:,1)*0.01,blank_correct(:,7),'wo','markerfacecolor','g','markersize',20);
%         hold on; 
% %         h4=plot(SSw,mp(:,2),'g','LineWidth',4);
% %         hold on; 
%         h5=plot(blank_correct(:,1)*0.01,blank_correct(:,6),'wd','markerfacecolor','b','markersize',20);
%         hold on;
% %         h6=plot(SSw,mp(:,1),'b','LineWidth',4);
% %         hold on;
%         set(gca,'FontSize',30);
%         set(gca,'xTick',[0 0.25 0.5 0.75 1]);
%         xlabel('Sw','fontsize',40);
%         ylabel('Mp(Pa)','fontsize',40);    
%         set(gca,'ylim',[14e9 25e9]);
%         
%         subplot(1,2,2)
%         set(gcf,'outerposition',get(0,'screensize'));
%         plot(SSw(2:11),vp_low(2:11,3),'r','LineWidth',4);
%         hold on;  
%         plot(SSw(2:11),vp_low(2:11,2),'g','LineWidth',4);
%         hold on;  
%         plot(SSw(2:11),vp_low(2:11,1),'b','LineWidth',4);
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
%         h1=plot(blank_correct(:,1)*0.01,blank_correct(:,5),'ws','markerfacecolor','r','markersize',20);
%         hold on;  
% %         h2=plot(SSw,vp(:,3),'r','LineWidth',4);
% %         hold on; 
%         h3=plot(blank_correct(:,1)*0.01,blank_correct(:,4),'wo','markerfacecolor','g','markersize',20);
%         hold on; 
% %         h4=plot(SSw,vp(:,2),'g','LineWidth',4);
% %         hold on; 
%         h5=plot(blank_correct(:,1)*0.01,blank_correct(:,3),'wd','markerfacecolor','b','markersize',20);
%         hold on;
% %         h6=plot(SSw,vp(:,1),'b','LineWidth',4);
% %         hold on;
%         set(gca,'FontSize',30);
%         set(gca,'yTick',[2700 2900 3100 3300 3500]);
%         set(gca,'xTick',[0 0.25 0.5 0.75 1]);
%         xlabel('Sw','fontsize',40);
%         ylabel('Vp(m/s)','fontsize',40);    
%         set(gca,'ylim',[2700 3500]);
        
    elseif Kdry_method==3 % Kdry from data at Sw=100 considering squirt process
        
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
%         h1=plot(blank_correct(:,1)*0.01,blank_correct(:,8),'ws','markerfacecolor','r','markersize',20);
%         hold on;  
% %         h2=plot(SSw,mp(:,3),'r','LineWidth',4);
% %         hold on; 
%         h3=plot(blank_correct(:,1)*0.01,blank_correct(:,7),'wo','markerfacecolor','g','markersize',20);
%         hold on; 
% %         h4=plot(SSw,mp(:,2),'g','LineWidth',4);
% %         hold on; 
%         h5=plot(blank_correct(:,1)*0.01,blank_correct(:,6),'wd','markerfacecolor','b','markersize',20);
%         hold on;
%         
%         plot(SSw(nSw:nSw),mp_final(nSw:nSw,3),'om','markersize',30,'LineWidth',5);
%         hold on;  
%         plot(SSw(nSw:nSw),mp_final(nSw:nSw,2),'om','markersize',30,'LineWidth',5);
%         hold on;  
%         plot(SSw(nSw:nSw),mp_final(nSw:nSw,1),'om','markersize',30,'LineWidth',5);
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
%         h1=plot(blank_correct(:,1)*0.01,blank_correct(:,5),'ws','markerfacecolor','r','markersize',20);
%         hold on;  
% %         h2=plot(SSw,vp(:,3),'r','LineWidth',4);
% %         hold on; 
%         h3=plot(blank_correct(:,1)*0.01,blank_correct(:,4),'wo','markerfacecolor','g','markersize',20);
%         hold on; 
% %         h4=plot(SSw,vp(:,2),'g','LineWidth',4);
% %         hold on; 
%         h5=plot(blank_correct(:,1)*0.01,blank_correct(:,3),'wd','markerfacecolor','b','markersize',20);
%         hold on;
% %         h6=plot(SSw,vp(:,1),'b','LineWidth',4);
% %         hold on;
%         plot(SSw(nSw:nSw),vp_final(nSw:nSw,3),'om','markersize',30,'LineWidth',5);
%         hold on;  
%         plot(SSw(nSw:nSw),vp_final(nSw:nSw,2),'om','markersize',30,'LineWidth',5);
%         hold on;  
%         plot(SSw(nSw:nSw),vp_final(nSw:nSw,1),'om','markersize',30,'LineWidth',5);
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
%         plot(SSw(1:4),mp_final(1:4,3),'m','LineWidth',4);
%         hold on;  
%         plot(SSw(1:4),mp_final(1:4,2),'m','LineWidth',4);
%         hold on;  
%         plot(SSw(1:4),mp_final(1:4,1),'m','LineWidth',4);
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
%         h1=plot(blank_correct(:,1)*0.01,blank_correct(:,8),'ws','markerfacecolor','r','markersize',20);
%         hold on;  
% %         h2=plot(SSw,mp(:,3),'r','LineWidth',4);
% %         hold on; 
%         h3=plot(blank_correct(:,1)*0.01,blank_correct(:,7),'wo','markerfacecolor','g','markersize',20);
%         hold on; 
% %         h4=plot(SSw,mp(:,2),'g','LineWidth',4);
% %         hold on; 
%         h5=plot(blank_correct(:,1)*0.01,blank_correct(:,6),'wd','markerfacecolor','b','markersize',20);
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
%         plot(SSw(1:4),vp_final(1:4,3),'m','LineWidth',4);
%         hold on;  
%         plot(SSw(1:4),vp_final(1:4,2),'m','LineWidth',4);
%         hold on;  
%         plot(SSw(1:4),vp_final(1:4,1),'m','LineWidth',4);
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
%         h1=plot(blank_correct(:,1)*0.01,blank_correct(:,5),'ws','markerfacecolor','r','markersize',20);
%         hold on;  
% %         h2=plot(SSw,vp(:,3),'r','LineWidth',4);
% %         hold on; 
%         h3=plot(blank_correct(:,1)*0.01,blank_correct(:,4),'wo','markerfacecolor','g','markersize',20);
%         hold on; 
% %         h4=plot(SSw,vp(:,2),'g','LineWidth',4);
% %         hold on; 
%         h5=plot(blank_correct(:,1)*0.01,blank_correct(:,3),'wd','markerfacecolor','b','markersize',20);
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
%         plot(SSw(4:10),mp_final(4:10,3),'m','LineWidth',4);
%         hold on;  
%         plot(SSw(4:10),mp_final(4:10,2),'m','LineWidth',4);
%         hold on;  
%         plot(SSw(4:10),mp_final(4:10,1),'m','LineWidth',4);
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
%         h1=plot(blank_correct(:,1)*0.01,blank_correct(:,8),'ws','markerfacecolor','r','markersize',20);
%         hold on;  
% %         h2=plot(SSw,mp(:,3),'r','LineWidth',4);
% %         hold on; 
%         h3=plot(blank_correct(:,1)*0.01,blank_correct(:,7),'wo','markerfacecolor','g','markersize',20);
%         hold on; 
% %         h4=plot(SSw,mp(:,2),'g','LineWidth',4);
% %         hold on; 
%         h5=plot(blank_correct(:,1)*0.01,blank_correct(:,6),'wd','markerfacecolor','b','markersize',20);
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
%         plot(SSw(4:10),vp_final(4:10,3),'m','LineWidth',4);
%         hold on;  
%         plot(SSw(4:10),vp_final(4:10,2),'m','LineWidth',4);
%         hold on;  
%         plot(SSw(4:10),vp_final(4:10,1),'m','LineWidth',4);
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
%         h1=plot(blank_correct(:,1)*0.01,blank_correct(:,5),'ws','markerfacecolor','r','markersize',20);
%         hold on;  
% %         h2=plot(SSw,vp(:,3),'r','LineWidth',4);
% %         hold on; 
%         h3=plot(blank_correct(:,1)*0.01,blank_correct(:,4),'wo','markerfacecolor','g','markersize',20);
%         hold on; 
% %         h4=plot(SSw,vp(:,2),'g','LineWidth',4);
% %         hold on; 
%         h5=plot(blank_correct(:,1)*0.01,blank_correct(:,3),'wd','markerfacecolor','b','markersize',20);
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
%         plot(SSw(11:nSw-1),mp_final(11:nSw-1,3),'m','LineWidth',4);
%         hold on;  
%         plot(SSw(11:nSw-1),mp_final(11:nSw-1,2),'m','LineWidth',4);
%         hold on;  
%         plot(SSw(11:nSw-1),mp_final(11:nSw-1,1),'m','LineWidth',4);
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
%         h1=plot(blank_correct(:,1)*0.01,blank_correct(:,8),'ws','markerfacecolor','r','markersize',20);
%         hold on;  
% %         h2=plot(SSw,mp(:,3),'r','LineWidth',4);
% %         hold on; 
%         h3=plot(blank_correct(:,1)*0.01,blank_correct(:,7),'wo','markerfacecolor','g','markersize',20);
%         hold on; 
% %         h4=plot(SSw,mp(:,2),'g','LineWidth',4);
% %         hold on; 
%         h5=plot(blank_correct(:,1)*0.01,blank_correct(:,6),'wd','markerfacecolor','b','markersize',20);
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
%         plot(SSw(11:nSw-1),vp_final(11:nSw-1,3),'m','LineWidth',4);
%         hold on;  
%         plot(SSw(11:nSw-1),vp_final(11:nSw-1,2),'m','LineWidth',4);
%         hold on;  
%         plot(SSw(11:nSw-1),vp_final(11:nSw-1,1),'m','LineWidth',4);
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
%         h1=plot(blank_correct(:,1)*0.01,blank_correct(:,5),'ws','markerfacecolor','r','markersize',20);
%         hold on;  
% %         h2=plot(SSw,vp(:,3),'r','LineWidth',4);
% %         hold on; 
%         h3=plot(blank_correct(:,1)*0.01,blank_correct(:,4),'wo','markerfacecolor','g','markersize',20);
%         hold on; 
% %         h4=plot(SSw,vp(:,2),'g','LineWidth',4);
% %         hold on; 
%         h5=plot(blank_correct(:,1)*0.01,blank_correct(:,3),'wd','markerfacecolor','b','markersize',20);
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
%         subplot(1,2,1)
        set(gcf,'outerposition',get(0,'screensize'));
        h1=plot(blank_correct(:,1)*0.01,blank_correct(:,8),'ws','markerfacecolor','r','markersize',20);
        hold on;          
        h2=plot(blank_correct(:,1)*0.01,blank_correct(:,7),'wo','markerfacecolor','g','markersize',20);
        hold on;         
        h3=plot(blank_correct(:,1)*0.01,blank_correct(:,6),'wd','markerfacecolor','b','markersize',20);
        hold on;
        h4=plot(SSw(1:20),mp_final(1:20,3),'m','LineWidth',4);
%         hold on;
%         h4=plot(SSw(1:20),mp_final(1:20,2),'m','LineWidth',4);
%         hold on; 
%         h6=plot(SSw(1:20),mp_final(1:20,1),'m','LineWidth',4);
%         hold on;
        h5=plot(SSw(20:50),mp_final(20:50,3),'m--','LineWidth',4);
%         hold on;
%         h4=plot(SSw(20:50),mp_final(20:50,2),'m--','LineWidth',4);
%         hold on; 
%         h6=plot(SSw(20:50),mp_final(20:50,1),'m--','LineWidth',4);
%         hold on;
        h6=plot(SSw(51:nSw-1),mp_final(51:nSw-1,3),'m.','markersize',25,'LineWidth',3);
%         hold on;
%         plot(SSw(51:nSw-1),mp_final(51:nSw-1,2),'m.','markersize',12,'LineWidth',3);
%         hold on; 
%         plot(SSw(51:nSw-1),mp_final(51:nSw-1,1),'m.','markersize',12,'LineWidth',3);
%         hold on;
        h7=plot(SSw(nSw:nSw),mp_final(nSw:nSw,3),'om','markersize',25,'LineWidth',3);
%         hold on;
%         plot(SSw(nSw:nSw),mp_final(nSw:nSw,2),'om','markersize',25,'LineWidth',3);
%         hold on; 
%        plot(SSw(nSw:nSw),mp_final(nSw:nSw,1),'om','markersize',25,'LineWidth',3);
        hold on;
        h8=plot(SSw,mp_low(:,3),'m');  
        hold on; 
        set(gca,'FontSize',30);
        set(gca,'xTick',[0 0.25 0.5 0.75 1]);
        xlabel('S_w','fontsize',40);
        ylabel('M_p(Pa)','fontsize',40);    
        set(gca,'ylim',[14e9 25e9]);       
        
        lgd1=legend([h1,h2,h3],'Data-90бу','Data-45бу','Data-0бу','orientation','vertical','location','south');
        set(lgd1,'FontName','Time New Roman','FontSize',30);
        legend boxoff;
        ah1=axes('position',get(gca,'position'),'visible','off');
        lgd2=legend(ah1,[h4,h5],'Squirt','Uniform','orientation','vertical','location','south');
        set(lgd2,'FontName','Time New Roman','FontSize',30);
        legend boxoff;
        ah2=axes('position',get(gca,'position'),'visible','off');
        lgd3=legend(ah2,[h6,h7],'Patchy','Sat','orientation','vertical','location','south');
        set(lgd3,'FontName','Time New Roman','FontSize',30);
        legend boxoff;
        
        figure()
%         subplot(1,2,2)
        set(gcf,'outerposition',get(0,'screensize'));
        h1=plot(blank_correct(:,1)*0.01,blank_correct(:,5),'ws','markerfacecolor','r','markersize',20);
        hold on;         
        h2=plot(blank_correct(:,1)*0.01,blank_correct(:,4),'wo','markerfacecolor','g','markersize',20);
        hold on;       
        h3=plot(blank_correct(:,1)*0.01,blank_correct(:,3),'wd','markerfacecolor','b','markersize',20);
        hold on;
      
        h4=plot(SSw(1:20),vp_final(1:20,3),'m','LineWidth',4);
%         hold on;
%         h5=plot(SSw(1:20),vp_final(1:20,2),'m','LineWidth',4);
%         hold on; 
%         h6=plot(SSw(1:20),vp_final(1:20,1),'m','LineWidth',4);
%         hold on;
        h5=plot(SSw(20:50),vp_final(20:50,3),'m--','LineWidth',4);
%         hold on;
%         h4=plot(SSw(20:50),vp_final(20:50,2),'m--','LineWidth',4);
%         hold on; 
%         h6=plot(SSw(20:50),vp_final(20:50,1),'m--','LineWidth',4);
%         hold on;
        h6=plot(SSw(51:nSw-1),vp_final(51:nSw-1,3),'m.','markersize',25,'LineWidth',3);
        hold on;
%         plot(SSw(51:nSw-1),vp_final(51:nSw-1,2),'m.','markersize',12,'LineWidth',3);
%         hold on; 
%         plot(SSw(51:nSw-1),vp_final(51:nSw-1,1),'m.','markersize',12,'LineWidth',3);
%         hold on;
        h7=plot(SSw(nSw:nSw),vp_final(nSw:nSw,3),'om','markersize',25,'LineWidth',3);
        hold on; 
        h8=plot(SSw,vp_low(:,3),'m');  
        hold on; 
%         plot(SSw(nSw:nSw),vp_final(nSw:nSw,2),'om','markersize',25,'LineWidth',3);
%         hold on; 
%        plot(SSw(nSw:nSw),vp_final(nSw:nSw,1),'om','markersize',25,'LineWidth',3);
        set(gca,'FontSize',30);
        set(gca,'yTick',[2700 2900 3100 3300 3500]);
        set(gca,'xTick',[0 0.25 0.5 0.75 1]);
        xlabel('S_w','fontsize',40);
        ylabel('V_p(m/s)','fontsize',40);    
        set(gca,'ylim',[2700 3500]);
        


    end
end

