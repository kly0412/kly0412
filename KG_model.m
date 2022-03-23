function [C11 C12 C13 C22 C23 C33 C44 C55 C66]=KG_model(hc,H,Kg,Ug,Kb,Ub,Kfb,Dfb,Pb,Fb,Kfc,dN,dT,w0)

    Lb=Kb+4./3*Ub;
    Ab=1-Kb/Kg; % Biot's parameter alpha
    Mb=1./((Ab-Fb)/Kg+Fb/Kfb); % Biot's parameter M
    Cb=Lb+Ab*Ab*Mb;
    %*********** for fracture  ******************

    Zn=dN/(Lb*(1-dN));% fracture compliances
    Zt=dT/(Ub*(1-dT));
    Lc=hc/Zn;

%%%%%%% Low frequency %%%%%%%%%%%%%%%%%%

    %%%%% C33_Low 
    F=Kfc/Lc;
    A=1/Cb+Zn/(1+F);
    N=(Ab*Mb/Cb-F/(1+F))^2;
    D=Mb*Lb/Cb+Kfc/hc/(1+F);
    C33_Low = 1/(A+N/D);

    %%%%% C13_Low
    A=1-2*Ub/Lb;
    N=2*Ub*Ab/Lb*(Ab/Lb+Zn);
    D=Cb/Mb/Lb+(1+F)/(Kfc/hc);
    C13_Low = C33_Low*(A+N/D);

    %%%%% C11_Low
    A=4*(Ub-Ub*Ub/Lb);
    N=4*(Ub*Ab/Lb)^2;
    D=Cb/Mb/Lb+(1+F)/(Kfc/hc);
    C11_Low = C13_Low^2/C33_Low + A + N/D;

    %%%%% C44_Low
    C44_Low = 1/(1/Ub+Zt);

    %%%%% C66_Low
    C66_Low = Ub;

    %%%%% C12_Low
    C12_Low = C11_Low-2*C66_Low;

    %%%%%%% High frequency %%%%%%%%%%%%%%%%%%

    %%%%% C33_High 
    A=1/Cb+Zn/(1+F);
    C33_High = 1/A;

    %%%%% C13_High
    A=1-2*Ub/Cb;
    C13_High = C33_High*A;

    %%%%% C11_High
    A=4*(Ub-Ub*Ub/Cb);
    C11_High = C13_High^2/C33_High + A;

    %%%%% C44_High
    C44_High = 1/(1/Ub+Zt);

    %%%%% C66_High
    C66_High = Ub;

    %%%%% C12_High
    C12_High = C11_High-2*C66_High;
    
%%%%%%%%%%%%%%%%%%%%% Compute scalar complex function R(w) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
        w=w0*H^2*Dfb*Mb/(4*Cb*Pb*Lb);
        A=1/Cb+Zn/(1+F);
        N=(Ab*Mb/Cb-F/(1+F))^2;
        Db1=Lb*sqrt(i*w);
        Db2=Cb/Mb*sqrt(i*w);
        Db=Db1* cot(Db2);
        Dc=Kfc/hc/(1+F);

        C330=1/(A+N/(Db+Dc));

        R = (C330-C33_High)/(C33_Low-C33_High);

        %%%%%%%%%%%%%%%%%%%%%% step 3:  compute effective stiffness tensor Cij %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        C11 = C11_High - R * (C11_High - C11_Low);
        C12 = C12_High - R * (C12_High - C12_Low);
        C13 = C13_High - R * (C13_High - C13_Low);
        C22 = C11;
        C23 = C13;
        C33 = C33_High - R * (C33_High - C33_Low);
        C44 = C44_High - R * (C44_High - C44_Low);
        C55 = C44;
        C66 = C66_High - R * (C66_High - C66_Low);
        
end