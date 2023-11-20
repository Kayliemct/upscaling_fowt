%% Upscale using thetaP limit for OC4

clear all; close all; clc

Power = [5; 10; 15; 20];%; 25; 30];  %rated power
Sp = 401;
speed_rated = 11.4;     % high specific power speed

nturbine = 4;   
for j=1:nturbine
    radius(j) = round((sqrt((Power(j).*10^6)./(pi().*Sp))));
    A(j) = pi().*radius(j).^2;
end


rhoA = 1.225;
rhoW = 1025;
rhoM = 7850;
g = 9.81;
Ct = 0.7;      %0.75    0.84;      %Cp = 0.49

%% Upscale 
alpha = [0:0.005:2]; %1.5];
alphaT = 2.2; %RNA assume geometric similarity upscaling volume   %1.5;
alphat = 2; %tower

for j=1:nturbine
m_t(j,1) = 249718./(radius(1)/radius(j)).^alphat;
m_rna(j,1) = 350000./(radius(1)/radius(j)).^alphaT; 
turbine(j,1) = m_t(j)+m_rna(j);
hub_ht(j,1) = radius(j)+30;  %radius(j).*1.43;
hub_ht(1,1) = 90;  %OC4 90m hub ht
Thrust(j,1) = 0.5.*rhoA.*A(j).*speed_rated.^2.*Ct;           %using efficient turbines for now
% Thrust(1,1) = 832769;       %Thrust at rated avg from FAST *ensure aero
CM_t(j,1) = hub_ht(j,1).*0.482;   %m   %this is approximate, OC4 CM is 48% hub ht

for i= 1:length(alpha)
    scaling(j,i) = (radius(1)/radius(j)).^alpha(i);
    draft(j,i) = 20;%./scaling(j,i);    %constant draft
    t(j,i) = 0.06;%./scaling(j,i);      %!!!changed from 6cm; constant wall thickness    %round((0.06./(scaling(j,i).^0.25)),2);       %approx 0.06 - 0.08m thickness
    rad_column(j,i) = 6./scaling(j,i);
    rad_center(j,i) = 3.25./scaling(j,i);
    dist_cc(j,i) = 50/scaling(j,i);           
    hp_thickness(j,i) = 6/scaling(j,i);
    hp_rad(j,i) = 12/scaling(j,i);
    freeboard(j,i) = 12/scaling(j,i);           %note: freeboard is 12, tower base is 10m
    cfreeb(j,i) = 10/scaling(j,i);              %note: tower base is 10m above SWL

    Vcol(j,i) = (pi()*(rad_column(j,i).^2-(rad_column(j,i)-t(j,i)).^2).*((draft(j,i)-hp_thickness(j,i))+freeboard(j,i)));
    Vhpwalls(j,i) = ((pi()*(hp_rad(j,i).^2-(hp_rad(j,i)-t(j,i)).^2)).*(hp_thickness(j,i)));
    Vhp(j,i) = (pi()*(hp_rad(j,i)-t(j,i)).^2.*t(j,i).*2);
    Vtop(j,i) = pi().*(rad_column(j,i) - t(j,i)).^2.*t(j,i);
    Vcenter(j,i) = pi().*(rad_center(j,i).^2 - (rad_center(j,i) - t(j,i)./2).^2).*(draft(j,i)+freeboard(j,i));
%     Vcbc(j,i) = pi().*(rad_center(j,i)-t(j,i)./2).^2.*(t(j,i)./2);
    Vplatform(j,i) = (3.*(Vcol(j,i)+Vhpwalls(j,i)+Vhp(j,i)+Vtop(j,i)))+Vcenter(j,i); %Volume 1% lower than OC4
    Mplatform(j,i) = Vplatform(j,i).*rhoM;
    V1 = pi().*(rad_column(j,i)).^2.*(draft(j,i) - hp_thickness(j,i));
    V2 = pi().*(hp_rad(j,i)).^2.*(hp_thickness(j,i));
    V3 = pi()*(rad_center(j,i)).^2.*draft(j,i);
    Vdisp(j,i) = (3.*(V1 + V2) + V3); %Volume 1% lower than OC4
    B(j,i) = Vdisp(j,i).*rhoW;
    Mtotal(j,i) = B(j,i) - turbine(j);
    Mballast(j,i) = Mtotal(j,i) - Mplatform(j,i);
    Vballast(j,i) = Mballast(j,i)./rhoW;
    b1(j,i) = (0.28.*Vballast(j,i))./(3.*pi().*(rad_column(j,i)-t(j,i)).^2);
    b2(j,i) = (0.72.*Vballast(j,i))./(3.*pi().*(hp_rad(j,i) - t(j,i)).^2) + t(j,i);
    CoB(j,i) = ((3.*V1)./Vdisp(j,i)).*((draft(j,i) - hp_thickness(j,i))./2) + ((3.*V2)./Vdisp(j,i)).*((draft(j,i) - hp_thickness(j,i)./2)) + (V3./Vdisp(j,i)).*(draft(j,i)./2); %OC4 13.15m
    CG1(j,i) = (((draft(j,i) - hp_thickness(j,i) + freeboard(j,i))./2) - freeboard(j,i));
    CG2(j,i) = (draft(j,i) - hp_thickness(j,i)./2);
    CG3(j,i) = (draft(j,i) - hp_thickness(j,i)+t(j,i)./2);
    CG4(j,i) = (draft(j,i)-t(j,i)./2);
    CG5(j,i) = (-freeboard(j,i)+t(j,i)./2);
    CG6(j,i) = (draft(j,i) - hp_thickness(j,i) - b1(j,i)./2);
    CG7(j,i) = (draft(j,i) - (b2(j,i) - t(j,i))./2 - t(j,i));
    CG8(j,i) = (((draft(j,i) + cfreeb(j,i))./2 - cfreeb(j,i)));
    CoG(j,i) = ((3.*Vcol(j,i).*rhoM)./Mtotal(j,i)).*CG1(j,i) + ((3.*Vhpwalls(j,i).*rhoM)./Mtotal(j,i)).*CG2(j,i) + ((3.*Vhp(j,i).*rhoM)./Mtotal(j,i)).*((CG3(j,i)+CG4(j,i))./2) + ((3.*Vtop(j,i).*rhoM)./Mtotal(j,i)).*CG5(j,i) + ((Vcenter(j,i).*rhoM)./Mtotal(j,i)).*CG8(j,i) + ((0.28.*Mballast(j,i))./Mtotal(j,i)).*CG6(j,i) + ((0.72.*Mballast(j,i))./Mtotal(j,i)).*CG7(j,i);% + ((Vcbc.*rhoM)./Mtotal(j,i)).*(draft(j,i)-t(j,i)./4);  %CoG OC4 - 13.46m
%     CoG(j,i) = 13.46/scaling(j,i);
%     CoB(j,i) = (draft(j,i)/2)*1.315;        %improve CoB estimate later, OC4 CoB 13.15m
%    
    CM(j,i) = (Mtotal(j,i)./B(j,i)).*(-CoG(j,i)) + (m_t(j,1)./B(j,i)).*CM_t(j,1) + (m_rna(j,1)./B(j,i)).*hub_ht(j,1);
    CM(1,1) = -9.9;
    Ap(j,i) = (pi()./4).*(rad_center(j,i)).^4 + 3*((pi()./4).*((rad_column(j,i))^4)) + 2.*(pi().*(rad_column(j,i))^2).*((dist_cc(j,i)/2).^2); %what about central column?

    KM(j,i) = CoB(j,i) + Ap(j,i)./Vdisp(j,i);

    Moment(j,i) = Thrust(j).*(hub_ht(j) - CM(j,i));
    Moment2(j,i) = Thrust(j).*(hub_ht(j) + draft(j,i) - KM(j,i));     
    Moment3(j,i) = Thrust(j).*hub_ht(j);
%     Moment3(j,i) = 2.*Thrust(j).*(hub_ht(j)+abs(CM(j,i))) + 2.*Thrust(j).*(hub_ht(j)+abs(CM(j,i))+2.*radius(j)+10);
    
    stiff1(j,i) = rhoW.*g.*Vdisp(j,i).*(-CoB(j,i) - (CM(j,i)));% + rhoW.*g.*Ap(j,i);
    stiff2(j,i) = rhoW.*g.*Ap(j,i);
%     Kp(j,i) = stiff2(j,i) + rhoW.*g.*Vdisp(j,i).*(-CoB(j,i) + (CoG(j,i)));
    Kequn(j,i) = (stiff1(j,i) + stiff2(j,i)); %
    thetaP(j,i) = ((Moment(j,i)./Kequn(j,i)).*(180/pi()));  %about CM
    thetaP2(j,i) = ((Moment2(j,i)./Kequn(j,i)).*(180/pi()));  %about metacenter  
    thetaP3(j,i) = (Moment3(j,i)./Kequn(j,i)).*(180./pi());     %about wl 
%     thetaP3(j,i) = ((Moment3(j,i)./Kequn(j,i)).*(180/pi()));  %round(,1)  
    ratio(j,i) = stiff1(j,i)./stiff2(j,i);
    
     
    i = i + 1;
end
end

% ratio = Mballast./(Mplatform + Mballast); 

%% MoI 

for m=1:nturbine
    % tower
    rhoTower = 8500;    %kg/m^3
    b = (radius(1)/radius(m)).^1.5;
    tower_r = 2.6./b;
    tower_t = 0.023./b;
    tower_ht = hub_ht(m) - cfreeb(m,1);
    for n = 1:length(alpha)   
    d1 = abs(CG1(m,n) + CM(m,n));    %upper walls
    d2 = abs(CG2(m,n) + CM(m,n));    %lower walls
    d3 = abs((CG3(m,n)+t(m,n)./2) + CM(m,n));    %upper hp
    d4 = abs((CG4(m,n)+t(m,n)./2) + CM(m,n));    %lower hp 
    d5 = abs((CG5(m,n)+t(m,n)./2) + CM(m,n));    %top cap 
    d6 = abs(CG6(m,n) + CM(m,n));    %upper ballast
    d7 = abs(CG7(m,n) + CM(m,n));    %lower ballast
    d8 = abs(CG8(m,n) + CM(m,n));    %center 
    % mass
    m1 = Vcol(m,n).*rhoM;       %kg  %including some cross brace mass
    m2 = Vhpwalls(m,n).*rhoM;
    m3 = (Vhp(m,n)*rhoM)./2;                          % heave plate upper & lower
    m4 = m3;
    m5 = Vtop(m,n).*rhoM;                        %top cap
    m6 = (0.28.*Mballast(m,n))./3;                       %kg %water ballast upper 
    m7 = (0.72.*Mballast(m,n))./3;                     %kg
    m8 = Vcenter(m,n).*rhoM;          %central column
    MM(m,n) = m1 + m2 + m3 + m4 + m5 + m6 + m7;         %mass of one column

    % Platform MoI 
    Iwalls1 = ((1/12).*m1.*(3.*(rad_column(m,n).^2 + (rad_column(m,n)-t(m,n)).^2) + (draft(m,n)-hp_thickness(m,n)+freeboard(m,n)).^2) + m1.*d1.^2);
    Iwalls2 = ((1/12).*m2.*(3.*(hp_rad(m,n).^2 + (hp_rad(m,n)-t(m,n)).^2) + (hp_thickness(m,n)).^2) + m2.*d2.^2);
    Icaps1 = ((1/12).*m3.*(3.*((hp_rad(m,n)-t(m,n)).^2) + 4.*(t(m,n).^2))) + m3.*(d3.^2); 
    Icaps2 = ((1/12).*m3.*(3.*((hp_rad(m,n)-t(m,n)).^2) + 4.*t(m,n).^2)) + m3.*(d4.^2);
    Icaps3 = ((1/12).*m5.*(3.*(rad_column(m,n)-t(m,n)).^2 + 4*t(m,n).^2)+m5.*(d5).^2);
    Iballast1 = ((1/12).*m6.*(3.*(rad_column(m,n)-t(m,n)).^2 + (b1(m,n)).^2) + m6.*(d6).^2);
    Iballast2 = ((1/12).*m7.*(3.*(hp_rad(m,n)-t(m,n)).^2 + (b2(m,n)-t(m,n)).^2) + m7.*d7.^2);
    Icenter = ((1/12).*m8.*(3.*(rad_center(m,n).^2 + (rad_center(m,n)-(t(m,n)./2)).^2) + (draft(m,n) + cfreeb(m,n)).^2)) + m8.*(d8).^2; %include middle column
    I_platform1 = 3.*(Iwalls1 + Iwalls2 + Icaps1 + Icaps2 + Icaps3 + Iballast1 + Iballast2) + Icenter;

    I_platform(m,n) = I_platform1 + 2.*MM(m,n).*(dist_cc(m,n)./2)^2; %
    % tower
    I_tower(m,n) = (((1/12).*m_t(m).*(3.*((tower_r.^2) + (tower_r-tower_t).^2) + 4.*(tower_ht).^2)) + m_t(m).*(cfreeb(m)+abs(CM(m,n))).^2);
    % RNA pt mass
    I_RNA(m,n) = m_rna(m).*(hub_ht(m)+abs(CM(m,n))).^2;
    % total MoI m_
    I_pitch(m,n) = I_platform(m,n) + I_tower(m,n) + I_RNA(m,n);
    % Wn
    Tn(m,n) = (sqrt(((I_pitch(m,n) + 0.63.*I_platform(m,n)))./Kequn(m,n))).*(2.*pi());
    wn(m,n) = (2.*pi())./Tn(m,n);       %rad/s
    n = n+1;
    end
    m = m+1;
end


%% plot 

figure(1) 
yline(thetaP3(1,1),'LineWidth',2);
hold on 
plot(alpha,thetaP3(2,:),'*','LineWidth',1);
plot(alpha,thetaP3(3,:),'*','LineWidth',1);
plot(alpha,thetaP3(4,:),'*','LineWidth',1);
xlabel('alpha')
ylabel('Pitch angle rated thrust (deg)')
axis([0 2 0 10])
legend('OC4','10MW','15MW','20MW') %'15MW','20MW')
set(gca,'FontSize',20)
hold off

figure(2)
yline(Tn(1,1),'LineWidth',2);
hold on 
plot(alpha,Tn(2,:),'*','LineWidth',1)
plot(alpha,Tn(3,:),'*','LineWidth',1)
plot(alpha,Tn(4,:),'*','LineWidth',1)
legend('OC4','10MW','15MW','20MW')
xlabel('alpha')
ylabel('Period (s)')
axis([0 1.5 0 60])
set(gca,'FontSize',20)
hold off 


% figure(3)
% plot(alpha,Mplatform,'*','LineWidth',2)
% legend('OC4','10MW','15MW','20MW')
% xlabel('alpha')
% ylabel('Steel mass (kg)')
% set(gca,'FontSize',20)

%%
%Mtotal 
%turbine 
%Mplatform 
%Mballast for the alpha 
