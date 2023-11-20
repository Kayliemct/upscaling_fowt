%% 15MW IEA Ref Turbine Upscale 

clear all; close all; clc


Power = [15; 20; 25; 30];  %rated power 
Sp = 332;   %W/m^2
nturbine = 4;

for j=1:nturbine
    radius(j) = round((sqrt((Power(j).*10^6)./(pi().*Sp))));
    radiusCHECK(j) = ((sqrt((Power(j).*10^6)./(pi().*Sp))));
    A(j) = pi().*radius(j).^2;
end
diameter = radius.*2;

rhoA = 1.225;
rhoW = 1025;
rhoS = 7850;
rhoORE = 4800;  %2500; %4300;                              %kg/m^3 
g = 9.81;

%% Turbine upscale 
alphaT = 2.2;%3;         %alpha for RNA
alphatow = 2;       %alpha for tower 
Ct = 0.6;      %0.8?     %thrust coefficient  
speed_rated = 10.59;     

for j=1:nturbine
    hub_ht1(j,1) = radius(j) + 30;  %1.25.*radius(j);  % updated to have 25% clearance, good for wind shear? 30m clearance from rotor to waterline
    m_rna(j,1) = (1017e3)./(radius(1)/radius(j)).^alphaT;       %15MW floating RNA mass
    m_t(j,1) = (1263e3)./(radius(1)/radius(j)).^alphatow;       %15MW floating tower mass
    turbine(j,1) = m_t(j)+m_rna(j);
    Thrust(j,1) = 0.5.*rhoA.*A(j).*speed_rated.^2.*Ct;          %N
    CM_t(j,1) = hub_ht1(j,1).*0.277;   %CM of tower is 28% hub height, this is approximate
end 

%% Platform Upscale 
alpha = [0:0.005:1.5];           %alpha for the platform

draft3 = 20;                                %m
t3 = 0.045;  %0.045                                 %m approx wall thickness 
r_columns3 = 6.25;                          %m
r_center3 = 5;                              %m
pontoon_approx_l3 = 42.85;                  %m approx length of pontoon
pontoon_width3 = 12.5;                      %m width of pontoon 
pontoon_height3 = 7;                        %m height of pontoon 
free3 = 15;                                 %m pontoon freeboard
col_spacing3 = 51.75;                       %m length from center to column 
angle = 60;                                 %deg this is the angle from center to pontoon about x axis 

% upscale parameters w/ alpha 
for j=1:nturbine 
for i= 1:length(alpha)
    scaling(j,i) = (radius(1)/radius(j)).^alpha(i);
    draft(j,i) = draft3./scaling(j,i);                         %m draft is changing
    t(j,i) = t3; %./scaling(j,i);   %t3                                          %m wall thickness is constant, increase based on hydrostatic pressure?              
    r_columns(j,i) = r_columns3./scaling(j,i);                          %m radius of outer columns
    r_center(j,i) = r_center3./scaling(j,i);                            %m radius of central column
    pontoon_approx_l(j,i) = pontoon_approx_l3./scaling(j,i);            %m
    pontoon_width(j,i) = pontoon_width3./scaling(j,i);                  %m
    pontoon_height(j,i) = pontoon_height3./scaling(j,i);                %m
    free(j,i) = free3./scaling(j,i);                                    %m
    h(j,i) = draft(j,i) + free(j,i);                                    %m total height of platform (draft + freeboard)
    col_spacing(j,i) = col_spacing3./scaling(j,i);                      %m
end 
end

for j=1:nturbine 
for i= 1:length(alpha)
    %volume 
    V1(j,i) = ((pi().*((r_columns(j,i)).^2 - (r_columns(j,i)-t(j,i)).^2)).*h(j,i));         % volume of column walls 
    V2(j,i) = (pi().*t(j,i).*((r_columns(j,i)-t(j,i)).^2));                                 % volume of column top & base 
    V3(j,i) = pi().*(r_center(j,i).^2 - (r_center(j,i) - t(j,i)).^2).*h(j,i);               % volume of center walls
    V4(j,i) = pi().*(r_center(j,i) - t(j,i)).^2.*t(j,i);                                    % volume of center base 
%     V5(j,i) = (pontoon_approx_l(j,i).*(pontoon_width(j,i) - 2.*t(j,i)).*t(j,i));            % volume of pontoon top & base
%     V6(j,i) = (pontoon_approx_l(j,i)).*pontoon_height(j,i).*t(j,i);       
    V7(j,i) = (pontoon_approx_l(j,i).*pontoon_height(j,i).*pontoon_width(j,i)) - pontoon_approx_l(j,i).*(pontoon_width(j,i) - 2.*t(j,i)).*(pontoon_height(j,i) - 2.*t(j,i));
    % volume of pontoon sides 
    V_columns(j,i) = 3.*(V1(j,i) + 2.*V2(j,i));
    V_center(j,i) = V3(j,i) + V4(j,i);
    V_pontoon(j,i) = 3.*(V7(j,i));
    Vplatform(j,i) = V_pontoon(j,i) + V_columns(j,i) + V_center(j,i);           % volume of the platform [m^3]
    Mplatform(j,i) = Vplatform(j,i).*rhoS;                                      % mass of platform [kg]
    %displaced volume 
    Vd_center(j,i) = pi().*r_center(j,i).^2.*(draft(j,i));              % displaced volume center
    Vd_pontoons(j,i) = 3.*(pontoon_approx_l(j,i).*pontoon_width(j,i).*pontoon_height(j,i)); % displaced volume pontoons 
    Vd_columns(j,i) = 3.*(pi().*r_columns(j,i).^2.*draft(j,i));                             % displaced volume outer columns 
    Vdisp(j,i) = Vd_center(j,i) + Vd_pontoons(j,i) + Vd_columns(j,i);                       % total displace folume 
    %buoyancy 
    B(j,i) = Vdisp(j,i).*rhoW;
    Mtotal(j,i) = B(j,i) - turbine(j);                      % total mass must match the buoyancy 
    Mballast(j,i) = Mtotal(j,i) - Mplatform(j,i);
%     Mtest(j,i) = B(j,i) - Mballast(j,i);
    %seawater ballast
    Wballast(j,i) = 3.*(pontoon_approx_l(j,i).*(pontoon_width(j,i)-2.*t(j,i)).*(pontoon_height(j,i)./2)).*rhoW;    %
    Wratio(j,i) = Wballast(j,i)./Mballast(j,i);
%     Wballast(j,i) = 0.815.*Mballast(j,i);                   % water ballast is 81.5% of total ballast for 15MW
%     h_w_b(j,i) = (Wballast(j,i)./(3.*rhoW))./(pontoon_approx_l(j,i).*(pontoon_width(j,i) - 2.*t(j,i)));
    %iron ore ballast
    Oballast(j,i) = Mballast(j,i) - Wballast(j,i);          % iron ore ballast is 18.5% of total ballast for 15MW
    %CoB
    CoB(j,i) = (draft(j,i)./2).*((Vd_center(j,i) + Vd_columns(j,i))./Vdisp(j,i)) + (draft(j,i) - pontoon_height(j,i)./2).*(Vd_pontoons(j,i)./Vdisp(j,i));
    %CoG
    CoG_columns(j,i) = (draft(j,i))./2;      %(draft - freeboard)./2?       %outer columns & central column 
    CoG_pontoons(j,i) = draft(j,i) - pontoon_height(j,i)./2; 
    CG_steel(j,i) = CoG_columns(j,i).*((V_columns(j,i) + V_center(j,i))./Vplatform(j,i)) + CoG_pontoons(j,i).*(V_pontoon(j,i)./Vplatform(j,i));
    CG_water(j,i) = draft(j,i) - pontoon_height(j,i)./4;      %draft(j,i) - t(j,i) - h_w_b(j,i)./2;
    h0(j,i) = Oballast(j,i)./(3.*pi().*(r_columns(j,i) - t(j,i)).^2.*rhoORE);       % height of iron ore ballast [m]
    CG_ore(j,i) = draft(j,i) - t(j,i) - h0(j,i)./2;
    CoG(j,i) = CG_steel(j,i).*(Mplatform(j,i)./Mtotal(j,i)) + CG_water(j,i).*(Wballast(j,i)./Mtotal(j,i)) + CG_ore(j,i).*(Oballast(j,i)./Mtotal(j,i));
    CM(j,i) = -CoG(j,i).*(Mtotal(j,i)./B(j,i)) + (hub_ht1(j)).*(m_rna(j)./B(j,i)) + (CM_t(j)).*(m_t(j)./B(j,i));    %CM_t + free?
%     CMtest(j,i) = -CG_steel(j,i).*(Mplatform(j,i)./Mtest(j,i))+(hub_ht1(j)).*(m_rna(j)./Mtest(j,i)) + (CM_t(j) + free(j,i)).*(m_t(j)./Mtest(j,i));
    %K
    A55(j,i) = (pi()./4).*(r_center(j,i)).^4 + 3.*(pi()./4).*(r_columns(j,i)).^4 + 2.*pi().*(r_columns(j,i)).^2.*(col_spacing(j,i).*sin(angle.*(pi()./180))).^2; %
    KM(j,i) = CoB(j,i) + A55(j,i)./Vdisp(j,i);
    K1(j,i) = rhoW.*g.*Vdisp(j,i).*(-CoB(j,i) - CM(j,i));
    K1P(j,i) = rhoW.*g.*Vdisp(j,i).*(-CoB(j,i) + CoG(j,i)); 
    K2(j,i) = rhoW.*g.*A55(j,i);
    ratio(j,i) = K1(j,i)./K2(j,i);  %overturning relative to stabilizing part 
    K_CM(j,i) = K1(j,i) + K2(j,i);
%     Kplatform(j,i) = K1P(j,i) + K2(j,i);
    %thetaP 
    M(j,i) = Thrust(j).*(hub_ht1(j) - CM(j,i));          %below the waterline is negative, add with hub height
    M2(j,i) = Thrust(j).*(hub_ht1(j) + draft(j,i) - KM(j,i));
    M3(j,i) = Thrust(j).*(hub_ht1(j));  %about waterline
    thetaP(j,i) = (M(j,i)./K_CM(j,i)).*(180./pi());     %about CM               % 'static pitch' - pitch angle from rated thrust 
    thetaP2(j,i) = (M2(j,i)./K_CM(j,i)).*(180./pi());    %about metacenter
    thetaP3(j,i) = (M3(j,i)./K_CM(j,i)).*(180./pi());
end 
end

%% simple math 

% Ksimple = K2(3,1);
% Msimple = Thrust(3).*(hub_ht1(3) - CM(3,1));
% pitch_simple = (Msimple./Ksimple).*(180./pi())

%%
figure(1) 
yline(thetaP3(1,1),'LineWidth',2);
hold on 
plot(alpha,thetaP3(2,:),'*','LineWidth',1)
plot(alpha,thetaP3(3,:),'*','LineWidth',1)
plot(alpha,thetaP3(4,:),'*','LineWidth',1)
axis([0 1.5 0 20])
legend('15MW IEA','20MW','25MW','30MW')
xlabel('alpha')
ylabel('Pitch angle rated thrust (deg)')
set(gca,'FontSize',20)
hold off


%% Moment of Inertia about CM 
% ***update MoI***

for j=1:nturbine 
for i= 1:length(alpha)
    d1(j,i) = ((-draft(j,i)+free(j,i))./2 - CM(j,i));           % center of columns relative to CM [m]
    d2(j,i) = ((free(j,i) - t(j,i)./2) - CM(j,i));              % center of top cap relative to CM [m]
    d3(j,i) = ((-draft(j,i) - t(j,i)./2) - CM(j,i));            % center of base cap relative to CM [m]
    I12(j,i) = (1/12).*3.*V1(j,i).*rhoS.*((3.*(r_columns(j,i).^2 + (r_columns(j,i) - t(j,i)).^2) + (h(j,i)).^2)) + 3.*V1(j,i).*rhoS.*(d1(j,i)).^2;  % MoI of outer columns
    I22(j,i) = (1/12).*3.*V2(j,i).*rhoS.*(3.*(r_columns(j,i)-t(j,i)).^2 + t(j,i).^2) + 3.*V2(j,i).*rhoS.*(d2(j,i)).^2;                              % MoI of top cap
    I32(j,i) = (1/12).*3.*V2(j,i).*rhoS.*(3.*(r_columns(j,i)-t(j,i)).^2 + t(j,i).^2) + 3.*V2(j,i).*rhoS.*(d3(j,i)).^2;                              % MoI of base cap 
    % I center 
    I42(j,i) = (1/12).*V3(j,i).*rhoS.*(3.*(r_center(j,i).^2 + (r_center(j,i) - t(j,i)).^2) + h(j,i).^2) + V3(j,i).*rhoS.*(d1(j,i)).^2;   %d4 = d1, MoI central column walls
    I52(j,i) = (1/12).*V4(j,i).*rhoS.*(3.*(r_center(j,i) - t(j,i)).^2 + t(j,i).^2) + V4(j,i).*rhoS.*(d3(j,i)).^2;                        %d5 = d3,  MoI base plate central column 
    % I ORE ballast 
    d6(j,i) = CG_ore(j,i) + CM(j,i);                % center of Ore ballast relative to CM
%     Oballast(3,1) = 2.54e6;
    I62(j,i) = (1/12).*Oballast(j,i).*(3.*(r_columns(j,i) - t(j,i)).^2 + h0(j,i).^2) + Oballast(j,i).*(d6(j,i)).^2;         % MoI Ore ballast
    % pontoon pt mass x4  
    Mp(j,i) = Wballast(j,i)./3 + V7(j,i).*rhoS;                  % mass of pontoon [kg]
    yp(j,i) = col_spacing(j,i).*sin(angle.*(pi()./180));                               % distance about y axis [m]
    CGp(j,i) = (CoG_pontoons(j,i)+CG_water(j,i))./2;                                   % CoG of combined pontoon walls with water ballast [m]
    Ip_cg(j,i) = 3.*Mp(j,i).*(CGp(j,i) - CM(j,i)).^2;                                  % MoI of pontoons
    Ip4(j,i) = (Mp(j,i)./4).*((yp(j,i)./8).^2 + ((3.*yp(j,i))./8).^2 + ((5.*yp(j,i))./8).^2 + ((7.*yp(j,i))./8).^2);                    % parallel axis thm pontoons (4 x point mass)
    % parallel axis thm
    p12(j,i) = (V1(j,i).*rhoS + V2(j,i).*2.*rhoS + Oballast(j,i)./3).*(col_spacing(j,i).*sin(angle.*(pi()./180))).^2;                   % parallel axis thm 2 outer columns 
    IplatformCM4(j,i) = (I12(j,i) + I22(j,i) + I32(j,i) + I42(j,i) + I52(j,i) + I62(j,i)) + 2.*p12(j,i) + Ip_cg(j,i) + 2.*Ip4(j,i);     % total MoI for platform 
    %MoI RNA & tower 
    Irna(j,i) = m_rna(j).*(hub_ht1(j) - CM(j,i)).^2;                                    % RNA MoI 
    Itower4(j,i) = ((m_t(j)./4)).*((CM_t(j)./4)+free(j,i)-CM(j,i)).^2 + (m_t(j)./4).*((3.*CM_t(j)./4)+free(j,i)-CM(j,i)).^2 + (m_t(j)./4).*((5.*CM_t(j)./4)+free(j,i)-CM(j,i)).^2 + (m_t(j)./4).*((7.*CM_t(j)./4)+free(j,i)-CM(j,i)).^2;        % tower MoI is four point masses 
    Itotal(j,i) = IplatformCM4(j,i) + Irna(j,i) + Itower4(j,i);% 
    % Natural period estimate 
    Tn(j,i) = (sqrt(((Itotal(j,i) + 0.63.*IplatformCM4(j,i)))./K_CM(j,i))).*(2.*pi()); % added mass term 0.63*MoI from report; %fn = 1./Tn; % wn = (2.*pi())./Tn;
end
end


%% plots 

figure(3)
yline(Tn(1,1),'LineWidth',2);
hold on 
plot(alpha,Tn(2,:),'*','LineWidth',1)
plot(alpha,Tn(3,:),'*','LineWidth',1)
plot(alpha,Tn(4,:),'*','LineWidth',1)
legend('15MW IEA','20MW','25MW','30MW')
xlabel('alpha')
ylabel('Period (s)')
axis([0 1.5 0 60])
set(gca,'FontSize',20)
hold off

% add plots for mass, Tn for four cases [20, 2x10MW, 3x6.7MW, 4x5MW]
