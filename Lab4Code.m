clc
clear
close all
%define constants
Patm = 30.22*3386.38;% [Pa]
Tatm = 297.0389;% [K]
R = 287.05;% [J/Kmol*K]
rho = Patm/(R*Tatm);% [Kg/m^3]
vel = 20;
pdyn = 0.5*rho*vel^2;
dp_0 = 0.965*249.08;% [Pa]
c = 6*0.0254;% [m]
confidence = 1.96;
alpha = [-3,0,3,6,9,12];

%load .mat files
datapath1 = 'C:\Users\jwmay\Documents\Semester 6\AE Lab\Lab4';
label = ["0","3","6","9","12","neg3"];
for i = 1:length(label)
    for j = 1:77
        filename = fullfile(datapath1,sprintf('AOA_%s_wake',label(i)),sprintf('AOA_%s_%d',label(i),j));
        temp_patm = cell2mat(struct2cell(load(filename,"pAtm")));
        temp_tatm = cell2mat(struct2cell(load(filename,"tAtm")));
        temp_zCurr = cell2mat(struct2cell(load(filename,"zCurr")));
        temp_dp = cell2mat(struct2cell(load(filename,"dp")));
        dp_avg = sum(temp_dp)/length(temp_dp);
        dp_er = confidence*std(temp_dp)/sqrt(length(temp_dp));
        zCurrData(j,i) = temp_zCurr.*0.0254;% z-position [m]
        dpData(j,i) = dp_avg.*249.08; % dp [Pa] @AOA: 0,3,6,9,12 & -3 deg
        dpEr(j,i) = dp_er.*249.08; % dp error [Pa]
    end
end

%load .csv files
for i = 1:length(label)
    temp_ptap = readmatrix(sprintf("airfoil_aoa%s.csv",label(i)),"Range","C2:T401");
    ptap(:,i) = 249.08.*sum(temp_ptap)./(size(temp_ptap,1));
    E_ptap(:,i) = 249.08*confidence*std(temp_ptap)/sqrt(size(temp_ptap,1));
end

%load panel code data
datapath2 = 'C:\Users\jwmay\Documents\Semester 6\AE Lab\Lab4\PotentialFlow';
panel_xl = cell2mat(struct2cell(load(fullfile(datapath2,'xl.mat'))));
panel_xu = cell2mat(struct2cell(load(fullfile(datapath2,'xu.mat'))));
for i = 1:length(label)
    panel_cl(i) = cell2mat(struct2cell(load(fullfile(datapath2,sprintf('Cl_%saoa.mat',label(i))))));
    panel_cmcby4(i) = cell2mat(struct2cell(load(fullfile(datapath2,sprintf('CmcBy4_%saoa.mat',label(i))))));
    panel_cmLE(i) = cell2mat(struct2cell(load(fullfile(datapath2,sprintf('CmLE_%saoa.mat',label(i))))));
    panel_cpl(i) = struct2cell(load(fullfile(datapath2,sprintf('Cpl_%saoa.mat',label(i)))));
    panel_cpu(i) = struct2cell(load(fullfile(datapath2,sprintf('Cpu_%saoa.mat',label(i)))));
end

%calculate Cp
cp = ptap./pdyn;
cp_upper = cp(1:10,:);
cp_lower = cat(1,cp(1,:),cp(11:18,:));

%ambient uncertainty
E_P = 0.005*249.08;
E_T = 0.05*(5/9);
E_q = 0.005*249.08;

%{
%propagated dynamic perssure uncertainty
syms P T
fun1(P,T) = 0.5*P*vel^2/(R*T);
delP = diff(fun1,P);
delT = diff(fun1,T);
E_q = double(sqrt((delP(Patm,Tatm)*E_P)^2 + (delT(Patm,Tatm)*E_T)^2));
%}

%propagated Cp uncertainty
syms q p
fun2(q,p) = p/(q);
delq = diff(fun2,q);
delp = diff(fun2,p);
for i = 1:size(ptap,2)
    for j = 1:size(ptap,1)
        E_cp(j,i) = double(sqrt((delq(pdyn,ptap(j,i))*E_q)^2 + (delp(pdyn,ptap(j,i))*E_ptap(j,i))^2));
    end    
end

E_cp_upper = E_cp(1:10,:);
E_cp_lower = cat(1,E_cp(1,:),E_cp(11:18,:));


%plot cp vs x/c for each AOA
x_upper = c*[0,0.15,0.3,0.6,1.2,1.7998,2.4001,3,4.2,5.4];
x_lower = c*[0,0.24,0.36,0.6,1.2,1.8,3.0,4.2,5.15];
x_target = 1;
te = (interp1(x_upper,cp_upper,x_target,'linear','extrap')+interp1(x_lower,cp_lower,x_target,'linear','extrap'))./2;
cp_upper = cat(1,cp_upper,te);
cp_lower = cat(1,cp_lower,te);
E_cp_upper = cat(1,E_cp_upper,zeros(1,6));
E_cp_lower = cat(1,E_cp_lower,zeros(1,6));
x_upper = cat(2,x_upper,x_target);
x_lower = cat(2,x_lower,x_target);

plot_label = ["0","3","6","9","12","-3"];
for i = 1:length(plot_label)
    figure
    grid on
    title(sprintf("%s = %s°",texlabel('alpha'),plot_label(i)))
    ylabel(texlabel('C_p'),'FontWeight','bold')
    xlabel(texlabel('x/c'),'FontWeight','bold')
    set(gca,'YDir','reverse')
    hold on
    errorbar(x_upper,cp_upper(:,i),E_cp_upper(:,i),'o-','linewidth',2,'color',[0.8500 0.3250 0.0980])
    errorbar(x_lower,cp_lower(:,i),E_cp_lower(:,i),'o-','linewidth',2,'color',[0 0.4470 0.7410])
    plot(panel_xu,panel_cpu{i},'--r',panel_xl,panel_cpl{i},'--b');
    legend('Upper Surface','Lower Surface','Upper Surface Potential Flow','Lower Surface Potential Flow')
    hold off
end
%calculate Cl
cl = trapz(x_lower,cp_lower)-trapz(x_upper,cp_upper);
cl = cat(2,cl(6),cl(1:5));%sort values -3:12
panel_cl = cat(2,panel_cl(6),panel_cl(1:5));%sort values -3:12

%calculate CmLE
for i = 1:length(alpha)
    CmLE(i) = -1*(trapz(x_lower',cp_lower(:,i).*x_lower')-trapz(x_upper',cp_upper(:,i).*x_upper'));
end

%calculate Cmcby4
for i = 1:length(alpha)
    Cmcby4(i) = -1*(trapz(x_lower',cp_lower(:,i).*(x_lower'-0.25))-trapz(x_upper',cp_upper(:,i).*(x_upper'-0.25)));
end

%plot Cl vs AoA
figure
hold on
grid on
plot(alpha,cl,'o-','linewidth',2,'color',[0.8500 0.3250 0.0980])
plot(alpha,panel_cl,'r--','linewidth',2)
xlabel(texlabel('alpha [deg]'),'FontWeight','bold')
ylabel(texlabel('C_l'),'FontWeight','bold')
title(texlabel('C_l vs alpha'))
legend('Experimental','Potential Flow')
hold off

%plot CmLE vs alpha
figure
hold on
grid on
plot(alpha,CmLE,'o-','linewidth',2,'color',[0.8500 0.3250 0.0980])
%plot(alpha,panel_cmLE,'r--','linewidth',2)
xlabel(texlabel('alpha [deg]'),'FontWeight','bold')
ylabel(texlabel('Cm_L_E'),'FontWeight','bold')
title(texlabel('Cm_L_E vs alpha'))
%legend('Experimental','Potential Flow')
hold off

%plot Cmcby4 vs alpha
figure
hold on
grid on
plot(alpha,Cmcby4,'o-','linewidth',2,'color',[0.8500 0.3250 0.0980])
%plot(alpha,panel_cmcby4,'r--','linewidth',2)
xlabel(texlabel('alpha [deg]'),'FontWeight','bold')
ylabel(texlabel('Cm_c_/_4'),'FontWeight','bold')
title(texlabel('Cm_c_/_4 vs alpha'))
%legend('Experimental','Potential Flow')
hold off

%calculate wake velocity @AOA: 0,3,6,9,12, & -3 deg
v_wake = sqrt(2.*dpData./rho);
parse_wake = double(isoutlier(v_wake));% this finds the location of wake datapoints using the isoutlier function
%plot(v_wake(:,1),zCurrData(:,1),'o')

%calculate drag and Cd values @AOA: 0,3,6,9,12, & -3 deg
for i = 1:size(v_wake,2)
     parse_vel = v_wake(:,i)-(parse_wake(:,i).*v_wake(:,i));% this parses the free stream velocities from the wake velocities
     v_inf(i) = sum(parse_vel)/length(nonzeros(parse_vel));% averages the free stream velocites
    for j = 1:size(v_wake,1)
        intd(j,i) = v_wake(j,i)*(v_inf(i)-v_wake(j,i));
    end
    d(i) = rho*trapz(zCurrData(:,i),(intd(:,i)));% this integrates the v_wake(v_inf-v_wake)) at each z-position
    cd(i) = d(i)./(.5*rho*c.*v_inf(i)^2);
end
cd = cat(2,cd(6),cd(1:5));

%plot Cd vs AoA
figure
hold on
grid on
plot(alpha,cd,'o','linewidth',2,'color',[0.8500 0.3250 0.0980])
plot(alpha,zeros(1,6),'r--','linewidth',2)
ylim([-0.01 0.04])
xlabel(texlabel('alpha [deg]'),'FontWeight','bold')
ylabel(texlabel('C_d'),'FontWeight','bold')
title(texlabel('C_d vs alpha'))
legend('Experimental','Potential Flow')
hold off

%plot Cd vs Cl
figure
hold on
grid on
plot(cd,cl,'o','linewidth',2,'color',[0.8500 0.3250 0.0980])
plot(zeros(1,6),panel_cl,'r--','linewidth',2)
xlim([-0.005,0.04])
xlabel(texlabel('C_d'),'FontWeight','bold')
ylabel(texlabel('C_l'),'FontWeight','bold')
title('Drag Polar','FontWeight','bold')
legend('Experimental','Potential Flow')
hold off

%plot L/D vs alpha (note that potential flow has an infinite L/D due to no drag)
figure
hold on
grid on
plot(alpha,cl./cd,'o','linewidth',2,'color',[0.8500 0.3250 0.0980])
xlabel(texlabel('alpha [deg]'),'FontWeight','bold')
ylabel(texlabel('L/D'),'FontWeight','bold')
title(texlabel('L/D vs alpha'))
%legend('Experimental')
hold off

%calculate zero lift angle of attack
alpha_zl = interp1(cl,alpha,0,'linear','extrap');
alpha_zl_panel = interp1(panel_cl,alpha,0,'linear','extrap');
fprintf('zero lift angle of attack: %.3f [deg]\n',alpha_zl)
fprintf('potential flow zero lift angle of attack: %.3f [deg]\n',alpha_zl_panel)

%{
 based on http://airfoiltools.com/airfoil/details?airfoil=naca4412-il the
 naca 4412 airfoil should have a zero lift AoA of ~4.9 deg.
%}

hold off