% conventional model with differents water scalar functions
%%
clear all;
close all;
format long e;
RandSeed=clock; %just for the rand change every time
rand('seed',RandSeed(6));
data =xlsread( 'your_path\forcing_data.xlsx');
GPP=data(:,7);
ER=data(:,6);
Tair=data(:,9);
SWC_all=data(:,12);
Tsoil=data(:,10);
RH=data(:,15);
%% initial data
x_ss  =[0  0  95.59]; % inital leaf, root, litter x_ss  =  [0 54.5 788.03],litter 100,litter 200   181
%soilC =  6614.236978; % 100cm soilC = 10076.97
soilC = 6458.964802;
%% parameters
para_A_Lit    = [0.25 0.1]; % litter to fast and slow soil C pools  
para_A0       = [0.33 0.0051 0.4185 0.01 0.50]; % f2s,f2p,s2f,s2p,p2f  
para_A0_max   =  [0.5 0.01 0.6 0.02 0.7];% 
para_A0_min   =  [0.1 0    0.1 0    0.3];% 
matrix_B      = [0.55*0.39/1.39 0.55/1.39 0 0 0 0];   % 
para_C_a      = [1/(0.75*365) 1/(0.6*365)   1/(3.1*365)]; % turnover rate of leaf, root and litter per day 
para_C_s0     = [1.44/365       0.2701/365    0.005548/365]; % turnover rate of soil C per day para_C_s0   
para_C_s0_max = [1.44*3/365     0.511/365      0.073/365];  
para_C_s0_min = [1.44/3/365     0.10001/365   0.0000365/365];  
para_opt0  = [para_A0     para_C_s0     ]; 
para_max   = [para_A0_max para_C_s0_max ];
para_min   = [para_A0_min para_C_s0_min ];
%% scalar
above_scalar=ones(365*9,1);
ab_scalar=ones(365*9,1);
Soil_scalar=ones(365*9,1);
for ii_yr = 2009:2017
    ii_d=(ii_yr-2008-1)*365+1;
    ii_d2=(ii_yr-2008)*365;
    Tair_yr=Tair(ii_d:ii_d2);
    RH_yr=RH(ii_d:ii_d2);
    SWC=SWC_all(ii_d:ii_d2);
    Tsoil_yr=Tsoil(ii_d:ii_d2);
    above_scalar_for       = (2.^((Tair_yr-10)./10)).*(RH_yr/max(RH_yr)); %leaf
    ab_scalar_for          = (2.^((Tsoil_yr-10)./10)).*(SWC/max(SWC));%litter and root
    Soil_scalar_for        = soil_scalar_TECO(Tsoil_yr,SWC);% soil
    above_scalar(ii_d:ii_d2)=above_scalar_for;
    ab_scalar(ii_d:ii_d2)=ab_scalar_for;
    Soil_scalar(ii_d:ii_d2)=Soil_scalar_for;   
end 
 scalar_all = [above_scalar';ab_scalar';ab_scalar';Soil_scalar';Soil_scalar';Soil_scalar'];   


%% paras for MCMC
para_old   = para_min+rand(1,8).*(para_max-para_min);
nsimu      = 100000;
Nt         = 365;
J_old      = 3000;
upgrade    = 0;
DJ         = std(ER);
%% start simulation

para_keep      = zeros(8,20000);
J_keep         = zeros(1,20000);
res_record_co2 = zeros(365*9,20000);
res_record_scf = zeros(365*9,20000);
res_record_scs = zeros(365*9,20000);
res_record_scp = zeros(365*9,20000);

for simu = 1:1000000  %500000
    para_new  = Generate_C(para_old,para_min,para_max);
    para_A    = [para_A_Lit para_new(1:5)]; % litter to soc, soc to soc
    para_C    = [para_C_a  para_new(6:8)];  % leaf,root,litter; soc fast,slow,passive
    soil_c_f  = 0.029915881*soilC; % *soilC
    soil_c_s  = 0.27*soilC;
    
    soil_c_p  = (1-0.029915881-0.27)*soilC;
    x0        = [x_ss soil_c_f soil_c_s soil_c_p];
    %% conventional soil carbon model
    simu_res  = solve_forward(scalar_all,GPP,matrix_B,para_A,para_C,x0,Nt); % 6Cpools,3co2,sumCO2
    %%
    simu_co2  = simu_res(12,:);
    %         J         = (norm(simu_co2(simu_train_num)-SRH(data_training,1)')).^2;
    J         = (norm(simu_co2-ER')).^2;
    J_new     = J/(2*DJ);
    delta_J    = J_new - J_old;
    if min(1,exp(-delta_J))>rand 
        upgrade=upgrade+1;
        para_keep(:,upgrade) = para_new;
        J_keep(upgrade)=J_new;
        res_record_co2(:,upgrade) = simu_co2;
        res_record_scf(:,upgrade) = simu_res(4,:);
        res_record_scs(:,upgrade) = simu_res(5,:);
        res_record_scp(:,upgrade) = simu_res(6,:);
        para_old=para_new;
        J_old=J_new;
    end
    
    if upgrade == 20000
        break;
    end
    disp(['simu=', num2str(simu), '; upgrade =', num2str(upgrade), '; J_new = ', num2str(J_new), '; J_old =', num2str(J_old)]);
end
para_keep_shuzhi=para_keep';
