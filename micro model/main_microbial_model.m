% microbial model modified from Allison et al.(2010)

clear all;
close all;
format long e;
RandSeed=clock; %just for the rand change every time
rand('seed',RandSeed(6));
data =xlsread( 'F:\Teco\micro\D01_normal_year.xlsx');
GPP=data(:,7);
ER=data(:,6);
Tair=data(:,9);
SWC_all=data(:,12);
Tsoil=data(:,10);
RH=data(:,15);
for ii_yr = 2009:2017
    ii_d=(ii_yr-2008-1)*365+1;
    ii_d2=(ii_yr-2008)*365;
    Tair_yr=Tair(ii_d:ii_d2);
    RH_yr=RH(ii_d:ii_d2);
    SWC=SWC_all(ii_d:ii_d2);
    Tsoil_yr=Tsoil(ii_d:ii_d2);
    Soil_scalar_for        = soil_scalar_TECO(Tsoil_yr,SWC);% soil
    Soil_scalar(ii_d:ii_d2)=Soil_scalar_for;   
end 
%% initial data                                     
x_ss  = [0 0 95.59]; % inital leaf, root, litter
soilC = 6458.964802; % 100cm  
mbc =115.5921758;%
doc = 17.09311853;
ezc= 49.29;
%% parameters
%            r_death ;      cue_mbc0; r_ebcProd;    r_ebcLoss   ;vmax_soc0; f_I2soc;vmax_mbcUdoc0
para_0     = [0.0002*24        0.65    0.00008*24      0.001*24      1e8*24    0.5 ];%                      para_0     = [0.0002*24       0.63       0.000005*24       0.001*24      1e8*24    0.5]
para_0_max = [0.0002*24*3      0.87    0.00008*24*3    0.001*24*3    1e8*24*3    1 ];% test from 3-10          para_0_max = [0.0002*24*3  0.87   0.000005*24*3       0.001*24*3     1e8*24*3  1]
para_0_min = [0                0.43    0.00008*24/3    0.001*24/3    1e8*24/3    0 ]; % test from 0.3 - 0.1  para_0_min = [0                 0.43        0                 0          1e8*24/3  0]
para_max = [para_0_max];
para_min = [para_0_min];
para_old = para_min+(para_max-para_min).*rand(1,6);
% nr=1;
%% paras for MCMC
Nt         = 365;
J_old      = 3000000;
upgrade    = 0;
DJ         = std(ER);
%% fixed parameters
R = 0.008314;  %j/molK
cue_mbc_s=0.016;     f_2soc=0.5;    vmax_mbcUdoc0=1e8*24; %vmax_mbcUdoc0=1e8*24
km_soc0=500000;      ea_up=47;      ea_soc=47;
k_mbcUdoc0=0.1*1000; km_soc_s=5000; k_mbcUdoc_slop=10;  %k_mbcUdoc_slop=10
simu_train_num = [1:365*9];
simu_train_end = 365*9;
para_keep      =zeros(6,20000);
J_keep         =zeros(1,20000);
res_record_co2 =zeros(365*9,20000);
res_record_mbc =zeros(5,20000);
res_record_ebc =zeros(5,20000);
res_record_doc =zeros(5,20000);
%% start simulation
for nrun = 1:10000000000
    para_new  = Generate_C(para_old,para_min,para_max);
    cpools    = [x_ss  mbc ezc doc soilC];
    
    predict_c = micro_process(Nt,GPP,Tair,RH,Tsoil,SWC_all,Soil_scalar,para_new,cpools);
    if predict_c ==0
        continue;
    else
        simu_co2  = predict_c(11,:);
        J         = (norm(simu_co2(simu_train_num)-ER')).^2;

        J_new     = J/(2*DJ);
        delta_J    = J_new - J_old;

        if  min(1,exp(-delta_J))>rand
            upgrade=upgrade+1;
            para_keep(:,upgrade) = para_new;
            J_keep(upgrade)=J_new;
            res_record_co2(:,upgrade) = simu_co2;
            res_record_mbc(:,upgrade) = predict_c(4,simu_train_end);
            res_record_ebc(:,upgrade) = predict_c(5,simu_train_end);
            res_record_doc(:,upgrade) = predict_c(6,simu_train_end);
            %             res_record_mbc(:,upgrade) = predict_c(4,simu_train_end);
            para_old=para_new;
            J_old=J_new;
            %             J_co_old = J_co_new;
        end
        if upgrade == 20000
            break;
        end
    end
    disp(['simu=', num2str(nrun), '; upgrade =', num2str(upgrade), '; J_new = ', num2str(J_new), '; J_old =', num2str(J_old)]);

end
para_keep_shuzhi=para_keep';
co2=res_record_co2';
mbc=res_record_mbc';
ebc=res_record_mbc';
doc=res_record_doc';


