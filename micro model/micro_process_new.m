%% microbial soil process
function func_output = micro_process(Nt,GPP,Ta,RH,Ts,SWC,scalar_tem_mosit,parameters,cpools)
%% fixed parameters
R = 0.008314;  %j/molK
cue_mbc_s=0.016;     f_2soc=0.5;    vmax_mbcUdoc0=1e8*24;
km_soc0=500000;      ea_up=47;      ea_soc=47;
k_mbcUdoc0=0.1*1000; km_soc_s=5000; k_mbcUdoc_slop=10;
%% read parameters
r_death   =parameters(1); cue_mbc0  =parameters(2);
r_ebcProd =parameters(3); r_ebcLoss =parameters(4);
vmax_soc0 =parameters(5); f_I2soc   =parameters(6);
%% output
leafc = 0;
rootc = cpools(2);
litc  = cpools(3);
mbc =cpools(4);
ebc =cpools(5);
doc =cpools(6);
soc =cpools(7);
res = zeros(11,Nt*9);
    n_step = 0;
    for i_yr = 2009:2017
    i_d=(i_yr-2008-1)*365+1;
    i_d2=(i_yr-2008)*365;
    for i=i_d:i_d2
        %% input
        out_leaf = leafc*(1/(0.75*365))*(2.^((Ta(i)-10)./10))*(RH(i)/max(RH));
        out_root = rootc*(1/(0.6*365))*(2.^((Ts(i)-10)./10))*(SWC(i)/max(SWC));
        out_lit  = litc *(1/(3.1*365))*(2.^((Ts(i)-10)./10))*(SWC(i)/max(SWC));
        ra=GPP(i)*0.45;    %
        leafc = leafc+GPP(i)*0.55*0.39/1.39 - out_leaf;  

        rootc = rootc+GPP(i)*0.55/1.39 - out_root;   
        litc  = litc+out_leaf+out_root-out_lit;
        if i== 302||i== 667||i== 1032||i== 1397||i== 1762||i== 2127||i== 2492||i== 2857||i== 3222||i== 3587||i== 3952||i== 4317||i== 4682
            litc  = litc + leafc;
            leafc = 0;
        end
        %% microbial process
   
        I_soc= f_I2soc*out_lit; % the portion of Input to soc
        I_doc= (1-f_I2soc)*out_lit; % the portion of Input to doc
        
        soc = soc + I_soc;
        doc = doc + I_doc;
        %% set temperature upper and ll
        if Ts(i)>31
            tem_vk = 31;
        elseif Ts(i)<0
            tem_vk =0;
        else
            tem_vk = Ts(i);
        end
        %% delete the wrong parameters
        cue_mbc=-cue_mbc_s*tem_vk+cue_mbc0;
        if cue_mbc>1
            cue_mbc=0.9999;
        end
        if cue_mbc <0 || cue_mbc ==0
            cue_mbc = 0.0001;
        end
        %% doc to mbc
        vmax_mbcUdoc=vmax_mbcUdoc0*exp(-ea_up/(R*(tem_vk+273.15)));
        k_mbcUdoc=k_mbcUdoc_slop*tem_vk+k_mbcUdoc0;
        
        f_mbcUdoc_r=scalar_tem_mosit(i)*vmax_mbcUdoc*mbc/(k_mbcUdoc+doc);
        if f_mbcUdoc_r>1
            f_mbcUdoc_r = 0.9999;
        end
        f_mbcUdoc=f_mbcUdoc_r*doc;
        
        mbc = mbc + f_mbcUdoc*cue_mbc;
        f_prod2e=r_ebcProd*mbc;
        ebc=ebc+f_prod2e;
        %% for soc
        vmax_soc=vmax_soc0*exp(-ea_soc/(R*(tem_vk+273.15)));
        k_soc=km_soc_s*tem_vk+km_soc0;
        f_decomp_r=scalar_tem_mosit(i)*vmax_soc*ebc/(k_soc+soc);
        if f_decomp_r>1
            f_decomp_r = 0.9999;
        end
        f_decomp=f_decomp_r*soc;
        %% for mbc and ebc loss
        f_death=r_death*mbc;
        ebc_loss=r_ebcLoss*ebc;
        mbc=mbc-f_death-f_prod2e;
        ebc=ebc-ebc_loss;
        %%
        soc=soc+f_death*f_2soc-f_decomp;
        doc=doc+f_death*(1-f_2soc)+f_decomp+ebc_loss-f_mbcUdoc;
        
        co2=f_mbcUdoc*(1-cue_mbc);
        co2_all=ra+co2;  %%youcuihai
        n_step= n_step+1;
        res(1,n_step) = leafc;
        res(2,n_step) = rootc;
        res(3,n_step) = litc;
        res(4,n_step) = mbc;
        res(5,n_step) = ebc;
        res(6,n_step) = doc;
        res(7,n_step) = soc;
        res(8,n_step) = f_decomp;
        res(9,n_step) = I_soc;
        res(10,n_step)= co2;
        res(11,n_step)= co2_all;
        %res(12,n_step)=soc;
    end
    end
	soc =  6458.964802;
func_output=res;
