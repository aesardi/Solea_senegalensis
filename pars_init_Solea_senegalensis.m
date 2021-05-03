function [par, metaPar, txtPar] = pars_init_Solea_senegalensis(metaData)

metaPar.model = 'abj'; 

%% reference parameter (not to be changed) 
par.T_ref = 293.15;   free.T_ref = 0;   units.T_ref = 'K';        label.T_ref = 'Reference temperature'; 

%% core primary parameters 
par.T_A = 1.652e+04	;  free.T_A   = 1;   units.T_A = 'K';          label.T_A = 'Arrhenius temperature'; 
par.z = 3.968;        free.z     = 1;   units.z = '-';            label.z = 'zoom factor'; 
par.F_m = 6.5;        free.F_m   = 0;   units.F_m = 'l/d.cm^2';   label.F_m = '{F_m}, max spec searching rate'; 
par.kap_X = 0.8;      free.kap_X = 0;   units.kap_X = '-';        label.kap_X = 'digestion efficiency of food to reserve'; 
par.kap_P = 0.1;      free.kap_P = 0;   units.kap_P = '-';        label.kap_P = 'faecation efficiency of food to faeces'; 
par.v = 0.0599;     free.v     = 1;   units.v = 'cm/d';         label.v = 'energy conductance'; 
par.kap = 0.7254;    free.kap   = 1;   units.kap = '-';          label.kap = 'allocation fraction to soma'; 
par.kap_R = 0.95;     free.kap_R = 0;   units.kap_R = '-';        label.kap_R = 'reproduction efficiency'; 
par.p_M = 26.68;    free.p_M   = 1;   units.p_M = 'J/d.cm^3';   label.p_M = '[p_M], vol-spec somatic maint'; 
par.p_T = 0;          free.p_T   = 0;   units.p_T = 'J/d.cm^2';   label.p_T = '{p_T}, surf-spec somatic maint'; 
par.k_J = 0.002;      free.k_J   = 0;   units.k_J = '1/d';        label.k_J = 'maturity maint rate coefficient'; 
par.E_G = 5169;  free.E_G   = 1;   units.E_G = 'J/cm^3';     label.E_G = '[E_G], spec cost for structure'; 
par.E_Hb = 0.2328; free.E_Hb  = 1;   units.E_Hb = 'J';         label.E_Hb = 'maturity at birth'; 
par.E_Hj = 4.81; free.E_Hj  = 1;   units.E_Hj = 'J';         label.E_Hj = 'maturity at metamorphosis'; 
par.E_Hp = 1.372e+06; free.E_Hp  = 1;   units.E_Hp = 'J';         label.E_Hp = 'maturity at puberty'; 
par.h_a = 1.958e-08;  free.h_a   = 1;   units.h_a = '1/d^2';      label.h_a = 'Weibull aging acceleration'; 
par.s_G = 0.001;      free.s_G   = 0;   units.s_G = '-';          label.s_G = 'Gompertz stress coefficient'; 

%% other parameters 
par.E_Hh = 0.04708; free.E_Hh  = 1;   units.E_Hh = 'J';         label.E_Hh = 'maturity at hatch'; 
par.T_AH = 100000;    free.T_AH  = 0;   units.T_AH = 'K';         label.T_AH = 'T_AH'; 
par.T_AL = 50000;     free.T_AL  = 0;   units.T_AL = 'K';         label.T_AL = 'T_AL'; 
par.T_H = 303;        free.T_H   = 0;   units.T_H = 'K';          label.T_H = 'T_H'; 
par.T_L = 285;        free.T_L   = 0;   units.T_L = 'K';          label.T_L = 'T_L'; 
par.del_M = 0.1921;  free.del_M = 1;   units.del_M = '-';        label.del_M = 'del_M, shape coefficient'; %from solea solea
par.del_Me = 0.2101;  free.del_Me = 1;   units.del_Me = '-';       label.del_Me = 'shape coefficient for embryo'; 
par.f = 0.7255;      free.f     = 1;   units.f = '-';            label.f = 'scaled functional response for YufeParr1999 data'; 
par.f_CanaFern1 = 0.4051;  free.f_CanaFern1 = 1;   units.f_CanaFern1 = '-';  label.f_CanaFern1 = 'scaled functional response for tWd_f1 data'; 
par.f_CanaFern2 = 0.4133;  free.f_CanaFern2 = 1;   units.f_CanaFern2 = '-';  label.f_CanaFern2 = 'scaled functional response for tWd_f2 data'; 
par.f_CanaFern3 = 0.546;  free.f_CanaFern3 = 1;   units.f_CanaFern3 = '-';  label.f_CanaFern3 = 'scaled functional response for tWd_f3 data'; 
par.f_CanaFern4 = 0.4543;  free.f_CanaFern4 = 1;   units.f_CanaFern4 = '-';  label.f_CanaFern4 = 'scaled functional response for tWd_f3 data'; 
par.f_TeixCabr = 0.7757;  free.f_TeixCabr = 1;   units.f_TeixCabr = '-';   label.f_TeixCabr = 'scaled functional response for tL3 data'; 
par.f_TeixCabr2 = 0.7264;  free.f_TeixCabr2 = 1;   units.f_TeixCabr2 = '-';  label.f_TeixCabr2 = 'scaled functional response for Tl4 data'; 
par.f_field = 0.9483;  free.f_field = 1;   units.f_field = '-';      label.f_field = 'scaled functional response for 0-var data fromthe field'; 
par.f_tL = 0.9034;   free.f_tL  = 1;   units.f_tL = '-';         label.f_tL = 'scaled functional response for tL data'; 
par.f_tL2 = 0.366;  free.f_tL2 = 1;   units.f_tL2 = '-';        label.f_tL2 = 'scaled functional response for tL2 data'; 
%par.f_Man = 1;  free.f_Man = 0;   units.f_Man = '-';        label.f_Man = 'scaled functional response for tL2 data'; 
%% set chemical parameters from Kooy2010 
[par, units, label, free] = addchem(par, units, label, free, metaData.phylum, metaData.class); 

%% Pack output: 
txtPar.units = units; txtPar.label = label; par.free = free; 
