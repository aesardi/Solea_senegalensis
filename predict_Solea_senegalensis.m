function [prdData, info] = predict_Solea_senegalensis(par, data, auxData)
  %modified the 10 of February by A. Sardi
  % unpack par, data, auxData
  cPar = parscomp_st(par); vars_pull(par); 
  vars_pull(cPar);  vars_pull(data);  vars_pull(auxData);

  
 % filters <-- you need the 'abj' specific filter (with s_M)
 filterChecks =   E_Hh > E_Hb || E_Hh <= 0 || f_field > 1 || f_tL > 1 || f_tL2 > 1 || ... % maturity at hatching has to be between 0 and Ehb
     f_tL < 0.1 || f_tL2 < 0.1 || ...
     ~reach_birth(g, k, v_Hb, f_tL) || ... % constraint required for reaching birth with that f
     ~reach_birth(g, k, v_Hb, f_tL2) ; %|| ...  % constraint required for reaching birth with that f
   
  %~reach_birth(g, k, v_Hb, f_field) || ...
 
  % k * v_Hp >= f_field^3 || ... % constraint constraint required for reaching puberty with f_field
%      ~reach_birth(g, k, v_Hb, f_TeixCabr) || ... % constraint required for reaching birth with that f
%      k * v_Hp >= f_TeixCabr^3 || ... % constraint constraint required for reaching puberty with f_TeixCabr
%      ~reach_birth(g, k, v_Hb, f_TeixCabr2) || ... % constraint required for reaching birth with that f
%      k * v_Hp >= f_TeixCabr2^3 || ... % constraint constraint required for reaching puberty with f_TeixCabr2
%      ~reach_birth(g, k, v_Hb, f_CanaFern1)|| ... % constraint required for reaching birth with that f
%      k * v_Hp >= f_CanaFern1^3 ||...  % constraint constraint required for reaching puberty with f_CanaFern1
%      ~reach_birth(g, k, v_Hb, f_CanaFern2)|| ... % constraint required for reaching birth with that f
%      k * v_Hp >= f_CanaFern2^3 ||...  % constraint constraint required for reaching puberty with f_CanaFern2
%      ~reach_birth(g, k, v_Hb, f_CanaFern3)|| ... % constraint required for reaching birth with that f
%      k * v_Hp >= f_CanaFern3^3 ||...  % constraint constraint required for reaching puberty with f_CanaFern3
%      ~reach_birth(g, k, v_Hb, f_CanaFern4)|| ... % constraint required for reaching birth with that f
%      k * v_Hp >= f_CanaFern4^3 ||...  % constraint constraint required for reaching puberty with f_CanaFern4
%      k * v_Hp >= f_DiniRibe^3 ||...     % constraint constraint required for reaching puberty with f_DineRibe
%     ~reach_birth(g, k, v_Hb, f_DiniRibe); % constraint required for reaching birth with that f
%       

 
  if filterChecks  
    info = 0;
    prdData = {};
    return;
  end  
 
% 0-var data 'ah';'ab';'tj';'aj';'tp';'am';'Lh';'Lb';'Lp';'Li';'Wwh';'Wwb';'Wdh';'Wdb';'Ri'
  % compute temperature correction factors
  T_pars=[T_A, T_L, T_H, T_AL, T_AH];
  TC_ah = tempcorr(temp.ah, T_ref, T_pars);
  TC_ab = tempcorr(temp.ab, T_ref, T_pars);
  TC_aj = tempcorr(temp.aj, T_ref, T_pars);
  TC_ap = tempcorr(temp.ap, T_ref, T_pars);
  TC_am = tempcorr(temp.am, T_ref, T_pars);
  TC_tj = tempcorr(temp.tj, T_ref, T_pars);
  TC_Ri = tempcorr(temp.Ri, T_ref, T_pars);
  % univariate data temp corrections
  %TC_Tah = tempcorr(temp.Tah, T_ref, T_pars);
  TC_tL = tempcorr(temp.tL, T_ref, T_pars);
  TC_tL2 = tempcorr(temp.tL2, T_ref, T_pars);
  TC_tL3 = tempcorr(temp.tL3, T_ref, T_pars);
  TC_tL4 = tempcorr(temp.tL4, T_ref, T_pars);
  TC_tWd = tempcorr(temp.tWd(1), T_ref, T_pars);
  TC_tWd2 = tempcorr(temp.tWd2(1), T_ref, T_pars);
  TC_tWd_f1 = tempcorr(temp.tWd_f1(1), T_ref, T_pars);
  TC_tWd_f2 = tempcorr(temp.tWd_f2(1), T_ref, T_pars);
  TC_tWd_f3 = tempcorr(temp.tWd_f3(1), T_ref, T_pars);
  TC_tWd_f4 = tempcorr(temp.tWd_f4(1), T_ref, T_pars);
 % TC_tM_N = tempcorr(temp.tM_N(1), T_ref,T_A); 
%   TC_tE = tempcorr(temp.tE(1), T_ref, T_A);
%   TC_tE2 = tempcorr(temp.tE2(1), T_ref, T_A);  
  
  %% % zero-variate data

  % life cycle
  pars_tj = [g; k; l_T; v_Hb; v_Hj; v_Hp];
  [tau_j, tau_p, tau_b, l_j, l_p, l_b, l_i, rho_j, rho_B, info] = get_tj(pars_tj, f);
  
  if info ~= 1 % numerical procedure failed
     fprintf('warning: invalid parameter value combination for get_tj \n')
  end
  
  % initial
  pars_UE0 = [V_Hb; g; k_J; k_M; v]; % compose parameter vector
  U_E0 = initial_scaled_reserve(f, pars_UE0); % d.cm^2, initial scaled reserve
  % ^-- U_E0 is underpredicted; this gives an E_0 of U_E0 * p_Am * TC_ah = 2.13e-04

  %EGG
   
  E_0 = U_E0 * p_Am ;          % J, energy in egg
  Wd_0 = E_0 * w_E/ mu_E;      % g, egg dry weight 
%   V0 = Wd_0/ d_E;             % cm^3, egg volume 
%   Lw_0 = (6 * V0/ pi)^(1/3);  % cm, egg diameter
%   
  
  
  % HATCH  
 [U_H, aUL] = ode45(@dget_aul, [0; U_Hh; U_Hb], [0 U_E0 1e-10], [], kap, v, k_J, g, L_m);
  aT_h = aUL(2,1)/ TC_ah;                   % d, age at hatch at f and T
  L_h = aUL(2,3);                           % cm, strucural length at hatch
  Lw_h = L_h/del_Me; 
  E_h = aUL(2,2) * p_Am *TC_ah;             % J, energy in reserves at hatch
  Wd_h = (d_V * L_h^3 + w_E/ mu_E * E_h) *1e6; % ug, dry weight at hatch; 

  % BIRTH
  L_b = L_m * l_b;                  % cm, structural length at birth at f
  Lw_b = L_b/ del_Me;                % cm, physical length at birth at f
  aT_b = tau_b/ k_M/ TC_ab;           % d, age at birth at f and T
% 
  
  Wd_b = d_V *L_b^3 * (1 + f * w) *1e6; % ug, dry weight at birth at f 

  %Ww_b = L_b^3 * (1 + f * w) * 1e6;       % ug, wet weight at birth at f 
  

  % metamorphosis decide if using start or end of metam (change del_M accordingly)
  L_j = L_m * l_j;                  % cm, structural length at metam
  Lw_j = L_j/ del_Me;                % cm, physical length at START of metam at f 
  
  tT_j = (tau_j - tau_b) / k_M/ TC_tj;  % d, time since birth at metam
  aT_j = tau_j/TC_aj/k_M;           % d, age at metamorphosis 
  
  %Ww_j = L_j^3 * (1 + f * w) * 1e6;     % ug, wet weight at metam 
  Wd_j = L_j^3 * d_V * (1 + f * w) * 1e6; % ug, dry weight at metam 
  
  
  [tau_j, tau_p, tau_b, l_j, l_p, l_b, l_i, rho_j, rho_B, info] = get_tj(pars_tj, f_field);
  
  % puberty 
  L_p = L_m * l_p;                  % cm, structural length at puberty at f
  aT_p = tau_p / k_M/ TC_ap;        %d, age at puberty at f and T
  Lw_p = L_p/ del_M;                % cm, physical length at puberty at f
  Ww_p = L_p^3 *(1 + f * w);        % g, wet weight at puberty 
 

  % ultimate
  L_i = L_m * l_i;                  % cm, ultimate structural length at f
  Lw_i = L_i/ del_M;                % cm, ultimate physical length at f
  Ww_i = L_i^3 * (1 + f_field * w);       % g, ultimate wet weight 
 
  % reproduction
  pars_R = [kap, kap_R, g, k_J, k_M, L_T, v, U_Hb, U_Hj, U_Hp];
  [R_i, UE0, Lb, Lj, Lp, info]  =  reprod_rate_j(L_i, f, pars_R);
  %RT_i = TC_Ri * R_i;% #/d, max reprod rate
  %RT_i = TC_Ri * reprod_rate_j(L_i, f, pars_R); --> 0
  
  
  %RT_i2 = TC_Ri * reprod_rate_j(Li / del_M, f, pars_R); %--> prediction
  %approachs actual value
  %RT_i3 = TC_Ri * reprod_rate_j(Lw_i, f, pars_R); %--> prediction = 94.7
  
 % L_b^3 = Ww_b / (1 + f * w) * 1e6  ; if i dont find the lenght but have
 % the weight for the Ri data, use this equation
  
  
  RT_i = TC_Ri * reprod_rate_j(Lw_i * del_M, f, pars_R); %--> prediction is 0

%   
   % life span
  pars_tm = [g; l_T; h_a/ k_M^2; s_G];  % compose parameter vector at T_ref
  t_m = get_tm_s(pars_tm, f_field, l_b);      % -, scaled mean life span at T_ref
  aT_m = t_m/ k_M/ TC_am;               % d, mean life span at T
  
  %pack to output
  prdData.ah = aT_h;
  prdData.ab = aT_b;
  prdData.aj = aT_j; prdData.aj2 = aT_j; % age at START and END of metam. 
  % Because metamorphosis is modeled as a 'discrete event' rather than something lasting more days, 
  % the predicted value should fall between these two observed values
  prdData.tj = tT_j;prdData.tj2 = tT_j;
  prdData.ap = aT_p;
  prdData.am = aT_m;
 % prdData.L0 = Lw_0;
  prdData.Lh = Lw_h;
  prdData.Lb = Lw_b;
  prdData.Lj = Lw_j;prdData.Lj2 = Lw_j;
  prdData.Lp_f = Lw_p; prdData.Lp_m = Lw_p; %lenght at puberty for females and males
  prdData.Li = Lw_i;
  prdData.Wd0 = Wd_0;
  %prdData.Wwb = Ww_b;
  %prdData.Wwj0 = Ww_j; prdData.Wwj = Ww_j; % wet weight at START and END of metam. 
  % Because metamorphosis is modeled as a 'discrete event' rather than something lasting more days, 
  % the predicted value should fall between these two observed values
  prdData.Wwp_f = Ww_p;prdData.Wwp_m = Ww_p; %wet weight at puberty for females and males
  prdData.Wwi = Ww_i;
  prdData.Wdh = Wd_h;
  prdData.Wdb = Wd_b;
  prdData.Wdj = Wd_j; prdData.Wdj2 = Wd_j; % dry weight at START and END of metam. 
  % Because metamorphosis is modeled as a 'discrete event' rather than something lasting more days, 
  % the predicted value should fall between these two observed values
  %prdData.Ri = RT_i;
  %prdData.Ri = RT_i2;
  %prdData.Ri = RT_i3;
  prdData.Ri = RT_i;
  prdData.E0 = E_0;
%   
  %% ------------- uni-variate data----------------
  % PARAMETERS for egg
%   pars_UE0 = [V_Hb; g; k_J; k_M; v];
%   [U_E0, L_b, info] = initial_scaled_reserve(f, pars_UE0);
%   [U_H aUL] = ode45(@dget_aul, [0; U_Hh; U_Hb], [0 U_E0 1e-10], [], kap, v, k_J, g, L_m);
%   
%   % Tah 
%   a_h = aUL(2,1);                             % d, age at birth
%   Eah = a_h*ones(length(Tah(:,1)),1) ./ TC_Tah; % d, age at birth temp corrected   
  
  %% time-length tL
  
  [tau_j, tau_p, tau_b, l_j, l_p, l_b, l_i, rho_j, rho_B, info] = get_tj(pars_tj, f_tL);
  
  %t-L RibeSara1999
  %time-length since hatching concatenated to juveniles (from day 12)
  kT_M = k_M * TC_tL;
  rT_B = rho_B * kT_M;  % 1/d, von Bert growth rate   
  rT_j = rho_j * kT_M;  % 1/d, exponential growth rate
  tT_j = (tau_j - tau_b)/ kT_M; % time since *birth* at metamorphosis


  L_j = l_j * L_m; 
  L_i = l_i * L_m;
  
  EL_bj = L_b * exp(tL((tL(:,1)<= tT_j),1) * rT_j/3); % exponential growth as V1-morph
%   EL_hj = L_h * exp (tL((tL(:,1)<= aT_j),1) * rT_j/3);
  EL_ji = L_i - (L_i - L_j) * exp( - rT_B * (tL((tL(:,1) >= tT_j),1) - tT_j)); % cm, expected length at time
  ELw = [EL_bj/del_Me; EL_ji/del_M]; % catenate lengths
  
  % tL2 RibeEngr2017 only juveniles
  [tau_j, tau_p, tau_b, l_j, l_p, l_b, l_i, rho_j, rho_B, info] = get_tj(pars_tj, f_tL2);
  kT_M2 = TC_tL2 * k_M;
  rT_B =  rho_B * kT_M2;  % 1/d, von Bert growth rate   
  rT_j =  rho_j * kT_M2;  % 1/d, exponential growth rate
  tT_j = (tau_j - tau_b)/ kT_M2; % time since birth at metamorphosis
%   aT_j = tau_j/ kT_M2 ; % time since *fertilization* at metamorphosis --> closer to hatching (which is where the dataset starts!)
  
  L_j = l_j * L_m; 
  L_i = l_i * L_m;
  EL_bj2 = Lw_b * exp(tL2((tL2(:,1)<= tT_j),1)  * rT_j/3); % exponential growth as V1-morph
  EL_ji2 = L_i - (L_i - L_j) * exp( - rT_B * (tL2((tL2(:,1) > tT_j),1)- tT_j)); % cm, expected length at time
  ELw2 = [EL_bj2; EL_ji2]/del_M; %TL2 is for juveniles all already metamorphosed
  
   %t-L3 TeixCabr2010 females
  [tau_j, tau_p, tau_b, l_j, l_p, l_b, l_i, rho_j, rho_B, info] = get_tj(pars_tj, f_TeixCabr);
  kT_M2 = TC_tL3 * k_M;
  rT_B =  rho_B * kT_M2;  % 1/d, von Bert growth rate   
  rT_j =  rho_j * kT_M2;  % 1/d, exponential growth rate
  tT_j = (tau_j - tau_b)/ kT_M2; % time since birth at metamorphosis

  L_j = l_j * L_m; 
  L_i = l_i * L_m;
  EL_bj3 = Lw_b * exp(tL3((tL3(:,1)<= tT_j),1)  * rT_j/3); % exponential growth as V1-morph
  EL_ji3 = L_i - (L_i - L_j) * exp( - rT_B * (tL3((tL3(:,1) > tT_j),1)- tT_j)); % cm, expected length at time
  ELw3 = [EL_bj3; EL_ji3]/del_M; 
  
   %t-L4 TeixCabr2010 males
  [tau_j, tau_p, tau_b, l_j, l_p, l_b, l_i, rho_j, rho_B, info] = get_tj(pars_tj, f_TeixCabr2);
  kT_M2 = TC_tL4 * k_M;
  rT_B =  rho_B * kT_M2;  % 1/d, von Bert growth rate   
  rT_j =  rho_j * kT_M2;  % 1/d, exponential growth rate
  tT_j = (tau_j - tau_b)/ kT_M2; % time since birth at metamorphosis

  L_j = l_j * L_m; 
  L_i = l_i * L_m;
  EL_bj4 = Lw_b * exp(tL4((tL4(:,1)<= tT_j),1)  * rT_j/3); % exponential growth as V1-morph
  EL_ji4 = L_i - (L_i - L_j) * exp( - rT_B * (tL4((tL4(:,1) > tT_j),1)- tT_j)); % cm, expected length at time
  ELw4 = [EL_bj4; EL_ji4]/del_M; %
  
  
    
 %% % length-weight
 %lenght wet weight fish at 190, 398 and 790 days of age (manchado data)
ELWw = (LWw(:,1) * del_M).^3 * (1 + f_field * w); 

% ELWw2 = (LWw2(:,1) * del_M).^3 * (1 + f_Man * w); %for females
% ELWw3 = (LWw3(:,1) * del_M).^3 * (1 + f_Man * w); %for males
 
 
%Length dry weight --> laboratory conditions assume ab libitum use f=1
%LWd (OrtiFune2019)
L1 = LWd(LWd(:,1)<data.Lj,1) * del_Me; % for data before metamorphosis
L2 = LWd(LWd(:,1)>=data.Lj,1) * del_Me; % for data after metamorphosis
ELWd1 = [L1; L2].^3 * d_V*(1 + f * w)*1e6; % ug, wet weight 

% here we assume that wga is the same before and after metamorphosis 

%LWd2 (RibeEngr2017) --> they are all metamorphosed 
%L3 = LWd2(LWd2(:,1)<data.Lj,1) * del_Me; %before metamorphosis
L4 = LWd2(LWd2(:,1)>data.Lj,1) * del_Me; %after metamorphosis
ELWd2 = L4.^3 * d_V* (1 + f * w)*1e6; % ug, dry weight 

%LWd3 (YufeParr1999)
L5 = LWd3(LWd3(:,1)<data.Lj,1) * del_Me; %before metamorphosis
L6 = LWd3(LWd3(:,1)>data.Lj,1) * del_M; %after metamorphosis
ELWd3 = [L5; L6].^3 * d_V* (1 + f * w)*1e6; % ug, dry weight 

%% % time-dry weight
% tWd (YufeParr1999)-->laboratory conditions assume ab libitum use f=1
  [tau_j, tau_p, tau_b, l_j, l_p, l_b, l_i, rho_j, rho_B, info] = get_tj(pars_tj, f);
  L_b = l_b * L_m; 
  L_j = l_j * L_m; 
  L_i = l_i * L_m;
  
    
  tT_j1 = (tau_j - tau_b)/(k_M * TC_tWd);    % d, time since birth at metamorphosis corrected at 19 degrees for dry weight data
  rT_j = rho_j * (k_M * TC_tWd);  
  rT_B = rho_B * (k_M * TC_tWd);  
  L_bj = L_b * exp(tWd((tWd(:,1) <= tT_j1),1) * rT_j/ 3);
  L_jm = L_i - (L_i - L_j) * exp( - rT_B * (tWd((tWd(:,1) > tT_j1),1) - tT_j1)); % cm, expected length at time
  EWd = [L_bj; L_jm].^3 * d_V * (1 + f * w) * 1e6 ;
  
 % tWd2 (ParrYufe2001) --> laboratory conditions assume ab libitum use f=1
  tT_j1 = (tau_j - tau_b)/(k_M * TC_tWd2);    % d, time since birth at metamorphosis corrected at 19 degrees for dry weight data
  rT_j = rho_j * (k_M * TC_tWd2);  
  rT_B = rho_B * (k_M * TC_tWd2);  
  L_bj = L_b * exp(tWd2((tWd2(:,1) <= tT_j1),1) * rT_j/ 3);
  L_jm = L_i - (L_i - L_j) * exp( - rT_B * (tWd2((tWd2(:,1) > tT_j1),1) - tT_j1)); % cm, expected length at time
  EWd2 = [L_bj; L_jm].^3 * d_V * (1 + f * w) * 1e6 ;
  
 %tWd_Feeding regimes (CañaFern1999)
 %1 L100
[tau_j, tau_p, tau_b, l_j, l_p, l_b, l_i, rho_j, rho_B] = get_tj(pars_tj, f_CanaFern1);
rT_j = rho_j * (k_M * TC_tWd_f1); 
rT_B = rho_B * (k_M * TC_tWd_f1); 
tT_j = (tau_j - tau_b)/ (k_M * TC_tWd_f1); 
L_bj = L_b * exp(tWd_f1(tWd_f1(:,1) < tT_j,1) * rT_j/ 3); % cm length and weight during V1-morph period
L_jm = L_i - (L_i - L_j) * exp( - rT_B * (tWd_f1(tWd_f1(:,1) >= tT_j,1) - tT_j));   % cm, length after V1-morph period
EWd_1 = 1e6 * [L_bj; L_jm].^3  * d_V * (1 + f_CanaFern1 * w); % ug, dry weight

% %2 L50
[tau_j, tau_p, tau_b, l_j, l_p, l_b, l_i, rho_j, rho_B] = get_tj(pars_tj, f_CanaFern2);
rT_j = rho_j * (k_M * TC_tWd_f2); 
rT_B = rho_B * (k_M * TC_tWd_f2); 
tT_j = (tau_j - tau_b)/ (k_M * TC_tWd_f2); 
L_bj = L_b * exp(tWd_f2(tWd_f2(:,1) < tT_j,1) * rT_j/ 3); % cm length and weight during V1-morph period
L_jm = L_i - (L_i - L_j) * exp( - rT_B * (tWd_f2(tWd_f2(:,1) >= tT_j,1) - tT_j));   % cm, length after V1-morph period
EWd_2 = 1e6 * [L_bj; L_jm].^3  * d_V * (1 + f_CanaFern2 * w); % ug, dry weight
  
%   
% %3 L100I50
[tau_j, tau_p, tau_b, l_j, l_p, l_b, l_i, rho_j, rho_B] = get_tj(pars_tj, f_CanaFern3);
rT_j = rho_j * (k_M * TC_tWd_f3); 
rT_B = rho_B * (k_M * TC_tWd_f3); 
tT_j = (tau_j - tau_b)/ (k_M * TC_tWd_f3); 
L_bj = L_b * exp(tWd_f3(tWd_f3(:,1) < tT_j,1) * rT_j/ 3); % cm length and weight during V1-morph period
L_jm = L_i - (L_i - L_j) * exp( - rT_B * (tWd_f3(tWd_f3(:,1) >= tT_j,1) - tT_j));   % cm, length after V1-morph period
EWd_3 = 1e6 * [L_bj; L_jm].^3  * d_V * (1 + f_CanaFern3 * w); % ug, dry weight
  
%4 L50I50
[tau_j, tau_p, tau_b, l_j, l_p, l_b, l_i, rho_j, rho_B] = get_tj(pars_tj, f_CanaFern4);
rT_j = rho_j * (k_M * TC_tWd_f4); 
rT_B = rho_B * (k_M * TC_tWd_f4); 
tT_j = (tau_j - tau_b)/ (k_M * TC_tWd_f4); 
L_bj = L_b * exp(tWd_f4(tWd_f4(:,1) < tT_j,1) * rT_j/ 3); % cm length and weight during V1-morph period
L_jm = L_i - (L_i - L_j) * exp( - rT_B * (tWd_f4(tWd_f4(:,1) >= tT_j,1) - tT_j));   % cm, length after V1-morph period
EWd_4 = 1e6 * [L_bj; L_jm].^3  * d_V * (1 + f_CanaFern4 * w); % ug, dry weight
  


%% %time energy content in larvae

%  kT_M = k_M * TC_tE;
%  rT_B = rho_B * kT_M;  % 1/d, von Bert growth rate   
%  rT_j = rho_j * kT_M;  % 1/d, exponential growth rate
%  tT_j = (tau_j - tau_b)/ kT_M; % time since *birth* at metamorphosis
%   
%   L_j = l_j * L_m; 
%   L_i = l_i * L_m;
%   
%   EE_bj = L_b * exp(tE((tE(:,1)<= tT_j),1) * rT_j/3); % exponential growth as V1-morph
%   EE_ji = L_i - (L_i - L_j) * exp( - rT_B * (tE((tE(:,1) >= tT_j),1) - tT_j)); % cm, expected length at time
%   E = [EE_bj/del_Me; EE_ji/del_M]; % catenate lengths
%   %E = [EE_bj; EE_ji]/del_M; %   
% 
%   EE_bj2 = L_b * exp(tE2((tE2(:,1)<= tT_j),1) * rT_j/3); % exponential growth as V1-morph
%   EE_ji2 = L_i - (L_i - L_j) * exp( - rT_B * (tE2((tE2(:,1) >= tT_j),1) - tT_j)); % cm, expected length at time
%   E2 = [EE_bj2/del_Me; EE_ji2/del_M]; % data are from individuals before metamorphosis
%    
%   EE = E.^3 * (M_V * mu_V + f_tL * E_m) * 1e-3;   %in J 
%   EE2 = E2.^3 * (M_V * mu_V + f_tL * E_m) * 1e-3; 

%  % N-data
%   rT_B = TC_tM_N * k_M/ 3/ (1 + f/ g);
%   L = L_i - (L_i - L_b) * exp( - rT_B * tM_N(:,1));
%  % EM_N = 1e6 * L.^3 * d_V * Y_N_W * (1 + f * w); %   ug N from
%  % Moniliformis dubius --> what is Y_N_W?
%    EM_N = 1e6 * L.^3 * d_V * (1 + f * w); %   ug N

  % pack to output
%   prdData.Tah = Eah;
  prdData.tL = ELw;
  prdData.tL2 =ELw2;
  prdData.tL3 =ELw3;
  prdData.tL4 =ELw4; 
  prdData.LWw = ELWw;
%   prdData.LWw2 = ELWw2;
%   prdData.LWw3 = ELWw3;
  prdData.LWd = ELWd1;
  prdData.LWd2 = ELWd2;
  prdData.LWd3 = ELWd3;
  prdData.tWd = EWd;
  prdData.tWd2 = EWd2;
  prdData.tWd_f1 = EWd_1;
  prdData.tWd_f2 = EWd_2;
  prdData.tWd_f3 = EWd_3;
  prdData.tWd_f4 = EWd_4;
  %prdData.tM_N = EM_N;
%   prdData.tE = EE;
%   prdData.tE2 = EE2;

  