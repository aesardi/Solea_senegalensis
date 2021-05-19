% A script to estimate T_A based on developmental rates; outside of DEBtool
% for Solea senegalensis; script by N Marn modified by (A. Sardi 06-may-2021)

%actualizado
close all

% T-ah data from Carballo et al 2018 
% laboratory experiment - reliable data
Tah = [... temperature (C), age at hatching (d) 
    16      2
    18      1.5
    20      1    ]; 
    
%units.Tah = {'deg C', 'd'}; label.Tab = {'temperature', 'age at hatching'};  
%bibkey.Tah = 'CarbFirm2018';
%comment.Tah = 'incubation time: from beginning of gastrula to hatching)';
%14683


% % T-ab data for 20 degrees from Manchado, 18 from RibeSara1999 and
% % 19.5 from YufeParr1999
% % laboratory experiment - reliable data
% Tab = [... temperature (C), age at birth (d)
%     18      6
%     19.5    5.5
%     20      4.5    ]; 
    
%units.Tab = {'deg C', 'd'}; label.Tab = {'temperature', 'age at birth, fertilization to birth'};  
%bibkey.Tab = 'Fond1979';
%comment.Tab = 'incubation time: from fertilisation to first feeding)';
%10574

%T-aj data for 16 and 20 degrees from Manchado, 18 from RibeSara1999 and
% 19.5 from YufeParr1999
% laboratory experiment - reliable data
Taj = [... % temperature (C), time since birth until END of metamorphosis (d) 
16      20
18      17
19      15
20      14];
% units.Taj = {'deg C', 'd'}; label.Taj = {'temperature', 'time since birth at endmetamorphosis'};  
% bibkey.Taj = 'Manchado, RibeSara1999 and YufeParr1999';
% comment.Taj = 'development time: from first feeding to metamorphosis)';

%7739
% calculating Arrhenius temperature (from Pecquerie2008 phD) according % to equation about temp affect on processes:  
% p(T) = exp (T_A/T1 - T_A/T) * p(Ti) 
% T is avg body temp, T1 a chosen reference temp, T_A the Arrh temp, p a
% physiological rate
noDsets = length( {'Tah','Taj'});
dset{1}.name = 'T-ah';dset{2}.name = 'T-aj';%dset{3}.name = 'T-ab'; 
Ta(1).data = Tah;Ta(2).data = Taj;% Ta(3).data = Tab; 

fig(1).sty = 'b*'; fig(2).sty = 'r*';%fig(3).sty = 'g*';
sim(1).sty = 'b-'; sim(2).sty = 'r-';%sim(3).sty = 'g-'; 

D(noDsets).data = 0; T(noDsets).data = 0; % initialize
for d = 1: noDsets
    D(d).data = Ta(d).data( :,2 ); % days, incubation duation
    T(d).data = 273.15 + Ta(d).data( :,1 ); % K, temp of incubation
end

% get coefficients and evaluate the function
for d = 1:noDsets
Ta(d).koef = polyfit(1./T(d).data, log(1./D(d).data), 1);
Ta(d).y = polyval(Ta(d).koef, 1./T(d).data);
end

figure

%choose datasets
dArray = [1 2]; %data array
sArray = [1 2]; % simulation array (for polyval)
for d = 1: length(dArray) % plot data
  hold on
    plot(1./T(dArray(d)).data, log(1./D(dArray(d)).data), fig(d).sty)
   txtd{d} = dset{dArray(d)}.name;
end

for s = 1: length(sArray) % plot f predictions
  hold on
   plot(1./T(sArray(s)).data, Ta(sArray(s)).y, sim(s).sty)
   txts{s} = dset{sArray(s)}.name;
   Arrh_T(s) = Ta(sArray(s)).koef(:,1);
end
 txt = [txtd  txts];
legend (txt)
xlabel('1/T, K')
ylabel(' ln (incub time), days')
title('Relationship between the incubation temperature and incubation duration')

fprintf('T_A for %s is %6.0f K,and for %s it is %6.0f K. \n', txtd{1}, Arrh_T(1), txtd{2}, Arrh_T(2))%,txtd{3}, Arrh_T(3) )
fprintf('Mean value for T_A is %6.0f K. \n', mean(Arrh_T) )


% check visual appearance of data
figure
hold on
plot(Tah(:,1), Tah(:,2), 'b*')
%plot(Tab(:,1), Tab(:,2), 'g*')
plot(Taj(:,1), Taj(:,2), 'r*')

legend({'ah','aj'});%, 'ab'})
xlabel('temperature (C)')
ylabel('age or time since fertilization (d)')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% data.T_dw = [... % temperature (C), dry weight at end of metamorphosis
%     16   994.7  
%     19   1281
%     20   2792];
% units.T_dw = {'deg C', 'ug'}; label.T_dw = {'temperature', 'dry weight at end of metamorphosis'};  
% bibkey.T_dw = 'Manchado and YufeParr1999';
% % comment.T_dw = ;
