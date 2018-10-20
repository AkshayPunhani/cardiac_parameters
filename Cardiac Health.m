MATLAB code for ECG and PPG
clear all
clc
close all
 
 
%% Importing signal into MATLAB
% ECG=importdata('C:\Users\Abhinav\Desktop\SEMESTER 6\IP\Files\Signal Bank\MIMIC databse- refined\Subject 20\ECG.txt');
% PPG_data=importdata('C:\Users\Abhinav\Desktop\SEMESTER 6\IP\Files\Signal Bank\MIMIC databse- refined\Subject 20\PPG.txt');
% HR=importdata('C:\Users\Abhinav\Desktop\SEMESTER 6\IP\Files\Signal Bank\MIMIC databse- refined\Subject 20\HR.txt');
% BP_Dia =importdata('C:\Users\Abhinav\Desktop\SEMESTER 6\IP\Files\Signal Bank\MIMIC databse- refined\Subject 20\BP_Dia.txt');
% BP_Sys=importdata('C:\Users\Abhinav\Desktop\SEMESTER 6\IP\Files\Signal Bank\MIMIC databse- refined\Subject 20\BP_Sys.txt');
% SPO2 = importdata('C:\Users\Abhinav\Desktop\SEMESTER 6\IP\Files\Signal Bank\MIMIC databse- refined\Subject 20\SPO2.txt');
 
% ECG = ECG.data;
% PPG_data = PPG_data.data;
% HR = HR.data;
% BP_Dia = BP_Dia.data;
% BP_Sys = BP_Sys.data;
% SPO2= SPO2.data;
 
samples = 800;
data =csvread('C:\Users\Abhinav\Desktop\SEMESTER 6\IP\Queensland\32.csv', 1,3);
HR = data(1:samples/10,1);
SPO2= data(1:samples/10,2);
BP_Sys = data(1:samples/10,3);
BP_Dia = data(1:samples/10,4);
ECG = data(1:samples,5);
PPG_data = data(1:samples,6);
 
 
fs = 100;
min_seperation = round(fs);
 
%% Signal pre-processing
%PPG_data = PPG_data(1:length(PPG_data));
PPG_data = PPG_data - mean (PPG_data);    % cancel DC conponents
PPG_data = PPG_data/ max( abs(PPG_data )); % normalize to one
[p,s,mu] = polyfit((1:numel(PPG_data))',PPG_data,6);
f_y = polyval(p,(1:numel(PPG_data))',[],mu);
 
PPG = PPG_data - f_y;        % Detrend data
PPG_inv = -PPG';
 
%ECG = ECG(1:1250);
ECG = ECG - mean (ECG);    % cancel DC conponents
ECG = ECG/ max( abs(ECG )); % normalize to one
 
ECG = ECG;
FLAG = 0;       %% Flag bit for wrong R-peak detection
FLAG1 = 0;       %% Flag bit for wrong PPG
 
Spacing_factor = 1;
 
 
%% Locating R-peaks, PCG foot and peak
 
[pks_ECG_R,locs_ECG_R] = findpeaks(ECG,'MinPeakDistance',min_seperation,'MinPeakHeight',0.5);
 
[pks_PPG,locs_PPG] = findpeaks(PPG,'MinPeakDistance',min_seperation,'MinPeakHeight',0.2);
[pks_PPG_inv,locs_PPG_inv] = findpeaks(PPG_inv,'MinPeakDistance',min_seperation,'MinPeakHeight',0.2);
 
 
%% R-peak verification
[pks_ECG_Rn,locs_ECG_Rn] = findpeaks(ECG,'MinPeakDistance',round(min_seperation/5),'MinPeakHeight',0.5);
 
for (j=2:4)
for i=1:length(locs_ECG_R)-1
    diff_ECG(i) = locs_ECG_R(i+1) - locs_ECG_R(i);
end
 
counter =0;
for i=1:length(diff_ECG)
    if ((mean(diff_ECG)-diff_ECG(i))/mean(diff_ECG)*100 > 10)
        counter = counter +1;    
    end
end
 
 
 if abs (length(pks_ECG_Rn)-(length(pks_ECG_R))) > 2
        [pks_ECG_R,locs_ECG_R] = findpeaks(ECG,'MinPeakDistance',round((min_seperation)/(j+1)),'MinPeakHeight',0.5);
        %Spacing_factor = j+1;
        display ('R-peak detection algorithm revised 2')
        FLAG = 1;
 
 elseif (counter/length(diff_ECG))*100 > 10
    [pks_ECG_R,locs_ECG_R] = findpeaks(ECG,'MinPeakDistance',round((min_seperation)/(j+1)),'MinPeakHeight',0.5);
    %Spacing_factor = j+1;
    display ('R-peak detection algorithm revised 1')
    FLAG = 1;
    
 else  
    FLAG = 11;
    Spacing_factor = j;
    display('Breaking out')
    break;
end
end
 
%% Peak vector position checking
 
if locs_ECG_R(1) > locs_PPG(1)
     locs_PPG_check = locs_PPG(2:length(locs_PPG));
end
 
   
for (i=1: min(length(locs_ECG_R),length(locs_PPG))-1)
    
    if (((locs_ECG_R(i) < locs_PPG(i)) && (locs_PPG(i) < locs_ECG_R(i+1))))
       continue;
        
    elseif (((locs_ECG_R(i) < locs_PPG(i)) && (locs_PPG(i) > locs_ECG_R(i+1))))
       FLAG1 = 12;
    else
        %FLAG1 = 12;
    end
end
  
if (FLAG1 == 12)
    %display ('Error in PTT: Location error of PTT foot/peak wrt R-peak')
    display('FLAG1 active')
end
        
%% PPG verification
 
locs_PPG_inv = locs_PPG_inv';
pks_PPG_inv = pks_PPG_inv';
for (i=2:6)
if (abs(length(locs_PPG) - length(locs_PPG_inv))) > 1
    %[pks_PPG,locs_PPG] = findpeaks(PPG,'MinPeakDistance',round((min_seperation)/(i+1)),'MinPeakHeight',0.1);
    %[pks_PPG_inv,locs_PPG_inv] = findpeaks(PPG_inv,'MinPeakDistance',round((min_seperation)/(i+1)),'MinPeakHeight',0.1);
    display ('PPG peak detection algorithm revised 1')
    %display ('Wrong PTT')
end    
if (length(pks_ECG_R) - length(pks_PPG)) > 1
    [pks_PPG,locs_PPG] = findpeaks(PPG,'MinPeakDistance',round((min_seperation)/(i-1)),'MinPeakHeight',0.1);
    %[pks_PPG_inv,locs_PPG_inv] = findpeaks(PPG_inv,'MinPeakDistance',round((min_seperation)/(i+1)),'MinPeakHeight',0.1);
    display ('PPG peak detection algorithm PPG -1')
    
 elseif (length(pks_PPG_inv)- length(pks_ECG_R)) > 1
   % [pks_PPG,locs_PPG] = findpeaks(PPG,'MinPeakDistance',round((min_seperation)/(i+1)),'MinPeakHeight',0.1);
    [pks_PPG_inv,locs_PPG_inv] = findpeaks(PPG_inv,'MinPeakDistance',round((min_seperation)/(i+1)),'MinPeakHeight',0.1);
    display ('PPG peak detection algorithm PPG_inv +1')
end
 
 if (length(pks_ECG_R) - length(pks_PPG_inv)) > 1
   % [pks_PPG,locs_PPG] = findpeaks(PPG,'MinPeakDistance',round((min_seperation)/(i+1)),'MinPeakHeight',0.1);
    [pks_PPG_inv,locs_PPG_inv] = findpeaks(PPG_inv,'MinPeakDistance',round((min_seperation)/(i-1)),'MinPeakHeight',0.1);
    display ('PPG peak detection algorithm PPG_inv -1')
    
 
elseif (length(pks_PPG)- length(pks_ECG_R)) > 1
   [pks_PPG,locs_PPG] = findpeaks(PPG,'MinPeakDistance',round((min_seperation)/(i+1)),'MinPeakHeight',0.1);
   % [pks_PPG_inv,locs_PPG_inv] = findpeaks(PPG_inv,'MinPeakDistance',round((min_seperation)/(i+1)),'MinPeakHeight',0.1);
    display ('PPG peak detection algorithm PPG +1')    
 
end
end
 
%% PTT verfication
if (abs(length(locs_ECG_R)- length(locs_PPG)) > 2)
    display ('Signal error: Wrong PTT calculations')
elseif (abs((length(locs_ECG_R)- length(locs_PPG_inv)) > 2))
    display ('Signal error: Wrong PTT calculations')
end
 
%% Plotting graphs
t_ECG=0:1/fs:(length(ECG)-1)/fs;
t = 1:length(ECG);
figure();
%hold on;
 
 
 
t_PPG=0:1/fs:(length(PPG)-1)/fs;
subplot(2,1,1)
hold on;
plot (t, PPG, '-k')
 
plot(locs_PPG,PPG(locs_PPG),'o','MarkerFaceColor','m');
plot(locs_PPG_inv,PPG(locs_PPG_inv),'o','MarkerFaceColor','c');
legend('PPG signal','PPG peak', 'PPG foot');
xlabel('# of Samples'); ylabel('Amplitude');
grid on;
hold off;
 
 
subplot(2,1,2)
hold on;
plot (t,ECG);
%plot(locs_ECG_Q,ECG(locs_ECG_Q),'rv','MarkerFaceColor','b');
plot(locs_ECG_R,ECG(locs_ECG_R),'<','MarkerFaceColor','r');
 
%plot(locs_ECG_S,ECG(locs_ECG_S),'>','MarkerFaceColor','g');
 
%plot(locs_ECG_T,ECG(locs_ECG_T),'^','MarkerFaceColor','y');
legend('ECG signal','R-peak');
xlabel('# of Samples'); ylabel('Amplitude');
title ('PPG and ECG features');
grid on
hold off;
pks_PPG_inv = -pks_PPG_inv;
 
%% Calculate PTT Max
if (locs_PPG(1)>locs_ECG_R(1))
    for i=1:min(length(locs_ECG_R),length(locs_PPG)) 
        ptt(i) = locs_PPG(i) - locs_ECG_R(i);
    end
    
elseif (locs_PPG(2)>locs_ECG_R(1))
    locs_PPG = locs_PPG(2:length(locs_PPG));
    for i=1:min(length(pks_ECG_R),length(pks_PPG))-1
        ptt(i) = locs_PPG(i) - locs_ECG_R(i);
    end
    
 elseif (locs_PPG(3)>locs_ECG_R(1))
    locs_PPG = locs_PPG(3:length(locs_PPG));
    for i=1:min(length(pks_ECG_R),length(pks_PPG))-1
        ptt(i) = locs_PPG(i) - locs_ECG_R(i);
    end
end
 
%% Caculate PTT Min
if (locs_PPG_inv(1)>locs_ECG_R(1))
    for i=1:min(length(locs_ECG_R),length(locs_PPG_inv)) 
        ptt_min(i) = locs_PPG_inv(i) - locs_ECG_R(i);
    end
    
elseif (locs_PPG_inv(2)>locs_ECG_R(1))
    locs_PPG_inv = locs_PPG_inv(2:length(locs_PPG_inv));
    for i=1:min(length(pks_ECG_R), length(pks_PPG_inv))-1
        ptt_min(i) = locs_PPG_inv(i) - locs_ECG_R(i);
    end
    
 elseif (locs_PPG_inv(3)>locs_ECG_R(1))
    locs_PPG_inv = locs_PPG_inv(3:length(locs_PPG_inv));
    for i=1:min(length(pks_ECG_R), length(pks_PPG_inv))-1
        ptt_min(i) = locs_PPG_inv(i) - locs_ECG_R(i);
    end
end
 
%% Calculate Heart Beat
Heart_Beat = 60 * (length(pks_ECG_R))/(length(ECG)/fs);
display(Heart_Beat);
ptt_mean = mean(ptt);
ptt_time_max = ptt_mean/fs
 
 
%% Calculate PTT
ptt_mean_min = mean(ptt_min);
ptt_time_min = ptt_mean_min/fs
 
%% Calculate R-R interval
for i=2:length(locs_ECG_R)
    RR_interval(i-1) = locs_ECG_R(i) - locs_ECG_R(i-1);
end
 
RR_interval_mean = mean(RR_interval)/fs
 
 
%% 3-pint cross check on caculated values
 
 
%% HR verification
HR_mean = mean(HR);
 
if ((abs(HR_mean-Heart_Beat)/HR_mean)*100 < 5)
    display ('Heart rate calculation correct')
    
elseif ((abs(HR_mean-Heart_Beat)/HR_mean)*100 < 10)
    display ('Heart rate calculation slightly off')
    %(abs(HR_mean-Heart_Beat)/HR_mean)*100
else
    display ('Heart rate calculation completely wrong')
end
 
%% R-peak verification
for i=1:length(locs_ECG_R)-1
    diff_ECG(i) = locs_ECG_R(i+1) - locs_ECG_R(i);
end
 
counter =0;
for i=1:length(diff_ECG)
    if ((mean(diff_ECG)-diff_ECG(i))/mean(diff_ECG)*100 > 5)
        counter = counter +1;    
    end
end
 
if ((counter/length(diff_ECG))*100 > 20 && FLAG == 1)
    display ('ECG signal detected irregular')
elseif (counter/length(diff_ECG))*100 > 20 
    display ('Possible R-peak detection error')
end
 
%% PTT verfication
if (length(locs_ECG_R)- length(locs_PPG) > 2)
    display ('Signal error: Wrong PTT calculations')
elseif ((length(locs_ECG_R)- length(locs_PPG_inv) > 2))
    display ('Signal error: Wrong PTT calculations')
end
 
if (length(locs_PPG) - length(locs_PPG_inv) > 1)
    display ('PCG signal error: Wrong PTT calculations')
end
 
%% HR, SpO2, BP from signals
heart = mean(HR)
SpO2 = mean(SPO2)
Bp_Dia = mean(BP_Dia)
Bp_Sys = mean(BP_Sys)
 
%% ST and DT calculation
if (locs_PPG_inv(1) < locs_PPG(1))
    for (i=1:min(length(locs_PPG_inv),length(locs_PPG))-1)
        ST(i) = locs_PPG_inv(i) - locs_PPG(i);
        DT(i) = locs_PPG_inv(i+1) - locs_PPG(i);
    end
elseif (locs_PPG_inv(1)> locs_PPG(1))
    for (i=1:min(length(locs_PPG_inv),length(locs_PPG))-1)
        ST(i) = locs_PPG_inv(i) - locs_PPG(i+1);
        DT(i) = locs_PPG_inv(i) - locs_PPG(i);
    end
end
ST = ST(1:length(ST)-1);
DT = DT(1:length(DT)-1);
Systolic_upstroke = abs(mean(ST))/fs
Diastolic_time = mean(DT)/fs
 
 
%% Calculate P-P interval
for i=2:length(locs_PPG)
    PP_interval(i) = locs_PPG(i) - locs_PPG(i-1);
end
 
PP_interval_mean = mean(PP_interval)/fs
