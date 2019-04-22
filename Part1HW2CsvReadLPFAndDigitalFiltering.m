#code adaption for HW#2 will be created herein

#written by Nathaniel P. Travanti on 4/10/19
#use of examples available for matlab and octave standard
#documentation used extensively

#needed packages
#below is command to inport package into octave 
#pkg load signal

#read in data
dos('del resulthw3.txt')
diary('resulthw2.txt')
clc
clear all
close all 
M = csvread("prot2_hw2_2108.csv");
#needs to not be the first line of file can be declared anywhere
#matlab requires at end of file

function out = rcfilt(fc,in,t)
  %Function for preforming the input RC filter
  %using good practice every function of the filter will be declared functionally
  %input variable vector of data points before filtering i[]
  %output LPF result signal o[] 
  %time vector from which to extract sample rate t[]
  %cut off frequency for ~0.707 = 1/sqrt(2) point of signal 
  RC = 1/(2* pi * fc);
  pval = 0; #helper var to hold previous value
  delt = t(2,1) - t(1,1);
  o = [];
  for k = 1:length(in)
    pval = (delt * in(k,1) + RC * pval)/(delt + RC); #recursive call with variable
    o = [o, pval];
    out = o;
  endfor
endfunction


time = M(:,1); #used for timing and native sample rate
Ia = M(:,5); #prefilter, primary side.
Ib = M(:,6);
Ic = M(:,7);
Va = M(:,2);
Vb = M(:,3);
Vc = M(:,4);

CTR = 800/5 ; #current transformer ratio primary to secondary
VTR = (500e3 / sqrt(3)) / 67; #voltage transformer ratio primary to secondary

IaSec = Ia / CTR; #convert to secondary current value
IbSec = Ib / CTR; 
IcSec = Ic / CTR;
VaSec = Va / VTR; #convert to secondary voltage value
VbSec = Vb / VTR;
VcSec = Vc / VTR;

#plot unfiltered results current
figure
plot(time,IaSec,time,IbSec,time,IcSec)
xlim([0 0.1])
grid on
title("Currents Ia, Ib, Ic no filter Secondary Value")
xlabel("time (seconds)")
ylabel("Current (amps)")
legend("Current Ia Secondary","Current Ib Secondary","Current Ic Secondary")

#plot unfiltered results voltage
figure
plot(time,VaSec,time,VbSec,time,VcSec)
xlim([0 0.1])
grid on
title("Volts Va, Vb, Vc no filter Secondary Value")
xlabel("time (seconds)")
ylabel("Volts (volts)")
legend("Volts Va Secondary","Volts Vb Secondary","Volts Vc Secondary")

fc = 600; #selected cut off frequency
IaLPF = rcfilt(fc,IaSec,time);
IbLPF = rcfilt(fc,IbSec,time);
IcLPF = rcfilt(fc,IcSec,time);
VaLPF = rcfilt(fc,VaSec,time);
VbLPF = rcfilt(fc,VbSec,time);
VcLPF = rcfilt(fc,VcSec,time);

figure
plot(time,IaLPF,time,IbLPF,time,IcLPF)
xlim([0 0.1])
grid on
title("Current Ia, Ib, Ic after LPF filtering")
xlabel("time (seconds)")
ylabel("Current (amps)")
legend("Current Ia LPF", "Current Ib LPF", "Current Ic LPF")

figure
plot(time,VaLPF,time,VbLPF,time,VcLPF)
xlim([0 0.1])
grid on
title("Voltage Va, Vb, Vc after LPF filtering")
xlabel("time (seconds)")
ylabel("Voltage (Volts)")
legend("Volts Va LPF", "Volts Vb LPF", "Volts Vc LPF")

#implementation of interpolation (resample)
#this makes our signal a digital signal from its virtual "analog signal"
#which is a very high sample rate digital... 16 samples per cycle at 60 cycl a second
#digital

newSamp = 16; #new number of samples per cycle of wave
fs = newSamp * 60; #convert to frequency in terms of time
newDelT = 1/fs;
timeDsmp = (0:newDelT:time(length(time)));
IaFr = interp1(time,IaLPF,timeDsmp); %resampling from original time, of source sig, to new time sample
IbFr = interp1(time,IbLPF,timeDsmp);
IcFr = interp1(time,IcLPF,timeDsmp);
VaFr = interp1(time,VaLPF,timeDsmp);
VbFr = interp1(time,VbLPF,timeDsmp);
VcFr = interp1(time,VcLPF,timeDsmp);

figure
plot(timeDsmp,IaFr,timeDsmp,IbFr,timeDsmp,IcFr)
xlim([0 0.1])
grid on;
title("Current Ia, Ib, Ic after postFIR filtering")
xlabel("time (seconds)")
ylabel("Current (amps)")
legend("Current Ia postFIR", "Current Ib postFIR", "Current Ic postFIR")

figure
plot(timeDsmp,VaFr,timeDsmp,VbFr,timeDsmp,VcFr)
xlim([0 0.1])
grid on;
title("Voltage Va, Vb, Vc after postFIR filtering")
xlabel("time (seconds)")
ylabel("Voltage (Volts)")
legend("Volts Va postFIR", "Volts Vb postFIR", "Volts Vc postFIR")

k = (1:newSamp); #to define our filter as a 1 cycle filter we force coefficent 
sinCoeff = (2/newSamp) * sin(2 * pi * k/newSamp);
cosCoeff = (2/newSamp) * cos(2 * pi * k/newSamp);
Iax = filter(cosCoeff,1,IaFr); #x component of the phasor
Iay = filter(sinCoeff,1,IaFr); #y component of the phasor

Ibx = filter(cosCoeff,1,IbFr); #etc.
Iby = filter(sinCoeff,1,IbFr); 

Icx = filter(cosCoeff,1,IcFr); #ic
Icy = filter(sinCoeff,1,IcFr);

Vax = filter(cosCoeff,1,VaFr); #Va
Vay = filter(sinCoeff,1,VaFr); 

Vbx = filter(cosCoeff,1,VbFr); #Vb
Vby = filter(sinCoeff,1,VbFr);

Vcx = filter(cosCoeff,1,VcFr); #Vc
Vcy = filter(sinCoeff,1,VcFr); 

#calculate out phasor components from x and ys given by filter
j= sqrt(-1); #imaginary
rmsConst = 1/ sqrt(2);
Iap = rmsConst * (Iax + j * Iay) .* exp(-j * 2 * pi * 60 * timeDsmp); #* is similar to convolving the signal in this case
Ibp = rmsConst * (Ibx + j * Iby) .* exp(-j * 2 * pi * 60 * timeDsmp); 
Icp = rmsConst * (Icx + j * Icy) .* exp(-j * 2 * pi * 60 * timeDsmp);
Vap = rmsConst * (Vax + j * Vay) .* exp(-j * 2 * pi * 60 * timeDsmp); #Va voltage
Vbp = rmsConst * (Vbx + j * Vby) .* exp(-j * 2 * pi * 60 * timeDsmp); 
Vcp = rmsConst * (Vcx + j * Vcy) .* exp(-j * 2 * pi * 60 * timeDsmp);


C_radDeg = 180/pi;

#plot the original secondary vs the final result phasor 
#Voltage Va magnitude vs its Original
figure;plot(time,VaSec),grid;hold on;plot(timeDsmp,abs(Vap),'r')
xlim([0 0.1])
title('Original instantaneous Phase-A voltage vs Magnitude of Phasor')
legend('Volts Va secondary','Volts Va phasor Mag')
ylabel('Voltage in volts')
xlabel('Time in seconds')
hold off;

#Current Ia magnitude vs its Original
figure;plot(time,IaSec),grid;hold on;plot(timeDsmp,abs(Iap),'r')
xlim([0 0.1])
title('Original instantaneous Phase-A current vs Magnitude of Phasor')
legend('Current Ia secondary','Current Ia phasor Mag')
ylabel('Current (amps)')
xlabel('Time (secs)')
hold off;

#magnitudes of all results (subplotted)
figure;
subplot(3,2,1)
plot(time,VaSec);hold on;plot(timeDsmp,abs(Vap),'r')
xlabel('Time in seconds')
ylabel('Va Secondary, Va Mag')
hold off;

subplot(3,2,2)
plot(time,VbSec);hold on;plot(timeDsmp,abs(Vbp),'r')
xlabel('Time in seconds')
ylabel('Vb Secondary, Vb Mag')
hold off;

subplot(3,2,3)
plot(time,VcSec);hold on;plot(timeDsmp,abs(Vcp),'r')
xlabel('Time in seconds')
ylabel('Vc Secondary, Vc Mag')
hold off;

subplot(3,2,4)
plot(time,IaSec);hold on;plot(timeDsmp,abs(Iap),'r')
xlabel('Time in seconds')
ylabel('Ia Secondary, Ia Mag')
hold off;

subplot(3,2,5)
plot(time,IbSec);hold on;plot(timeDsmp,abs(Ibp),'r')
xlabel('Time in seconds')
ylabel('Ib Secondary, Ib Mag')
hold off;

subplot(3,2,6)
plot(time,IcSec);hold on;plot(timeDsmp,abs(Icp),'r')
xlabel('Time in seconds')
ylabel('Ic Secondary, Ic Mag')
hold off;


#angles of all results (subplotted)
figure;
subplot(3,2,1)
plot(timeDsmp,angle(Vap)*C_radDeg,'r') #no primary value to compare to
xlabel('Time in seconds')
ylabel('Va angle in deg')

subplot(3,2,2)
plot(timeDsmp,angle(Vbp)*C_radDeg,'r')
xlabel('Time in seconds')
ylabel('Vb angle in deg')

subplot(3,2,3)
plot(timeDsmp,angle(Vcp)*C_radDeg,'r')
xlabel('Time in seconds')
ylabel('Vc angle in deg')

subplot(3,2,4)
plot(timeDsmp,angle(Iap)*C_radDeg,'r')
xlabel('Time in seconds')
ylabel('Ia angle in deg')

subplot(3,2,5)
plot(timeDsmp,angle(Ibp)*C_radDeg,'r')
xlabel('Time in seconds')
ylabel('Ib angle in deg')

subplot(3,2,6)
plot(timeDsmp,angle(Icp)*C_radDeg,'r')
xlabel('Time in seconds')
ylabel('Ic angle in deg')

#print out results to the command window
disp('Homework 2 results:')
disp('')
disp('Voltage and current phasors as read by the relay are the following (at steady steate):')
disp(['Va = ' num2str(abs(Vap(length(Vap)-1))) ' V /angle ' num2str(angle(Vap(length(Vap)-1)) *C_radDeg) ' degrees'])
disp(['Vb = ' num2str(abs(Vbp(length(Vbp)-1))) ' V /angle ' num2str(angle(Vbp(length(Vbp)-1)) *C_radDeg) ' degrees'])
disp(['Vc = ' num2str(abs(Vcp(length(Vcp)-1))) ' V /angle ' num2str(angle(Vcp(length(Vcp)-1)) *C_radDeg) ' degrees'])
disp(['Ia = ' num2str(abs(Iap(length(Iap)-1))) ' A /angle ' num2str(angle(Iap(length(Iap)-1)) *C_radDeg) ' degrees'])
disp(['Ib = ' num2str(abs(Ibp(length(Ibp)-1))) ' A /angle ' num2str(angle(Ibp(length(Ibp)-1)) *C_radDeg) ' degrees'])
disp(['Ic = ' num2str(abs(Icp(length(Icp)-1))) ' A /angle ' num2str(angle(Icp(length(Icp)-1)) *C_radDeg) ' degrees'])


#traditional fault calculation for comparison (without capacitance)
# 500kv source, 60Hz (american)
# Zs1 = Zs2 = Zs0 = 2+30j ohms (sequence subs)
# ZL1 = ZL2 = 0.073 = j 0.8 ohm/mi [@20 miles of TL0
# ZL0 = 0.1 + j2.6 ohm/mi [0 sequence]
# neglecting capacitance for model between TLs
#find line voltage
E = 500e3/sqrt(3);

#find sequence impedances of source, for relay side
Zs1 = 2+30*j;
Zs0 = Zs1; #sources impedances are equal in this case
Zs2 = Zs1;
Zl1 = (0.073 + j * 0.8) * 20;
Zl2 = Zl1;
Zl0 = (0.1+j*2.6) * 20;
If1 = E/(Zs1+Zl1+Zs2+Zl2+Zs0+Zl0); #single line to ground fault impedance network.
If2 = If1; If0 = If1; #in sequence network all sequences are in series so current is the same
Vq1 = E - If1 * Zs1;  #source only has effect on adding to postivie sequence
Vq2 = -Zs2 * If2; #neg is return fault voltage
Vq0 = -Zs0 * If0; #zero also
a = 1* exp(j*120/C_radDeg); #by def, radian values for octave cartesian form
IaNoC = (If0 + If1 + If2)/CTR; #a is first forward sequence phasor
IbNoC = (If0 + (a^2) * If1 + a * If2) / CTR;
IcNoC = (If0 + a * If1 + (a^2) * If2) / CTR; #reminder 0 sequence has no phase angle

VaNoC = (Vq0 + Vq1 + Vq2) / VTR;
VbNoC = (Vq0 + (a^2) * Vq1 + a * Vq2) / VTR;
VcNoC = (Vq0 + a * Vq1 + (a^2) * Vq2) / VTR;

disp('')
disp('Compare to a traditional fault calculation, neglecting capacitance:')
disp(['Va = ' num2str(abs(VaNoC)) ' V /angle ' num2str(angle(VaNoC)*C_radDeg) ' degrees'])
disp(['Vb = ' num2str(abs(VbNoC)) ' V /angle ' num2str(angle(VbNoC)*C_radDeg) ' degrees'])
disp(['Vc = ' num2str(abs(VcNoC)) ' V /angle ' num2str(angle(VcNoC)*C_radDeg) ' degrees'])
disp(['Ia = ' num2str(abs(IaNoC)) ' A /angle ' num2str(angle(IaNoC)*C_radDeg) ' degrees'])
disp(['Ib = ' num2str(abs(IbNoC)) ' A /angle ' num2str(angle(IbNoC)*C_radDeg) ' degrees'])
disp(['Ic = ' num2str(abs(IcNoC)) ' A /angle ' num2str(angle(IcNoC)*C_radDeg) ' degrees'])
diary off