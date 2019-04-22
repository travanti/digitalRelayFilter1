#HW #1 this is the first step in simulating a digital relay. At this point the
#program will be able to read a 'no_header' .csv file produced with a simulated
#fault condition from an ATP simulation and will then graph the wavefroms after
#filtering the data from a low pass filter
#written by Nathaniel P. Travanti on 4/7/19
#use of examples available for matlab and octave standard
#documentation used extensively

#first steap read the .csv from the file, note the expectation is that the file
#file name choosen will be avaialable in the file path that this session of octave
#or matlab is running in.

#preamble to clear up run enviornment, user's preference
#clear
clc
format compact

M = csvread('example02a_nohead.csv');

#then the colums for each data value will be read, this assumes the following data
#order from the atp simulation: time | Ia current | Ib current | Ic current |
# Va volts | Vb volts | Vc volts |
time = M(:,1); #first column (A)
Ia = M(:,2);
Ib = M(:,3);
Ic = M(:,4);
Va = M(:,5);
Vb = M(:,6);
Vc = M(:,7); #last column (G)

#Calculate Secondary values for the relay after the PT has reduced the value
VTR = 4308; #Set the VTR value, to be changed for different VTs.
VaSec = Va / VTR;
VbSec = Vb / VTR;
VcSec = Vc / VTR;

#plot for the unfiltered results
figure
plot(time, Ia, time, Ib, time, Ic)
title("Currents Ia, Ib, Ic no filter simulated value")
xlabel("time (seconds)")
ylabel("Current (amps)")
legend("Current Ia","Current Ib","Current Ic")

figure
plot(time, Va, time, Vb, time, Vc)
title("Volts Va, Vb, Vc no filter simulated value")
xlabel("time (seconds)")
ylabel("Volts (volts)")
legend("Volts Va","Volts Vb","Volts Vc")

#implement filter RC for values
pval = 0;
VaLPF = []; #octave requires predeclaration
delt = time(2,1) - time(1,1); #assumes constant sample rate on data set
fc = 500;
RC = 1/(2* pi * fc);
for i = 1:length(VaSec)
  pval = (delt * VaSec(i,1) + RC * pval)/(delt + RC);
  VaLPF = [VaLPF, pval]; 
endfor

pval = 0;
VbLPF = []; #octave requires predeclaration
for i = 1:length(VbSec)
  pval = (delt * VbSec(i,1) + RC * pval) / (delt + RC);
  VbLPF = [VbLPF, pval];
endfor

pval = 0;
VcLPF = []; #again predeclaration
for i = 1:length(VcSec)
  pval = (delt * VcSec(i,1) + RC * pval) / (delt + RC);
  VcLPF = [VcLPF, pval];
endfor

figure
plot(time,VaLPF,time,VbLPF,time,VcLPF)
title("Volts Va, Vb, Vc LPF filter Applied")
xlabel("time (seconds)")
ylabel("Volts (volts)")
legend("Volts VaSec","Volts VbSec","Volts VcSec")

