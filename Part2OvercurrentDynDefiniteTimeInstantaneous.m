#code for HW#2 part 2 created here

#written by Nathaniel P. Travanti 4-21-19
#use of examples available for matlab and octave standard
#documentation used extensively
#relies on .csv values of fault from ATP program, provided as part of class.

#read in data
M1 = csvread("case01.csv");
M2 = csvread("case02.csv");
M3 = csvread("case03.csv");
newSamp = 2; #throughout this exercise 16 samples a cycle is sufficent based on input
#fidelity. This would be all dependant on the quality of input signal, faults,
#bus word fidelity in practice.

#funciton to handle the virtual analog to digital conversion of this signal
function out = AtoD(in,newSamp,timeOg)
  fs = newSamp * 60;
  newDelT = 1/fs;
  newTimeRange = (0:newDelT:timeOg(length(timeOg)));
  out = interp1(timeOg,in,newTimeRange);
  
endfunction
#end
#function to handle the sin filter portion
function out = sinFilt(in,newSamp)
  #Reminder only one cycle of signal is needed for filter range
  k = (1:newSamp);
  sinCoeff = (2/newSamp) * sin(2 * pi * k/newSamp);
  out = filter(sinCoeff,1,in);
  
endfunction 
#end
#function to handle the cos filter portion
function out = cosFilt(in,newSamp)
  k = (1:newSamp);
  cosCoeff = (2/newSamp) * cos(2 * pi * k/newSamp);
  out = filter(cosCoeff,1,in); #for why one is used see feedback loop on filter
  
endfunction
#end
#function to read case data for dynamic element
function [tp,tStrt,tStp,Imax] = dynamicCalc(curve,Ipud,TDS,newTimeI1,I1NoHfreq2)
Imax = 0;
tStrt = newTimeI1(1,length(newTimeI1));
tStp = 0;
for k = 1:length(I1NoHfreq2)
  #find inital start time and updated if excceded
  if(I1NoHfreq2(1,k) >= Ipud && I1NoHfreq2(1,k) >= Imax)
  Imax = I1NoHfreq2(1,k); #remove semi-colon for debugging
  #disp('I1max is: '),disp(I1max) #debugging
  endif


  #update start time
  if(newTimeI1(1,k) < tStrt && I1NoHfreq2(1,k) >= Ipud)
  tStrt = newTimeI1(1,k); #remove semi-colon for debugging
  #disp('tStrt: '),disp(tStrt) #debugging
  endif
  
  if(newTimeI1(1,k) > tStp && I1NoHfreq2(1,k) >= Ipud)
    tStp = newTimeI1(1,k); #remove semi-colon for debugging
    #disp('tStp: '),disp(tStp) #debugging
  endif
  
  #reset start time because element has been energized and reenergized
  if(k-1 > 0) #check for history
    if(I1NoHfreq2(1,k-1) < Ipud && I1NoHfreq2(1,k) > Ipud) #reset start value index
    Imax = 0;
    tStrt = newTimeI1(1,k);
    tStp = newTimeI1(1,k);
    #disp('reset occured.') #debugging
    endif
  endif 
endfor

M = Imax/Ipud;
if curve == 1
    tp = TDS * (0.0226 + (0.0104)/((M)^0.02-1));
endif

if curve == 2
  tp = TDS * (0.180 + (5.95)/((M)^2 - 1));
endif

if curve == 3
  tp = TDS * (0.0963 + (3.88)/((M)^2-1));
endif

if curve == 4
  tp = TDS * (0.0352+(5.67)/((M)^2-1));
endif

endfunction
#end
  
time1 = M1(:,1); #time (column one) case one
time2 = M2(:,1);
time3 = M3(:,1);
I1 = M1(:,2);
I2 = M2(:,2);
I3 = M3(:,2);

rmsCst = 1 / sqrt(2);
j = sqrt(-1);

#do processing stepwise... our collective processing lines are not giving satisfactory
#results
newTimeI1 = (0:1/(newSamp * 60):time1(length(time1)));
newTimeI2 = (0:1/(newSamp * 60):time2(length(time2)));
newTimeI3 = (0:1/(newSamp * 60):time3(length(time3)));

I1Dsamp = AtoD(I1,newSamp,time1);
I1x = cosFilt(I1Dsamp,newSamp);
I1y = sinFilt(I1Dsamp,newSamp);

I2Dsamp = AtoD(I2,newSamp,time2);
I2x = cosFilt(I2Dsamp,newSamp);
I2y = sinFilt(I2Dsamp,newSamp);

I3Dsamp = AtoD(I3,newSamp,time3);
I3x = cosFilt(I3Dsamp,newSamp);
I3y = sinFilt(I3Dsamp,newSamp);

I1NoHfreq2 = rmsCst * abs(I1x + j *I1y .* exp(-j * 2 * pi * 60 * newTimeI1));
I2NoHfreq2 = rmsCst * abs(I2x + j *I2y .* exp(-j * 2 * pi * 60 * newTimeI2));
I3NoHfreq2 = rmsCst * abs(I3x + j *I3y .* exp(-j * 2 * pi * 60 * newTimeI3));

#I1 current plotted

figure
subplot(2,2,1)
plot(time1,I1,'r')
grid on
#title('Fault Case #1 plotted')
xlabel('time in sec');
ylabel('I1 (amps) OG');

#figure
subplot(2,2,2)
plot(newTimeI1,I1Dsamp);grid;
xlabel('time in sec');
ylabel('I1 dwn samp');

#figure
subplot(2,2,3)
plot(newTimeI1,I1x);grid;
xlabel('time in sec');
ylabel('I1 cos comp')

#figure
subplot(2,2,4)
plot(newTimeI1,I1NoHfreq2);grid;
ylabel('FIR filt I1');
xlabel('time in sec');

#I2 current plotted
figure
subplot(2,2,1)
plot(time2,I2,'g')
grid on;
xlabel('time in sec');
ylabel('12 (amps) OG');

#figure
subplot(2,2,2)
plot(newTimeI2,I2Dsamp);grid;
xlabel('time in sec');
ylabel('I2 dwn samp');

#figure
subplot(2,2,3)
plot(newTimeI2,I2x);grid;
xlabel('time in sec');
ylabel('I2 cos comp');

#figure
subplot(2,2,4)
plot(newTimeI2,I2NoHfreq2);grid;
ylabel('FIR filt I2');
xlabel('time in sec');

#I3 current plotted

figure
subplot(2,2,1)
plot(time3,I3,'c')
grid on
xlabel('time in sec');
ylabel('13 (amps) OG');

#figure
subplot(2,2,2)
plot(newTimeI3,I3Dsamp);grid;
xlabel('time in sec');
ylabel('I3 dwn samp');

#figure
subplot(2,2,3)
plot(newTimeI3,I3x);grid;
xlabel('time in sec');
ylabel('I3 cos comp');

#figure
subplot(2,2,4)
plot(newTimeI3,I3NoHfreq2);grid;
ylabel('FIR filt I2');
xlabel('time in sec');

#predeclare output collection inputs.
Ipud = 0;
TDS = 0;
IpuDef = 0;
defT = 0;
IpuInst = 0;
curve = 0; #to be confined to IEEE US curves

#gather settings from the user, only needed after I magnitude is found.
#only magnitude of Irms is needed, all of these settings are 'electromechanical'
#or classic relay functions
#dynamic curve eleement settings selection
Ipud = input('Choose a dynamic curve pick up value: ');
disp('dynamic curve pickup is: '), disp(Ipud);
disp('');
while curve <= 0 || curve > 4
  curve = input('Chose a dynamic curve type (1 for U1, 2 for U2, 3 for U3, 4 for U4)')
  switch curve
    case 1
    disp('Curve chosen was: ')
    disp('U1.')
    case 2
    disp('Curve chosen was: ')
    disp('U2.')
    case 3
    disp('Curve chosen was: ')
    disp('U3.')
    case 4
    disp('Curve chosen was: ')
    disp('U4.')
    otherwise
    disp('plese choose a valid value.')  
  
  endswitch 
  disp('');
endwhile

TDS = input('Choose a time dial setting for the curve. Values from 1 to 15 recommended: ');
disp('TDS setting is: '), disp(TDS);
disp('');

IpuDef = input('Choose a definite time current pick up value: ');
disp('definite time pickup is: '), disp(IpuDef);
disp('');

defT = input('Choose a definite time delay for definite element: ');
disp('Definite time chosen: '), disp(defT);
disp('');

IpuInst = input('Choose an instantanteous pickup time (instant) pick up value: ');
disp('Instant time pickup is:'), disp(IpuInst);
disp('');

I1max = 0;
tStrt = newTimeI1(1,length(newTimeI1));
tStp = 0;

#calculate the dynamic curves using homemade predefined function.
#We will include the max current value found because this will be used for
#other elements
[tp1,tStrt1,tStp1,Imax1] = dynamicCalc(curve,Ipud,TDS,newTimeI1,I1NoHfreq2);
[tp2,tStrt2,tStp2,Imax2] = dynamicCalc(curve,Ipud,TDS,newTimeI2,I2NoHfreq2);
[tp3,tStrt3,tStp3,Imax3] = dynamicCalc(curve,Ipud,TDS,newTimeI3,I3NoHfreq2);

#report the operational results to the console...
if tp1 <= (tStp1 - tStrt1)
disp('Bang! The dynamic element has tripped the breaker... for Case 1');
disp('');
else
disp('No dynamic element trip... for Case 1');
disp('');
endif

if tp2 <= (tStp2 - tStrt2)
disp('Bang! The dynamic element has tripped the breaker... for Case 2');
disp('');
else
disp('No dynamic element trip... for Case 2');
disp('');
endif

if tp3 <= (tStp3 - tStrt3)
disp('Bang! The dynamic element has tripped the breaker... for Case 3');
disp('');
else
disp('No dynamic element trip... for Case 3');
disp('');
endif

if  (tStp1 - tStrt1 >= defT && IpuDef <= Imax1)
  disp('Bang! The definite element has tripped the breaker... for Case 1');
  disp('');
else
  disp('No definite element trip... for Case 1');
  disp('');
endif

if  (tStp2 - tStrt2 >= defT && IpuDef <= Imax2)
  disp('Bang! The definite element has tripped the breaker... for Case 2');
  disp('');
else
  disp('No definite element trip... for Case 2');
  disp('');
endif
if  (tStp3 - tStrt3 >= defT && IpuDef <= Imax3)
  disp('Bang! The definite element has tripped the breaker... for Case 3');
  disp('');
else
  disp('No definite element trip... for Case 3');
  disp('');
endif

if  (Imax1 > IpuInst)
  disp('Bang! The instantaneous element has tripped the breaker... for Case 1');
  disp('')
else
  disp('No Instantaneous element trip... for Case 1')
  disp('')
endif

if  (Imax2 > IpuInst)
  disp('Bang! The instantaneous element has tripped the breaker... for Case 2');
  disp('')
else
  disp('No Instantaneous element trip... for Case 2')
  disp('')
endif

if  (Imax3 > IpuInst)
  disp('Bang! The instantaneous element has tripped the breaker... for Case 3');
  disp('')
else
  disp('No Instantaneous element trip... for Case 3')
  disp('')
endif
#fin