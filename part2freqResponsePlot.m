#Used to plot frequency response of RC filter of HW1 part2
#limited by the annotation abilities of Octave, was finished in matlab to use
#the data cursor there

%Used to plot frequency response of RC filter of HW1 part2
RC = 1/(2*pi*500) %uses fc=500hz

wp = linspace(0,100000,10000);

resp = 1./sqrt(1+(wp.*RC).^2);

figure

semilogy(wp, resp)
grid;
title("RC Filter Frequency Response @ 500 Hz")
xlabel("Frequency (Rads)")
ylabel("Log Magnitude of Attenuation")