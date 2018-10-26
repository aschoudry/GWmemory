% USAGE: psi4toh
%
% AUTHOR: Sean McWilliams
% DATE: 12/01/2011
%
% DESCRIPTION: Converts psi4 to strain.  Has multiple options, and file
% names are hardcoded, which is why this is a script instead of a function
% - it is far from a black box.  It will either filter and integrate in
% frequency space, or integrate in time (using 1st or 2nd order, with 1st
% MUCH faster for now), and filter in frequency.

format long

hintmeth = 3;    %if 1, integrate psi4 in time to 1st order in delta and filter (fastest), if 2, same but integrate to higher order, if 3, use FFI (best)

tstart = 100;    %time in M to start, ie length of time at beginning to cut out
tstop = 50;      %time in M to stop, measured from end
filename='rMPsi4_Sz1_m0p9_Sz2_0p9_q1p5dataClean';
save_location='/home/ashok/gravitational_wave_memory_project/data/Spinning_Sz_antialigned/Memory_data/';
filename_and_loc=strcat(save_location, filename);
wave = load(strcat(strcat('/home/ashok/gravitational_wave_memory_project/data/Spinning_Sz_antialigned/',filename),'.txt'));
timeNR = wave(:,1);
psi_plus = wave(:,2);
psi_cross = wave(:,3);

%idx_ini= dsearchn(timeNR,100);
%idx_fin= dsearchn(timeNR,3700);
%timeNR=timeNR(idx_ini:idx_fin);
%psi_plus=psi_plus(idx_ini:idx_fin);
%psi_cross=psi_cross(idx_ini:idx_fin);

psi_tot = complex(psi_plus,psi_cross);

plot(timeNR,psi_plus)
hold all
plot(timeNR,psi_cross)
clear wave

phase = unwrap(angle(psi_tot));
instfreq = zeros(length(timeNR),1);
instfreq(1)=(phase(2)-phase(1))/(timeNR(2)-timeNR(1));
for i=2:length(timeNR)-1
    instfreq(i) = (phase(i+1)-phase(i-1))/(timeNR(i+1)-timeNR(i-1));
end
instfreq(end)=(phase(end)-phase(end-1))/(timeNR(end)-timeNR(end-1));    %this is GW angular freq (NOT orbital)
figure
plot(timeNR,instfreq)
hold all
clear wave

if hintmeth == 1
    hdot_tot = (timeNR(2)-timeNR(1))*cumsum(psi_tot);
    h_tot = (timeNR(2)-timeNR(1))*cumsum(hdot_tot);
    fminnorm = min(instfreq)*2*(timeNR(2)-timeNR(1));
%    if wfswitch == 3; fminnorm = -fminnorm; end
    [b,a] = butter(3,fminnorm,'high');
    hdot_tot = filter(b,a,hdot_tot);
    hdot_plus = real(hdot_tot);
    hdot_cross = imag(hdot_tot);     
    h_tot = filter(b,a,h_tot);
    h_plus = real(h_tot);
    h_cross = imag(h_tot);
end

if hintmeth == 2
    psi4int = @(t) interp1(timeNR,psi_tot,t,'spline');
    for i=1:length(timeNR)
        hdot_tot(i) = quad(psi4int,timeNR(1),timeNR(i));
    end
    hdotint = @(t) interp1(timeNR,hdot_tot,t,'spline');
    for i=1:length(timeNR)
        h_tot(i) = quad(hdotint,timeNR(1),timeNR(i));
    end
    fminnorm = min(instfreq)*2*(timeNR(2)-timeNR(1));
    [b,a] = butter(3,fminnorm,'high');
    hdot_tot = filter(b,a,hdot_tot);
    hdot_plus = real(hdot_tot);
    hdot_cross = imag(hdot_tot);
    h_tot = filter(b,a,h_tot);
    h_plus = real(h_tot);
    h_cross = imag(h_tot);
end

if hintmeth == 3
    fftsize = 2^(nextpow2(length(timeNR)));
    hfreq = ((1/(2*(timeNR(2)-timeNR(1))))*(0:fftsize/2+1)/(fftsize/2+1))';
    hfreq(fftsize/2+2:fftsize) = -hfreq(fftsize/2:-1:2);
    hfreq(1) = 1e-30; %to avoid NaN
    if min(instfreq) > 0
        fo = dsearchn(hfreq,min(instfreq)/5);
        fe = dsearchn(hfreq,-min(instfreq)/5);
    else
        fo = 1;
        fe = fftsize;
    end
    psi_fft = (timeNR(2)-timeNR(1))*fft(psi_tot,fftsize);
    hdotfft = psi_fft/(timeNR(2)-timeNR(1));
    hdotfft(1:fo) = -hdotfft(1:fo)./hfreq(fo)/2/pi;
    hdotfft(fe:end) = -hdotfft(fe:end)./hfreq(fo)/2/pi;
    hdotfft(fo+1:fe-1) = -hdotfft(fo+1:fe-1)./hfreq(fo+1:fe-1)/2/pi;
    hdot_tot = ifft(hdotfft,fftsize);
    hdot_tot = hdot_tot(1:length(timeNR));
    hdot_plus = real(hdot_tot);
    hdot_cross = imag(hdot_tot);
    hfft = psi_fft/(timeNR(2)-timeNR(1));
    hfft(1:fo) = -hfft(1:fo)./hfreq(fo)./hfreq(fo)/4/pi/pi;
    hfft(fe:end) = -hfft(fe:end)./hfreq(fo)./hfreq(fo)/4/pi/pi;
    hfft(fo+1:fe-1) = -hfft(fo+1:fe-1)./hfreq(fo+1:fe-1)./hfreq(fo+1:fe-1)/4/pi/pi;
    h_tot = ifft(hfft,fftsize);
    h_tot = h_tot(1:length(timeNR));
    h_plus = real(h_tot);
    h_cross = imag(h_tot);

end


maghdot = abs(hdot_tot);
phasehdot = unwrap(angle(hdot_tot));

[peak1, index1] = max(maghdot);
indexleft = dsearchn(maghdot(1:index1),0.8*peak1);
indexright = index1 + dsearchn(maghdot(index1+1:length(maghdot)),0.8*peak1);
% indexleft = index1-1;
% indexright = index1+1;
pcoeff = polyfit(timeNR(indexleft:indexright),maghdot(indexleft:indexright),2);
tpeak = -pcoeff(2)/2/pcoeff(1);
%hpeak = -pcoeff(2)^2/4/pcoeff(1) + pcoeff(3);
phpeak = spline(timeNR(indexleft:indexright),phasehdot(indexleft:indexright),tpeak);

timeNR = timeNR-tpeak;
phasehdot = phasehdot-phpeak;
hdot_plus = maghdot.*cos(phasehdot);
hdot_cross = maghdot.*sin(phasehdot);
hdot_tot = complex(hdot_plus,hdot_cross);



figure
plot(timeNR,hdot_plus)
hold all
plot(timeNR,hdot_cross)
A=[timeNR h_plus h_cross];
%save /home/ashok/gravitational_wave_memory_project/data/hdotNR.dat A -ASCII -double
file_hNR=strcat(strcat(filename_and_loc,'_hdotNR'),'.dat');
save(file_hNR,'A','-ASCII','-double')

magh = abs(h_tot);
phaseh = unwrap(angle(h_tot));

[peak1, index1] = max(magh);
indexleft = dsearchn(magh(1:index1),0.8*peak1);
indexright = index1 + dsearchn(magh(index1+1:length(magh)),0.8*peak1);
% indexleft = index1-1;
% indexright = index1+1;
pcoeff = polyfit(timeNR(indexleft:indexright),magh(indexleft:indexright),2);
tpeak = -pcoeff(2)/2/pcoeff(1);
%hpeak = -pcoeff(2)^2/4/pcoeff(1) + pcoeff(3);
phpeak = spline(timeNR(indexleft:indexright),phaseh(indexleft:indexright),tpeak);

timeNR = timeNR-tpeak;
phaseh = phaseh-phpeak;
h_plus = magh.*cos(phaseh);
h_cross = magh.*sin(phaseh);
h_tot = complex(h_plus,h_cross);

figure
plot(timeNR,h_plus)
hold all
plot(timeNR,h_cross)
A=[timeNR h_plus h_cross];
%save /home/ashok/gravitational_wave_memory_project/data/hNR.dat A -ASCII -double
file_hNR=strcat(strcat(filename_and_loc,'_hNR'),'.dat');
save(file_hNR,'A','-ASCII','-double')

h_plus_mem =cumtrapz(abs(hdot_plus.*hdot_plus))*(timeNR(2)-timeNR(1))/(192*pi);
h_plus_total=h_plus+h_plus_mem;
figure
plot(timeNR,h_plus_mem)
hold all

A=[timeNR h_plus_mem, h_plus_total];
%save /home/ashok/gravitational_wave_memory_project/data/hMemNR.dat A -ASCII -double
file_hNR=strcat(strcat(filename_and_loc,'_hMemNR'),'.dat');
save(file_hNR,'A','-ASCII','-double')