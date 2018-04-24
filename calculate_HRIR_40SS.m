%% calculate_HRIR_xx.m
%
%
%
%% create_sweep_xx.m
% Generate exponention sweep. 
% S. G. Tanyer
% 180207
%%
% Analyze and design cross-correlation 
% properties for differing setting parameters.
% TX and RX implementation of Farina's formula.
% S. Gokhun Tanyer
% 171231
%% HISTORY
%
%
%
%
%
%               
%               SS  Fig.1 is the autocorr for TX signal.
%                   Fig.2 short cct signal's cross-correlation
%                   Fig.3 is the coupling between wires.
%                   resolution is found to be 0.00005 seconds = 1.7 cm.
%...40  180227  ..  Check the shorted data for xcorr.
% calculate_HRIR_40.m

%...40
%               NOTE! LEFT and RIGHT TX signals should have the same H(w).
%                   Saved to Data_180225.wav
%                   10 LEFT/UP, then 10 RIGHT/DOWN chirp pulses completed.
%...39          SS  Left is left, right is right. 
%               OK  data saved. Everything works just left and right needs to be switched.
%...38          ..  Now save in stereo. 
%               ..  saved waveform into audio
%               ..  OK I checked if st ts flips at t=103.6th second.
%...37          ..  PW=10sec. Psilence=0.03sec (30ms)(10m)(1440 samples)
%                   Now create the audio.
%               SS  UP and DOWN chirps are ready! Cross-correlations are very low SUPER!
%...36          ..  Create UP/DOWN codes and check the cross-corelation.
%                   This audio will be played once on Laptop for all.
%...35  180225  ..  Will create 10X2 pulses NLFM Up/Down for LEFT/R.
%...34              ..  I will try new recording
%                       Kirks-03-Loudspeaker-3cut-180221.wav
%                       I need fcor of length: 12761
%                   ..  now calculate H(w) with the same length as Y(w).
%                   ..  HRTF's show the spectrum of the sweep. Should I
%                       normalize out(w) = Y(w)/|H(w)|^2 where |H(w)|^2
%                       is the FFT of auto-correlation h[n]
%                   OK  Find a window wide and clean. HRTF's calculated.
%...33  180214      ..  Calculate HTRF for T=10.
%                   ..  Now find HRTF, the frequency spectrum.
%...32              ..  T=2 case.
%                   SS  T=10 case.This is comparison between T = 2 and 10 seconds.
%                       Reading both: rawdataT_10_pn_10_48K4.wav and
%...31              OK  Data: rawdataT_2_pn_10.wav
%                       Now try a much shorter pulse, say T=2.
%...30              SS  sweep recording shows high background.
%                   ..  kinda impulse response observed. But there is a typo
%...29              ..  Reading rawdataT_10_pn_10_48K4.wav. It is 48000!
%...28              OK  Analyze the phase of stfilter.
%                   This file calculates the Sweep function, 
%                   analyzes it using MYFFT, MYSPECTROGRAM and
%                   calculates the cross-correlation betwee st and stfilter.
%                   Correlation can be done in time or frequency domain.
%
%...27              SS  Super. It is SOO fast to convolve in frequency. 
%                       Now you have a tool to compare computation time.
%...26              OK  Now use functions for simplification of the code.
%...                ..  Now go for myautocor and myxcor, tomorrow.
%...25              OK  myspectrogram works too. dB norm etc all.
%...24              ..  mysweep and myfft is used.
%                   ..  Now use my functions for simplification.
%...23              ..  I am checking again it is 48K.
%                   ..  Works but data is long T=1 takes more than 6
%                       minutes. and FS is NOT 48000, it is 44100!
%...22              ..  Study file rawdataT_10_pn_10_48K4.wav for T=10.
%                   ..  So before the second recording find S and A.
%                   ..  I played and recorded. But the right speaker gain
%                       was so low. I need to calibrate them to be equal??
%                       Recording is completed. Now, play and record.
%...21              SS  Left channel is recorded using this version.
%...20              OK  Saving MONO. Now go to LEFT channel only. 
%                       Need pulse train of pulses now. 
%                   OK  saving to Filename dataNLFM_mono_T_1_5000_20000.wav
%                       Save an example file.
%...19              ..  Now save waveform and the corresponding filter.
%                       for f1 = 5e3; f2 = 20e3; and for
%                       T = 1, 2, 3, 4, 5, 6, 7, 8, 9, 10
%                       This is the version to determine  PW and f1/f2!
%                   SS  Check cross-correlation. I need to have largest gain.
%...18              ..  Now I need to save st and stfilter both. 
%                       (st in wav and stfilter in .mat format.)
%                       Frequency values and duration are CORRECTED.
%...17              OK  DC component of st and stfilter cleared.
%...16              OK  stfilter, and spectrogram is working back again.
%                       Need to get rid of powers of 2.
%                   OK  Now OK. I corrupted the file. Now corrected back.
%                       Sweep not OK. Passes beyond f2.
%...15              ..  sweep OK. Check spectrogram.
%                   OK  Now we start again. Need to stop forcing powers of
%...14                  2 for the length since it is going to be analog.
%create_sweep_14.m  ..  Copy of ExpoSweep_14.m
%%%
%                   OK  stfilter and Voice is saved. Now deconvolution.
%                       Why Voice is shorter than stfilter??
%....14             ..  This is branch of 13. Will save stfilter and Voice
%                       to be deconvoluted on the next version.
%                   ..  will do deconvolution for IR
%....13             OK  Voice is windowed and shortened. Then convolved
%                       with stfilter.
%                   SS  DC corrected, and phase plot added version.
%                   OK  A1) Recorded sound file has DC component. Cleared.
%                       Did not fix it. 
%                   xx  stfilter now has no DC component.
%                   ..  Kevin's comments. Q1) I have a DC component.
%                       A1) I checked stfilter and found DC. Cleared it.
%....12             OK  Clean version. No file writing. Just analysis.
%....11             OK  Write ExpoSweepLEFT.wav and ExpoSweepLEFT.wav files
%                   SS  Saving roomIR.
%....10             ..  Save IR to a file to be used in convolution of music recordings.
%....09             OK  Capture IR of the room.%                   OK  This recording worked. Now we are ready for the
%                       real recording.
%....08             ..  Now use a closer recording Voice_008.wav on touchpad
%                   OK  Away about 1 meters below PC level VoiceVoice_005.wav
%....07     180101  ..  Find IR using the recorded exposweep file.
%....06             SS  File play/write.
%                   SS  Saving.
%....05             OK  Worked. You are done for now.
%                   ..  Correlation is not good. Flip st filter and try again.
%....04             ..  Now check cross-correlation.
%                   OK  Previous formulas had typos which I could not
%                       figure out. But Farina's did work. This is a working
%                       version.
%....03             ..  Try Farina's.
%                   ..  Same problem. Generated freqs are very small.
%....02             ..  try wikipedia.
%....01     171231  xx  Carson's formula does not give correct freqs.
%ExpoSweep_01.m  


%%  INTRO: Initialization 
    clear;clc;
    addpath_MyFcns
    mytic;
    
    load NewColorMapPDF;colormap(NewColorMapPDF);
    figure(1), clf; 
    drawnow;
    
    
%% CREATE SWEEP FUNCTION
    Fs = 48e3; Ts = 1/Fs;
    T  = 10; %seconds
    f1 = 1e2;        
    f2 = 20e3;
    NoiseLev = 0;
    
    fprintf('CREATE: sweep\n');  
    % st is the up chirp, ts is the down chirp.
    [st,ts,time,N] = mysweep_01(Fs, f1, f2, T, NoiseLev);
    %fprintf('CREATE: sweep DONE\n');
    %normalize st and ts to 1.00
    dum = max(max([st,-st,ts,-ts]));
    st = 0.8 .* st ./ dum;  %not to overload.
    ts = 0.8 .* ts ./ dum;
    
%% CALCULATE FFT OF SWEEP
    norm = 1;
    dB = 1;
    
    [fft_st, freq, N] = myfft_01(Fs,length(st),norm,dB,st);
    [fft_ts, freq, N] = myfft_01(Fs,length(ts),norm,dB,ts);
    
    phase_st = unwrap(atan2(imag(fft_st), real(fft_st)));
    phase_ts = unwrap(atan2(imag(fft_ts), real(fft_ts)));
    
    
%% CALCULATE SPECTROGRAM
    nfft = 1024*8;
    overlapr = 0.5;
    noverlap = nfft*overlapr;
    window = nfft;
    dB = 1;
    dBmax = 0;
    dBmin = -30;
    norm = 1;
    
    [Ssig_st,xt,yf] = myspectrogram_01(st,Fs,nfft,norm,dB,dBmin,dBmax);
    [Ssig_ts,xt,yf] = myspectrogram_01(ts,Fs,nfft,norm,dB,dBmin,dBmax);
    
    
%% CALCULATE AUTO- AND CROSS-CORRELATION
    fdom = 1;
    dB = 0;
    dBmax = 0;
    dBmin = -30;
    
    [autocor,fcor,cortime,corfreq,N] = myxcorr_01(st,ts,Fs,fdom,dB,dBmin,dBmax);
    [croscor,fcor,cortime,corfreq,N] = myxcorr_01(st,st,Fs,fdom,dB,dBmin,dBmax);
    
    %here we need a shorter FFT of cor. So here, we need to time window cor.
    %shorter cor will yield shorter fcor, length of  12761 is necessary.
    Ncap = 12761;
    %474000 474000+Ncap-1 
    N1 = 474000;
    N2 = 474000+Ncap;
    corcap([1:N2-N1]) = autocor([N1:N2-1]);
    Hw2 = abs( fft(corcap).*fft(corcap) );  % this will normalize output
    
    timecorcap = [0:length(corcap)-1]./Fs;
    freqcorcap = [0:length(corcap)-1]./length(corcap).*Fs;
    
    
%% PLOTTING SECTION
%1    
    figure(1), clf, hold off;
    subplot(211),
    plot(time,st,       'r-'), hold on
    plot(time,ts, 'b-'), hold on
    V=axis;
    axis([0 70/f2 V(3) V(4)]); %axis([T-70/f2 T V(3) V(4)]);
    title('Exponential UP/DOWN NLFM Sweep (100 Hz - 20 KHz)');
    xlabel('Time (seconds)');
    
    subplot(212),
    semilogx(freqcorcap, 20.*log10(abs(Hw2)));
    %semilogx(freq, 20.*log10(abs(fftst)));
    V=axis; axis([f1-100 Fs/2 V(3) V(4)]);
    %V=axis; axis([1e3 Fs/2 -30 5]);
    %V=axis; axis([f1-100 Fs/2 -30 25]); %V(3) V(4)]);
    title('This is the Spectrum of |H(w)|^2');
    %title('Spectrum of Exponential Sweep (10 Hz to 20 KHz)');
    xlabel('Frequency (Hertz)');
    %pause(1)
    
%2    
    figure(2),clf, hold off;
    subplot(211),
    mesh(xt, yf, (Ssig_st))
    shading('interp');axis('tight');view( [0 90]); 
    V = axis; axis([V(1) V(2) 0 Fs/2]);
    view( [0 90])
    xlabel('Time (seconds)'); ylabel('Frequency (Hz)');
    title('Spectrogram of TX-1 signal'); colorbar;
    colormap(NewColorMapPDF);colorbar;
    clear Ssig_st; 
    subplot(212),
    mesh(xt, yf, (Ssig_ts))
    shading('interp');axis('tight');view( [0 90]); 
    V = axis; axis([V(1) V(2) 0 Fs/2]);
    view( [0 90])
    xlabel('Time (seconds)'); ylabel('Frequency (Hz)');
    title('Spectrogram of TX-2 signal'); colorbar;
    colormap(NewColorMapPDF);colorbar;
    clear Ssig_ts; 
    %drawnow; pause(2);
    
%3
    figure(1), clf, hold off;
    subplot(211),
    %plot(cortime, croscor), hold on;
    plot(cortime, autocor);grid on;axis('tight');
    delT = 0.05;
    V=axis; axis([10-delT 10+delT V(3) V(4)]);
    %V=axis; axis([10-delT 10+delT 2*V(3) 1.3*V(4)]);
    
    %plot(croscor), hold on;
    %plot(autocor);grid on;axis('tight');
    %plot(corcap);grid on;axis('tight');
    %plot(cortime, cor);grid on;axis('tight');
    %V=axis; axis([474000 474000+Ncap-1 2*V(3) 1.3*V(4)]);
    %V=axis; axis([0 2*T V(3) V(4)]);
    title('Auto-Correlation')
    xlabel('Time (seconds)');

    dum = 20 .* log10( abs(autocor) );
    dum = dum - max(dum);
    
    subplot(212),
     plot(cortime, dum), hold on;
    V=axis; axis([10-delT 10+delT -50 V(4)]); grid on
    %plot(cortime, , hold on;
    %plot(cortime, 20 .* log10( autocor));grid on;axis('tight');
    %delT = 0.05;
    %V=axis; axis([10-delT 10+delT 45  V(4)]);
    %V=axis; axis([10-delT 10+delT 45  V(4)]);
    %V=axis; axis([10-delT 10+delT 2*V(3) 1.3*V(4)]);
    
    %plot(croscor), hold on;
    %plot(autocor);grid on;axis('tight');
    %plot(corcap);grid on;axis('tight');
    %plot(cortime, cor);grid on;axis('tight');
    %V=axis; axis([474000 474000+Ncap-1 2*V(3) 1.3*V(4)]);
    %V=axis; axis([0 2*T V(3) V(4)]);
    title('Auto-Correlation')
    xlabel('Time (seconds)');

    
    
    
    
    %subplot(212),
    %plot(corfreq, 20 .* log10( abs(fcor)));grid on;
    %plot(corfreq, 20 .* log10( abs(fcor)));grid on;
    %title('Spectrum of Auto-Correlation')
    %xlabel('Frequency (Hertz)');
    %V=axis; axis([0 Fs/48 V(3) V(4)]);
    
    
    
    
    
    
%% ANALYZE DATA FILE:  
     filename = 'Record_180227_3.wav';
    %filename = 'rawdataT_2_pn_10.wav';
    %filename = 'rawdataT_10_pn_10_48K4.wav';
    [RoomIR, Fs]=audioread(filename);
    RoomIR_left  = RoomIR (:,1)';    
    RoomIR_right = RoomIR (:,2)';
    
    %figure(1), clf
    %plot(RoomIR_right)
    
    fRoomIR_left  = fft(RoomIR_left);
    fRoomIR_right = fft(RoomIR_right);
    
    
    lenIR = length(RoomIR_left);
    lenst = length(st);
    lendif = lenIR - lenst;
    dum = zeros(1,lendif);
    
    st = [st, dum];
    fftst = fft(st);

    ts = [ts, dum];
    fftts = fft(ts);
        
    fftS = fRoomIR_left  .* conj(fftts);
    fftA = fRoomIR_right .* conj(fftts);
    %fftS = fRoomIR_left  .* conj(fftst);
    %fftA = fRoomIR_right .* conj(fftst);
    
    fftSdB = 20 .* log10( abs(fftS));
    fftAdB = 20 .* log10( abs(fftA));
    Maxx = max(max(fftSdB,fftAdB));
    fftSdB = fftSdB - Maxx;
    fftAdB = fftAdB - Maxx;
    
    
    S = ifft(fftS);
    A = ifft(fftA);
    lenS = length(S);
    timeS = [0:lenS-1]./Fs;
    freqS = [0:lenS-1]./lenS .* Fs;
    
    
%% CALCULATE SPECTROGRAM OF S AND A
    nfft = 1024*8;
    overlapr = 0.5;
    noverlap = nfft*overlapr;
    window = nfft;
    dB = 1;
    dBmax = 0;
    dBmin = -40;
    norm = 1;
    
    [SpecIR_S,xt,yf] = myspectrogram_01(RoomIR_left,Fs,nfft,norm,dB,dBmin,dBmax);
    [SpecIR_A,xt,yf] = myspectrogram_01(RoomIR_right,Fs,nfft,norm,dB,dBmin,dBmax);
    
    
%4
    figure(2), clf, hold off;
    subplot(211),
    plot(timeS,abs(S), 'r');grid on;axis('tight'); hold on
    plot(timeS,abs(S), 'r+');grid on;axis('tight'); hold on
    
    t0 = 36.4633; dt = 0.0005; %2/344; %This is for T=10;
    %t0 = 3.872; dt = 20/344; %This is for T=2;
    %t0 = 20.306; dt = 2/344; %This is for T=10;
    V=axis; axis([t0 t0+dt V(3) 3.5e4]); %This is for T=10;
    title('Cross-correlation Short CCT Test')
    xlabel('Time (seconds)');
    
    
    dum = 20 .* log10(abs(S));
    dum = dum - max(dum);
    
    subplot(212),
    plot(timeS,dum, 'r');grid on;axis('tight'); hold on
    plot(timeS,dum, 'r+');grid on;axis('tight'); hold on
    V=axis; axis([t0 t0+dt -50 0]); %This is for T=10;
    title('Cross-correlation Short CCT Test (dB)')
    xlabel('Time (seconds)');
    
    
    
    
    
    
    
    
    figure(3),clf
    %subplot(212),
    t0 = 36.45; dt = 0.05; %2/344; %This is for T=10;
    plot(timeS,abs(A), 'b');grid on;axis('tight');hold on
    plot(timeS,abs(A), 'b+');grid on;axis('tight');
    V=axis; axis([t0 t0+dt V(3) 5000]); %This is for T=10;
    %V=axis; axis([t0 t0+dt -800 800]); %This is for T=2;
    %V=axis; axis([t0-1.5 t0+dt+1.5 -3500 3500]); %This is for T=10;
    %V=axis; axis([t0 t0+dt -35000 35000]); %This is for T=10;
    title('The OPPOSITE-SIDE channel response: HRIR_A')
    xlabel('Time (seconds)');

    
    
    
    
    
    
    dur
%2    
    figure(2),clf, hold off;
    mesh(xt, yf, SpecIR_S)
    shading('interp');axis('tight');view( [0 90]); 
    %t0 = 20.306; dt = T; %This is for T=10;
    %V = axis; axis([t0-0.5 t0+dt+0.5 0 Fs/2]);
    %V = axis; axis([V(1) V(2) 0 Fs/2]);
    view( [0 90])
    xlabel('Time (seconds)'); ylabel('Frequency (Hz)');
    title('Spectrogram of the signal return'); colorbar;
    colormap(NewColorMapPDF);colorbar;
    clear Ssig; 
    drawnow; pause(2);
    
%% CAPTURE HRIR in window then calculate HRTF

t1  = t0-0.01;
t2  = t0+dt+0.25;
it1 = floor(t1*Fs);
it2 = ceil (t2*Fs);

Scap(1,[1:it2-it1+1]) = S(1,[it1:it2]);
Acap(1,[1:it2-it1+1]) = A(1,[it1:it2]);

fScap = fft(Scap);
fAcap = fft(Acap);

% this did not worked. Normalization did not give any meaningful result??
fScap = fScap ./Hw2;
fAcap = fAcap ./Hw2;
% this did not worked


fScap = 20 .*log10(abs(fScap));
fAcap = 20 .*log10(abs(fAcap));

lenS = length(fScap);
timeS = [0:lenS-1]./Fs;
freqS = [0:lenS-1]./lenS.*Fs;
%1
    figure(2), clf, hold off;
    subplot(211),
    plot(timeS, Scap, 'r');grid on;axis('tight');
    %t0 = 3.872; dt = 20/344; %This is for T=2;
    %t0 = 20.306; dt = 2/344; %This is for T=10;
    %V=axis; axis([t0 t0+dt -80 80]); %This is for T=2;
    %V=axis; axis([t0-1.5 t0+dt+1.5 -3500 3500]); %This is for T=10;
    title('The time windowed: HRIR_S')
    xlabel('Time (seconds)');
    
    subplot(212),
    plot(timeS, Acap, 'b');grid on;axis('tight');
    %V=axis; axis([t0 t0+dt -800 800]); %This is for T=2;
    %V=axis; axis([t0-1.5 t0+dt+1.5 -3500 3500]); %This is for T=10;
    %V=axis; axis([t0 t0+dt -35000 35000]); %This is for T=10;
    title('The time windowed: HRIR_A')
    xlabel('Time (seconds)');

%2
    figure(3), clf, hold off;
    subplot(311),
    plot(freqS, fScap, 'r');grid on;axis('tight');
    V=axis; axis([1000 Fs/2 V(3) V(4)]); %This is for T=2;
    title ('Normalized: fScap/|H(w)|^2')
    xlabel('Frequency (Hertz)');
    
    subplot(312),
    plot(freqS, fAcap, 'r');grid on;axis('tight');
    V=axis; axis([1000 Fs/2 V(3) V(4)]); %This is for T=2;
    title ('Normalized: fAcap/|H(w)|^2')
    xlabel('Frequency (Hertz)');








    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    durr
   
%% CREATE THE TRANSMITTER SIGNAL 
    %stereo signal!
    
    Nsilence = 14400;
    ss = zeros(1,Nsilence);
    s2 = zeros(2,Nsilence);
    
    pL (2,:) = [st ss];
    pL (1,:) = zeros(1, length([st ss]));
    
    pR (1,:) = [st ss];
    pR (2,:) = zeros(1, length([st ss]));
    
    %pR (1,:) = [ts ss];
    %pR (2,:) = zeros(1, length([ts ss]));
    
    
    sLEF = [pL pL pL pL pL pL pL pL pL pL];
    sRIG = [pR pR pR pR pR pR pR pR pR pR];
    
    data = [s2, sLEF, s2, sRIG, s2];
    
    lendata = length(data)
    lendata/Fs/60
    
    timed = [0:lendata-1]./Fs;
    Filnam = 'Data_180226.wav';
    audiowrite(Filnam,data',Fs);
    fprintf('AUDIOWRITE: to file\n');  
    
    figure(4), clf
    subplot(211), 
    plot(timed, data(1,:), 'k-'), hold on
    plot(timed, data(1,:), 'k+'), axis('tight');
    %axis([103-0.005 103 -1 1]);
    
    subplot(212), 
    plot(timed, data(2,:), 'k-'), hold on
    plot(timed, data(2,:), 'k+'), axis('tight');
    %axis([103.6 103.6+0.005 -1 1]);
    
    
    
mytoc; 



%% THE END  THE END  THE END  THE END  THE END  THE END  THE END  THE END 