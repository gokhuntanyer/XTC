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
%...43          ..
%               ..  You should be sure Fs=48000 and need to wait between
%                   pulses.
%...42          ..  Create audio but this time make record it 48000!
%               ..  I measured T=1 and it is for some reason worst!
%               XX  Waveform T=10 gives a very wide IR so should be much
%                   less. 
%...41          ..  Check/compare the GAIN factor for T= 10 and 1 seconds.
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
    st = st ./ dum;
    ts = ts ./ dum;
    
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
%    Ncap = 12761;
    %474000 474000+Ncap-1 
%    N1 = 474000;
%    N2 = 474000+Ncap;
%    corcap([1:N2-N1]) = autocor([N1:N2-1]);
%    Hw2 = abs( fft(corcap).*fft(corcap) );  % this will normalize output
    
%    timecorcap = [0:length(corcap)-1]./Fs;
%    freqcorcap = [0:length(corcap)-1]./length(corcap).*Fs;
    
    
%% PLOTTING SECTION
%1    
    figure(1), clf, hold off;
    subplot(211),
    plot(time,st,       'r-'), hold on
    plot(time,ts, 'b-'), hold on
    V=axis;
    %axis([0 70/f2 V(3) V(4)]); %axis([T-70/f2 T V(3) V(4)]);
    title('Exponential UP/DOWN NLFM Sweep (100 Hz - 20 KHz)');
    xlabel('Time (seconds)');
    
%    subplot(212),
%    semilogx(freqcorcap, 20.*log10(abs(Hw2)));
    %semilogx(freq, 20.*log10(abs(fftst)));
%    V=axis; axis([f1-100 Fs/2 V(3) V(4)]);
    %V=axis; axis([1e3 Fs/2 -30 5]);
    %V=axis; axis([f1-100 Fs/2 -30 25]); %V(3) V(4)]);
%    title('This is the Spectrum of |H(w)|^2');
    %title('Spectrum of Exponential Sweep (10 Hz to 20 KHz)');
%    xlabel('Frequency (Hertz)');
    %pause(1)
    
%2    
    figure(2),clf, hold off;
    subplot(211),
    mesh(xt, yf, (Ssig_st))
    shading('interp');axis('tight');view( [0 90]); 
    %V = axis; axis([V(1) V(2) 0 Fs/2]);
    view( [0 90])
    xlabel('Time (seconds)'); ylabel('Frequency (Hz)');
    title('Spectrogram of TX-1 signal'); colorbar;
    colormap(NewColorMapPDF);colorbar;
    clear Ssig_st; 
    subplot(212),
    mesh(xt, yf, (Ssig_ts))
    shading('interp');axis('tight');view( [0 90]); 
    %V = axis; axis([V(1) V(2) 0 Fs/2]);
    view( [0 90])
    xlabel('Time (seconds)'); ylabel('Frequency (Hz)');
    title('Spectrogram of TX-2 signal'); colorbar;
    colormap(NewColorMapPDF);colorbar;
    clear Ssig_ts; 
    %drawnow; pause(2);
    
%3
    figure(3), clf, hold off;
    subplot(211),
    autocor = autocor./max(autocor);
    %plot(cortime, croscor), hold on;
    plot(cortime, autocor);grid on;axis('tight');
    delT = 0.05;
    V=axis; axis([T-delT T+delT 2*V(3) 1.1*V(4)]);
    
    %plot(croscor), hold on;
    %plot(autocor);grid on;axis('tight');
    %plot(corcap);grid on;axis('tight');
    %plot(cortime, cor);grid on;axis('tight');
    %V=axis; axis([474000 474000+Ncap-1 2*V(3) 1.3*V(4)]);
    %V=axis; axis([0 2*T V(3) V(4)]);
    title('Auto-Correlation')
    xlabel('Time (seconds)');

    subplot(212),
    
    dum  = 20 .* log10(abs(autocor));
    
    plot(cortime, dum);grid on;
    V=axis; axis([T-delT T+delT -50 1 ]);
    title('Auto-Correlation (dB)')
    
    %plot(corfreq, 20 .* log10( abs(fcor)));grid on;
    %plot(corfreq, 20 .* log10( abs(fcor)));grid on;
    %title('Spectrum of Auto-Correlation')
    %xlabel('Frequency (Hertz)');
    %V=axis; axis([0 Fs/48 V(3) V(4)]);
    
    
   
%% CREATE THE TRANSMITTER SIGNAL 
    %stereo signal!
    
    %Nsilence = 14400;
    Nsilence = Fs;
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
    
    data = 0.8 .* [s2, sLEF, s2, sRIG, s2];
    
    lendata = length(data)
    lendata/Fs/60
    
    timed = [0:lendata-1]./Fs;
    Filnam = 'Data_180228_8_T10sec.wav';
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