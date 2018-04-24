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
    clear;
    addpath_MyFcns
    mytic;
    
    load NewColorMapPDF;colormap(NewColorMapPDF);
    figure(1), clf; 
    drawnow;
    
    
%% CREATE SWEEP FUNCTION
    Fs = 48e3; Ts = 1/Fs;
    T  = 10; %seconds
    f1 = 5e3;        
    f2 = 20e3;
    NoiseLev = 0;
    
    fprintf('CREATE: sweep\n');  
    [st,stfilter,time,N] = mysweep_01(Fs, f1, f2, T, NoiseLev);
    fprintf('CREATE: sweep DONE\n');

    
%% CALCULATE FFT OF SWEEP
    norm = 1;
    dB = 1;
    
    fprintf('CALCULATE: fft of sweep\n');
    [fftst, freq, N] = myfft_01(Fs,length(st),norm,dB,st);
    fprintf('CALCULATE: fft of sweep DONE\n');
    phasedist = unwrap(atan2(imag(fftst), real(fftst)));
    
    
%% CALCULATE SPECTROGRAM
    nfft = 1024*8;
    overlapr = 0.5;
    noverlap = nfft*overlapr;
    window = nfft;
    dB = 1;
    dBmax = 0;
    dBmin = -30;
    norm = 1;
    
    fprintf('CALCULATE: spectrogram\n');
    [Ssig,xt,yf] = myspectrogram_01(st,Fs,nfft,norm,dB,dBmin,dBmax);
    fprintf('CALCULATE: spectrogram DONE\n');
    
    
%% CALCULATE AUTO- AND CROSS-CORRELATION
    fdom = 1;
    dB = 1;
    dBmax = 0;
    dBmin = -30;
    
    fprintf('CALCULATE: correlation\n');
    [cor,fcor,cortime,corfreq,N] = myxcorr_01(st,stfilter,Fs,fdom,dB,dBmin,dBmax);
    fprintf('CALCULATE: correlation DONE\n');
    
%% PLOTTING SECTION
%1    
    figure(1), clf, hold off;
    subplot(311),
    plot(time,st,       'r-'), hold on
    plot(time,stfilter, 'b-'), hold on
    V=axis;
    axis([T-70/f2 T V(3) V(4)]);
    title('Exponential Sweep (10 Hz to 20 KHz)');
    xlabel('Time (seconds)');
    
    subplot(312),
    plot(freq, (fftst));
    V=axis; axis([1e3 Fs/2 -30 5]);
    %V=axis; axis([f1-100 Fs/2 -30 0]);
    title('Spectrum of Exponential Sweep (10 Hz to 20 KHz)');
    xlabel('Frequency (Hertz)');

    subplot(313),
    semilogx(freq, phasedist);
    V=axis; axis([f1-100 Fs/2 V(3) V(4)]);
    title('Phase Distribution');
    xlabel('Frequency (Hertz)');
    drawnow; pause(2);
    

%2    
    figure(1),clf, hold off;
    mesh(xt, yf, (Ssig))
    shading('interp');axis('tight');view( [0 90]); 
    V = axis; axis([V(1) V(2) 0 Fs/2]);
    view( [0 90])
    xlabel('Time (seconds)'); ylabel('Frequency (Hz)');
    title('Spectrogram of test signal'); colorbar;
    colormap(NewColorMapPDF);colorbar;
    clear Ssig; 
    drawnow; pause(2);
    
%3
    figure(1), clf, hold off;
    subplot(211),
    plot(cortime, cor);grid on;axis('tight');
    %V=axis; axis([0 2*T V(3) V(4)]);
    title('Auto-Correlation')
    xlabel('Time (seconds)');

    subplot(212),
    plot(corfreq, fcor);grid on;
    title('Spectrum of Auto-Correlation')
    xlabel('Frequency (Hertz)');
    V=axis; axis([0 Fs/2 V(3) V(4)]);
    
    
%% ANALYZE DATA FILE: rawdataT_10_pn_10_48K4.wav    
    filename = 'rawdataT_10_pn_10_48K4.wav';
    [RoomIR, Fs]=audioread(filename);
    RoomIR_left  = RoomIR (:,1)';    
    RoomIR_right = RoomIR (:,2)';
    
    fRoomIR_left  = fft(RoomIR_left);
    fRoomIR_right = fft(RoomIR_right);
    
    lenIR = length(RoomIR_left);
    lenst = length(st);
    lendif = lenIR - lenst;
    dum = zeros(1,lendif);
    
    st = [st, dum];
    fftst = fft(st);
    
    fftS = fRoomIR_left  .* conj(fftst);
    fftA = fRoomIR_right .* conj(fftst);
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
    
    %fftS = fRoomIR_left  .* (fftst);
    %fftA = fRoomIR_right .* (fftst);
    %Maxx = max(fftS,fftA);
    %fftS = fftS./Maxx;
    %fftA = fftA./Maxx;
    %fftst = fftst ./ max(fftst);
    
    
%% CALCULATE SPECTROGRAM OF S AND A
    nfft = 1024*8;
    overlapr = 0.5;
    noverlap = nfft*overlapr;
    window = nfft;
    dB = 1;
    dBmax = 0;
    dBmin = -40;
    norm = 1;
    
    fprintf('CALCULATE: spectrogram\n');
    [SpecIR_S,xt,yf] = myspectrogram_01(RoomIR_left,Fs,nfft,norm,dB,dBmin,dBmax);
    [SpecIR_A,xt,yf] = myspectrogram_01(RoomIR_right,Fs,nfft,norm,dB,dBmin,dBmax);
    fprintf('CALCULATE: spectrogram DONE\n');    
    
    
    
%4
    figure(1), clf, hold off;
    subplot(211),
    plot(timeS,S, 'r');grid on;axis('tight');
    t0 = 20.3;
    dt = 0.03;
    V=axis; axis([t0 t0+dt -5000 5000]);
    %V=axis; axis([t0 t0+dt V(3) V(4)]);
    xlabel('Time (seconds)');
    
    subplot(212),
    plot(timeS,A, 'b');grid on;axis('tight');
    V=axis; axis([t0 t0+dt -5000 5000]);
    %V=axis; axis([t0 t0+dt V(3) V(4)]);
    xlabel('Time (seconds)');
    
    %subplot(413),
    %plot(freqS, fftSdB,'r');grid on;axis('tight');
    %V=axis; axis([1e3 Fs/2 -35 5]);
    %title('HRTF - S')
    %xlabel('Frequency (Hertz)');

    %subplot(414),
    %plot(freqS, fftAdB,'r');grid on;axis('tight');
    %V=axis; axis([1e3 Fs/2 -35 5]);
    %title('HRTF - A')
    %xlabel('Frequency (Hertz)');
    
    figure(2),clf, hold off;
    mesh(xt, yf, SpecIR_S)
    shading('interp');axis('tight');view( [0 90]); 
    V = axis; axis([V(1) V(2) 0 Fs/2]);
    view( [0 90])
    xlabel('Time (seconds)'); ylabel('Frequency (Hz)');
    title('Spectrogram of test signal'); colorbar;
    colormap(NewColorMapPDF);colorbar;
    clear Ssig; 
    drawnow; pause(2);
    
    
mytoc; 



%% THE END  THE END  THE END  THE END  THE END  THE END  THE END  THE END 