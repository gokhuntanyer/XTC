%% ExpoSweep_xx.m
% Generate exponention sweep. Analyze and design cross-correlation 
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
%
%
%
%
%....
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
    clear;clc;tic;
fprintf('\n\nRunning:______________\n');
    load NewColorMapPDF;colormap(NewColorMapPDF);
    figure(1), clf; figure(2), clf; %figure(3), clf;
    %figure(4), clf; figure(5), clf; figure(3), clf; 
    drawnow;
%%  DEFINE:  ANALYSIS GRID
    Fs = 48e3; Ts = 1/Fs;
    vels = 334; %velocity of sound
    T = 15; %seconds
    Lextend = vels * T; %10 meters of analysis window
    N = T*Fs; n = floor( log(T*Fs)/log(2)); N = 2 * 2^n;
    T = N/Fs;
    time = [0:N-1].*Ts;
    %T = 43.6907 seconds
    f1 = 10;    f2 = 20e3;
    w1 = 2*pi*f1; w2 = 2*pi*f2;
    Term1 = log(w2/w1);
    Term2 = exp(time ./T .* Term1)-1;
    arg = w1 .* T .* Term2 ./ Term1;
    st = sin(arg);
    NoisLev = 0;
    st = st + NoisLev .* randn(size(st));
    stfilter([1:N]) = st([N:-1:1]);
    stfilter = stfilter -sum(stfilter)./length(stfilter);
fprintf('Sweep and filter ready:\n');
%1    
    figure(1), clf, hold off;
    subplot(211),
    plot(time,st), V=axis;
    axis([0 30/f1 V(3) V(4)]);
    title('Exponential Sweep (10 Hz to 20 KHz)');
    xlabel('Time (seconds)');
    subplot(212),
    fftst = 20 .* log10(abs(fft(st)));
    fftst = fftst - max(max(fftst));
    freqst = [0:length(st)-1];
    plot(freqst, fftst);
    V=axis; axis([0 Fs/2 -30 0]);
    title('Spectrum of Exponential Sweep (10 Hz to 20 KHz)');
    xlabel('Frequency (Hertz)');

    %% GENERATE SPECTROGRAMS
fprintf('Spectrogram ready:\n');
    nfft = 1024*4;
    overlapr = 0.5;
    noverlap = nfft*overlapr;
    window = nfft;
    %calculate the spectrogram
    Ssig = ( abs( spectrogram(st,window,noverlap,nfft,Fs) ));
    [sizeF, sizeT] = size(Ssig);   %
    yf = [0:sizeF-1] .* Fs ./ nfft;
    xt = [0:sizeT-1] .* (nfft - noverlap) .* Ts;
%2    
    figure(1),clf, hold off;
    mesh(xt, yf, Ssig)
    shading('interp');axis('tight');view( [0 90]); 
    V = axis; axis([V(1) V(2) 0 f2]);
    view( [0 90])
    xlabel('Time (seconds)'); ylabel('Frequency (Hz)');
    title('Spectrogram of test signal'); colorbar;
    colormap(NewColorMapPDF);colorbar;
    clear Ssig; 
fprintf('Spectrogram ready:READY\n');
%% Check auto-correlation of the exponential sweep waveform..
fprintf('Auto-correlation:\n');
    dum = conv(st,stfilter,'full');
    dumtime = toc/60;  %runtime in minutes
    fprintf('\n Auto-correlation calculated at the %4.0f minutes.\n', dumtime);
    lendum = length(dum);
    Nshift = 115200;; %delete the earlier part of the recording
    ACorr = dum;
    %ACorr([1: lendum-Nshift+1]) = dum([Nshift:lendum],1);
    ACorrdb = 10 .* log10(abs(ACorr));
    ACorrdb = ACorrdb - max(max(ACorrdb));
    Acorrtime = [0:length(ACorr)-1].* Ts;
    fftACorr = abs(fft(ACorr));
    fftACorrdb = 10 .* log10(fftACorr);
    fftACorrdb = fftACorrdb - max(max(fftACorrdb));
    freqACorr = [0:length(fftACorr)-1]./length(fftACorr).*Fs;
%3
    figure(1), clf, hold off;
    subplot(311),
    %plot(ACorr);grid on;
    plot(Acorrtime, ACorr);grid on;
    title('Auto-Correlation')
    xlabel('Time (seconds)');
    center = 21.85; %floor(max(Acorrtime)/2); %1048600;
    V=axis;axis([center-2 center+2 -1e5 6e5]);
    %V=axis;axis([V(1) V(2) -1e5 6e5])
    %V=axis;axis([center-0.3*fs center+0.3*fs V(3) V(4)])
    %clear ACorr;
    subplot(312),
    plot(Acorrtime, ACorrdb);grid on;
    title('Auto-Correlation (dB)')
    xlabel('Time (seconds)');
    V=axis; axis([0 V(2) -70 0]);
    %clear ACorrdb;
    subplot(313),
    plot(freqACorr, fftACorrdb);grid on;
    title('Spectrum of the Auto-Correlation Function (dB)')
    xlabel('Frequency (Hertz)');
    V=axis; axis([0 Fs/2 -80 0]);
    %clear fftACorrdb; clear fftACorr;clear freqAorr;
fprintf('Auto-correlation:READY\n');
    
%% Output files
    %sound(st,Fs);
    %p = audioplayer(st, Fs);
    %play(p, [1 N]);
    %play(p, [1 (get(p, 'SampleRate') * 3)]);
    %filename = 'ExpoSweep_01.wav';
    %audiowrite(filename,st,Fs);
    
    %ExpoSweepLEFT  = zeros(2, length(st));
    %ExpoSweepRIGHT = zeros(2, length(st));
    
    %ExpoSweepLEFT (1,:) = st(1,:);
    %ExpoSweepRIGHT(2,:) = st(1,:);
    
    %audiowrite('ExpoSweepLEFT.wav', ExpoSweepLEFT', Fs);
    %audiowrite('ExpoSweepRIGHT.wav',ExpoSweepRIGHT',Fs);
    
    %toc/60  %runtime in minutes (4 minutes at my hope PC.)
%% Check cross-correlation at the receiver: 
% Cross-correlation between input (sound emitted) and output (recordings)
% It is the IR of the system defined by the recording setup
fprintf('Reading file\n');
    [dum,fs] = audioread('Voice_008.wav');
    %[dum,fs] = audioread('Voice_005.wav');
    %lendum = length(dum);   %1276159
    Rec(:) = dum(:,1); clear dum;
    Rec = Rec - sum(Rec)./length(Rec);
fprintf('IR HT (Cross-correlation):\n');
    dum = conv(stfilter',Rec,'full');
    dumtime = toc/60;  %runtime in minutes
fprintf('\n Cross-correlation calculated at the %4.0f minutes.\n', dumtime);
    lendum = length(dum);
    Npeak = 1152193; %keep data in {Npeak-57, Npeak+1437}
    %XCorr = dum;
    XCorr(1,[1: 57+1437+1]) = dum(1,[Npeak-57:Npeak+1437]);
    XCorr = XCorr - sum(XCorr)./length(XCorr);
    roomIR = XCorr; 
    %roomIR = roomIR - sum(roomIR)./length(roomIR);
    save('roomIR','roomIR');
    %   vels=334; lngt = 0.40; timedelay = lngt/vels; %1.2 miliseconds
    %   no_of_samples = floor(timedelay*Fs);  %=57 samples.
    %   no_of_samples_per_meter = 1/vels*Fs; %144 samples/meter
    %TIME:      3 miliseconds/meter  1/vels    = 0.003 seconds
    %SAMPLES: 144 samples/meter      1/vels*Fs = 143.7126 samples
    %          57 samples/40cm     0.4/vels*Fs = 57.4850 samples
    %        1437 samples/10 m      10/vels*Fs = 1.4371e+03 samples
    
    XCorrdb = 10 .* log10( abs(XCorr));
    %XCorrdb = XCorrdb - max(max(XCorrdb));
    lenCorr = length(XCorr);
    corrtime = [0:lenCorr-1].* Ts;
    %
    Hw =                   fft(XCorr);
    Hwabs =           (abs(fft(XCorr)));
    Hwdb = 10 .* log10(Hwabs);
    
    Hwdum = atan2(imag(Hw),real(Hw));
    Hwphase = unwrap(Hwdum)./pi.* 180;
    
    %Hwdb = Hwdb - max(max(Hwdb));
    freq = [0:lenCorr-1]/lenCorr .* Fs;
%2
    figure(2), clf, hold off;
    subplot(211),
    %plot(XCorr);grid on; cntr=1152193; %peak at 1152193
    %axis('tight'); V=axis; axis([cntr-57 cntr+1437 1.1*V(3) 1.1*V(4)]);
    %axis('tight'); %V=axis; axis([1152000 1153000 1.1*V(3) 1.1*V(4)]);
    plot(corrtime, XCorr);grid on;
    title('Impulse Response')
    xlabel('Time (seconds)');
    %V=axis; axis([23.99 24.12 -2.8e4 2.8e4 ]);
    %V=axis; axis([1152000 length(XCorr) -2.8e4 2.8e4 ]);
    axis('tight');

    subplot(212),
    plot(corrtime, XCorrdb);grid on;
    title('Impulse Response (dB)')
    xlabel('Time (seconds)');
    axis('tight');
    %V=axis; axis([0 V(2) -70 0]);
    clear Corrdb;
%3
    figure(3), clf, hold off;
    subplot(311),
    plot(freq, Hwabs);grid on;
    ylabel('Room TF abs(H(w))')
    xlabel('Frequency (Hertz)');
    V=axis; axis([0 Fs/2 0 15e5]);  %real valued
    
    subplot(312),
    plot(freq, Hwdb);grid on;
    ylabel('Room TF abs(H(w)) in dB')
    xlabel('Frequency (Hertz)');
    V=axis; axis([0 Fs/2 0  80]);  %desibels
    
    subplot(313),
    plot(freq, Hwphase);grid on;
    ylabel('Room TF Phase (Degrees)')
    xlabel('Frequency (Hertz)');
    axis('tight'); V=axis; axis([0 Fs/2 V(3) V(4)]);%-180 180 ]);
    
    
    
    
    
    fprintf('IR HT (Cross-correlation):READY\n');
       

    
%%  THE END   THE END   THE END   THE END   THE END   THE END   THE END 
    dumtime = toc/60;  %runtime in minutes (12.43 minutes at my hope PC.)
    fprintf('\n Completed in%4.0f minutes.\n', dumtime);
    fprintf('_______________________\n');


    




