%% calculate_HRIR_xx.m
%
% Load recordings, compress pulse returns using the AmpESS sweep signal.
% Find S and A for the LEFT ear. Find HRIR/HRTF's
% Get ready to work for the inverse filters.
% S.G. Tanyer, 180307 Victoria
%

%% HISTORY
%
%
%               ..  Now use time window function for capturing HRIRs. Select 
%                   the left side to be sharp, and right side hanning etc.
%...77          ..  Kinda inprovement. 1) I multiplied by the linear term.
%                   But it shifted IR to the left and warped from the right.
%               ..  Then subtract this term.
%...76          OK  I corrected the code. Now back to linear phase term.
%               ..  This is again Kevin's semilogx phase plots. Really ugly looking plots.
%                   Now I am correcting them.
%...75          ..  This is Kevin's plots. semilogx for phase looks ugly.
%...74  180312  ..  (1) Zero the phase of shorted cable
%                   (2) Plot phase responses in terms of time delay, too.
%                   (3) Plan the required measurement setup/scenario items. 
%                   Have zero phase response for the shorted cable. 
%               SS  labels are updated. Saved version.
%...73          ..  I am providing information on the plots to answer Kevin's questions.
%...72          OK  You don't do funny things to HRTF. It is just FFT of
%                   the HRIR, so don't worry, you are right.
%               ..  OK seems super, did you try to see what HRIR is when it is IFFTed.
%               SS  Super!! S and A are compared with the shorted data.
%                   S and A are normalized using the same constant.
%                   ALL pulse measurements are normalized using the same constant!
%                   Find ALL the peak values of Shorted, then use them to compare
%...71          ..  Check if it is really true.
%...70          !!  S and A is exact. NO DIFFERENCE?? Too bad.
%...69          OK  Ohhh. Now we read new Shorted data and we have A=S.
%               ..  Do shorted recording and come back.
%               xx  Instead of using record of shorted, I directly used ESS
%...67          ..  Re-record shorted and see if timing is changing or constant?
%               xx  I don't have shortedA pulse to sync with A. Need to re-record shorted case. 
%...66          OK  A = S is completed. Now shift A to really capture A.
%               OK  Notation is changed to S and ShortedS, now do A
%...65          ..  Calculate A now. 
%...64          OK  Now direct (S) channel is calculated and normalized!.
%                   Use the same number to normalize the side (A) channel correlation.
%               OK  Finally kinda works. Uses myxcorr_03test.m for
%                   amplitude normalization at the correlation step.
%               ..  HRIR and shortIR are shifted to t=0. You need to
%                   multiply HRIR by a phase factor due to delay.
%               ..  Write a code that reads every data file and processes
%                   with no further index number adjustments.
%...63  180308  ..  I hope to solve HRIR today. I might need to rewrite the
%                   correction function to calculate fft()/fft() instead of fft()*conj(fft())
%...62          SoSo time delay OK, cable's TF is OK, but HRIR does not look good!
%               OK  Captured. Now chect FFT's.
%...61          ..  Capture HRIR using the rough time information.
%               ..  Isolate both in windows, then do whatever you want to do. 
%               xx  Not a good idea since there are a train of pulses to work with!!!
%...60          ..  Can you try to cross-correlate shorted and Kirk to find
%                   the actual delay more accurately??
%               ..  Now time window HRIR finally!!
%               SS  Peak autocorr @ t= 82.6999582 sec, Kirk:82.6936,
%                   Short:82.6925. This makes a lead of 7.5 msec.
%                   Kirk head actual time delay is then; 1.1 msec.
%...59          OK  Now cross correlations are ready; (1) played audio,
%                   (2) recorded shorted audio, (3)recorded by Kirk head.
%                   Matched filter uses single time reversed AmpESS pulse.
%               OK  This shows that shorted signal is earlier as expected.
%...58          ..  Do the syncing. Check the corelator outputs.
%...57  180307  OK  Shorted cable gives very high noise level. Not
%                   balanced, and no power is given on Scarlett.
%               OK  Works but there is no t=0 sync.
%...56          ..  Now separate/window filter HRIR from the 
%               OK  AmpESS and REC data gives very good correlation.
%...55          ..  New AmpESS waveform, recorded Kirk's LEFT ear.
%                   Calculate HRIR and HRTF
%...54          SS  Now you have the PHASE plots too, and they are also correct.
%                   Phase distribution is NOT calculated.
%...53          SS  This is super! Delay measured. New convolution method
%                   tested to work. The shorted cable TF is almost allpass as expected.
%               OK  I used the function myMatchedFilt_01.m here to test
%                   what is the impulse response of st and ts. Frequency domain
%                   is almost ONE and the time domain is almost asymmetric IMPULSE.
%...52              Correct the function myxcorr_01.m
%...51          SS  Did the convolution in frequency NOT BY FFT(x).* conj(FFT(y))
%                   but do this convolution by; FFT(x)./ FFT(y) !!!
%               SS  SUUUUUPEEER. Finally I solved the equalization problem.
%               OK  This is SINGLE CHANNEL (analysis-wise)short circuit 
%                   transfer function analysis code. 
%...50          ..  Examine: Record_180302_10_Shorted_T10sec.wav
%               ..  I lost this version. 49 lost its input. New recording is coming.
%               ..  Shorted cable tested: time delay 4 msec = 137.6 cm But
%                   there is a BIG ERROR. Recording is earlier then PLAY !!??
%...49          ..  Nothing done. Copy of 48SS. Please come back to this
%                   version after checking the short cct delay and TF.
%                   now time to capture S and A for HRIR /then HRTF.
%...48          SS  Measured time delay; direct:9msec/310cm, cross:10msec/344cm
%...47          OK  I calculated the time delay of 0.0005sec 0.5msec, 17cm.
%                   Now cut and paste for S and A, then integrate to check SNR
%...46          OK  This is a similar version with spectrogram of recording.
%...45          OK  Pulse compression of the recorded signal
%                   Data: record_180228_8_T5sec.wav and  record_180228_8_T10sec.wave
%calculate_HRIR_45  A new fresh start with the time synchronized data.

%               OK  This is calculating S and A somehow normalized.
%...44          ..  New data; Record_180227_5_T5.wav  T=5
%...43          ..  Going back to T=10 sec. This is not good. 
%                   resolution is now bad; 6.88 cm. for T=1 sec.
%...41          ..  Now calculate again for record T=1 sec.               
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
    %Add paths
    myaddpath_02;
    mytic;
    load NewColorMapPDF;colormap(NewColorMapPDF);
    %Cleaning up
    figure(1), clf; 
    figure(2), clf; 
    figure(3), clf; 
    figure(4), clf; 
    figure(5), clf; 
    drawnow; pause(1)
    
    
%% CREATE SWEEP FUNCTION
    Fs = 48e3; Ts = 1/Fs;
    T  = 10; %seconds
    f1 = 1e2;        
    f2 = 20e3;
    
    %% Load input files for TX and Rec data.    
    %TXfilename = 'LEFTChAllPulses-T10-1e2TO20e3-180308.wav';
    TXfilename = 'AmpESS-T10-1e2TO20e3-180306.wav';
    [dum,Fs]=audioread(TXfilename);
    ts(:) = dum(:,1);
    timets = ([1:length(ts)]./Fs)';

    %% Analysis of the TX audio data is not necessary any more.
    %TXfilename = 'TXData-T10-1e2To20e3-180306.wav';
    %[Tx,Fs]=audioread(TXfilename);
    %Tx_left(:)  = Tx(:,1);
    %timeTx = [1:length(Tx_left)]./Fs;
    
     Recfilename = 'RecData-T10-1e2To20e3-180308-Shorted.wav';
    %Recfilename = 'LEFTChAllPulses-T10-1e2TO20e3-180308.wav';
    %Recfilename = 'RecData-T10-1e2To20e3-180307-Shorted.wav';
    [Shorted,Fs]=audioread(Recfilename);
    shorted_left(:)  = Shorted(:,1);
    %shorted_right(:) = Shorted(:,2);  %not used
    timeshr = [1:length(shorted_left)]./Fs;

    Recfilename = 'RecData-T10-1e2To20e3-180306.wav';
    [RoomIR,Fs]=audioread(Recfilename);
    sig_left(:)  = RoomIR (:,1);
    sig_right(:) = RoomIR (:,2);
    timerec = [1:length(sig_left)]./Fs;

%1  FIG TIME: TS  SHORTED_LEFT   SIG_LEFT
    figure(1)
        subplot(312),  plot(timeshr, shorted_left) 
    title('Recorded Left Channel (Shorted)');
    xlabel('Time (seconds)');
    axis('tight'); V= axis; 
        subplot(313),  plot(timerec, sig_left)   %this is recording of A then S
        %subplot(313),  plot(sig_right)  %this is garbage channel
    title('Recorded Left Channel (Kirk head)');
    xlabel('Time (seconds)');
    axis('tight');
        subplot(311),  plot(timets, ts)         %this is the excitation AmpESS signal
    title('AmpESS Matched Filter Signal (100 Hz - 20 KHz)');
    xlabel('Time (seconds)');
    axis([V(1) V(2) -1.1 1.1]);


    
%% CALCULATE AUTO- AND CROSS-CORRELATION
    fdom = 1;
    dB = 0;
    dBmax = 0;
    dBmin = -30;
    
    [shortedcor,shortedfcor,shortedcortime,shortedcorfreq,N] = ...
        myxcorr_03test(shorted_left,ts,Fs,fdom,dB,dBmin,dBmax);
    [croscor,fcor,cortime,corfreq,N] = ...
        myxcorr_03(sig_left,ts,Fs,fdom,dB,dBmin,dBmax);
    
%1  FIND THE DELAY TIMES TO CAPTURE S AND A    
    %manual search for peaks
    figure(1), clf, hold off;
        subplot(211), 
    plot(shortedcortime, shortedcor)
    title('Received signal (compressed) - Cable ');
    xlabel('Time (seconds)');
    legend('First 5 is A, and next 5 is S');
        subplot(212), 
    plot(cortime, croscor)
    title('Received signal (compressed) - Loudspeaker');
    xlabel('Time (seconds)');
    legend('First 5 is A, and next 5 is S');
    
    
    %figure(2), clf, hold off;
    %    subplot(211), 
    %plot(shortedcortime, shortedcor)
    %    subplot(212), 
    %plot(cortime, croscor)
    
    
    
    
    
    %% CAPTURE DIRECT PATH (S) capture data using time windowing
    %Shorted pulse times
    % first  5: for A: 10.2954, 20.59495, 30.8946, 41.1942, 51.49385
    % second 5: for S: 62.09348,72.3931,  82.6927, 92.9921, 103.293

    % HRIR: this is (S)
    Tdelay = -22.07e-3; %this is for nfft=2048
    %Tdelay = -10.07e-3; %this is for nfft=512
    %Tdelay = -10.05e-3;
    %Tdelay = 0.00121 - 10.05e-3;
    %Tdelay = 0.00121;

    T0short = 62.09348 - Tdelay;         %delT = 0.003;
    T0short = 72.3931 - Tdelay;          %delT = 0.003;
    T0short = 82.6937 - Tdelay;          %delT = 0.003;
    T0short = 92.9921 - Tdelay;          %delT = 0.003;
    T0short = 103.293 - Tdelay;          %delT = 0.003;
    Nshort  = floor( T0short * Fs);
    Nfft    = 2048; %255; %delT * Fs
    Ndel    = Nfft; 
    tcor = [-Ndel:Ndel]./Fs;
    %
    ShortedS([1: 2*Ndel+1]) = shortedcor([Nshort-Ndel:Nshort+Ndel]);
    S([1: 2*Ndel+1])        = croscor   ([Nshort-Ndel:Nshort+Ndel]);
    
    
%% NOW CAPTURE SIDE PATH (A) capture data using time windowing
    % A: 
    %Tdelay = 0.00121;
    %Shorted A
    %10.2954, 20.59495, 30.8946, 41.1942, 51.49385
    T0short = 10.29524 - Tdelay;  
    T0short = 20.59495 - Tdelay;  
    T0short = 30.8946 - Tdelay;  
    T0short = 41.1942 - Tdelay;  
    T0short = 51.49385 - Tdelay;  
    NshortA  = floor( T0short * Fs);
    
    ShortedA([1: 2*Ndel+1]) = shortedcor([Nshort-Ndel:Nshort+Ndel]);
    A([1: 2*Ndel+1])        = croscor   ([Nshort-Ndel:Nshort+Ndel]);    
    
    
 
    
    %% CALCULATE IN FREQUENCY DOMAIN: HRTR
    fS      = fft(S);  %this is the HRTF
    fA      = fft(A);  %this is the HRTF
    fshortS  = fft(ShortedS);
    fshortmax = 1; %0.7476; %max(abs(fshortS))
    fshortS  = fshortS ./ fshortmax;
        
    fshortA  = fft(ShortedA);
    fshortA  = fshortA ./ fshortmax;
    
    fSmax   = 2.0618e+04; %max(abs(fS))   %Now use this factor also to normalize A!
    fS      = fS ./ fSmax;
    fA      = fA ./ fSmax;
    
    fSnorm = fS./fshortS;
    fAnorm = fA./fshortA;
    
    
    
    
    %% LINEAR PHASE TERM: Calculate the to-be-subtracted linear phase term
    %linear_phase_deg is -1.424586663e5 at 20 KHz.
    %linear_phase_rad is -2.4864e+03    at 20 KHz
    freq = [0:length(fshortS)-1]./(length(fshortS)-1).*Fs;
    
    lineartimefactor = freq./20e3.*2.4864e3;
    linearphaseterm = exp(i.*lineartimefactor);
    %linearphaseterm = exp(i.*(2.*pi.*freq.*lineartimefactor));
    %linearphaseterm = 1; 
    
    fS      = fS      .* linearphaseterm;
    fA      = fA      .* linearphaseterm;
    fshortS = fshortS .* linearphaseterm;
    fshortA = fshortA .* linearphaseterm;
    fSnorm  = fSnorm  .* linearphaseterm;
    fAnorm  = fAnorm  .* linearphaseterm;
    
    S = ifft(fS);
    A = ifft(fA);
    
    ShortedS = ifft(fshortS);
    ShortedA = ifft(fshortA);
    
    
    % calculate phase responses
    phase_S         = unwrap(atan2(imag(fS),      real(fS)));
    phase_A         = unwrap(atan2(imag(fA),      real(fA)));
    phase_shortS    = unwrap(atan2(imag(fshortS), real(fshortS)));
    phase_shortA    = unwrap(atan2(imag(fshortA), real(fshortA)));
    phase_Snorm     = unwrap(atan2(imag(fSnorm),  real(fSnorm)));
    phase_Anorm     = unwrap(atan2(imag(fAnorm),  real(fAnorm)));
    %
    phase_S         = phase_S /pi * 180;
    phase_A         = phase_A /pi * 180;
    phase_shortS    = phase_shortS /pi * 180;
    phase_shortA    = phase_shortA /pi * 180;
    phase_Snorm     = phase_Snorm /pi * 180;
    phase_Anorm     = phase_Anorm /pi * 180;
    
    
    
    % calculate time-delay responses
    delay_S         = 1000 .* phase_S ./ 360 ./ freq;
    delay_A         = 1000 .* phase_A ./ 360 ./ freq;
    delay_shortS    = 1000 .* phase_shortS ./ 360 ./ freq;
    delay_shortA    = 1000 .* phase_shortA ./ 360 ./ freq;
    delay_Snorm     = 1000 .* phase_Snorm ./ 360 ./ freq;
    delay_Anorm     = 1000 .* phase_Anorm ./ 360 ./ freq;
    
    
    Fmn = 20;
    Fmx = 20e3;
    
    
    
    
    
       
%2  PLOT CORRELATION    
    figure(2), clf
        subplot(211), 
    plot(tcor, real(ShortedS), 'r-'), hold on;
    plot(tcor, imag(ShortedS), 'm-'), hold on;
    plot(tcor, real(ShortedA)-0.1, 'b-')
    plot(tcor, imag(ShortedA)-0.1, 'm-')
    axis('tight'); grid on;
    title('HRTFs (red)S (blue)A-0.1 (Shorted cable)')
    xlabel('Time (seconds)')
        subplot(212), 
    plot(tcor, real(S), 'r-'); hold on;
    plot(tcor, imag(S), 'm-'); hold on;
    plot(tcor, real(A)-0.1, 'b-');
    plot(tcor, imag(A)-0.1, 'm-');
    axis('tight'); grid on;
    title('HRTFs (red)S (blue)A-0.1 (Kirkhead-Loudspeaker-PC)')
    xlabel('Time (seconds)')
    
    
    
    
    figure(3), clf
        subplot(311),
    semilogx(freq, 20 .* log10(abs(fshortS)), 'r-'); hold on;
    semilogx(freq, 20 .* log10(abs(fshortA))-10, 'b-');
    V=axis; axis([Fmn Fmx -80 20]); 
    set(gca,'xtick',[Fmn  100 1000 10000 Fmx]); grid on;
    title ('HRTFs in dB (red)S (blue)(A-10dB) (Cable transmission)');
    xlabel('Frequency (Hertz)');
    ylabel('Desibels');
        subplot(312),
    semilogx(freq, 20 .* log10(abs(fS)), 'r-'); hold on;
    semilogx(freq, 20 .* log10(abs(fA))-10, 'b-');
    V=axis; axis([Fmn Fmx  -80 20]); grid on;
    set(gca,'xtick',[Fmn  100 1000 10000 Fmx]); grid on;
    title ('HRTFs in dB (red)S (blue)(A-10dB) (Short distance loudspeaker)');
    xlabel('Frequency (Hertz)')
    ylabel('Desibels');
        subplot(313),
    semilogx(freq, 20 .* log10(abs(fSnorm)), 'r-'); hold on;
    semilogx(freq, 20 .* log10(abs(fAnorm))-10, 'b-');
    V=axis; axis([Fmn Fmx  -80 20]); grid on;
    set(gca,'xtick',[Fmn  100 1000 10000 Fmx]); grid on;
    title ('Divide TF of loudspeaker / TF of cable');
    xlabel('Frequency (Hertz)')
    ylabel('Desibels');
    
    figure(4),
        subplot(311), 
    plot(freq, phase_shortS,'r-'); hold on;
    plot(freq, phase_shortA-5,'b-'); grid on;
    %axis('tight'); 
    V=axis; axis([Fmn Fmx  -15e3 0.1e4]); grid on;
    %V=axis; axis([Fmn Fmx  V(3)/2 V(4)]); grid on;
    set(gca,'xtick',[Fmn  5000 10000 15000 Fmx]); grid on;
    %set(gca,'xtick',[Fmn  100 1000 10000 Fmx]); grid on;
    title ('Phase of (red) S and (blue:A-5) (Cable transmission)');
    xlabel('Frequency (Hertz)');
    ylabel('Phase (degrees)');
    
        subplot(312), 
    plot(freq, phase_S,'r-'); hold on; 
    plot(freq, phase_A-5,'b-'); grid on;
    axis('tight'); 
    V=axis; axis([Fmn Fmx  -15e3 0.1e4]); grid on;
    %V=axis; axis([Fmn Fmx  V(3)/1.5 V(4)]); grid on;
    set(gca,'xtick',[Fmn  5000 10000 15000 Fmx]); grid on;
    %set(gca,'xtick',[Fmn  100 1000 10000 Fmx]); grid on;
    title('Phase of (red) S and (blue: A-5)(Short distance loudspeaker)');        
    xlabel('Frequency (Hertz)');
    ylabel('Phase (degrees)');
    
        subplot(313), 
    plot(freq, phase_Snorm,'r-'); hold on;
    plot(freq, phase_Anorm-500,'b-'); grid on;
    axis('tight'); 
    V=axis; axis([Fmn Fmx  V(3) V(4)/2]); grid on;
    set(gca,'xtick',[Fmn  5000 10000 15000 Fmx]); grid on;
    %set(gca,'xtick',[Fmn 100 1000 10000 Fmx]); grid on;
    title('Phase difference of two S (red) A-500 (blue) measurements');        
    xlabel('Frequency (Hertz)');
    ylabel('Phase (degrees)');
    
    
 
       figure(5),clf
        subplot(311), 
    plot(freq, delay_shortS,'r-'); hold on;
    plot(freq, delay_shortA-1,'b-'); grid on;
    axis('tight'); 
    %V=axis; axis([Fmn Fmx  -25 -15]); grid on;
    %V=axis; axis([Fmn Fmx  V(3) V(4)]); grid on;
    set(gca,'xtick',[Fmn  5000 10000 15000 Fmx]); grid on;
    %set(gca,'xtick',[Fmn  100 1000 10000 Fmx]); grid on;
    V=axis; axis([0 Fs/2 -5 2]); grid on;
    title('Time delay response (red) S and (blue:A-1) (Cable transmission)');
    xlabel('Frequency (Hertz)');
    ylabel('Time delay (mseconds)');
    
        subplot(312), 
    plot(freq, delay_S,'r-'); hold on; 
    plot(freq, delay_A-1,'b-'); grid on;
    axis('tight'); 
    V=axis; axis([Fmn Fmx  -5 2]); grid on;
    %V=axis; axis([Fmn Fmx  V(3) V(4)]); grid on;
    set(gca,'xtick',[Fmn  5000 10000 15000 Fmx]); grid on;
    %set(gca,'xtick',[Fmn  100 1000 10000 Fmx]); grid on;
    title('Time delay response of (red) S and (blue: A-1)(Short distance loudspeaker)');        
    xlabel('Frequency (Hertz)');
    ylabel('Time delay (mseconds)');
    
        subplot(313), 
    plot(freq, delay_Snorm,'r-'); hold on;
    plot(freq, delay_Anorm-1,'b-'); grid on;
    axis('tight'); 
    V=axis; axis([Fmn Fmx  10 25]); grid on;
    %V=axis; axis([Fmn Fmx  V(3) V(4)]); grid on;
    set(gca,'xtick',[Fmn  5000 10000 15000 Fmx]); grid on;
    %set(gca,'xtick',[Fmn  100 1000 10000 Fmx]); grid on;
    title('Time delay difference of two S (red) and A (blue)-1 measurements');        
    xlabel('Frequency (Hertz)');
    ylabel('Time delay (mseconds)');
    
 
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
mytoc; 



%% THE END  THE END  THE END  THE END  THE END  THE END  THE END  THE END 