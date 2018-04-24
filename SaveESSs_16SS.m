%% SaveESSs_xx.m
% x copy of Amplitude Modulated Sweep design
% Save wav files for ESS which have amplitude modulation for 
% unity spectrum distribution (check the phase)
%
% S. Gokhun Tanyer
% 180306

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
%...16          OK  All are recorded once with new naming convention. T=20
%               ..  I still miss things. All should be once. AmpESS, TXdata
%                   and the TXdata for the shorted test.
%...15          SS  This is for T=20.
%               SS  This is for T=10.
%...14  180327  ..  Will now save AmpESS too as a standard.
%...13SS
%               SS  Data ready.
%                   2.5 seconds of initial silence. 0.5 seconds between.
%...13  180326  OK  Total of 422 seconds, 7.03 minutes.
%               ..  There is not enough silence at the beginning.
%       180326  OK  I am changing the save data filename.
%...12  180323  ..  Micro mics are almost ready. Now generate data for
%                   them. Alternate one by one for 20 pulses each.
%                   SUPER!! keep this version.
%...11_FS       SS  This is the SS version for free-space transmission.
%...11          OK  (2+2) pulses. Left and Right channels corrected.
%...10          OK  Better version. 5+5 pulses. Now do 2+2 pulses.
%                   T=20 and frequency band {20,20K} are updated
%...09  180314  ..  This is version 07SS for {20,20K} (2+2)

%REVIEW OF VERSIONS: 07SS and 08SS's are in 100, 20K Hertz band.
%                   They use 5+5 pulses for left and right.
%                   07SS is for HRIR, 08SS is for shorted audio tests.
%                   now adjust all for the new tests (20, 20K) 2+2 pulse

%...08  180308  ..  Create signal file for the short circuited recordings.
%                   where I need all pulses on one channel, do both.
%               SS  This is the version to generate matched filter and TX
%                   signals.
%                   Do a full analysis and save data.
%...07  180306  ..  Goooo and save T = 1, 5, 10, 15 and 20 second wav data.
%SaveESSs_07.m

%...07  180306  ..  Record ESS files of various lengths to test.
%...06          ..  Go from here. 
%...05          OK  This version uses mysweep_02.m and shows s, sinv and
%                   the cross-correlation function in time and frequency. 
%                   They look beautiful!
%               ..  Now copy this to mysweep_02.m
%...04          SS  Superrr!   You found a much better ESS with almost flat
%                   response.
%                   Now I will try directly multiplying m(t) of (7)
%...03          OK  This is Carson's second inverse filter and still not
%                   good and also it can be calculated due to large matrix sizes.
%...02          OK  This shows the s * sinv of Carson is not allpass.
%                   That also means sinv is not a good inverse filter.
%               OK  Carson's sweep/time reversed sweep and EQUAL with mines.
%...01  180305  ..  OK goooo! Use Peter's students; M.Carson, H.Giesbrecht and T.Perry.
%AmpModSweep_01.m
%new branched version.

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
    myaddpath_02;
    %addpath_MyFcns
    mytic;
    
    load NewColorMapPDF;colormap(NewColorMapPDF);
    figure(1), clf; 
    figure(2), clf; 
    figure(3), clf; 
    drawnow;
    
    
%% CREATE SWEEP FUNCTION
    Fs = 48e3; Ts = 1/Fs;
    T  = 20; %seconds
    f1 = 100; %1e2;        
    f2 = 20e3;
    
    fprintf('CREATE: sweep\n');  
    % st is the up chirp, ts is the down chirp.
    [st,ts,time,N] = mysweep_02(Fs, f1, f2, T); %this is amplitude equalized pulse
    %fprintf('CREATE: sweep DONE\n');
    %normalize st and ts to 1.00
    dum = max(max([st,-st,ts,-ts]));
    st = st ./ dum;  %not to overload.
    ts = ts ./ dum;
    
    
%% CALCULATE FFT OF SWEEP
    norm = 1;
    dB = 1;
    
    [fft_st,  freq, N] = myfft_01(Fs,length(st  ),norm,dB,st  );
    [fft_ts,  freq, N] = myfft_01(Fs,length(ts  ),norm,dB,ts  );
    
    phase_st   = unwrap(atan2(imag(fft_st  ), real(fft_st  )));
    phase_ts   = unwrap(atan2(imag(fft_ts  ), real(fft_ts  )));

    
%% CALCULATE AUTO- AND CROSS-CORRELATION
    fdom = 1;
    dB = 1;
    dBmax =  0;
    dBmin = -90;
    
    %function [cor,fcor,time,freq,N] = myxcorr_02(data1,data2,Fs,fdom, dB,dBmin,dBmax) 
    [autocor,fautocor,cortime,corfreq,N] = myxcorr_03(st,st,Fs,fdom,dB,dBmin,dBmax);
    [croscor,fcroscor,cortime,corfreq,N] = myxcorr_03(st,ts,Fs,fdom,dB,dBmin,dBmax);
    
    
    
%% PLOTTING SECTION
%1    
    figure(1), clf, hold off;
    subplot(211),
    plot(time, st, 'r-'), hold on, grid on
    title('My Exponential Sine Sweep (ESS) st');
    xlabel('Time (seconds)');

    subplot(212),
    plot(time, ts, 'b-'), hold on, grid on
    title('My Exponential Inverse Sine Sweep (inv-ESS) ts');
    xlabel('Time (seconds)');
%2    
    figure(2), clf, hold off;
    subplot(311),
     plot((cortime-T).*1000, 20 .* log10( abs(croscor))...
         -max(20 .* log10( abs(croscor))),'r-'), hold on;
    V=axis;grid on;
    axis([V(1) V(2) -150 0]);
    title('Cross-Correlation Function (st and ts)');
    xlabel('Time (seconds)');
    
    subplot(312),
    %plot(cortime, croscor,'r-'), hold on;
     plot((cortime-T).*1000, 20 .* log10( abs(croscor))...
         -max(20 .* log10( abs(croscor))),'r-'), hold on;
    V=axis;grid on;
    axis([-4 4 -60 0]);
    title('Cross-Correlation Function (st and ts)');
    xlabel('Time (microseconds)');
    
    subplot(313),
    semilogx(corfreq, fautocor, 'r');hold on;
    semilogx(corfreq, fcroscor, 'b');hold on;
    V=axis; 
    axis([20 Fs/2 -50 5]);grid on;
    %axis([f1-10 Fs/2 -50 5]);grid on;
    title('Spectrum of the Cross correlation Function');
    xlabel('Frequency (Hertz)');
    pause(1)
        
   
%%  CREATE THE TRANSMITTER SIGNAL   / also the SHORTED SIGNAL
    %TX:    stereo signal with single active channel!
    %SHORT: stereo signal with both same active channels!
    
    Nsilence = Fs/2; %14400;
    Nlong    = 2*Fs;
    ss    = zeros(1,Nsilence);
    s2    = zeros(2,Nsilence);
    slong = zeros(2,Nlong);
    
    pL (1,:) = [st ss];
    pL (2,:) = zeros(1, length([st ss]));
      
    pR (2,:) = [st ss];
    pR (1,:) = zeros(1, length([st ss]));
    
    pSH(1,:) = [st ss];
    pSH(2,:) = [st ss];
    
    sLEF  = [s2 pL ];
    sRIG  = [s2 pR ];
    sSHR  = [s2 pSH];
    %sLEF = [s2 pL s2 pL];
    %sRIG = [s2 pR s2 pR];
    
    %now have alternating 10 pulses per channel.
    data = [slong sLEF, sRIG,  sLEF, sRIG,  sLEF, sRIG,  sLEF, sRIG,  sLEF, sRIG, ...
                  sLEF, sRIG,  sLEF, sRIG,  sLEF, sRIG,  sLEF, sRIG,  sLEF, sRIG, ];
    %data = [slong, sLEF, slong, sRIG, s2];
  dataSH = [slong sSHR, sSHR,  sSHR, sSHR,  sSHR, sSHR,  sSHR, sSHR,  sSHR, sSHR, ...
                  sSHR, sSHR,  sSHR, sSHR,  sSHR, sSHR,  sSHR, sSHR,  sSHR, sSHR, ];
    
    lendata = length(data);
    
    timed = [0:lendata-1]./Fs;
    
    figure(3), clf
    subplot(211), 
    plot(timed, dataSH(1,:), 'r-'), hold on, axis('tight');
    title('Left Channel')
    
    subplot(212), 
    plot(timed, dataSH(2,:), 'b-'), hold on
    %plot(timed, data(2,:), 'k+'), 
    axis('tight');
    %axis([103.6 103.6+0.005 -1 1]);
    title('Right Channel')
    
    
    
    %% WRITE DATA OUT
    
    %p = audioplayer(data', Fs);
    %play(p, [1 Fs*20*3]);
 
    %% write the original time warped pulse signal (st)
    filename = '180327-AmpESSts-T20-1e2To20e3.wav';
    audiowrite(filename,ts',Fs);
   
    %% write the transmitted stereo test signal FOR THE SHORTED CABLE (dataSH)
    filename = '180327-TXDataSH-T20-1e2To20e3.wav';
    audiowrite(filename,dataSH',Fs);
   
    %% write the transmitted stereo test signal FOR THE FREE SPACE TEST (data)
    filename = '180327-TXDataFS-T20-1e2To20e3.wav';
    audiowrite(filename,data',Fs);
   
mytoc; 




















%% THE END  THE END  THE END  THE END  THE END  THE END  THE END  THE END 