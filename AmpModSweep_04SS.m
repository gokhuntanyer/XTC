%% Amplitude Modulated Sweep design
%
% S. Gokhun Tanyer
% 180302

%% HISTORY
%
%
%
%
%
%
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
    addpath_MyFcns
    mytic;
    
    load NewColorMapPDF;colormap(NewColorMapPDF);
    figure(1), clf; 
    figure(2), clf; 
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
    
    
    
    %% Carson's sweep is here
    % Create the swept sine tone
    w1 = 2*pi*f1;
    w2 = 2*pi*f2;
    K = T*w1/log(w2/w1); %parameter for creating the sweep
    L = T/log(w2/w1); %%parameter for creating the sweep
    t = linspace(0,T-1/Fs,Fs*T); %time vector
    s = sin(K*(exp(t/L) - 1)); %sweep vector
    LT = length(s);
    
    % Let us multiply my sweep with Carson's m(t) given in (7)
    arg = -time./L;
    mt = 1 ./ (K./L .* exp(arg./2));
    s = mt .* s;   % this is a SUPER flat-spectrum waveform.
    
    s = s -sum(s)./length(s);
    dum = max(max(s));
    s = 0.8 .* s ./ dum;
    invs([1:LT]) = s([LT:-1:1]);
    
    
    %% Carson's practical but not very good inverse sweep
    % there is no explanation why we have the power .^3.
    %x = linspace(f1,f2,length(s));  %This inverse filter is not perfect
    %slope = fliplr(x).^3;           %but is sufficient to prove the method works
    %sinv = slope.*fliplr(s);        %Create the inverse filter
    %sinv = sinv./max(sinv);         %normalize the inverse filter
    %fsinv = fft(sinv);
    
    % Carson's 2: Different approach to inverse filter (based on eq (7) report)
    %f_range=linspace(f1,f2,length(s)); %This inverse filter is not perfect
    %A = f1;
    %m = f1./f_range; % modulating function (8)
    %s_rev = fliplr(s); %time reveresal of s
    %sinv2 = m'*s_rev;
    %sinv2 = sinv2./max(sinv); %normalize the inverse filter
    %SINV2 = fft(sinv2);
    
    
    
%% CALCULATE FFT OF SWEEP
    norm = 1;
    dB = 1;
    
    [fft_st,  freq, N] = myfft_01(Fs,length(st  ),norm,dB,st  );
    [fft_ts,  freq, N] = myfft_01(Fs,length(ts  ),norm,dB,ts  );
    [fft_s,   freq, N] = myfft_01(Fs,length(s   ),norm,dB,s   );
    %[fft_sinv,freq, N] = myfft_01(Fs,length(sinv),norm,dB,sinv);
    %fft_sCsinv = fft_s + fft_sinv;      %check how s and sinv cancels
    %fft_sCsinv = fft_sCsinv - max(fft_sCsinv);
    
    phase_st   = unwrap(atan2(imag(fft_st  ), real(fft_st  )));
    phase_ts   = unwrap(atan2(imag(fft_ts  ), real(fft_ts  )));
    phase_s    = unwrap(atan2(imag(fft_s   ), real(fft_s   )));
    %phase_sinv = unwrap(atan2(imag(fft_sinv), real(fft_sinv)));
    
    
%% CALCULATE SPECTROGRAM
    %nfft = 1024*8;
    %overlapr = 0.5;
    %noverlap = nfft*overlapr;
    %window = nfft;
    %dB = 1;
    %dBmax = 0;
    %dBmin = -30;
    %norm = 1;
    %[Ssig_st,xt,yf] = myspectrogram_01(st,Fs,nfft,norm,dB,dBmin,dBmax);
    %[Ssig_ts,xt,yf] = myspectrogram_01(ts,Fs,nfft,norm,dB,dBmin,dBmax);
    
%% CALCULATE AUTO- AND CROSS-CORRELATION
    fdom = 1;
    dB = 0;
    dBmax =  0;
    dBmin = -30;
    
    [autocor ,fcor ,cortime,corfreq,N] = myxcorr_01(st,ts  ,Fs,fdom,dB,dBmin,dBmax);
    [croscor ,fcor ,cortime,corfreq,N] = myxcorr_01(st,st  ,Fs,fdom,dB,dBmin,dBmax);
    [croscor2,fcor2,cortime,corfreq,N] = myxcorr_01(s ,invs,Fs,fdom,dB,dBmin,dBmax);
    croscor2 = croscor2 ./ max(croscor2);
    fcor2= fcor2 ./ max(fcor2) ;
    
    
    
    
%% PLOTTING SECTION
%1    
    figure(1), clf, hold off;
    subplot(211),
    %plot(time,st,       'r-'), hold on
    %plot(time,ts, 'b-'), hold on
    plot(time, s, 'r-'), hold on
    %plot(time,invs, 'c-'), hold on
    %V=axis;
    %axis([9 11 V(3) V(4)]); %axis([T-70/f2 T V(3) V(4)]);
    title('My Exponential Sine Sweep (ESS)');
    xlabel('Time (seconds)');
    
    
    
    
    %figure(2), clf, hold off;
     subplot(212),
    %semilogx(freq, fft_st  , 'r-');hold on
    %semilogx(freq, fft_ts  , 'b-');hold on
    semilogx(freq, fft_s + 0.71   , 'r-');hold on
    %semilogx(freq, fft_sinv, 'b-');hold on
    %semilogx(freq, fft_sCsinv, 'k-');hold on
    legend('Spectrum of my ESS', 'Carson-invESS', 's times sinv');
    V=axis; axis([50 Fs/2 -15 2]), grid on;
    xlabel('Frequency (Hertz)');
    title('Spectrum of My Exponential Sine Sweep (ESS)');
    
    figure(2), clf, hold off;
    
    subplot(211),
    plot(cortime, 20 .* log10( abs(croscor2)),'r-'), hold on;
    V=axis;
    axis([9.99 10.01 -80 0]);
    %axis([9.99 10.01 V(3) V(4)]); %axis([T-70/f2 T V(3) V(4)]);
    title('Cross-Correlation Function (s and sinv)');
    xlabel('Time (seconds)');
    
    subplot(212),
    %semilogx(corfreq, abs(fcor2));
     semilogx(corfreq, 20 .* log10( abs(fcor2)) + 2.67);
    %loglog(corfreq, abs(fcor));
    %semilogx(freq, 20.*log10(abs(fftst)));
     V=axis; axis([f1-10 Fs/2 -50 5]);
    title('Spectrum of the Correlation Function (s and sinv)');
    xlabel('Frequency (Hertz)');
    
    
    
    dur
    pause(1)
%2    
    figure(1),clf, hold off;
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
    drawnow; pause(2);
    %dur
%1
    figure(1), clf, hold off;
    subplot(211),
    %plot(cortime, croscor), hold on;
    plot(cortime, autocor);grid on;axis('tight');
    delT = 0.05;
    V=axis; axis([10-delT 10+delT V(3) V(4)]);
    title('Auto-Correlation')
    xlabel('Time (seconds)');

    dum = 20 .* log10( abs(autocor) );
    dum = dum - max(dum);
    
    subplot(212),
    plot(cortime, dum), hold on;
    V=axis; axis([10-delT 10+delT -50 V(4)]); grid on
    title('Auto-Correlation')
    xlabel('Time (seconds)');
    drawnow; pause(1);
    
    
    
    
    
    dur   
    
%% ANALYZE DATA FILE:  
%   1. load recording; RoomIR (recorded) and the Sweep (played)
%   2. take FFT of both sides.
%   3. multiply them with FFT of the original pulse.
%   4. take IFFT and compare RoomIR and Sweep

%   Recorded signal in the room
    filename = 'Record_180302_10_Shorted_T10sec.wav';
    %filename = 'Record_180228_9_Shorted_T10sec.wav';
    %filename = 'record_180228_8_T10sec.wav';
    [RoomIR, Fs]=audioread(filename);
    %this is just the recorded signal from the room
    RoomIR_left  = RoomIR (:,2)';    
    RoomIR_right = RoomIR (:,1)'; %changing left and right places
    %RoomIR_left  = RoomIR (:,1)';    
    %RoomIR_right = RoomIR (:,2)'; %RIGHT CHANNEL IS GARBAGE!
    %
    fRoomIR_left  = fft(RoomIR_left);
    %fRoomIR_right = fft(RoomIR_right); GARBAGE

%   Played signal in the room    
    filename2 = 'Data_180228_8_T10sec.wav';
    [Sweep, Fs]=audioread(filename2);
    Sweep_left  = Sweep (:,1)';
    Sweep_right = Sweep (:,2)';
    %
    fSweep_left  = fft(Sweep_left);
    fSweep_right = fft(Sweep_right);

    % pulse compression of the received signal
    lenIR = length(RoomIR_left);
    lenst = length(st);
    lendif = lenIR - lenst;
    dum = zeros(1,lendif);
    
    st_long = [st, dum];
    fftst = fft(st_long);
        
    fftIRleft = fRoomIR_left  .* conj(fftst);
    IRleft = ifft(fftIRleft);
    
    
    % pulse compression of the played recording
    lenIR = length(Sweep_left);
    lenst = length(st);
    lendif = lenIR - lenst;
    dum = zeros(1,lendif);
    
    st_long = [st, dum];
    fftst = fft(st_long);
        
    fftSweepleft  = fSweep_left  .* conj(fftst);
    fftSweepright = fSweep_right .* conj(fftst);
    Sweepleft  = ifft(fftSweepleft);
    Sweepright = ifft(fftSweepright);
    
%1    
    figure(1),clf,
    subplot(311),     plot(RoomIR_left,  'r');
    subplot(312),     plot(RoomIR_right, 'b');

    subplot(313),     plot(Sweep_left,   'r'); hold on
                      plot(Sweep_right,  'b'); hold on
    drawnow; pause(1);
    %dur   
%1    
    figure(1),clf,
    len = length(IRleft);
    time1 = [0:len-1]./Fs;
    
    subplot(211),     plot(time1, IRleft, 'r');
    
    No = 5.3755e6; dN =  1400;
    %V=axis; axis([No./Fs (No+dN)./Fs V(3) V(4)]);
    %V=axis; axis([5.3755e6 5.3768e6 V(3) V(4)]);                  
    title('Direct Path (S) in Time Domain');
    %
    len = length(Sweepleft);
    time2 = [0:len-1]./Fs;
    
    subplot(212),     plot(time2, Sweepleft, 'r'); hold on;
                      plot(time2, Sweepright, 'b');
    %V=axis; axis([No./Fs (No+dN)./Fs V(3) V(4)]);                  
    %V=axis; axis([5.3755e6 5.3768e6 V(3) V(4)]);                  
    title('Direct Path: Synchronization signal in Time Domain');
    dur   
%2    
    figure(2),clf,
    
    subplot(211),     plot(time1, IRleft, 'b');
    No = 4.79945e6; dN =  1400;
    V=axis; axis([No/Fs (No+dN)/Fs V(3) V(4)]);                  
    title('Cross Path (A) in Time Domain');
    %
    subplot(212),     plot(time2, Sweepleft, 'r'); hold on;
                      plot(time2, Sweepright, 'b');
    V=axis; axis([No/Fs (No+dN)/Fs V(3) V(4)]);                  
    %V=axis; axis([4.79945e6 4.8009e6 V(3) V(4)]);                  
    title('Cross Path: Synchronization signal in Time Domain');
    
    drawnow; pause(1)
    
    
     
    
%2    
    figure(2),
    subplot(211), axis('tight');
    title('Direct and Cross Paths (S) and (A) in Time Domain');
    
    subplot(212), axis('tight');
    title('Direct and Cross Paths: Synchronization signal in Time Domain');
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
                      
                      dur
    
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
%% CALCULATE SPECTROGRAM OF S AND A
    nfft = 1024*8;
    overlapr = 0.5;
    noverlap = nfft*overlapr;
    window = nfft;
    dB = 1;
    dBmax = 0;
    dBmin = -40;
    norm = 1;
    
    [SpecIRleft,xt,yf] = myspectrogram_01(RoomIR_left,Fs,nfft,norm,dB,dBmin,dBmax);
    %[SpecIR_A,xt,yf] = myspectrogram_01(RoomIR_right,Fs,nfft,norm,dB,dBmin,dBmax);
%2
    figure(2), clf, hold off;
    %subplot(212),
    mesh(xt, yf, SpecIRleft)
    shading('interp');axis('tight');view( [0 90]); 
    V = axis; axis([V(1) V(2) 0 Fs/2]);
    view( [0 90])
    xlabel('Time (seconds)'); ylabel('Frequency (Hz)');
    title('Spectrogram of TX-2 signal'); colorbar;
    colormap(NewColorMapPDF);colorbar;
    drawnow; pause(1)
    
    
%1
    figure(1), clf, hold off;
    Mnorm = max(abs(IRleft));
    IRleft = IRleft ./ Mnorm;
    
    subplot(311),    
    plot(timeS,IRleft, 'r');grid on;axis('tight'); hold on
    %V=axis; axis([109.5 112.5 V(3) 1.1*V(4)]); %This is for T=10;
    title('Cross-correlation Short CCT Test')
    xlabel('Time (seconds)');
    
    subplot(312),    
    plot(timeS,IRleft, 'r');grid on;axis('tight'); hold on
    V=axis; axis([112.0085 112.02 V(3) 1.1*V(4)]); %This is for T=10;
    title('Cross-correlation Short CCT Test')
    xlabel('Time (seconds)');

    
    IRleftdB = 20 .* log10(abs(IRleft));
    IRleftdB = IRleftdB - max(IRleftdB);
    
    subplot(313),
    plot(timeS,IRleftdB, 'r');grid on;axis('tight'); hold on
    plot(timeS,IRleftdB, 'r+');grid on;axis('tight'); hold on
    V=axis; axis([109.5 112.5 -50 1.1*V(4)]); %This is for T=10;
    %V=axis; axis([t0 t0+dt -50 0]); %This is for T=10;
    title('Cross-correlation Short CCT Test (dB)')
    xlabel('Time (seconds)');
 
    
    drrr
    
    figure(3), clf, hold off;
    subplot(211),
    %Mnorm = max(abs(S));
    A = A./ Mnorm;
    
    plot(timeS,(A), 'r');grid on;axis('tight'); hold on
    %plot(timeS,(S), 'r+');grid on;axis('tight'); hold on
    
    t0 = 70.726; dt = 2/344; %This is for T=10;
    %t0 = 36.4633; dt = 0.0005; %2/344; %This is for T=10;
    %t0 = 3.872; dt = 20/344; %This is for T=2;
    %t0 = 20.306; dt = 2/344; %This is for T=10;
    V=axis; axis([t0 t0+dt V(3) 1.1*V(4)]); %This is for T=10;
    title('Cross-correlation Short CCT Test')
    xlabel('Time (seconds)');
    
    
    dum = 20 .* log10(abs(A));
    dum = dum - max(dum);
    
    subplot(212),
    plot(timeS,dum, 'r');grid on;axis('tight'); hold on
    plot(timeS,dum, 'r+');grid on;axis('tight'); hold on
    V=axis; axis([t0 t0+dt -50 0]); %This is for T=10;
    title('Cross-correlation Short CCT Test (dB)')
    xlabel('Time (seconds)');
 
    
    
    
    
    
    
    
%    figure(3),clf
%    t0 = 36.45; dt = 0.05; %2/344; %This is for T=10;
%    plot(timeS,abs(A), 'b');grid on;axis('tight');hold on
%    plot(timeS,abs(A), 'b+');grid on;axis('tight');
%    V=axis; axis([t0 t0+dt V(3) 5000]); %This is for T=10;
%    %V=axis; axis([t0 t0+dt -800 800]); %This is for T=2;
%    %V=axis; axis([t0-1.5 t0+dt+1.5 -3500 3500]); %This is for T=10;
%    %V=axis; axis([t0 t0+dt -35000 35000]); %This is for T=10;
%    title('The OPPOSITE-SIDE channel response: HRIR_A')
%    xlabel('Time (seconds)');

    
    
    
    
    
   
    
%2    
%    figure(2),clf, hold off;
%    mesh(xt, yf, SpecIR_S)
%    shading('interp');axis('tight');view( [0 90]); 
%    %t0 = 20.306; dt = T; %This is for T=10;
%    %V = axis; axis([t0-0.5 t0+dt+0.5 0 Fs/2]);
%    %V = axis; axis([V(1) V(2) 0 Fs/2]);
%    view( [0 90])
%    xlabel('Time (seconds)'); ylabel('Frequency (Hz)');
%    title('Spectrogram of the signal return'); colorbar;
%    colormap(NewColorMapPDF);colorbar;
%    clear Ssig; 
%    drawnow; pause(2);
    
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
%fScap = fScap ./Hw2;
%fAcap = fAcap ./Hw2;
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
    %subplot(311),
    semilogx(freqS, fScap-36.3, 'r');grid on;axis('tight'); hold on
    V=axis; axis([1000 Fs/2 -70 V(4)]); %This is for T=2;
    title ('Normalized: fScap')
    xlabel('Frequency (Hertz)');
    
    %subplot(312),
    semilogx(freqS, fAcap-36.3, 'b');grid on;axis('tight');
    V=axis; axis([1000 Fs/2 -70 V(4)]); %This is for T=2;
    %title ('Normalized: fAcap/|H(w)|^2')
    xlabel('Frequency (Hertz)');



    
    
    
    
    
    
    
mytoc; 



%% THE END  THE END  THE END  THE END  THE END  THE END  THE END  THE END 