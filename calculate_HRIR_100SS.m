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
%..101  180329  ..  Time window should be 2^N and need to calculate HRTF of
%                   that length not more. Then calculate also the phase plots.
%               SS  I believe this time we got S and A. Report it then plan
%                   the next step INVERSE FILTERING.
%..100          ..  New data. In-ear Knowles!! Good luck!
%                   Overlapping CENTER shows that it is consistent above f>400 Hz.
%...99  180328  SS  Will overlap plots of the HRTF-CENTER for each pulse
%                   Now captured IRs are ready for FFT.
%               OK  Super! I found that Shorted cable is sampled slower
%                   such that samples are relatively faster; order 20/Fs!!
%...98          ..  Shorted cable 20/Fs faster per pulse. Check if TRUE.
%...97          OK  Signals are captured and time-windowed. Now FFT them.
%...96          xx  For some reason PRI of the shorted cable is shorter??
%                   Compression OK, and capture values are found.
%...95          OK  180327-RecDataSH-T10-1e2To20e3.wav is OK. 
%               OK  Compressed signals OK. Now capture HRIR.
%               ..  all lines are using myxcorr_03 not the test version.
%...95  180327  ..  I now have the center recording with micro mics.
%...94          ..  This version is out of control. Wrap it up with the 
%                   new micro-mic recordings. 
%               SS  Now try to normalize with time windowed CENTER signal!
%               ..  Center ch S and A are now equal amplitude in TF.
%...93          ..  Nooo. Center was NOT IDENTICAL.
%               OK  Center is NOT a problem. I deleted the DC term to find
%                   the left and right channels of CENTER to be identical!
%               ..  Center is still a problem. Mics have directivity and TF
%                   variation as a function of angle! and I believe the TFs
%                   of the speakers vary a bit.
%                   You should multiply ALL the TFs and then design an
%                   inverse filter for ALL those TFs. Do it then.
%...92          OK  Shorter window does NOT help. It is an uncausal filter!
%...91          xx  Check Fig.1 time warping destroys HRTFs.Filtered a wider range.
%...90          SS  This is the first working time-windowed version.  
%...89          OK  Synced them. Now time window them.
%...88          OK  Check Fig.3 to see that S and A are not sync in time.
%               ..  I will shift all captured signals to center, then
%                   will short them using time windowing.
%                   You need to time window the signal to narrow it down.
%...87          OK  Time domain looks ugly (as expected) since there is a 
%                   long convolution of both signals; FS and center. 
%                   Note: I normalized using CENTER signal. Check the time
%                   domain signal if it is still an impulse??
%               SS  This is the version of this morning. HRTF's looks good.
%                   Data looks better when have enough amplitude. It was weak before.
%...86  180316  ..  Early morning recordings. Hope center works this time.
%               Kinda. Loads laptop data. HRTF is OK but CENTER is not good
%                   looking. Couldn't normalize successfully.
%...85          ..  Load laptop recordings.
%                   This was the PC loudspeaker signals. Try laptops!
%...84          xx  I believe center channel has a problem. 
%                   I will try to add S and A for Center to average it.
%...83          Kinda. Center signal has low compressed signal amplitude?
%...82          OK  Capture SHORTED CENTER SIG are plotted in time.
%...81          ..  CENTER channel very noisy for In Ear recordings. Go PC
%                   new data load. not sure about the channels. Check later.
%...80          ..  Here is what I need to do: (1) load new data files for
%                   the in-ear microphone recordings. (2) Load also the
%                   center free-space. (3)Finding HRTFs: First use cable 
%                   response to find t=0 instants, then normalize S and A 
%                   HRTFs with the center response.
%...79          ..  TX wave is changed to f1:20 f2:20K.
%...78  180315  OK  Working copy of 77. Now will process in-ear recordings for PC and laptop
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
    myaddpath_03;
    mytic;
    load NewColorMapPDF;colormap(NewColorMapPDF);
    %Cleaning up
    figure(1), clf; 
    figure(2), clf; 
    %figure(3), clf; 
    %figure(4), clf; 
    %figure(5), clf; 
    drawnow; pause(1)
    
    
%% CREATE SWEEP FUNCTION
    %Fs = 48e3; Ts = 1/Fs;
    %T  = 10; %seconds
    %f1 = 1e2;        
    %f2 = 20e3;
    
%%NOW LOAD: ts,  cntr_left,  shorted_left,  sig_left/sig_right
    
    %% 1. Load single TX pulse  
    filename = '180327-AmpESSts-T10-1e2To20e3.wav';
    [dum,Fs] = audioread(filename);
    ts(:) = dum(:,1);
    timets = ([1:length(ts)]./Fs)';

    %% 2. Load recorded SHORTED reference signal
    filename = '180327-RecDataSH-T10-1e2To20e3.wav';
    [dum,Fs] = audioread(filename);
    shorted_left(:)  = dum(:,1);
    timeshr = [1:length(shorted_left)]./Fs;

    %% 3. Load recorded CENTER reference signal
    filename = '180328-RecDataCEN-T10-1e2To20e3.wav';
    %filename = '180327-RecDataFS-T10-1e2To20e3.wav';
    [dumcntr,Fs] = audioread(filename);
    Chan = 1; %should be 1, tested!
    cntr_left(:)  = dumcntr(:,Chan);   %you need to check if this is the right channel
    timecntr = [1:length(cntr_left)]./Fs;
    
        %% 4. Load recorded FS (free-space) signal
    filename = '180328-RecDataRIG-T10-1e2To20e3.wav';
    %filename = '180327-RecDataFS-T10-1e2To20e3.wav';
    %filename = 'RecData2-T20-20TO20K-180316-FS-Laptop-InEar-FJ.wav';
    %filename = 'RecData2-T20-20TO20K-180314-FS-Laptop-InEar-nohand.wav';
    %filename = 'RecData2-T20-20TO20K-180314-FS-Loudsp-wohands.wav';
    %filename = 'RecData2-T20-20TO20K-180314-FS-Laptop-InEar-nohand.wav';
    [RoomIR,Fs] = audioread(filename);
    sig_left (:) = RoomIR (:,1);
    %sig_right(:) = RoomIR (:,2);  %this channel is garbage
    timerec = [1:length(sig_left)]./Fs;


%1  
    figure(1), clf, hold off;
        subplot(412), 
    plot(timeshr, shorted_left)
    axis('tight'); V=axis; grid on;
    axis([0 V(2) -1 1]);
    title('Received signal (uncompressed) - Shorted');
        subplot(413), 
    plot(timecntr, cntr_left)
    axis('tight'); V=axis; grid on;
    %axis([0 V(2) -1 1]);
    title('Received signal (uncompressed) - Center');
        subplot(414), 
    plot(timerec, sig_left)
    axis('tight'); V=axis; grid on;
    %axis([0 V(2) -1 1]);
    title('Received signal (uncompressed) - Free-space');
        subplot(411), 
    plot(timets, ts)
    axis('tight'); V=axis; grid on;
    %axis([0 V(2) -1 1]);
    title('AmpESS waveform');
    
    
    
    
    
%% CALCULATE AUTO- AND CROSS-CORRELATION
    fdom = 1;
    dB = 0;
    dBmax = 0;
    dBmin = -30;
    
    %OUTPUT DATA: shortedcor cntrcor sigLcor
    [shortedcor,shortedfcor,shortedcortime,shortedcorfreq,N] = ...
        myxcorr_03(shorted_left,ts,Fs,fdom,dB,dBmin,dBmax);
    [cntrcor,cntrfcor,cntrcortime,cntrcorfreq,N] = ...
        myxcorr_03(cntr_left,ts,Fs,fdom,dB,dBmin,dBmax);
    [sigLcor,sigLfcor,sigLcortime,sigLcorfreq,N] = ...
        myxcorr_03(sig_left,ts,Fs,fdom,dB,dBmin,dBmax);
    
    %% FIND THE DELAY TIMES TO CAPTURE S AND A    
    % manual search for peaks
    figure(2), clf, hold off;
    
    Nfft    = 1024; %512; % 2048; %255; %delT * Fs
    Ndel    = Nfft; 
    tdel    = Nfft/Fs;
    tpls    = 11; % next pulse
    Npls    = tpls * Fs;
    
    Ndm = 20; %for some reason PRI of shorted cable is shorter??
    Nplshort = Npls - Ndm;
    
    
    Nshrt = 599793 -Nfft/4 +18 .*Nplshort; Tshrt = Nshrt/Fs;
    Ncntr = 600247 -Nfft/4 +18 .*Npls    ; Tcntr = Ncntr/Fs;
    Nsign = 600247 -Nfft/4 +18 .*Npls    ; Tsign = Ncntr/Fs;
    
    %There are 11 seconds time-length between S   and A pulses.
    %There are 22 seconds time-length between S/A and next S/A pulses.
    
        subplot(311), 
    %plot(shortedcor,'r-')
    %V=axis;%axis([Nshrt Nshrt+Ndel V(3) V(4)]);
    plot(shortedcortime, shortedcor)
    axis('tight'), V=axis; 
    title('Compressed signal - Shorted');
    xlabel('Time (seconds)');grid on;
    %legend('Alternating S and A','location','Northwest');
        subplot(312), 
    %plot(cntrcor,'r-')
    %V=axis;%axis([Ncntr Ncntr+Ndel V(3) V(4)]);
    plot(cntrcortime, cntrcor)
    axis('tight')
    V2=axis; axis([ V(1) V(2) V2(3) V2(4) ]);
    %V=axis;axis([xmin xmin+xdel V(3) V(4)]);
    title('Compressed signal - Center');
    xlabel('Time (seconds)');grid on;
    %legend('Alternating S and A','location','Northwest');
        subplot(313), 
    %plot(sigLcor,'r-');
    %V=axis;%axis([Ncntr Ncntr+Ndel V(3) V(4)]);
    %axis([Nsign Nsign+Ndel V(3) -V(3)]);
    %V=axis;axis([Nshrt Nshrt+Ndel V(3) V(4)]);
    plot(sigLcortime, sigLcor,'r-');
    %axis('tight')
    %V=axis;axis([xmin xmin+xdel V(3) V(4)]);
    title('Compressed signal - Signal');
    xlabel('Time (seconds)');grid on;
    %legend('Alternating S and A','location','Northwest');

    
    
    %% CAPTURE signals using rectangular time-windows.
    %Signals; ts,  cntr_left,  shorted_left,  sig_left/sig_right
    
    %Nfft    = 1024; %512; % 2048; %255; %delT * Fs
    %Ndel    = Nfft; 
    %tdel    = Nfft/Fs;
    %tpls    = 11; % next pulse
    %Npls    = tpls * Fs;
    
    
    Ncor = [1:Nfft];
    tcor = (Ncor-1)./Fs;
    
    shrS(Ncor) = shortedcor([Nshrt         : Nshrt+Nfft-1]);
    shrA(Ncor) = shortedcor([Nshrt+Nplshort: Nshrt+Nplshort+Nfft-1]);
    
    cntrS(Ncor) = cntrcor([Ncntr      : Ncntr+Nfft-1]);
    cntrA(Ncor) = cntrcor([Ncntr+Npls : Ncntr+Npls+Nfft-1]);
    
    %we only have cntrl so copy that. Later use the recorded signal.
    signS(Ncor) = sigLcor([Ncntr      : Ncntr+Nfft-1]);
    signA(Ncor) = sigLcor([Ncntr+Npls : Ncntr+Npls+Nfft-1]);
    %signS(Ncor) = cntrcor([Ncntr      : Ncntr+Nfft-1]);
    %signA(Ncor) = cntrcor([Ncntr+Npls : Ncntr+Npls+Nfft-1]);
    
    
     %% now take the FFT after zero padding. 
    twin = zeros(1,Nfft);
    twin([200:800]) = ones(1,601);
    twin([801:1000]) = [200:-1:1]./200;  %[49:-1:0]./49;
    % time window is tested to be correct.
    
    shrS = shrS .* twin;
    shrA = shrA .* twin;
    
    cntrS = cntrS .* twin;
    cntrA = cntrA .* twin;
    
    signS = signS .* twin;
    signA = signA .* twin;
    
    
    %2   
    figure(3), %clf, hold off;
        subplot(411), 
    plot(tcor, shrS,'r-'), hold on
    plot(tcor, shrA,'b-'), axis('tight');
    %V=axis;axis([0 500 V(3) V(4)]);
    %V=axis;axis([0 xdel V(3) V(4)]);
    title('Captured signal - SHORTED');
    xlabel('Time (seconds)');
        subplot(412), 
    plot(tcor, cntrS,'r-'), hold on
    plot(tcor, cntrA,'b-'), axis('tight');
    %V=axis;axis([0 500 V(3) V(4)]);
    %V=axis;axis([0 xdel V(3) V(4)]);
    title('Captured signal - CENTER');
    xlabel('Time (seconds)');
        subplot(413), 
    plot(tcor, signS, 'r-'), hold on;
    plot(tcor, signA, 'b-'), axis('tight');
%    %V=axis;axis([0 500 V(3) V(4)]);
%    %V=axis;axis([0 xdel V(3) V(4)]);
%    title('Compressed signal - SIGNAL');
%    xlabel('Time (seconds)');
        subplot(414), 
    plot(tcor, twin, 'r-'), axis('tight');
    %axis([0 500 0 1])
    title('Time windowing function');
    xlabel('Time (seconds)');
    grid on;    
    
    
    
    
    
    %% CALCULATE IN FREQUENCY DOMAIN: HRTR
    %Calculate:  fshrS, fshrA   fcntrS, fcntrA   fsignS, fsignA   
    
    fshrS = abs( fft(shrS));
    fshrA = abs( fft(shrA));
    
    fcntrS = abs( fft(cntrS));
    fcntrA = abs( fft(cntrA));
    
    fsignS = abs( fft(signS));
    fsignA = abs( fft(signA));
    
    
    freq = [0:length(fsignS)-1]./(length(fsignS)-1).*Fs;
    
    fSHmax = max(max(abs(fshrS ), abs(fshrA )));
    fSCmax = max(max(abs(fcntrS), abs(fcntrA)));
    fSGmax = max(max(abs(fsignS), abs(fcntrA)));
    
    %fshrS    = fshrS ./ fSHmax;
    %fshrA    = fshrA ./ fSHmax;
    
    %fcntrS    = fcntrS ./ fSCmax;
    %fcntrA    = fcntrA ./ fSCmax;
    
    %fsignS    = fsignS ./ fSGmax;
    %fsignA    = fsignA ./ fSGmax;
    
    %fSnorm = fsignS./fshrS;
    %fAnorm = fsignA./fshrA;
    
    fSnorm = fsignS./fcntrS;
    fAnorm = fsignA./fcntrA;
    
    %SIGSnorm = ifft(fSnorm);
    %SIGAnorm = ifft(fAnorm);
    
    
    %figure(1),clf
    %plot(SIGSnorm, 'r-'), hold on;
    %plot(SIGAnorm, 'b-'), hold on;
    %axis('tight');
    
    
    
    
    %% LINEAR PHASE TERM: Calculate the to-be-subtracted linear phase term
    %linear_phase_deg is -1.424586663e5 at 20 KHz.
    %linear_phase_rad is -2.4864e+03    at 20 KHz
    %freq = [0:length(fshortS)-1]./(length(fshortS)-1).*Fs;
    %
    %lineartimefactor = freq./20e3.*2.4864e3;
    %linearphaseterm = exp(i.*lineartimefactor);
    %
    %fS      = fS      .* linearphaseterm;
    %fA      = fA      .* linearphaseterm;
    %fshortS = fshortS .* linearphaseterm;
    %fshortA = fshortA .* linearphaseterm;
    %fSnorm  = fSnorm  .* linearphaseterm;
    %fAnorm  = fAnorm  .* linearphaseterm;
    %
    %S = ifft(fS);
    %A = ifft(fA);
    %
    %ShortedS = ifft(fshortS);
    %ShortedA = ifft(fshortA);
    
    % calculate phase responses
    %phase_SSSS      = unwrap(atan2(imag(fS),      real(fS)));
    %phase_AAAA     = unwrap(atan2(imag(fA),      real(fA)));
    %phase_SHRT    = unwrap(atan2(imag(fshortS), real(fshortS)));
    %phase_shortA    = unwrap(atan2(imag(fshortA), real(fshortA)));
    %phase_Snorm     = unwrap(atan2(imag(fSnorm),  real(fSnorm)));
    %phase_Anorm     = unwrap(atan2(imag(fAnorm),  real(fAnorm)));
    %
    %phase_S         = phase_S /pi * 180;
    %phase_A         = phase_A /pi * 180;
    %phase_shortS    = phase_shortS /pi * 180;
    %phase_shortA    = phase_shortA /pi * 180;
    %phase_Snorm     = phase_Snorm /pi * 180;
    %phase_Anorm     = phase_Anorm /pi * 180;
    %
    % calculate time-delay responses
    %delay_S         = 1000 .* phase_S ./ 360 ./ freq;
    %delay_A         = 1000 .* phase_A ./ 360 ./ freq;
    %delay_shortS    = 1000 .* phase_shortS ./ 360 ./ freq;
    %delay_shortA    = 1000 .* phase_shortA ./ 360 ./ freq;
    %delay_Snorm     = 1000 .* phase_Snorm ./ 360 ./ freq;
    %delay_Anorm     = 1000 .* phase_Anorm ./ 360 ./ freq;
    
    Fmn = 100;
    Fmx = 20e3;
    
       
%2  PLOT CORRELATION    
    %figure(2), clf
    %    subplot(211), 
    %plot(tcor, real(ShortedS), 'r-'), hold on;
    %plot(tcor, imag(ShortedS), 'm-'), hold on;
    %plot(tcor, real(ShortedA)-0.1, 'b-')
    %plot(tcor, imag(ShortedA)-0.1, 'm-')
    %axis('tight'); grid on;
    %title('HRTFs (red)S (blue)A-0.1 (Shorted cable)')
    %xlabel('Time (seconds)')
    %    subplot(212), 
    %plot(tcor, real(S), 'r-'); hold on;
    %plot(tcor, imag(S), 'm-'); hold on;
    %plot(tcor, real(A)-0.1, 'b-');
    %plot(tcor, imag(A)-0.1, 'm-');
    %axis('tight'); grid on;
    %title('HRTFs (red)S (blue)A-0.1 (Kirkhead-Loudspeaker-PC)')
    %xlabel('Time (seconds)')
    
    %Plot: fSIGS, fSIGA   fSHRS, fSHRA     fCENS, fCENA
    figure(4), %clf, hold off
        subplot(411),
    semilogx(freq, 20 .* log10(abs(fshrS)), 'r-'); hold on;
    semilogx(freq, 20 .* log10(abs(fshrA)), 'b-');
    V=axis; axis([Fmn Fmx 25 95]); 
    %V=axis; axis([Fmn Fmx -80 20]); 
    set(gca,'xtick',[Fmn  1000 10000 Fmx]); grid on;
    title ('HRTF Shorted Cable');
    legend('S cable','A cable','location','Southwest');
    xlabel('Frequency (Hertz)');
    ylabel('Desibels');
    
        subplot(412),
    semilogx(freq, 20 .* log10(abs(fcntrS)), 'r-'); hold on;
    semilogx(freq, 20 .* log10(abs(fcntrA)), 'b-');
    set(gca,'xtick',[Fmn  1000 10000 Fmx]); grid on;
    V=axis; axis([Fmn Fmx  25 95]); grid on;
    %V=axis; axis([Fmn Fmx  -60 10]); grid on;
    title ('HRTF  Freespace Center');
    legend('S cable','A cable','location','Southwest');
    xlabel('Frequency (Hertz)')
    ylabel('Desibels');
    
        subplot(413),
    semilogx(freq, 20 .* log10(abs(fsignS)), 'r-'); hold on;
    semilogx(freq, 20 .* log10(abs(fsignA)), 'b-');
    set(gca,'xtick',[Fmn  1000 10000 Fmx]); grid on;
    V=axis; axis([Fmn Fmx  25 95]); grid on;
    %    V=axis; axis([Fmn Fmx  -60 10]); grid on;
%    set(gca,'xtick',[Fmn  1000 10000 Fmx]); grid on;
%    legend('HRTF S/A FS','location','Southwest');
%    xlabel('Frequency (Hertz)')
%    ylabel('Desibels');
    
        subplot(414),
    semilogx(freq, 20 .* log10(abs(fSnorm)), 'r-'); hold on;
    semilogx(freq, 20 .* log10(abs(fAnorm)), 'b-');
    set(gca,'xtick',[Fmn  1000 10000 Fmx]); grid on;
    V=axis; axis([Fmn Fmx  -40 30]); grid on;
%    V=axis; axis([Fmn Fmx  -60 10]); grid on;
%    set(gca,'xtick',[Fmn  1000 10000 Fmx]); grid on;
%    legend('HRTF(FreeSpace)/HRTF(ShortedCable)','location','South');
%    xlabel('Frequency (Hertz)')
%    ylabel('Desibels');
    
    
    
    
    durrr
    
    
    
    figure(4),
        subplot(311), 
    plot(freq, phase_shortS,'r-'); hold on;
    plot(freq, phase_shortA-5,'b-'); grid on;
    %axis('tight'); 
    %V=axis; axis([Fmn Fmx  -15e3 0.1e4]); grid on;
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