%% process_HRIR_xx.m
%
% Load recordings, all HRTFs
% divide by system response using long DFT solution.
% First trial for the inverse filters.
% S.G. Tanyer, 180331 Victoria
%

%% HISTORY
%
%
%
%
%...08  
%                   Now use MONO audio to play on only one side...
%               SS  I believe REC1 worked!!!!!!!
%               xx  REC1 and REC2 does not seem to work. Not even stereo.
%...07          ..  Use the TX audio file for testing crosstalk.
%...06          SS  S, A, invS, invA, 1/(S+A) and 1/(S-A) is ready. USE IT!
%                   S and A seems to not multiplied by twin!??
%...05          OK  Time warped inverse invS and invA ready to be processed.
%               Keep this version in case you need upsampling and filtering
%               (1) Try the inverse by time warping. 
%               xx  Don't do this. Try something else.
%...04          ..  Now (2) Low pass filter HRTFs for smoothing.
%               OK  This is;(1)zeros between each sample of HRIRs. 
%...03          ..  Try something now. 
%               (1) Place zeros between each sample of HRIRs. 
%               (2) Low pass filter HRTF. 
%               (3) Then ifft low pass filter the HRTFnorm to smooth high freqs.
%...02          OK  This shows that HRTFnorm should be filtered to limit
%                   the peaks to avoid time-warping for Snorm and Anorm.
%...01  180331  OK  Now HRTFs are loaded and plotted to be correct data.
%
%started with the version calculate_HRIR_105SS.m
%       180331  SS  Saves all HRTFs for different pulse numbers.
%..105  180330      Now is the time to record HRTFs for the INVERSE studies.


%%  INTRO: Initialization 
    clear;clc;
    %Add paths
    myaddpath_04;
    mytic;
    
    %Cleaning up
    figure(1), clf; 
    figure(2), clf; 
    figure(3), clf; 
    %figure(4), clf; 
    %figure(5), clf; 
    drawnow; pause(1)
   
    
    
%% NOW LOAD HTRFs: 
    %shorted S and A
    %filename =  '180331-HRTF-shrS-T10-1e2To20e3-p16.mat';
    %load(filename,         'shrS');
    %filename =  '180331-HRTF-shrA-T10-1e2To20e3-p16.mat';
    %load(filename,         'shrA');
    
    %center channel
    filename = '180331-HRTF-cntrS-T10-1e2To20e3-p16.mat';
    load(filename,         'cntrS');
    filename = '180331-HRTF-cntrA-T10-1e2To20e3-p16.mat';
    load(filename,         'cntrA');
    
    %in-ear micro-mic Knowles recording for the right ear
    filename = '180331-HRTF-signS-T10-1e2To20e3-p16.mat';
    load(filename,         'signS');
    filename = '180331-HRTF-signA-T10-1e2To20e3-p16.mat';
    load(filename,         'signA');
    
    
    %%
    Fs = 48e3;
    Nfft = length(cntrS);
    Nfft2 = 2*Nfft;
    time = [0:Nfft-1]./Fs;
    freq = [0:2*Nfft-1]./(2*Nfft-1).*Fs;
    Fmn = 100;
    Fmx = 20e3;
   
    
    
 %% Inverse filter by time warping
    invS([1:Nfft]) = signS([Nfft:-1:1]);
    invA([1:Nfft]) = signA([Nfft:-1:1]);
       
    invSpA = invS + invA;
    invSmA = invS - invA;
    
    
    
    %% Calculate raw HRTF's
    %fshrS = ( fft(shrS, Nfft));
    %fshrA = ( fft(shrA, Nfft));
    
    fcntrS = fft(cntrS, Nfft2);
    fcntrA = fft(cntrA, Nfft2);
    
    fsignS = fft(signS, Nfft2);
    fsignA = fft(signA, Nfft2);
    
    finvS = fft(invS, Nfft2);
    finvA = fft(invA, Nfft2);
    
    finvSpA = fft(invSpA, Nfft2);
    finvSmA = fft(invSmA, Nfft2);
    
    
    %fSnorm = fsignS./fcntrS;
    %fAnorm = fsignA./fcntrA;
    
%% Fill zeros between each sample
    %Nfft2 = 2 * Nfft;
    %freq2 = [0:Nfft2-1]./(Nfft2-1).*Fs;
    %time2 = [0:Nfft2-1]./Fs;
    
    %cntrS2 = zeros(1,Nfft2);
    %cntrA2 = cntrS2;
    %signS2 = cntrS2;
    %signA2 = cntrS2;
    
    %cntrS2([1:2:Nfft2]) = cntrS([1:Nfft]);
    %cntrA2([1:2:Nfft2]) = cntrA([1:Nfft]);
    %cntrS2([1:2:Nfft2]) = cntrS([1:Nfft]);
    %cntrA2([1:2:Nfft2]) = cntrA([1:Nfft]);
    
    %signS2([1:2:Nfft2]) = signS([1:Nfft]);
    %signA2([1:2:Nfft2]) = signA([1:Nfft]);
    %signS2([1:2:Nfft2]) = signS([1:Nfft]);
    %signA2([1:2:Nfft2]) = signA([1:Nfft]);
    
    
%% Inverse filtering to compansate for system TF.
    % Use an Nfft as long as 'necessary'
    
    %fcntrS2 = ( fft(cntrS2, Nfft2));
    %fcntrA2 = ( fft(cntrA2, Nfft2));
    
    %fsignS2 = ( fft(signS2, Nfft2));
    %fsignA2 = ( fft(signA2, Nfft2));
    
    %fSnorm2 = fsignS2./fcntrS2;
    %fAnorm2 = fsignA2./fcntrA2;

%% Calculate raw HRTF phase responses
    %phase_shrS  = rd2deg(unwrap(atan2(imag(fshrS),  real(fshrS))));
    %phase_shrA  = rd2deg(unwrap(atan2(imag(fshrA),  real(fshrA))));
    
    %phase_cntrS = rd2deg(unwrap(atan2(imag(fcntrS2), real(fcntrS2))));
    %phase_cntrA = rd2deg(unwrap(atan2(imag(fcntrA2), real(fcntrA2))));
    
    phase_signS = rd2deg(unwrap(atan2(imag(fsignS), real(fsignS))));
    phase_signA = rd2deg(unwrap(atan2(imag(fsignA), real(fsignA))));
    
    phase_invS = rd2deg(unwrap(atan2(imag(finvS), real(finvS))));
    phase_invA = rd2deg(unwrap(atan2(imag(finvA), real(finvA))));
    
    phase_invSpA = rd2deg(unwrap(atan2(imag(finvSpA), real(finvSpA))));
    phase_invSmA = rd2deg(unwrap(atan2(imag(finvSmA), real(finvSmA))));
    
    %phase_normS = rd2deg(unwrap(atan2(imag(fSnorm2), real(fSnorm2))));
    %phase_normA = rd2deg(unwrap(atan2(imag(fAnorm2), real(fAnorm2))));

    %% Calculate normalized HRIR
    %Snorm = ifft(fSnorm2, Nfft2);
    %Anorm = ifft(fAnorm2, Nfft2);
    
    
    %% PLOTTING
%1
    figure(1), %clf, hold off;
    %    subplot(311), 
    %plot(time, shrS,'r-'), hold on
    %plot(time, shrA,'b-'), axis('tight');
    %title('Captured signal - SHORTED');
    %xlabel('Time (seconds)');
        subplot(311), 
    plot(time, cntrS,'r-'), hold on
    plot(time, cntrA,'b-'), axis('tight');
    title('Captured signal - CENTER');
    xlabel('Time (seconds)');
        subplot(312), 
    plot(time, signS, 'r-'), hold on;
    plot(time, signA, 'b-'), axis('tight');
    title('Compressed signal - SIGNAL');
    xlabel('Time (seconds)');
        subplot(313), 
    plot(time, invSpA, 'r-'), hold on;
    plot(time, invSmA, 'b-'), axis('tight');
    title('Compressed signal - SIGNAL');
    xlabel('Time (seconds)');
%2    
    figure(2), %clf, hold off
    %    subplot(411),
    %semilogx(freq2, 20 .* log10(abs(fshrS))-80, 'r-'); hold on;
    %semilogx(freq2, 20 .* log10(abs(fshrA))-80, 'b-');
    %V=axis; axis([Fmn Fmx 25-80 95-80]); 
    %set(gca,'xtick',[Fmn  1000 10000 Fmx]); grid on;
    %title ('HRTF Shorted Cable');
    %legend('A cable','S cable','location','Southwest');
    %xlabel('Frequency (Hertz)');
    %ylabel('Desibels');
    %    subplot(411),
    %plot(freq, 20 .* log10(abs(fcntrS))-80, 'r-'); hold on;
    %plot(freq, 20 .* log10(abs(fcntrA))-80, 'b-');
    %set(gca,'xtick',[Fmn  1000 10000 Fmx]); grid on;
    %V=axis; axis([Fmn Fmx  25-80 95-80]); grid on;
    %title ('HRTF  Freespace - CENTER');
    %legend('A cable','S cable','location','Southwest');
    %xlabel('Frequency (Hertz)')
    %ylabel('Desibels');
        subplot(311),
    plot(freq, 20 .* log10(abs(fsignS))-80, 'r-'); hold on;
    plot(freq, 20 .* log10(abs(fsignA))-80, 'b-');
    set(gca,'xtick',[Fmn  1000 10000 Fmx]); grid on;
    V=axis; axis([Fmn Fmx  25-80 95-80]); grid on;
    title ('HRTF  Freespace - RIGHT');
    legend('A cable','S cable','location','Southwest');
    xlabel('Frequency (Hertz)')
    ylabel('Desibels');
        subplot(312),
    plot(freq, 20 .* log10(abs(finvS))-80, 'r-'); hold on;
    plot(freq, 20 .* log10(abs(finvA))-80, 'b-');
    set(gca,'xtick',[Fmn  1000 10000 Fmx]); grid on;
    V=axis; axis([Fmn Fmx  25-80 95-80]); grid on;
    title ('HRTF  Freespace - inverse S and A');
    legend('A cable','S cable','location','Southwest');
    xlabel('Frequency (Hertz)')
    ylabel('Desibels');
        subplot(313),
    plot(freq, 20 .* log10(abs(finvSpA))-80, 'r-'); hold on;
    plot(freq, 20 .* log10(abs(finvSmA))-80, 'b-');
    set(gca,'xtick',[Fmn  1000 10000 Fmx]); grid on;
    V=axis; axis([Fmn Fmx  25-80 95-80]); grid on;
    title ('HRTF  Freespace - inverse SpA and SmA');
    legend('A cable','S cable','location','Southwest');
    xlabel('Frequency (Hertz)')
    ylabel('Desibels');
    %    subplot(414),
    %plot(freq2, 20 .* log10(abs(fSnorm2)), 'r-'); hold on;
    %plot(freq2, 20 .* log10(abs(fAnorm2)), 'b-');
    %set(gca,'xtick',[Fmn  1000 10000 Fmx]); grid on;
    %V=axis; %axis([Fmn Fmx  -40 30]); grid on;
    %title ('HRTF  Freespace - RIGHT/CENTER');
    %legend('A cable','S cable','location','Southwest');
    %xlabel('Frequency (Hertz)')
    %ylabel('Desibels');
%3 
    figure(3), clf; hold off;
    %    subplot(411), 
    %semilogx(freq2, phase_shrS,'r-'); hold on;
    %semilogx(freq2, phase_shrA,'b-'); grid on;
    %axis('tight'); 
    %V=axis; axis([Fmn Fmx  -3000 0]); grid on;
    %title ('Phase HRTF Shorted Cable');
    %xlabel('Frequency (Hertz)');
    %ylabel('Phase (degrees)');
    %    subplot(412), 
    %semilogx(freq2, phase_cntrS,'r-'); hold on; 
    %semilogx(freq2, phase_cntrA,'b-'); grid on;
    %axis('tight'); 
    %V=axis; axis([Fmn Fmx  -3000 0]); grid on;
    %title ('Phase HRTF  Freespace - CENTER');
    %xlabel('Frequency (Hertz)');
    %ylabel('Phase (degrees)');
        subplot(311), 
    semilogx(freq, phase_signS,'r-'); hold on;
    semilogx(freq, phase_signA,'b-'); grid on;
    axis('tight'); 
    V=axis; axis([Fmn Fmx  -3000 0]); grid on;
    title ('Phase HRTF  Freespace - RIGHT');
    xlabel('Frequency (Hertz)');
    ylabel('Phase (degrees)');
        subplot(312), 
    semilogx(freq, phase_invS,'r-'); hold on;
    semilogx(freq, phase_invA,'b-'); grid on;
    axis('tight'); 
    V=axis; axis([Fmn Fmx  V(3) V(4)]); grid on;
    %V=axis; axis([Fmn Fmx  -3000 0]); grid on;
    title ('Phase HRTF  Freespace - RIGHT');
    xlabel('Frequency (Hertz)');
    ylabel('Phase (degrees)');
        subplot(313), 
    semilogx(freq, phase_invSpA,'r-'); hold on;
    semilogx(freq, phase_invSmA,'b-'); grid on;
    axis('tight'); 
    V=axis; axis([Fmn Fmx  V(3) V(4)]); grid on;
    %V=axis; axis([Fmn Fmx  -3000 0]); grid on;
    title ('Phase HRTF  Freespace - RIGHT');
    xlabel('Frequency (Hertz)');
    ylabel('Phase (degrees)');
    %    subplot(414), 
    %semilogx(freq2, phase_normS,'r-'); hold on;
    %semilogx(freq2, phase_normA,'b-'); grid on;
    %axis('tight'); 
    %V=axis; axis([Fmn Fmx  -2000 1000]); grid on;
    %title ('Phase HRTF  Freespace - RIGHT/CENTER');
    %xlabel('Frequency (Hertz)');
    %ylabel('Phase (degrees)');
%4 
    %figure(4), clf; hold off;
    %    plot(time2, Snorm, 'r-'); hold on, grid on;
    %    plot(time2, Anorm, 'b-');
        
    


    %% 1. Load single TX pulse  
    filename = '180327-AmpESSts-T10-1e2To20e3.wav';
    [dum,Fs] = audioread(filename);
    ts(:) = dum(:,1);
    timets = ([1:length(ts)]./Fs)';

    %% 2. Load recorded SHORTED reference signal
    filename = '180327-RecDataSH-T10-1e2To20e3.wav';
    [dum,Fs] = audioread(filename);
    shorted(:)  = dum(:,1);
    timeshr = [1:length(shorted)]./Fs;

    %% 3. Load recorded CENTER reference signal
    filename = '180328-RecDataCEN-T10-1e2To20e3.wav';
    [dumcntr,Fs] = audioread(filename);
    Chan = 1; 
    cntr(:)  = dumcntr(:,Chan);   
    timecntr = [1:length(cntr)]./Fs;
    
    %% 4. Load recorded FS (free-space) signal
    filename = '180328-RecDataRIG-T10-1e2To20e3.wav';
    [RoomIR,Fs] = audioread(filename);
    sig(:) = RoomIR (:,1);
    timerec = [1:length(sig)]./Fs;
%1  
    figure(1), clf, hold off;
        subplot(412), 
    plot(timeshr, shorted)
    axis('tight'); V=axis; grid on;
    axis([0 V(2) -1 1]);
    title('Received signal (uncompressed) - Shorted');
        subplot(413), 
    plot(timecntr, cntr)
    axis('tight'); V=axis; grid on;
    %axis([0 V(2) -1 1]);
    title('Received signal (uncompressed) - Center');
        subplot(414), 
    plot(timerec, sig)
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
        myxcorr_03(shorted,ts,Fs,fdom,dB,dBmin,dBmax);
    [cntrcor,cntrfcor,cntrcortime,cntrcorfreq,N] = ...
        myxcorr_03(cntr,ts,Fs,fdom,dB,dBmin,dBmax);
    [sigLcor,sigLfcor,sigLcortime,sigLcorfreq,N] = ...
        myxcorr_03(sig,ts,Fs,fdom,dB,dBmin,dBmax);
    
    %% FIND THE DELAY TIMES TO CAPTURE S AND A    
    % manual search for peaks
    
    %select wide analysis window to be 1024. Capture 512.
    Nfft    = 1024; %512; % 2048; %255; %delT * Fs
    Ndel    = Nfft; 
    tdel    = Nfft/Fs;
    tpls    = 11; % next pulse
    Npls    = tpls * Fs;
    
    Ndm = 20; %for some reason PRI of shorted cable is shorter??
    Nplshort = Npls - Ndm;
    
    Npulse = 18; %select an even number in {0,18}
    Nshrt = 250 + 599793 -Nfft/4 +Npulse .*Nplshort; Tshrt = Nshrt/Fs;
    Ncntr = 250 + 600247 -Nfft/4 +Npulse .*Npls    ; Tcntr = Ncntr/Fs;
    Nsign = 250 + 600247 -Nfft/4 +Npulse .*Npls    ; Tsign = Ncntr/Fs;
    
    %There are 11 seconds time-length between S   and A pulses.
    %There are 22 seconds time-length between S/A and next S/A pulses.
%2
    figure(2), clf, hold off;
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
    
    Nff2 = Nfft/2;
    Ncor = [1:Nff2];
    tcor = (Ncor-1)./Fs;
    
    shrS(Ncor) = shortedcor([Nshrt         : Nshrt+Nff2-1]);
    shrA(Ncor) = shortedcor([Nshrt+Nplshort: Nshrt+Nplshort+Nff2-1]);
    
    cntrS(Ncor) = cntrcor([Ncntr      : Ncntr+Nff2-1]);
    cntrA(Ncor) = cntrcor([Ncntr+Npls : Ncntr+Npls+Nff2-1]);
    
    %we only have cntrl so copy that. Later use the recorded signal.
    signS(Ncor) = sigLcor([Ncntr      : Ncntr+Nff2-1]);
    signA(Ncor) = sigLcor([Ncntr+Npls : Ncntr+Npls+Nff2-1]);
    
%%  SAVE HRTFs.
    %shorted S and A
    %filename =  '180331-HRTF-shrS-T10-1e2To20e3-p6.mat';
    %save(filename,         'shrS');
    %filename =  '180331-HRTF-shrA-T10-1e2To20e3-p6.mat';
    %save(filename,         'shrA');

    %filename = '180331-HRTF-cntrS-T10-1e2To20e3-p6.mat';
    %save(filename,         'cntrS');
    %filename = '180331-HRTF-cntrA-T10-1e2To20e3-p6.mat';
    %save(filename,         'cntrA');

    %filename = '180331-HRTF-signS-T10-1e2To20e3-p6.mat';
    %save(filename,         'signS');
    %filename = '180331-HRTF-signA-T10-1e2To20e3-p6.mat';
    %save(filename,         'signA');
    

    
    %% now take the FFT after zero padding. 
    twin = zeros(1,Nff2);
    twin([1:Nff2/2]) = ones(1,Nff2/2);
    twin([Nff2/2+1:Nff2-Nff2/4]) = 4 .*[Nff2/4:-1:1]./Nff2;
    
    shrS = shrS .* twin;
    shrA = shrA .* twin;
    
    cntrS = cntrS .* twin;
    cntrA = cntrA .* twin;
    
    signS = signS .* twin;
    signA = signA .* twin;
 %3   
    figure(3), %clf, hold off;
        subplot(411), 
    plot(tcor, shrS,'r-'), hold on
    plot(tcor, shrA,'b-'), axis('tight');
    title('Captured signal - SHORTED');
    xlabel('Time (seconds)');
        subplot(412), 
    plot(tcor, cntrS,'r-'), hold on
    plot(tcor, cntrA,'b-'), axis('tight');
    title('Captured signal - CENTER');
    xlabel('Time (seconds)');
        subplot(413), 
    plot(tcor, signS, 'r-'), hold on;
    plot(tcor, signA, 'b-'), axis('tight');
    title('Compressed signal - SIGNAL');
    xlabel('Time (seconds)');
        subplot(414), 
    plot(tcor, twin, 'r-'), axis('tight');
    title('Time windowing function');
    xlabel('Time (seconds)');
    grid on;    
    
    
    %% CALCULATE IN FREQUENCY DOMAIN: HRTR
    %Calculate:  fshrS, fshrA   fcntrS, fcntrA   fsignS, fsignA   
    
    fshrS = ( fft(shrS, Nfft));
    fshrA = ( fft(shrA, Nfft));
    
    fcntrS = ( fft(cntrS, Nfft));
    fcntrA = ( fft(cntrA, Nfft));
    
    fsignS = ( fft(signS, Nfft));
    fsignA = ( fft(signA, Nfft));
    
    freq = [0:length(fsignS)-1]./(length(fsignS)-1).*Fs;
    
    %fSHmax = max(max(abs(fshrS ), abs(fshrA )));
    %fSCmax = max(max(abs(fcntrS), abs(fcntrA)));
    %fSGmax = max(max(abs(fsignS), abs(fcntrA)));
    %
    %fshrS    = fshrS ./ fSHmax;
    %fshrA    = fshrA ./ fSHmax;
    %
    %fcntrS    = fcntrS ./ fSCmax;
    %fcntrA    = fcntrA ./ fSCmax;
    %
    %fsignS    = fsignS ./ fSGmax;
    %fsignA    = fsignA ./ fSGmax;
    %
    %fSnorm = fsignS./fshrS;
    %fAnorm = fsignA./fshrA;
    
    fSnorm = fsignS./fcntrS;
    fAnorm = fsignA./fcntrA;
%4    
    %Plot: fSIGS, fSIGA   fSHRS, fSHRA     fCENS, fCENA
    figure(4), %clf, hold off
        subplot(411),
    semilogx(freq, 20 .* log10(abs(fshrS))-80, 'r-'); hold on;
    semilogx(freq, 20 .* log10(abs(fshrA))-80, 'b-');
    V=axis; axis([Fmn Fmx 25-80 95-80]); 
    set(gca,'xtick',[Fmn  1000 10000 Fmx]); grid on;
    title ('HRTF Shorted Cable');
    legend('A cable','S cable','location','Southwest');
    xlabel('Frequency (Hertz)');
    ylabel('Desibels');
        subplot(412),
    semilogx(freq, 20 .* log10(abs(fcntrS))-80, 'r-'); hold on;
    semilogx(freq, 20 .* log10(abs(fcntrA))-80, 'b-');
    set(gca,'xtick',[Fmn  1000 10000 Fmx]); grid on;
    V=axis; axis([Fmn Fmx  25-80 95-80]); grid on;
    title ('HRTF  Freespace - CENTER');
    legend('A cable','S cable','location','Southwest');
    xlabel('Frequency (Hertz)')
    ylabel('Desibels');
        subplot(413),
    semilogx(freq, 20 .* log10(abs(fsignS))-80, 'r-'); hold on;
    semilogx(freq, 20 .* log10(abs(fsignA))-80, 'b-');
    set(gca,'xtick',[Fmn  1000 10000 Fmx]); grid on;
    V=axis; axis([Fmn Fmx  25-80 95-80]); grid on;
    title ('HRTF  Freespace - RIGHT');
    legend('A cable','S cable','location','Southwest');
    xlabel('Frequency (Hertz)')
    ylabel('Desibels');
        subplot(414),
    semilogx(freq, 20 .* log10(abs(fSnorm)), 'r-'); hold on;
    semilogx(freq, 20 .* log10(abs(fAnorm)), 'b-');
    set(gca,'xtick',[Fmn  1000 10000 Fmx]); grid on;
    V=axis; axis([Fmn Fmx  -40 30]); grid on;
    title ('HRTF  Freespace - RIGHT/CENTER');
    legend('A cable','S cable','location','Southwest');
    xlabel('Frequency (Hertz)')
    ylabel('Desibels');

    
  
%% Calculate PHASE RESPONSES
    
    fshrS = ( fft(shrS, Nfft));
    fshrA = ( fft(shrA, Nfft));
    
    fcntrS = ( fft(cntrS, Nfft));
    fcntrA = ( fft(cntrA, Nfft));
    
    fsignS = ( fft(signS, Nfft));
    fsignA = ( fft(signA, Nfft));
    
    
    phase_shrS  = rd2deg(unwrap(atan2(imag(fshrS),  real(fshrS))));
    phase_shrA  = rd2deg(unwrap(atan2(imag(fshrA),  real(fshrA))));
    
    phase_cntrS = rd2deg(unwrap(atan2(imag(fcntrS), real(fcntrS))));
    phase_cntrA = rd2deg(unwrap(atan2(imag(fcntrA), real(fcntrA))));
    
    phase_signS = rd2deg(unwrap(atan2(imag(fsignS), real(fsignS))));
    phase_signA = rd2deg(unwrap(atan2(imag(fsignA), real(fsignA))));
    
    phase_normS = rd2deg(unwrap(atan2(imag(fSnorm), real(fSnorm))));
    phase_normA = rd2deg(unwrap(atan2(imag(fAnorm), real(fAnorm))));
%1    
    figure(1), clf; hold off;
        subplot(411), 
    semilogx(freq, phase_shrS,'r-'); hold on;
    semilogx(freq, phase_shrA,'b-'); grid on;
    axis('tight'); 
    V=axis; axis([Fmn Fmx  -3000 0]); grid on;
    title ('Phase HRTF Shorted Cable');
    xlabel('Frequency (Hertz)');
    ylabel('Phase (degrees)');
        subplot(412), 
    semilogx(freq, phase_cntrS,'r-'); hold on; 
    semilogx(freq, phase_cntrA,'b-'); grid on;
    axis('tight'); 
    V=axis; axis([Fmn Fmx  -3000 0]); grid on;
    title ('Phase HRTF  Freespace - CENTER');
    xlabel('Frequency (Hertz)');
    ylabel('Phase (degrees)');
        subplot(413), 
    semilogx(freq, phase_signS,'r-'); hold on;
    semilogx(freq, phase_signA,'b-'); grid on;
    axis('tight'); 
    V=axis; axis([Fmn Fmx  -3000 0]); grid on;
    title ('Phase HRTF  Freespace - RIGHT');
    xlabel('Frequency (Hertz)');
    ylabel('Phase (degrees)');
        subplot(414), 
    semilogx(freq, phase_normS,'r-'); hold on;
    semilogx(freq, phase_normA,'b-'); grid on;
    axis('tight'); 
    V=axis; axis([Fmn Fmx  -2000 1000]); grid on;
    title ('Phase HRTF  Freespace - RIGHT/CENTER');
    xlabel('Frequency (Hertz)');
    ylabel('Phase (degrees)');
    
    
    
%%  Now load an audio file and make it MONO. Then, play on one side only.
    filename = 'CubanSalsa.mp3';
    [dum,Fs] = audioread(filename);
    dur
    
    
%%  Now we have all the filters;  S, A, invS, invA, 1/(S+A) and 1/(S-A)  
    %180327-TXDataFS-T10-1e2To20e3
    %filename = '180327-TXDataFS-T10-1e2To20e3.wav';
    %[dum,Fs] = audioread(filename);
    
    
    
    
    
    len = length(dum);
    L(:,1) = dum(:,1);
    R(:,1) = dum(:,2);
    
    Lprime = conv( (L-R) , invSmA) + conv( (L+R) , invSpA);
    Rprime = conv( (L+R) , invSpA) - conv( (L-R) , invSmA);
    
    REC1(:,1) =  Lprime;
    REC1(:,2) = -Rprime;
    RECmx = max(max(REC1));
    REC1 = REC1 ./ RECmx;
    
    %normalize REC
    REC2(:,2) = Lprime;
    REC2(:,1) = Rprime;
    RECmx = max(max(REC2));
    REC2 = REC2 ./ RECmx;
    
    %% write audio files 
    filename = '180402-REC1-T10-1e2To20e3.wav';
    audiowrite(filename,REC1,Fs);
    
    %filename = '180402-REC2-T10-1e2To20e3.wav';
    %audiowrite(filename,REC2,Fs);
    
    
mytoc; 



%% THE END  THE END  THE END  THE END  THE END  THE END  THE END  THE END 