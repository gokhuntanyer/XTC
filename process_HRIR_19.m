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
%...19  180403  ..  Check the XTC function implementation. Clean up.
%               ..  Checked S.conv.invS and A.conv.invA etc. found nothing
%...18          ..  Work on Audio again. CubanSalsa LEFT.
%...17          ..  testREC4: Now RIGHT channel only. Does not seem to work.
%...16          ..  testREC3: Now this is the LEFT channel only.
%...15          ..  All dummy lines are converted to comment. Still seems
%                   like it is working similar to 07SS.
%               OK  Seem to be working like 07SS. Now comment the
%                   unnecessary lines following new part which is forgetten to
%                   be deleted.
%...14 180402   ..  This is the copy of 07SS.
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
    %filename = '180331-HRTF-cntrS-T10-1e2To20e3-p16.mat';
    %load(filename,         'cntrS');
    %filename = '180331-HRTF-cntrA-T10-1e2To20e3-p16.mat';
    %load(filename,         'cntrA');
    
    %in-ear micro-mic Knowles recording for the right ear
    filename = '180331-HRTF-signS-T10-1e2To20e3-p16.mat';
    load(filename,         'signS');
    filename = '180331-HRTF-signA-T10-1e2To20e3-p16.mat';
    load(filename,         'signA');
    
    
    %%
    Fs = 48e3;
    Nfft = length(signS);
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
    
    %fcntrS = fft(cntrS, Nfft2);
    %fcntrA = fft(cntrA, Nfft2);
    
    fsignS = fft(signS, Nfft2);
    fsignA = fft(signA, Nfft2);
    
    finvS = fft(invS, Nfft2);
    finvA = fft(invA, Nfft2);
    
    finvSpA = fft(invSpA, Nfft2);
    finvSmA = fft(invSmA, Nfft2);
    
    
    
%% Calculate raw HRTF phase responses
    phase_signS = rd2deg(unwrap(atan2(imag(fsignS), real(fsignS))));
    phase_signA = rd2deg(unwrap(atan2(imag(fsignA), real(fsignA))));
    
    phase_invS = rd2deg(unwrap(atan2(imag(finvS), real(finvS))));
    phase_invA = rd2deg(unwrap(atan2(imag(finvA), real(finvA))));
    
    phase_invSpA = rd2deg(unwrap(atan2(imag(finvSpA), real(finvSpA))));
    phase_invSmA = rd2deg(unwrap(atan2(imag(finvSmA), real(finvSmA))));
    
    
    %% PLOTTING
%1
    figure(1), %clf, hold off;
        subplot(311), 
    plot(time, signS, 'r-'), hold on;
    plot(time, signA, 'b-'), axis('tight');
    title('Compressed signal - SIGNAL');
    xlabel('Time (seconds)');
        subplot(312), 
    plot(time, invS, 'r-'), hold on;
    plot(time, invA, 'b-'), axis('tight');
    title('Compressed signal - SIGNAL');
    xlabel('Time (seconds)');
        subplot(313), 
    plot(time, invSpA, 'r-'), hold on;
    plot(time, invSmA, 'b-'), axis('tight');
    title('Compressed signal - SIGNAL');
    xlabel('Time (seconds)');

%%  ANALYZE S * invS 
    test1 = conv( (signS+signA), invSpA );
    test2 = conv( (signS-signA), invSmA );
    %test2 = conv( (signS+signA), invSmA);
    
    
    
    
    
    
    %psdSpA = conv(signS+signA, invSmA);
    %psdSmA = conv(signS-signA, invSpA);
    %lenpsd = length(psdS);
    %tpsd = [0:lenpsd-1]./Fs;
    %
    figure (2)
    %    subplot(211)
    plot(test1, 'r-'), hold on; grid on;
    plot(test2, 'b-')
    %plot(tpsd-0.01065, psdS, 'r-'), hold on;
    %plot(tpsd-0.01065, psdA, 'b-'), grid on;
    %    subplot(212)    
    %plot(tpsd-0.01065, psdSpA, 'r-'), hold on;
    %plot(tpsd-0.01065, psdSmA, 'b-'), grid on;
    
    dur
    
    %2    
    figure(2), %clf, hold off
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
%3 
    figure(3), clf; hold off;
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
    %close all;
    
%%  Now load an audio file and make it MONO. Then, play on one side only.
    %filename = 'Cuban55-48K.wav';
    %[dum,Fs] = audioread(filename);
    
    
%%  Now we have all the filters;  S, A, invS, invA, 1/(S+A) and 1/(S-A)  
    %180327-TXDataFS-T10-1e2To20e3
    filename = '180327-TXDataFS-T10-1e2To20e3.wav';
    [dum,Fs] = audioread(filename);
    
    
    len = length(dum);
    %L(:,1) = dum(:,1);
    %R(:,1) = zeros(size(L)); %dum(:,2);
    R(:,1) = dum(:,2);
    L(:,1) = zeros(size(R)); 
    clear dum;
   
   
    Lprime = conv( (L+R) , invSpA) + conv( (L-R) , invSmA);
    Rprime = conv( (R+L) , invSpA) + conv( (R-L) , invSmA);
    
    Ltest  = conv( Lprime, signS)  + conv( Rprime, signA);
    Rtest  = conv( Rprime, signS)  + conv( Lprime, signA);
    
    figure(2), clf, hold off;
        subplot(211), 
    plot(Ltest,'r-'), hold on; grid on;
    axis('tight'); V=axis;
        subplot(212), 
    plot(Rtest,'b-'), grid on;
    axis([V]); V;
    
    REC1(:,1) =  Lprime;
    REC1(:,2) = -Rprime;
    RECmx = max(max(REC1));
    REC1 = REC1 ./ RECmx;
    
    %normalize REC
    %REC2(:,2) = L;
    %REC2(:,1) = R;
    %RECmx = max(max(REC2));
    %REC2 = REC2 ./ RECmx;
    
    %% write audio files 
    %filename = '180402-testREC5-T10-1e2To20e3.wav';
    %audiowrite(filename,REC1,Fs);
    
    %filename = '180402-REC2-T10-1e2To20e3.wav';
    %audiowrite(filename,REC2,Fs);
    
    
mytoc; 



%% THE END  THE END  THE END  THE END  THE END  THE END  THE END  THE END 