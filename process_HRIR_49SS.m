%% process_HRIR_xx.m
%
% Load recordings, all HRTFs
% divide by system response using long DFT solution.
% Second trial for the inverse filters.
% S.G. Tanyer, 180331 Victoria
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
%...50          ..  To the wild thing and play the audio.
%               SS  Farina's network works for the simulated case.
%...49          ..  Implement Farina's network.
%               OK  I am equating Nfft2 for ALL vectors.
%               ..  I am canceling CENTER calculations. 
%                   I will implement Farina's XLR solution. Gooooo...
%..48   180410  ..  47 was mixed up. I am jumping to 48 from 46SS.
%                   This is the version for the Report.
%             46SS  AmpESS, plots everything and save the audio signal
%...46          ..  Check AmpESS Lprime3. 
%               SS  AmpESS earL earR check version. Looking good.
%...45          OK  I now use AmpESS and earL earR is correct. Check Lprime3
%               OK  44 tests CubanSalsa to be working.
%...44          ..  Things got mixed up. This is 42SS and checked earL earR
%                   once again to find it is working for CubanSalsa.
%               xx  Audio files did not work. But angle of the monitor and
%                   the position of arms are very effective, not good.
%               ..  Going home to listen audio test files on my laptop!!
%             42SS  I also included LeftSwitched and RightSwitched.
%             42SS  All test audio files are recorded. LEFT, RIGHT, STEREO
%                   and also the StereoShifted where S and A are shifted in 
%                   case.
%...42          ..  Write different files and go home and listen on your laptop
%...41        41SS  PLOTS FOR THE REPORT.
%                   Input and Output: time and frequency domain. SUPER!
%               OK  Super. CubanSalsa single channel is XTCed to be ZEROO!
%...40          ..  Now try the audio file CubanSalsa
%             39SS  SUPERRRRR! It works with AmpESS too. Now try music!
%...39          ..  I will try with real audio. I am really excited to my
%                   bones!
%               SS  SUPER!! this version is better than 38SS. Here we have
%                   the input and the output on the same plot. Simply
%                   beautiful and even sexy when you put your full two
%                   weeks for a 65 hours/week tempo...
%               ..  Check the time domain signal.
%...38          SS  SUPERRR! I love it. I got some XTC working for the
%                   first time. I thank to Nick to just listen to me. I 
%                   figured it out when I was telling to him. 
%               ..  S A, CENTER read and FFTed.
%...37  180406  ..  Do every calculation in the frequency domain!
%...36          xx  (S+A)/(S+A) and (S-A)/(S-A) do not cancel each other in
%                   time domain. This is due to taking inverse two consecutive 
%                   times. But I hope that doing everything in the frequency 
%                   domain can solve this problem.
%               xx  It seems like inverses of 1/(S+A) and 1/(S-A) are incorrect.
%...35          ..  XTC completed. Results are WRONG! Lprime = Rprime??
%...34          OK  Now we are reading AmpESS.       
%...33          OK  HRTFs are artificially normalized by a constant.
%                   Since HRTFs are not normalized, their inverses are very
%                   very small. So, artificially divide S and A for norm.
%               ??  Note: the power of the inverse signal is 1/334.5 !!
%...32          ..  Now design the XTC filters.
%...31              So I return back to unnormalized plots.
%               xx  I believe NORMALIZATION to CENTER did not work. The
%               time domain invS invA are spread in time due to instability
%...30          ..  Now try to normalize S and A using CENTER recordings.
%...29          SS  This is the INVERSE FILTER DESIGN using all my FUNCTIONs.
%               OK  now delete the old myinvir lines.
%                   So write a separate FUNCTION file for MYINVIR().
%               OK  This version compares the Function myinvir to be
%                   correctly giving the same result as before.
%...28          ..  Now use FUNCTION for the inverse filtering.
%               !!  I found an error on Test1 and Test2. Phase of the
%                   inverse filter was not the one corresponding to the
%                   minimum phase calculations. Now, we observe EXCESS
%                   PHASE TERM as the error.
%...27          OK  Function MYMINPHASEIR() works.
%...26          OK  Function MYCLIPDB() is also working now. 
%...25          OK  Function MYFOLD() is now working.
%               OK  All pulse returns give good invS and invA results.
%...24  180405  ..  Now massage your code with writing all reusable functions.
%               SS  SUPER!! invS invA are calculated to be correct!!
%               ..  Check S*invS in frequency domain.
%               ..  You need an inverse filter FUNCTION.
%               ..  Test1 and Test2 means what??
%...23          ..  Now calculate the XTC filters. So exciting, what ifs!
%               OK  minS/A and invS/A are found. S*invS tested, not bad??
%...22          ..  Check minS.conv.invS
%               OK  Now you have the inverse IR and TFs. Calculate 1/(S-A)
%               kinda OK. Last FFT is 2*Nfft. Correct it to Nfft.
%...21          ..  Now do the wild thing and do a minimum phase inversion
%               OK  I have HRTFs and their inverses on the plots. Ready.
%...20  180404  ..  Gooooo
%Fresh start for inverse filtering based on Edgar J. Berdahl's work.

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

%NICE CODE LINES TO BE RECYCLED============================================

% Equate the lengths of two vectors.
%   [b,a] = eqtflength(b,a);  

% Get the transfer function from the poles and zeros. 
%   [b, a] = zp2tf(z_minP,p,km);

% Z plot of a filter
%   zplane(b,a);

% Frequency response of a filter
%   freqz(b,a);

% Given vector r, put/flip it into (1,n) vector form
%   [m,n] = size(r);
%   if m*n ~= m+n-1
%     error('fold.m: input must be a vector'); 
%   end
%   flipped = 0;
%   if (m > n)
%     n = m;
%     r = r.';
%     flipped = 1; %vector of size (m,1) is flipped to size(1,m)
%   end

% Find the smallest of the larger 2^n
%function [n,Nout] = myfindnfft_01(Nin)
%    n = ceil( log(Nin)/log(2));
%    Nout = 2^n;






%%  INTRO: Initialization 
    pause(3);
    clear;clc;
    myaddpath_04;
    mytic;
    figure(1), clf; 
    figure(2), clf; 
    figure(3), clf; 
    drawnow; pause(1)
    %
    opengl('save', 'software');
    
    
    
%%
    Fs = 48e3;
    %Nfft = length(signS);
    Nfft2 =  524288;
    %Nfft2 = 480000;  %this is for AmpESS
    %Nfft2 = 2640000; %this is for CubanSalsa
    %time  = [0:Nfft-1]./Fs;
    time2 = [0:Nfft2-1]./Fs;
    freq  = [0:Nfft2-1]./(Nfft2-1).*Fs;
    freq2 = [0:Nfft2-1]./Nfft2.*Fs;
    Fmn = 100;
    Fmx = 20e3;
    
%%  NOW LOAD HTRFs: 
%   in-ear micro-mic Knowles recording for the right ear
    filename = '180331-HRTF-signS-T10-1e2To20e3-p16.mat';
    load(filename,         'signS');
    filename = '180331-HRTF-signA-T10-1e2To20e3-p16.mat';
    load(filename,         'signA');
    %
    %filename = '180331-HRTF-cntrS-T10-1e2To20e3-p16.mat';
    %load(filename,         'cntrS');
    %filename = '180331-HRTF-cntrA-T10-1e2To20e3-p16.mat';
    %load(filename,         'cntrA');
    
    [signS,time2] = eqtflength(signS,time2);  
    [signA,time2] = eqtflength(signA,time2);  
    
    %switch S and A here in case it was already switched in error.
    %dum = signA;
    %signA = signS;
    %signS = dum;
    
    signS = signS ./ 100;
    signA = signA ./ 100;
    %cntrS = cntrS ./ 100;
    %cntrA = cntrA ./ 100;
    
    
    
%%  Calculate raw HRTF's; These are just the compressed raw'sign'al returns.
    fsignS = fft(signS, Nfft2);
    fsignA = fft(signA, Nfft2);
    %fcntrS = fft(cntrS, Nfft2);
    %fcntrA = fft(cntrA, Nfft2);
%% Calculate raw HRTF phase responses
    phase_fsignS = rd2deg(unwrap(atan2(imag(fsignS), real(fsignS))));
    phase_fsignA = rd2deg(unwrap(atan2(imag(fsignA), real(fsignA))));
    %phase_fcntrS = rd2deg(unwrap(atan2(imag(fcntrS), real(fcntrS))));
    %phase_fcntrA = rd2deg(unwrap(atan2(imag(fcntrA), real(fcntrA))));
%%  1
    figure(1), clf, hold off;
        subplot(311), 
    plot(time2.*1e3, signS, 'r-'), hold on; grid on;
    plot(time2.*1e3, signA, 'b-'), axis('tight'); 
    V=axis; axis([0 5 1.1*V(3) 1.1*V(4)]);
    %plot(time2, cntrS-20, 'r-'), hold on; grid on;
    %plot(time2, cntrA-20, 'b-'), axis('tight');
    title('Raw HRTF, S and A');
    xlabel('Time (msecs)');
        subplot(312),
    semilogx(freq2, 20 .* log10(abs(fsignS)), 'r-'); hold on;
    semilogx(freq2, 20 .* log10(abs(fsignA)), 'b-');
    %semilogx(freq2, 20 .* log10(abs(fcntrS))-40, 'r-'); hold on;
    %semilogx(freq2, 20 .* log10(abs(fcntrA))-40, 'b-');
    set(gca,'xtick',[Fmn  1000 10000 Fmx]); grid on;
    axis('tight'); V=axis; axis([Fmn Fmx  V(3) V(4)]); grid on;
    title ('Raw HRTF, S and A');
    xlabel('Frequency (Hertz)')
    ylabel('Desibels');
        subplot(313), 
    semilogx(freq2, phase_fsignS,'r-'); hold on;
    semilogx(freq2, phase_fsignA,'b-'); grid on;
    %semilogx(freq2, phase_fcntrS-1000,'r-'); hold on;
    %semilogx(freq2, phase_fcntrA-1000,'b-'); grid on;
    axis('tight'); 
    V=axis; axis([Fmn Fmx  -4000 0]); grid on;
    title ('Phase Raw HRTF, S and A)');
    xlabel('Frequency (Hertz)');
    ylabel('Phase (degrees)');

    
%%  Now load an audio file and make it MONO. Then, play on one side only.
    %filename = 'Cuban55-48K.wav';
    %[dum,Fs] = audioread(filename);
    %dum = dum';

    %write STEREO channel only
    %inpL = dum(1, :);
    %inpR = dum(2, :);
    %finpL = fft(inpL);
    %finpR = fft(inpR);
    
    %write RIGHT channel only
    %inpR = dum(2, :);
    %inpL = zeros(size(inpR));
    %finpR = fft(inpR);
    %finpL = zeros(size(finpR));
    
    %write LEFT channel only
    %inpL = dum(1, :);
    %inpR = zeros(size(inpL));
    %finpL = fft(inpL);
    %finpR = zeros(size(finpL));
    
    %then load your test signal      
    filename = 'AmpESS-180327-T10-1e2To20e3.wav';
    [AmpESS, Fs]=audioread(filename);
    AmpESS = AmpESS';
    
    [AmpESS,time2] = eqtflength(AmpESS,time2);  
    
    
    %Left=0
    %inpR = AmpESS;
    %inpL = zeros(size(inpR));
    %finpR = fft(inpR);
    %finpL = fft(inpL);
    
    %Right=0
    
    inpL = AmpESS;
    inpR = zeros(size(inpL));
    finpL = fft(inpL, Nfft2);
    finpR = fft(inpR, Nfft2);

%%  FARINA'S CORRECTED NETWORK

    %first calculate your filters
    SpA = signS + signA;
    SmA = signS - signA;

    fSpA = fsignS + fsignA;
    fSmA = fsignS - fsignA;

    fC = fSpA .* fSmA;
    eps = 0.0001;
    finvD = conj(fC) ./ (conj(fC) .* fC + eps);

    phase_fC    = rd2deg(unwrap(atan2(imag(fC   ), real(fC   ))));
    phase_finvD = rd2deg(unwrap(atan2(imag(finvD), real(finvD))));
    

    
%%  2
    figure(2), clf, hold off;
        subplot(311)
    plot(time2.*1e3, SpA, 'r-'); hold on; grid on;
    plot(time2.*1e3, SmA, 'b-'); hold on;    
    V=axis; axis([0 5 1.1*V(3) 1.1*V(4)]);
    title('Raw HRTF, S+A and S-A');
    xlabel('Time (msecs)');
    set(gca,'xtick',[Fmn  1000 10000 Fmx]); grid on;
        subplot(312),
    semilogx(freq2, 20 .* log10(abs(fC)), 'r-'); hold on;
    semilogx(freq2, 20 .* log10(abs(finvD)), 'b-'); hold on;
    set(gca,'xtick',[Fmn  1000 10000 Fmx]); grid on;
    axis('tight'); V=axis; axis([Fmn Fmx  V(3) V(4)]); grid on;
    set(gca,'xtick',[Fmn  1000 10000 Fmx]); grid on;
    title ('Raw fC and finvD');
    xlabel('Frequency (Hertz)')
    ylabel('Desibels');
        subplot(313), 
    semilogx(freq2, phase_fC   ,'r-'); hold on; grid on;
    semilogx(freq2, phase_finvD,'b-'); hold on;
    set(gca,'xtick',[Fmn  1000 10000 Fmx]); grid on;
    axis('tight'); V=axis; axis([Fmn Fmx  V(3) V(4)]); grid on;
    title ('Phase Raw HRTF, fC and finvD)');
    xlabel('Frequency (Hertz)');
    ylabel('Phase (degrees)');

        
    
    %%
    fSbar = fsignS .* finvD;
    fAbar = fsignA .* finvD;
    
    fLprime = fSbar .* finpL - fAbar .* finpR;
    fRprime = fSbar .* finpR - fAbar .* finpL;
    
    fLprime2 = fsignS .* fLprime + fsignA .* fRprime;
    fRprime2 = fsignS .* fRprime + fsignA .* fLprime;
    
    fearL = fLprime2;
    fearR = fRprime2;
    
    earL = real(ifft(fearL));
    earR = real(ifft(fearR));
    
    
%%  KEVIN'S METHOD:
%   Design XTF filters!
%   CORRECTION: none of the inversion and minimum phase operation will be
%   done here. Every calculation will be done in the frequency domain. Time
%   domain signal will be obtained at the very end!
%
%   Up to here we have S, A, Cntr and their FFTs.
%   Now calculate the OVERALL filter for the right and left audio channels
%   for the complete line of process; audio file to in-ear point!
%   I really hope that this time it works...
%
%    fLprime1 = finpL - finpR;
%    fRprime1 = finpL + finpR;
%        
%    fLprime2 = fLprime1 ./ fSmA;
%    fRprime2 = fRprime1 ./ fSpA;
%    
%    %this is going to be the output of speakers!
%    fLprime3 =  fLprime2 + fRprime2;
%    fRprime3 = -fLprime2 + fRprime2;
%        
%    %this part is the simulation of propagation of sound into ears!
%    fearL = fLprime3 .* fsignS + fRprime3 .* fsignA;
%    fearR = fLprime3 .* fsignA + fRprime3 .* fsignS;
%    fearL = fearL ./ 2;
%    fearR = fearR ./ 2;
%    earL = ifft(fearL);
%    earR = ifft(fearR);

%%  3    
    figure(3), clf, hold off;
        subplot(411)
    plot(time2, inpL, 'r-'); hold on;        
    plot(time2, inpR, 'b-'); hold on;        
    axis('tight');V=axis; grid on; 
    title('Input Audio File (RIGHTinp=0)');
    xlabel('Time (seconds)');
        subplot(412)
    semilogx(freq2, 20 .* log10(abs(finpL)), 'r-'); hold on;        
    semilogx(freq2, 20 .* log10(abs(finpR))-1, 'b-'); hold on;        
    axis('tight');V=axis; axis([Fmn Fmx  40 50]); grid on;
    %axis('tight');V=axis; axis([Fmn Fmx  V(3) V(4)]); grid on;
    title('Input Audio File (RIGHTinp=0)');
    xlabel('Frequency (Hertz)')
    ylabel('Desibels');
        subplot(413)
    semilogx(freq2, 20 .* log10(abs(fearL)), 'r-'); hold on;        
    semilogx(freq2, 20 .* log10(abs(fearR))-1, 'b-'); hold on;        
    axis('tight');V=axis; axis([Fmn Fmx  -300 85]); grid on;
    %axis('tight');V=axis; axis([Fmn Fmx  V(3) V(4)]); grid on;
    set(gca, 'Ytick', [-300 -247.6 -200 -100 0 47.39 ]);
    title('In-Ear Sound Signal (RIGHTinp=0)');
    xlabel('Frequency (Hertz)')
    ylabel('Desibels');
        subplot(414)
    plot(time2, earL, 'r-'); hold on;        
    plot(time2, earR, 'b-'); hold on;        
    axis('tight');V=axis; %axis([Fmn Fmx  V(3) V(4)]); grid on;
    title('In-Ear Sound Signal (RIGHTinp=0)');
    xlabel('Time (seconds)'); grid on; 
    
%%  THIS IS THE REAL THING, NOT SIMULATION OF SOUND PROPAGATION
    
    %this is the speaker output audio for the listening check
    %Lprime3 = ifft(fLprime3);
    %Rprime3 = ifft(fRprime3);
    
    %Lprime3 = Lprime3(1, [96000:380e3]); %Nfft2]);
    %Rprime3 = Rprime3(1, [96000:380e3]); %Nfft2]);
    
    %Lprime3 = 1 .* Lprime3 ./ max(max(Lprime3));
    %Rprime3 = 1 .* Rprime3 ./ max(max(Rprime3));
    
    
    
    
    %thresh = 0.2;
    %clipped = output;
    %as = abs(output);
    %mas = max(as(:));
    %toobig = find(as > thresh);
    %clipped(toobig) = thresh;
    %output = clipped;
    
    
    
    %output = zeros(2, length(Lprime3));
    %output(1,:) = real( earL(1,:) );
    %output(2,:) = real( earR(1,:) );
    %output = 0.999 .* output./max(max(abs(output)));
    %output = zeros(2, length(Lprime3));
    
    %Output is generated using correct left and right order.
    %output(1,:) = real( Lprime3(1,:) );
    %output(2,:) = real( Rprime3(1,:) );
    %output = 0.999 .* output./max(max(abs(output)));
    
    %Output is generated using interchanged left and right order.
    %output(2,:) = real( Lprime3(1,:) );
    %output(1,:) = real( Rprime3(1,:) );
    %output = 0.999 .* output./max(max(abs(output)));
    
    %p = audioplayer(output, Fs);
    %play(p, [1 6*Fs]);
%%  3
    %figure(3), clf, hold off;
    %    subplot(211)
    %plot(real(Lprime3), 'r-'); hold on;        
    %plot(imag(Lprime3), 'b-'); hold on;        
    %axis('tight');V=axis;grid on;
    %title('Output at the Left Speaker Terminal (RIGHTinp=0)');
    %xlabel('Time (seconds)');
    %legend('real','imag');
    %xlabel('Time (seconds)');
    %    subplot(212)
    %plot( real(Rprime3), 'r-'); hold on;        
    %plot( imag(Rprime3), 'b-'); hold on;        
    %title('Output at the Right Speaker Terminal (RIGHTinp=0)');
    %axis(V);grid on;
    %xlabel('Time (seconds)');
    %legend('real','imag');
    
    
    %% write audio files 
    %pwd
    %filename = 'AmpESS-XTC-RightIsZero2.wav';
    %audiowrite(filename,output',Fs);
    
    %filename = 'CubanSalsa-XTC-RightShifted.wav';
    %filename = 'CubanSalsa-XTC-StereoShifted.wav';
    %filename = 'CubanSalsa-XTC-Right.wav';
    %filename = 'CubanSalsa-XTC-Left.wav';
    %audiowrite(filename,output',Fs);
    
    
    



    
    
mytoc; 




%% THE END  THE END  THE END  THE END  THE END  THE END  THE END  THE END 