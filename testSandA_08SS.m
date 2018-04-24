%% testSandA_xx.m
%
% test S and A by convolving an audio signal and
% listening using headphones. 
% S.G. Tanyer, 180412 Victoria
%

%% HISTORY
%
%   
%
%
%               ..  Use their HRIR's here.
%             08SS  This is ircam's testnoise and our HRIR's. Did not work
%                   or we cannot understand clearly.
%...08  180417  ..  Use ircam's testnoise_02.m function's test signal.
%...07SS            Plotting version. Input S/A all mixed up for better 
%                   hearing, but there is some mistake to be found.              
%00006SS   kindaOK  Colter found one mistake. I was clipping the signal. It
%                   had a gain of 100! Now spectrum is good but no I have
%                   direction.
%             05SS  Tried both, S and A sucks big time! check PULSE COMPRESSION
%...05          ..  Try CubanSalsa to be sure.
%...04          ..  S and A are WRONG! I don't hear the direction. There
%                   are lots of artifacts, and why?? Simple convolution??
%...03          ..  earL/R and convL/R are calculated. now LISTEN!!
%...02          OK  I now have S A inpL and inpR for AmpESS. now CONVOLVE!!
%...01  180412  ..  Get S and A and some audio to start with.
%started using version process_HRIR_60

%...60          OK  Back to AmpESS to compare with CubanSalsa. XTC is bad
%                   for low and high pitch sounds. AmpESS sounds better
%                   since the input is frequency chirp. Salsa has so many
%                   components that it mushes in the time domain. So??
%...59          ..  Processing CubanSalsa-20secs
%...58        58SS  Writing to audio: 3D-AmpESS-FARINA-180411-1.wav
%               OK  Once again tried changing S/A and inpL/inpR, didnot help
%...57  180411  ..  Function does not sound?? Now play with plus and minuses.

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
    opengl('save', 'software');
    
%%
    Fs = 48e3;
    %Nfft = length(signS);
    %Nfft2 =  8 * 524288;
    Nfft2 = 4 * 480000;  %this is for AmpESS
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
    
    
    %signS = signA;
    
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
    %
    signS = signS ./ 70000;
    signA = signA ./ 70000;
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
    V=axis; axis([V(1) V(2) 1.1*V(3) 1.1*V(4)]);
    %V=axis; axis([0 5 1.1*V(3) 1.1*V(4)]);
    %plot(time2, cntrS-20, 'r-'), hold on; grid on;
    %plot(time2, cntrA-20, 'b-'), axis('tight');
    title('Raw HRTF, S and A');
    xlabel('Time (msecs)');
        subplot(312),
    plot(freq2, 20 .* log10(abs(fsignS)), 'r-'); hold on;
    plot(freq2, 20 .* log10(abs(fsignA)), 'b-');
    %semilogx(freq2, 20 .* log10(abs(fcntrS))-40, 'r-'); hold on;
    %semilogx(freq2, 20 .* log10(abs(fcntrA))-40, 'b-');
    set(gca,'xtick',[Fmn  1000 10000 Fmx]); grid on;
    axis('tight'); V=axis; axis([Fmn Fmx  V(3) V(4)]); grid on;
    %axis('tight'); V=axis; axis([Fmn Fmx  V(3) V(4)]); grid on;
    title ('Raw HRTF, S and A');
    xlabel('Frequency (Hertz)')
    ylabel('Desibels');
        subplot(313), 
    plot(freq2, phase_fsignS,'r-'); hold on;
    plot(freq2, phase_fsignA,'b-'); grid on;
    %semilogx(freq2, phase_fcntrS-1000,'r-'); hold on;
    %semilogx(freq2, phase_fcntrA-1000,'b-'); grid on;
    axis('tight'); 
    V=axis; axis([Fmn Fmx  -4000 0]); grid on;
    title ('Phase Raw HRTF, S and A)');
    xlabel('Frequency (Hertz)');
    ylabel('Phase (degrees)');

    
    
    
%%  HERE IS ircam's testnoise_02.m will be used
    %filename = 'Cuban-ORIG-Stereo-10.wav';
    Fs = 48e3; duration_sec = 2;
    [dum] = test_noise_02(duration_sec,Fs);
    dum = dum';
    AmpESS = dum; clear dum;

    
    
%%  HERE IS CUBANSALSA as your Input signal      
%   Now load an audio file and make it MONO. Then, play on one side only.
    %filename = 'Cuban-ORIG-Stereo-10.wav';
    %[dum,Fs] = audioread(filename);
    %dum = dum';
    %AmpESS = dum; clear dum;

%   Tdur = 10;
%    output = dum(:, [25*Fs : (25+Tdur)*Fs]);
    %filename = 'Cuban-ORIG-Stereo-10.wav';
    %audiowrite(filename,output',Fs);
    
    %output(1,:) = zeros(1, length(output));
    %filename = 'Cuban-ORIG-Left-10.wav';
    %audiowrite(filename,output',Fs);
    
    %write STEREO channel only
%    inpL = output(1, :);  %THIS IS ALL ZEROS.
%    inpR = output(2, :);
%    [inpL, time2] = eqtflength(inpL,time2);  
%    [inpR, time2] = eqtflength(inpR,time2);  
    
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
    
    %% HERE IS AmpESS as your Input signal      
    %filename = 'AmpESS-180327-T10-1e2To20e3.wav';
    %[AmpESS, Fs]=audioread(filename);
    %AmpESS = AmpESS';
    
    
    %
    %AmpESS(1,:) = AmpESS(1,:) + AmpESS(2,:);
    %AmpESS = AmpESS(1,:);
    [AmpESS,time2] = eqtflength(AmpESS,time2);  
    
    
    %
    %Right=0
    %inpL = AmpESS;
    %inpR = zeros(size(inpL));
    %Left=0
    inpR = AmpESS;
    inpL = zeros(size(inpR));
    
    
    
    [inpL, time2] = eqtflength(inpL,time2);  
    [inpR, time2] = eqtflength(inpR,time2);  
    %finpR = fft(inpR);
    %finpL = fft(inpL);
    %
    %Right=0
    %inpL = AmpESS;
    %inpR = zeros(size(inpL));
    finpL = fft(inpL, Nfft2);
    finpR = fft(inpR, Nfft2);
%2
%    figure(2), clf, hold off;
%        subplot(211)
%    plot(time2, inpR, 'r-'); hold on;        
%    plot(time2, inpL, 'b-'); hold on;        
%    axis('tight');V=axis; grid on; 
%    title('Input Audio File');
%    xlabel('Time (seconds)');
%        subplot(212)
%    semilogx(freq2, 20 .* log10(abs(finpR)), 'r-'); hold on;        
%    semilogx(freq2, 20 .* log10(abs(finpL))-50, 'b-'); hold on;        
%    %axis('tight');V=axis; axis([Fmn Fmx  40 50]); grid on;
%    axis('tight');V=axis; axis([Fmn Fmx  V(3) V(4)]); grid on;
%    title('Input Audio File');
%    xlabel('Frequency (Hertz)')
%    ylabel('Desibels');
    
    
    
%%  Calculate In-ear Signals
    fearL = fsignS .* finpL   + fsignA .* finpR;
    fearR = fsignS .* finpR   + fsignA .* finpL;
    
    earL = ifft(fearL);
    earR = ifft(fearR);
        
    convL = conv(signS, inpL)   + conv(signA, inpR);
    convR = conv(signA, inpL)   + conv(signS, inpR);
    
    convL = convL([1:Nfft2]);
    convR = convR([1:Nfft2]);
    
    fconvL = fft(convL);
    fconvR = fft(convR);
    
    output(2,:) = earL(1,:);
    output(1,:) = earR(1,:);
    
    
    Tdur = 10;
    Ndur = Tdur*Fs;
    output = 0.9 .* output ./ max(max(output));
    output = output(:,[1:Ndur]);
    %output(1,:) = earL(1,:);
    %output(2,:) = earR(1,:);
    %output(1,:) = inpL(1,:);
    %output(2,:) = inpR(1,:);
%%  2
    figure(2), clf, hold off;
        subplot(311)
    semilogx(freq2, 20 .* log10(abs(finpL)), 'r-'); hold on;        
    semilogx(freq2, 20 .* log10(abs(finpR)), 'b-'); hold on;grid on;     
    %axis('tight');V=axis; axis([Fmn Fmx  V(3) V(4)]); V=axis;grid on;
    title('Input Audio File');
    xlabel('Frequency (Hertz)')
    ylabel('Desibels');
    axis('tight');axis([Fmn Fmx -50 100]);
        subplot(312)
    semilogx(freq2, 20 .* log10(abs(fearL)), 'r-'); hold on;grid on;
    semilogx(freq2, 20 .* log10(abs(fearR)), 'b-'); hold on;grid on;
    xlabel('Frequency (Hertz)');
    title('Headphone signal (convolution-in-freq with S and A)');
    axis('tight');axis([Fmn Fmx -50 100]);
        subplot(313)
    semilogx(freq2, 20 .* log10(abs(finpR)), 'y-'); hold on;        
    semilogx(freq2, 20 .* log10(abs(fconvL)), 'r-'); hold on;grid on;
    semilogx(freq2, 20 .* log10(abs(fconvR)), 'b-'); hold on;grid on;
    xlabel('Frequency (Hertz)');
    title('Input(yellow) and The Headphone signal');
    axis('tight');axis([Fmn Fmx -50 100]);
    
    
    
    
%        subplot(311)
%    plot(output(1,:), 'r-'); hold on;        %Left
%    plot(output(2,:), 'b-'); hold on;        %Right
%    %axis('tight');V=axis;grid on;
%    %title('Output at the Left Speaker Terminal (RIGHTinp=0)');
%    xlabel('Time (seconds)');
%        subplot(312)
%    plot(convL, 'r-'); hold on;        
%    plot(convR, 'b-'); hold on;        
%    %axis('tight');V=axis;grid on;
%    %title('Output at the Left Speaker Terminal (RIGHTinp=0)');
%    xlabel('Time (seconds)');
%        subplot(313)
%    plot(freq2, 20 .* log10(abs(fearL)), 'r-'); hold on;grid on;
%    plot(freq2, 20 .* log10(abs(fearR)), 'b-'); hold on;grid on;
%    axis('tight');V=axis; axis([Fmn Fmx -10 100]);
%    xlabel('Frequency (Hertz)');

    
    
    p = audioplayer(output, Fs);
    play(p, [1, Ndur]);
    
    %% write audio files 
    pwd
    %filename = 'testSandA-180412-2.wav';
    %audiowrite(filename,output',Fs);
    
    %filename = '3D-AmpESS-FARINA-180411-2.wav';
    %filename = 'CubanSalsa-XTC-RightShifted.wav';
    %filename = 'CubanSalsa-XTC-StereoShifted.wav';
    %filename = 'CubanSalsa-XTC-Right.wav';
    %filename = 'CubanSalsa-XTC-Left.wav';
    %audiowrite(filename,output',Fs);
    
    
    
    
    
    
    
    
    
    
    
    dur
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
%%  FARINA'S CORRECTED SOLUTION  - CALCULATE FILTER NETWORK

    %first calculate your filters
    %SpA = signS + signA;
    %SmA = signS - signA;

    fSpA = fsignS + fsignA;
    fSmA = fsignS - fsignA;
  
    
    %determine eps
    delf = Fs/Nfft2;  %that is 0.0916 Hertz!
    %given Fmn and Fmx
    Nmn = floor(Fmn/delf);
    Nmx = ceil(Fmx/delf);
    
    small = 1e-16;
    eps = small .* ones(1, Nfft2);
    
    %supress the DC region of eps.
    Ndecayorder = 2;
    for iN = Nmn:-1:2
        delN = Nmn - iN;
        eps(iN) = eps(iN) + 1e-10 .* delN .^Ndecayorder;
        
        idum = Nfft2-Nmn +(Nmn-iN);
        eps(idum) = eps(idum) + 1e-10 .* delN .^Ndecayorder;
    end
    eps(1) = max(eps);

    %supress the High-Pass region of eps.
    Ndecayorder = 2;
    for iN = Nmx:Nfft2/2
        delN = iN - Nmx;
        eps(iN) = eps(iN) + 1e-13 .* delN .^Ndecayorder;

        dum = Nfft2/2 - Nmx;
        
        idum = (Nfft2/2+dum) -(iN-Nmx);
        eps(idum) = eps(idum) + 1e-13 .* delN .^Ndecayorder;
    end
    %eps = small .* ones(1, Nfft2);
    
    %figure(2), clf, hold off;
    %plot(freq2, eps, 'r+')
    %set(gca, 'xtick', [0 Fmn Fmx Fs/2 Fs]);
    %axis('tight'); %V=axis; axis([Fmn Fmx  V(3) V(4)]); grid on;
    %dur
    
    fC = fSpA .* fSmA;
    finvD = conj(fC) ./ (conj(fC) .* fC + eps);

    phase_fC    = rd2deg(unwrap(atan2(imag(fC   ), real(fC   ))));
    phase_finvD = rd2deg(unwrap(atan2(imag(finvD), real(finvD))));
    
    
    
%%  2   PLOT FARINA'S  C  invD
    figure(2), clf, hold off;
        subplot(311)
    %plot(time2.*1e3, SpA, 'r-'); hold on; grid on;
    %plot(time2.*1e3, SmA, 'b-'); hold on;    
    V=axis; axis([0 5 1.1*V(3) 1.1*V(4)]);
    title('Raw HRTF, S+A and S-A');
    xlabel('Time (msecs)');grid on
    set(gca,'xtick',[Fmn  1000 10000 Fmx]); grid on;
        subplot(312),
    %plot(freq2, 20 .* log10(abs(fC)), 'r-'); hold on;
    %plot(freq2, 20 .* log10(abs(finvD)), 'b-'); hold on;
    %semilogx(freq2, 20 .* log10(abs(conj(fC) .* fC)), 'r-'); hold on;
    %semilogx(freq2, 20 .* log10(abs(eps)), 'b-'); hold on;
    plot(freq2, 20 .* log10(abs(fSpA)), 'r-'); hold on;
    plot(freq2, 20 .* log10(abs(fSmA)), 'b-'); hold on;
    set(gca,'xtick',[10  Fmn 1000 10000 Fmx Fs]); grid on;
    axis('tight'); V=axis; axis([10 Fs  -380 120]); grid on;
    set(gca,'xtick',[Fmn  1000 10000 Fmx Fs/2 Fs]); grid on;
    title ('Raw fC and finvD');
    xlabel('Frequency (Hertz)')
    ylabel('Desibels');
%        subplot(313), 
%    semilogx(freq2, phase_fC   ,'r-'); hold on; grid on;
%    semilogx(freq2, phase_finvD,'b-'); hold on;
%    set(gca,'xtick',[Fmn  1000 10000 Fmx]); grid on;
%    axis('tight'); V=axis; axis([Fmn Fmx  V(3) V(4)]); grid on;
%    title ('Phase Raw HRTF, fC and finvD)');
%    xlabel('Frequency (Hertz)');
%    ylabel('Phase (degrees)');
    %drawnow; pause(1)
       
    
    
 
    
    dur
    
    
        
    %%  THIS IS THE REAL THING, NOT SIMULATION OF SOUND PROPAGATION
    
    %% FARINA'S CORRECTED METHOD
    fSbar = fsignS .* finvD;
    fAbar = fsignA .* finvD;
    
    fLprime = fSbar .* finpL - fAbar .* finpR;
    fRprime = fSbar .* finpR - fAbar .* finpL;
    
    fLprime2 = fsignS .* fLprime + fsignA .* fRprime;
    fRprime2 = fsignS .* fRprime + fsignA .* fLprime;
    
    fearL = fLprime2;
    fearR = fRprime2;
    
    %earL = (ifft(fearL));   %NOT anymore if eps is used
    %earR = (ifft(fearR));   %NOT anymore if eps is used
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

%%  3    PLOT  inpR/L  R/Lprime
    figure(3), clf, hold off;
        subplot(411)
    plot(time2, inpR, 'r-'); hold on;        
    plot(time2, inpL, 'b-'); hold on;        
    axis('tight');V=axis; grid on; 
    title('Input Audio File');
    xlabel('Time (seconds)');
        subplot(412)
    semilogx(freq2, 20 .* log10(abs(finpR)), 'r-'); hold on;        
    semilogx(freq2, 20 .* log10(abs(finpL))-50, 'b-'); hold on;        
    axis('tight');V=axis; axis([Fmn Fmx  40 50]); grid on;
    %axis('tight');V=axis; axis([Fmn Fmx  V(3) V(4)]); grid on;
    title('Input Audio File');
    xlabel('Frequency (Hertz)')
    ylabel('Desibels');
        subplot(413)
    semilogx(freq2, 20 .* log10(abs(fRprime)), 'r-'); hold on;        
    semilogx(freq2, 20 .* log10(abs(fLprime))-50, 'b-'); hold on;        
    %semilogx(freq2, 20 .* log10(abs(fearL)), 'r-'); hold on;        
    %semilogx(freq2, 20 .* log10(abs(fearR))-1, 'b-'); hold on;        
     axis('tight');V=axis; axis([Fmn Fmx  -300 85]); grid on;
    %axis('tight');V=axis; axis([Fmn Fmx  V(3) V(4)]); grid on;
    set(gca, 'Ytick', [-300 -247.6 -200 -100 0 47.39 ]);
    title('Emitted Sound Signal');
    xlabel('Frequency (Hertz)')
    ylabel('Desibels');
        subplot(414)
    plot(time2, earR, 'r-'); hold on;        
    plot(time2, earL, 'b-'); hold on;        
    axis('tight');V=axis; %axis([Fmn Fmx  V(3) V(4)]); grid on;
    title('In-Ear Sound Signal');
    xlabel('Time (seconds)'); grid on; 
    
   
    
    
    %% FARINA'S CORRECTED METHOD - SPEAKER OUTPUT CALCULATIONS
    %this is the speaker output audio for the listening check
    %Lprime = (ifft(fLprime));
    %Rprime = (ifft(fRprime));
    Lprime = real(ifft(fLprime));
    Rprime = real(ifft(fRprime));
    
    %Lprime3 = conv( Lprime, [1 0 -1]); Low and High reject filter.
    %Lprime(:) = Lprime3([2:Nfft2+1]);
    %dur
    
    %sm = 0.2; 
    sm = 15e-3;  %this is the expected output level. The larger is coming from the inverse of zero.
    [jj, kk] = find( abs(Lprime) > sm); Lprime(jj,kk) = 0;
    [jj, kk] = find( abs(Rprime) > sm); Rprime(jj,kk) = 0;
    
    Maxx = max(max(Lprime), max(Rprime));
    Lprime = 0.99 .* Lprime ./ Maxx;
    Rprime = 0.99 .* Rprime ./ Maxx;
%2    
    figure(2), clf, hold off;
    plot(time2, Lprime, 'r'), hold on;grid on;
    plot(time2, Rprime, 'b'), hold on;grid on;
    %plot(time2, imag(Lprime), 'r'), hold on;grid on;
    %plot(time2, imag(Rprime), 'b'), hold on;grid on;
    %plot(time2, imag(Lprime)./abs(Lprime), 'r'), hold on;grid on;
    %plot(time2, imag(Rprime)./abs(Rprime), 'b'), hold on;grid on;
    %axis('tight');V=axis; axis([2.5 7 -10e-3 10e-3]);
     axis('tight');V=axis; axis([0 10 V(3) V(4)]);
    title('Transmitted Sound Signal');
    xlabel('Time (seconds)'); grid on; 
    
    %DC and HighPass reject filtering
    %Lprime = real (conv(Lprime, [1 0 1]));
    %Rprime = real (conv(Rprime, [1 0 1]));
    %Lprime = Lprime(1,[2:Nfft2+1]);
    %Rprime = Rprime(1,[2:Nfft2+1]);
    
    %Lprime = real (conv(Lprime, [1 0 1]));
    %Rprime = real (conv(Rprime, [1 0 1]));
    %Lprime = Lprime(1,[2:Nfft2+1]);
    %Rprime = Rprime(1,[2:Nfft2+1]);

    %Lprime = real (conv(Lprime, [1 0 1]));
    %Rprime = real (conv(Rprime, [1 0 1]));
    %Lprime = Lprime(1,[2:Nfft2+1]);
    %Rprime = Rprime(1,[2:Nfft2+1]);

    
    
    output = zeros(2, Nfft2);
    output(2,:) = real( Lprime(1,:) );
    output(1,:) = real( Rprime(1,:) );
    output = 0.8 .* output./max(max(abs(output)));
    
    %DC and HighPass reject filtering
    %output = conv(output, [1 0 1])
    
    
    %Tdur = 20;
    output = output(:,[1:Tdur*Fs]);
 
    %KEVIN'S METHOD
    %this is the speaker output audio for the listening check
    %Lprime3 = ifft(fLprime3);
    %Rprime3 = ifft(fRprime3);
    %
    %Lprime3 = Lprime3(1, [96000:380e3]); %Nfft2]);
    %Rprime3 = Rprime3(1, [96000:380e3]); %Nfft2]);
    %
    %Lprime3 = 1 .* Lprime3 ./ max(max(Lprime3));
    %Rprime3 = 1 .* Rprime3 ./ max(max(Rprime3));
    %
    %output = zeros(2, length(Lprime3));
    %output(1,:) = real( earL(1,:) );
    %output(2,:) = real( earR(1,:) );
    %output = 0.999 .* output./max(max(abs(output)));
    %output = zeros(2, length(Lprime3));
    %
    %Output is generated using correct left and right order.
    %output(1,:) = real( Lprime3(1,:) );
    %output(2,:) = real( Rprime3(1,:) );
    %output = 0.999 .* output./max(max(abs(output)));
    %
    %Output is generated using interchanged left and right order.
    %output(2,:) = real( Lprime3(1,:) );
    %output(1,:) = real( Rprime3(1,:) );
    %output = 0.999 .* output./max(max(abs(output)));
    
    
    Ndur = Tdur*Fs;
    p = audioplayer(output, Fs);
    play(p, [1, Ndur]);
    %
    %p = myplayer1(output, Fs, Nfft2);
    
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
    %filename = '3D-Cuban-10-FARINA-180411-2.wav';
    %audiowrite(filename,output',Fs);
    
    %filename = '3D-AmpESS-FARINA-180411-2.wav';
    %filename = 'CubanSalsa-XTC-RightShifted.wav';
    %filename = 'CubanSalsa-XTC-StereoShifted.wav';
    %filename = 'CubanSalsa-XTC-Right.wav';
    %filename = 'CubanSalsa-XTC-Left.wav';
    %audiowrite(filename,output',Fs);
    
    



mytoc; 




%% THE END  THE END  THE END  THE END  THE END  THE END  THE END  THE END 