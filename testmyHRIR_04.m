%% testmyHRIR_xx.m
%
% Load my HRIRs plot them with HRTF
% Then divide HRTF with the CENTER measurement and IFFT
% see if HRIRs improve or not. 
% S.G. Tanyer, 180419, Victoria

%% HISTORY
%
%
%               ..  I couldn't do it. My files and his struct don't match up.
%..3            ..  Use ircam's minphase/excessphase routines.                   
%               >>  Not working. CENTER normalization using FARINA is unstable?? 
%                   Should use minimum/excess phase method of ircam
%..2            ..  Divide by CENTER using Farina's method.               
%..1    180419  OK  Plots signA signS ffts and the CENTER channel.

%% NICE CODE LINES TO BE RECYCLED============================================

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
    %figure(3), clf; 
    drawnow; pause(1)
    %
    opengl('save', 'software');
    

    
%%  NOW LOAD HTRFs: 
%   in-ear micro-mic Knowles recording for the right ear
    filename = '180331-HRTF-signS-T10-1e2To20e3-p16.mat';
    load(filename,         'signS');
    filename = '180331-HRTF-signA-T10-1e2To20e3-p16.mat';
    load(filename,         'signA');
    %
    filename = '180331-HRTF-cntrS-T10-1e2To20e3-p16.mat';
    load(filename,         'cntrS');
    filename = '180331-HRTF-cntrA-T10-1e2To20e3-p16.mat';
    load(filename,         'cntrA');
    
    
    
%%
    Fs = 48e3;
    Nfft = length(signS);
    Nfft2 =  2 * 512;
    %Nfft2 =  8 * 524288;
    time  = [0:Nfft-1]./Fs;
    time2 = [0:Nfft2-1]./Fs;
    freq  = [0:Nfft-1]./Nfft .*Fs;
    freq2 = [0:Nfft2-1]./Nfft2.*Fs;
    Fmn = 100;
    Fmx = 20e3; 
    
   
%%   
    [signS,time2] = eqtflength(signS,time2);  
    [signA,time2] = eqtflength(signA,time2);  
    [cntrS,time2] = eqtflength(cntrS,time2);  
    [cntrA,time2] = eqtflength(cntrA,time2);  
    
    
    %switch S and A here in case it was already switched in error.
    %dum = signA;
    %signA = signS;
    %signS = dum;
    
    signS = signS ./ 1000;
    signA = signA ./ 1000;
    cntrS = cntrS ./ 1000;
    cntrA = cntrA ./ 1000;
   

    
%%  Calculate raw HRTF's; These are just the compressed raw'sign'al returns.
    fsignS = fft(signS, Nfft2);
    fsignA = fft(signA, Nfft2);
    fcntrS = fft(cntrS, Nfft2);
    fcntrA = fft(cntrA, Nfft2);
%% Calculate raw HRTF phase responses
    phase_fsignS = rd2deg(unwrap(atan2(imag(fsignS), real(fsignS))));
    phase_fsignA = rd2deg(unwrap(atan2(imag(fsignA), real(fsignA))));
    phase_fcntrS = rd2deg(unwrap(atan2(imag(fcntrS), real(fcntrS))));
    phase_fcntrA = rd2deg(unwrap(atan2(imag(fcntrA), real(fcntrA))));
%%  1
    figure(1), clf, hold off;
        subplot(311), 
    plot(time2, signS, 'r-'), hold on; grid on;
    plot(time2, signA, 'b-'), axis('tight'); 
    plot(time2, cntrS-20, 'r-'), hold on; grid on;
    plot(time2, cntrA-20, 'b-'), axis('tight');
    %V=axis; axis([0 5e-3 1.1*V(3) 1.1*V(4)]);
    title('Raw HRTF, S and A');
    xlabel('Time (msecs)');
        subplot(312),
    semilogx(freq2, 20 .* log10(abs(fsignS)), 'r-'); hold on;
    semilogx(freq2, 20 .* log10(abs(fsignA)), 'b-');
    semilogx(freq2, 20 .* log10(abs(fcntrS))-40, 'r-'); hold on;
    semilogx(freq2, 20 .* log10(abs(fcntrA))-40, 'b-');
    set(gca,'xtick',[Fmn  1000 10000 Fmx]); grid on;
    axis('tight'); V=axis; axis([Fmn Fmx  V(3) V(4)]); grid on;
    title ('Raw HRTF, S and A');
    xlabel('Frequency (Hertz)')
    ylabel('Desibels');
        subplot(313), 
    semilogx(freq2, phase_fsignS,'r-'); hold on;
    semilogx(freq2, phase_fsignA,'b-'); grid on;
    semilogx(freq2, phase_fcntrS-1000,'r-'); hold on;
    semilogx(freq2, phase_fcntrA-1000,'b-'); grid on;
    axis('tight'); 
    V=axis; axis([Fmn Fmx  -4000 0]); grid on;
    title ('Phase Raw HRTF, S and A)');
    xlabel('Frequency (Hertz)');
    ylabel('Phase (degrees)');

    
    
%% Equalization with CENTER signal using FARINA's method.
    meps = eps; %seem to work with uglgy 1e-1;
    finvcntrS = conj(fcntrS) ./ ( conj(fcntrS) .* fcntrS + meps);
    finvcntrA = conj(fcntrA) ./ ( conj(fcntrA) .* fcntrA + meps);
    
    %use ircam's routines here.
    %[magnitudes_S,excess_phases_S] = phase_decomposition(transfer_functions_S);
    
    [magS,xssS] = phase_decomposition1(finvcntrS);
    [magA,xssA] = phase_decomposition1(finvcntrA);
    
    dur
    
    
    fnormS = fsignS .* finvcntrS;
    fnormA = fsignA .* finvcntrA;
    
    normS = real(ifft(fnormS));
    normA = real(ifft(fnormA));
    
    phase_fnormS = rd2deg(unwrap(angle(fnormS))); %rd2deg(unwrap(atan2(imag(fnormS), real(fnormS))));
    phase_fnormA = rd2deg(unwrap(angle(fnormA))); %rd2deg(unwrap(atan2(imag(fnormA), real(fnormA))));

    phase_finvcntrS = rd2deg(unwrap(angle(finvcntrS)));
    phase_finvcntrA = rd2deg(unwrap(angle(finvcntrA)));
    %phase_finvcntrS = rd2deg(unwrap(atan2(imag(finvcntrS), real(finvcntrS))));
    %phase_finvcntrA = rd2deg(unwrap(atan2(imag(finvcntrA), real(finvcntrA))));
% 2
    figure(2), clf, hold off;
        subplot(311), 
    plot(time2, normS, 'r-'), hold on; grid on;
    plot(time2, normA, 'b-'), axis('tight'); 
    %V=axis; axis([0 5e-3 1.1*V(3) 1.1*V(4)]);
    title('Norm HRIR, S and A');
    xlabel('Time (msecs)');
        subplot(312),
    %semilogx(freq2, 20 .* log10(abs(finvcntrS)), 'r-'); hold on;
    %semilogx(freq2, 20 .* log10(abs(finvcntrA)), 'b-');
    semilogx(freq2, 20 .* log10(abs(fnormS)), 'r-'); hold on;
    semilogx(freq2, 20 .* log10(abs(fnormA)), 'b-');
    set(gca,'xtick',[Fmn  1000 10000 Fmx]); grid on;
    axis('tight'); V=axis; axis([Fmn Fmx  V(3) V(4)]); grid on;
    title ('Norm HRTF, S and A');
    xlabel('Frequency (Hertz)')
    ylabel('Desibels');
        subplot(313), 
    semilogx(freq2, phase_fsignS,'g-'); hold on;
    semilogx(freq2, phase_fsignA,'g-'); grid on;
    semilogx(freq2, phase_fcntrS,'r-'); hold on;
    semilogx(freq2, phase_fcntrA,'b-'); grid on;
    %semilogx(freq2, phase_fnormS,'r-'); hold on;
    %semilogx(freq2, phase_fnormA,'b-'); grid on;
    semilogx(freq2, phase_finvcntrS,'k-'); hold on;
    semilogx(freq2, phase_finvcntrA,'k-'); grid on;
    
    axis('tight'); 
    %V=axis; axis([Fmn Fmx  -4000 0]); grid on;
    title ('Phase Norm HRTF, S and A)');
    xlabel('Frequency (Hertz)');
    ylabel('Phase (degrees)');
    
    
    
    
    dur
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    dur
    
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
    %plot(freq2, 20 .* log10(abs(fC)), 'r-'); hold on;
    %plot(freq2, 20 .* log10(abs(finvD)), 'b-'); hold on;
    %semilogx(freq2, 20 .* log10(abs(conj(fC) .* fC)), 'r-'); hold on;
    %semilogx(freq2, 20 .* log10(abs(eps)), 'b-'); hold on;
    semilogx(freq2, 20 .* log10(abs(fC)), 'r-'); hold on;
    semilogx(freq2, 20 .* log10(abs(finvD)), 'b-'); hold on;
    set(gca,'xtick',[10  Fmn 1000 10000 Fmx Fs]); grid on;
    axis('tight'); V=axis; axis([10 Fs  -380 120]); grid on;
    set(gca,'xtick',[Fmn  1000 10000 Fmx Fs/2 Fs]); grid on;
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
    %drawnow; pause(1)
        
    
    %%
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
    semilogx(freq2, 20 .* log10(abs(fLprime)), 'r-'); hold on;        
    semilogx(freq2, 20 .* log10(abs(fRprime))-1, 'b-'); hold on;        
    %semilogx(freq2, 20 .* log10(abs(fearL)), 'r-'); hold on;        
    %semilogx(freq2, 20 .* log10(abs(fearR))-1, 'b-'); hold on;        
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
    
    %FARINA'S CORRECTED METHOD
    %this is the speaker output audio for the listening check
    %Lprime = (ifft(fLprime));
    %Rprime = (ifft(fRprime));
    Lprime = real(ifft(fLprime));
    Rprime = real(ifft(fRprime));
    
    %Lprime3 = conv( Lprime, [1 0 -1]); Low and High reject filter.
    %Lprime(:) = Lprime3([2:Nfft2+1]);
    %dur
    
    sm = 15e-3;
    [jj, kk] = find( abs(Lprime) > sm); Lprime(jj,kk) = 0;
    [jj, kk] = find( abs(Rprime) > sm); Rprime(jj,kk) = 0;
    
    Maxx = max(max(Lprime), max(Rprime));
    Lprime = 0.99 .* Lprime ./ Maxx;
    Rprime = 0.99 .* Rprime ./ Maxx;
    
    figure(2), clf, hold off;
    plot(time2, Lprime, 'r'), hold on;grid on;
    plot(time2, Rprime, 'b'), hold on;grid on;
    %plot(time2, imag(Lprime), 'r'), hold on;grid on;
    %plot(time2, imag(Rprime), 'b'), hold on;grid on;
    %plot(time2, imag(Lprime)./abs(Lprime), 'r'), hold on;grid on;
    %plot(time2, imag(Rprime)./abs(Rprime), 'b'), hold on;grid on;
    %axis('tight');V=axis; axis([2.5 7 -10e-3 10e-3]);
     axis('tight');V=axis; axis([0 10 V(3) V(4)]);
    
    
    
    output = zeros(2, Nfft2);
    output(2,:) = real( Lprime(1,:) );
    output(1,:) = real( Rprime(1,:) );
    %output = 0.999 .* output./max(max(abs(output)));
 
    dum = myplay(output,Fs,Nfft2);
    
    


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



%function [dum] = myplay(output,Fs,Nfft2)
    p = audioplayer(output, Fs);
    play(p, [1 Nfft2]);
    play(p, [1 10*Fs]);
    dum=0;
%end
    



%% THE END  THE END  THE END  THE END  THE END  THE END  THE END  THE END 