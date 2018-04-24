%% Amplitude Modulated Sweep design
%
% S. Gokhun Tanyer
% 180302

%% HISTORY
%
%               SS  This is a much better version. Comparing two
%                   alternative waveforms; magn and phase responses.
%...09  180413  ..  Check the phase of fcor.
%...08  180314  SS  Clean up the code. Change to f1:20 f2:20K
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
    fprintf('CREATE: sweep DONE\n');
    %normalize st and ts to 0.9
    dum = max(max([st,-st,ts,-ts]));
    st = 0.999999 .* st ./ dum;  %not to overload.
    ts = 0.999999 .* ts ./ dum;
    
    
    
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
    
    %s = s -sum(s)./length(s);
    dum = max(max(s));
    s = s ./ dum;
    invs([1:LT]) = s([LT:-1:1]);
    len = length(s);
    
%Find the smallest of the larger 2^n
%function [n,Nout] = myfindnfft_01(Nin)
%    n = ceil( log(Nin)/log(2));
%    Nout = 2^n;
    [n,Nout] = myfindnfft_02(len);
    N = 2*Nout;
    time = [0:N-1]./Fs;
    freq = [0:N-1]./N.*Fs;
    
% Equate the lengths of two vectors.
%   [b,a] = eqtflength(b,a);      
    [ts  ,time] = eqtflength(ts  ,time);  
    [st  ,time] = eqtflength(st  ,time);  
    [s   ,time] = eqtflength(s   ,time);  
    [invs,time] = eqtflength(invs,time);  
    
%% CALCULATE FFT OF SWEEP
    fft_ts      = fft(ts);
    fft_st      = fft(st);
    fft_s       = fft(s );
    fft_invs    = fft(invs);
    
    phase_st   = 180./pi .* unwrap(atan2(imag(fft_st  ), real(fft_st  )));
    phase_ts   = 180./pi .* unwrap(atan2(imag(fft_ts  ), real(fft_ts  )));
    phase_s    = 180./pi .* unwrap(atan2(imag(fft_s   ), real(fft_s   )));
    phase_invs = 180./pi .* unwrap(atan2(imag(fft_invs), real(fft_invs)));
    
    fcor  = fft_ts .* fft_st;
    fcor2 = fft_s  .* fft_invs;
    
    cor  = ifft(fcor);
    cor2 = ifft(fcor2);
    
    Max  = max(cor);  cor  = cor ./Max; fcor  = fft(cor);
    Max2 = max(cor2); cor2 = cor2./Max2; fcor2 = fft(cor2);
    
    phase_fcor  = 180./pi .* unwrap(atan2(imag(fcor ), real(fcor )));
    phase_fcor2 = 180./pi .* unwrap(atan2(imag(fcor2), real(fcor2)));
    
    
    
%% PLOTTING SECTION
% 1    
    figure(1), clf, hold off;
    plot(time,st    , 'r-'), hold on;grid on;
    plot(time,ts-3  , 'b-'), hold on;grid on;
    plot(time,s-6   , 'k-'), hold on,grid on
    plot(time,invs-9, 'g-'), hold on;grid on;
    axis('tight'); 
    title('My Exponential Sine Sweep (ESS)');
    xlabel('Time (seconds)');
% 2    
    figure(2), clf, hold off;
        subplot(311),
    plot(time-10, 20 .* log10( abs(cor)),'r-'), hold on;grid on;
    plot(time-10, 20 .* log10( abs(cor)),'ro'), hold on;grid on;
    plot(time-10, 20 .* log10( abs(cor2))+14.57,'b-'), hold on;grid on;
    plot(time-10, 20 .* log10( abs(cor2))+14.57,'bo'), hold on;grid on;
    plot(time-10, zeros(size(cor2)), 'k-');
    axis('tight'); V=axis; axis([-0.001 0.001 -25 20]);
    title('Cross-Correlation Functions (st*ts and s*sinv)');
    ylabel('Desibels(dB)');
    xlabel('Time (seconds)');
    ylabel('Desibels(dB)');
        subplot(312),
    semilogx(freq, 20 .* log10( abs(fcor ))-13.11 ,'r');hold on;grid on;
    semilogx(freq, 20 .* log10( abs(fcor2)) -1.63  ,'b');hold on;grid on;
    V=axis; axis([f1-10 Fs/2 -50 45]);
    title('Spectrum of the Correlation Functions (st*ts and s*sinv)');
    xlabel('Frequency (Hertz)');
    ylabel('Desibels(dB)');
        subplot(313),
    plot(freq, phase_fcor ,'b-');hold on;grid on;
    plot(freq(1:1e5:N), phase_fcor(1:1e5:N) ,'bo');hold on;grid on;
    plot(freq, phase_fcor2,'r-');hold on;grid on;
    plot(freq(1:1e5:N), phase_fcor2(1:1e5:N),'r+');hold on;grid on;
    title('Phase Response of the Correlation Function (s and sinv)');
    xlabel('Frequency (Hertz)');axis('tight');grid on;
    ylabel('Degrees');
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
mytoc; 



%% THE END  THE END  THE END  THE END  THE END  THE END  THE END  THE END 