function Dissonance_Bowling_interval(cpn, dur, rmp, spl, nrep, note_shift, BMF)
%% Entire "experiment" is run with this function    (adapted from ARO 2019)
% Note that "Bowling_Fig1_data.mat" file is needed for plot that compares model to data.

%% Prameters
% cpn: number of components; cpn = 6 was used in ARO 2019;
% dur: stimulus duration in sec; to match Tufts et al. dur = 0.750; (was
%       0.100 s for ARO 2019);
% rmp: ramp time; to match Tufts et al. rmp = 0.05; (was 0.01 for ARO 2019)
% spl: sound level in dB SPL; match this to Tufts et al, who used 83 dB
%       SPL; (was 70 dB PSL for ARO 2019);
% nrep: number of replication for the IHC and Synapse function
% note_shift: note shift (modulation) of F0s from A3; E.g., note_shift=2
% means root note shifts to B3;
% BMF: Best modulation frequency for IC model cell; (was 200 Hz for ARO
%       2019);
% 
%% Stimulus Generation
Fs = 1e5; % sampling rate for AN model  (100 kHz)

%Root A3 F0=220 Hz with 6 components 0.5s

F0_A3 = 220*2^(note_shift/12);
F0s = 220*2^(note_shift/12)*2.^([1:12]/12); % Bowling freqs

filter_type = 0; % filter to be applied to complex tone (0 = no filter)
Wn_freq = 0; %  Note: NA if filter_type = 0
include_fundmntl = 1; % include the fundamental freq
[pin_A3] = complex_tone(dur,rmp,F0_A3,cpn,0,0,1,spl,Fs);  % pressure wave (pressure 'input') for standard complex tone, A3


% For Interval
for ifreq = 1:length(F0s)
   pin = pin_A3 + complex_tone(dur,rmp,F0s(ifreq),cpn,filter_type,Wn_freq,include_fundmntl,spl,Fs); % compute 2nd complex tone and add to A3
   
   pin = 20e-6 * 10.^(spl/20) * pin/rms(pin);
   [mean_ic_sout_BE,std_ic_sout_BE,mean_ic_sout_BS,std_ic_sout_BS,CFs] = model_IC_LHC(pin, BMF, Fs, nrep);  % this function is included in this file, below.
%% IC model
figure
plot(CFs,std_ic_sout_BE,'b','linewidth',1.5); % blue line for BE
hold on
plot(CFs,std_ic_sout_BS,'r','linewidth',1.5); % red line for BS

std_mean_interval_BE(ifreq) = std(mean_ic_sout_BE);
std_std_interval_BE(ifreq) = std(std_ic_sout_BE);
std_mean_interval_BS(ifreq) = std(mean_ic_sout_BS);
std_std_interval_BS(ifreq) = std(std_ic_sout_BS);
%max_minus_min(ifreq) = max(mean_ic_sout) - min(mean_ic_sout);

%peak_interval(ifreq) = max(mean_ic_sout);
%m880(ifreq) = mean_ic_sout(17); % mean rate for CF = 880 Hz

%std_interval=[sig_std_Un sig_std_m2 sig_std_M2 sig_std_m3 sig_std_M3 sig_std_P4...
%    sig_std_TT sig_std_P5 sig_std_m6 sig_std_M6 sig_std_m7 sig_std_M7 sig_std_oct];
%std_chord=[sig_std_major sig_std_minor sig_std_diminished sig_std_augmented];

end
legend()

figure
subplot(2,1,1)
%intervals=0:12; % for B&H plot
intervals = 1:length(F0s); % for T&L plot
plot(intervals,std_std_interval_BE,'-*b','linewidth',1.5);
hold on
plot(intervals,std_std_interval_BS,'-*r','linewidth',1.5);
%chord_2=1:3:12;
%plot(chord_2,std_chord,'-*','linewidth',1.5);
%hold off
ylabel('Standard Deviation');
ax = gca;
ax.FontSize = 18;

std_interval_BE_norm = std_mean_interval_BE/max(std_mean_interval_BE);
std_interval_BS_norm = std_mean_interval_BS/max(std_mean_interval_BS);
%std_interval_BE_norm = (std_interval_BE-min(std_interval_BE))/(max(std_interval_BE)-min(std_interval_BE));  % normalize the std values to the maximum value
%std_interval_BS_norm = (std_interval_BS-min(std_interval_BS))/(max(std_interval_BS)-min(std_interval_BS));  % normalize the std values to the maximum value
%subplot(2,1,2)
figure
plot(intervals,std_interval_BE_norm - mean(std_interval_BE_norm) + 0.5,'-*b','linewidth',1.5); % plot normalized std values, referenced to the mean std
%plot(intervals,std_interval_BE_norm,'-*b','linewidth',1.5);
hold on
plot(intervals,std_interval_BS_norm - mean(std_interval_BS_norm) + 0.5,'-*r','linewidth',1.5); % plot normalized std values, referenced to the mean std
%plot(intervals, std_interval_BS_norm,'-*r','linewidth',1.5);
load('Bowling_Fig1_data.mat','x','y','y_norm','h');
plot(y_norm,'o-k','linewidth',0.5); % could plot Freq ratios (stored in loot.x)
plot([1 length(intervals)],[0.5, 0.5],':k')
%chord_2=1:3:12;
%plot(chord_2,std_chord,'-*','linewidth',1.5);
%hold off
ylabel('Conssonance Score (norm)');
xlabel('Dyad No.')
title('Consonance: std(Model IC Rate Profile)','fontsize',18)
% Bottom label: Freq ratio.
ax = gca;
ax.FontSize = 18;
legend('Model BE','Model BS','Bowling')


% Compute correlations between Tufts data and model std's 
[diss_corr_BE_bh, p_BE_bh] = corrcoef(y,std_std_interval_BE);
[diss_corr_BS_bh, p_BS_bh] = corrcoef(y,std_std_interval_BS);
[diss_corr_BE_hs, p_BE_hs] = corrcoef(h,std_std_interval_BE);
[diss_corr_BS_hs, p_BS_hs] = corrcoef(h,std_std_interval_BS);


end

function  [mean_ic_sout_BE,std_ic_sout_BE,mean_ic_sout_BS,std_ic_sout_BS,CFs] = model_IC_LHC(pin, BMF, Fs, nrep)
% modified from Tong's model_IC.m for ARO 2019
% BMF, Hz    % Might try lower BMF, say 100 or 150 Hz??

CFs=logspace(log10(430),log10(1600),30); % this range was reduced from Tong's wider CF range - for ARO 2019
% CFs=logspace(log10(430),log10(1600),60);

%nrep = 1; % was 10; ony 1 rep needed for longer durations
onset = 0.050*Fs; % number of points to skip at onset
dur = length(pin)/Fs; % sec
sim_dur = dur + 0.05; % simulation must be longer than stimulus duration
cohc = 1; % for healthy-ear simulation
cihc = 1; % ditto
species = 2; % 1=cat; 2=human AN model parameters (with Shera tuning sharpness)
fiberType= 3; %   % AN fiber type. (1 = low SR, 2 = medium SR, 3 = high SR)
noiseType = 1; % 0 for fixed fGn (1 for variable fGn) - this is the 'noise' associated with spontaneous activity of AN fibers - see Zilany et al., 2009. "0" lets you "freeze" it.
implnt = 0; % 0 = approximate model, 1=exact powerlaw implementation(See Zilany etal., 2009)

for i=1:length(CFs)     
    vihc = model_IHC(pin,CFs(i),nrep,1/Fs,sim_dur,cohc,cihc,species);   % "0.11" = Duration of 0.1 + 0.01 for simulation
    [an_sout,~,~] = model_Synapse(vihc,CFs(i),nrep,1/Fs,fiberType,noiseType,implnt); % an_sout = auditory-nerve synapse output

%    sigOpt = SFIE_BE_BMF(meanrate, bmf, Fs); % old version
    [ic_sout_BE,ic_sout_BS,~] = SFIE_BE_BS_BMF(an_sout, BMF, Fs);

    mean_ic_sout_BE(i) = mean(ic_sout_BE(onset:end)); % onset response excluded 
    std_ic_sout_BE(i) = std(ic_sout_BE(onset:end));
    mean_ic_sout_BS(i) = mean(ic_sout_BS(onset:end)); % onset response excluded 
    std_ic_sout_BS(i) = std(ic_sout_BS(onset:end));
end
end

function [ic_sout_BE,ic_sout_BS,cn_sout] = SFIE_BE_BS_BMF(an_sout, BMF, fs)
% Expanded the original SFIE model (Nelson & Carney, 2004 JASA) to include
% Band-Suppressed MTF (i.e. Low-Pass/Notch or High-pass) by adding a cell that is excited by
% the CN input and inhibited by the Band-Enhanced Cell (i.e. Bandpass cell) - see eNeuro, Carney et al., 2015.
% Adjustable BMF, per Ken Henry's model
% Adjustable BMF BE and BS model is Described in Carney & McDonough, 2019
%% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CN MODEL PARAMETERS:
tau_ex_cn = 0.5e-3;         % CN exc time constant
tau_inh_cn = 2.0e-3;         % CN inh time constant
cn_delay = 1.0e-3;        % "disynaptic inhibition delay" (all ANFs excitatory)
inh_str_cn = 0.6;       % re: excitatory strength == 1
afamp_cn = 1.5;         % alpha function area --> changes RATE of output cell

% IC MODEL PARAMETERS:
% BMF-dependent SFIE parameters
tau_ex_ic = 1/(10*BMF); %[0.00025 0.0005 0.001 0.002];		% Time constant excitation in seconds
tau_inh_ic = tau_ex_ic*1.5;								% Time constant inhibition in seconds
ic_delay_inh = tau_ex_ic*2;% Delay of inhibition in seconds
afamp_ic = 1;             % alpha function area --> changes RATE of output IC BE cell 
inh_str_ic = 0.9; % inhibitory strength

% BS parameters
inh_str_bs = 4;  
tau_inh_bs = tau_inh_ic; %1.0e-3; % relatively long inhibition, from BE to BS
ic_delay_bs = 1.0e-3;  % Delay from BE to BS cell (local)
Aex = 0.5; %0.3; % Rate Scalar for BS cell; note that this is effectively multiplied by afamp_ic (for Table in eNeuro)

% CN model:
% Generate frequency-domain equivalent of alpha functions
[B1, A1] = get_alpha_norm(tau_ex_cn, fs, 1);
[B2, A2] = get_alpha_norm(tau_inh_cn, fs, 1);
cn_ex = [afamp_cn*(1/fs)*(filter(B1, A1, [an_sout])) zeros(1,fs*cn_delay)];
cn_inh = [zeros(1,fs*cn_delay) afamp_cn*inh_str_cn*(1/fs)*(filter(B2, A2,[an_sout]))];

% final CN model response:
cn_sout = ((cn_ex-cn_inh) + abs(cn_ex-cn_inh))/2;   % subtract inhibition from excitation and half-wave-rectify
%cn_t = [0:(length(cn_sout)-1)]/fs;        % time vector for plotting CN responses

% IC Model #1: (SFIE; Bandpass MRF)
% Generate alpha functions for BP IC model (same as CN model, but with different taus) See Nelson & Carney 2004
[B3, A3] = get_alpha_norm(tau_ex_ic, fs, 1);
[B4, A4] = get_alpha_norm(tau_inh_ic, fs, 1);
ic_lp_ex1 = [afamp_ic*(1/fs)*(filter(B3, A3, [cn_sout])) zeros(1,floor(fs*ic_delay_inh))];
ic_lp_inh1 = [zeros(1,floor(fs*ic_delay_inh)) afamp_ic*inh_str_ic*(1/fs)*(filter(B4, A4, [cn_sout]))];
ic_sout_BE = ((ic_lp_ex1-ic_lp_inh1) + abs(ic_lp_ex1-ic_lp_inh1))/2; % half-wave rectified; standard SFIE model

%  Band-suppressed cell (see Carney et al., 2015)
[B5, A5] = get_alpha_norm(tau_inh_bs, fs, 1);
ic_bs_ex = Aex * [ic_lp_ex1 zeros(1,floor(fs*ic_delay_bs))]; % add zeros at end to match lengths  
ic_bs_inh = [zeros(1,floor(fs*ic_delay_bs)) Aex*inh_str_bs*(1/fs)*(filter(B5, A5,[ic_sout_BE]))];
ic_sout_BS = ((ic_bs_ex-ic_bs_inh) + abs(ic_bs_ex-ic_bs_inh))/2; % half-wave rectified
end

function [pin]=complex_tone(dur,rampdur,f0,ncomponents,filter_type,Wn_freq,include_fundmntl,stimdB,Fs)
%% Parameters - Complex tone (harmonic complex)  Version from UR_Ear (could be pruned!)

%dur is total duration of stimulus (s) (not same as Moore's dur)
%rampdur is duration of on or off ramp (s)
%f0 is fundamental frequency
%ncomponents is the number of components (including the fundamental, whether or not it is ultimately present)
%filter_type: 0=none,1=lowpass,2=highpass,3=bandpass,4=bandreject
%Wn_freq is the cutoff frequency or frequencies of the filter specified.
%include_fundmntl: 1 for include, 0 for omit
%stimdB is the level of the stimulus
%Fs is the sampling freq in Hz

%Other parameters for generating tones
compnt_step = 1; %Ratio between included components: 1 for every harmonic, 2 for every odd harmonic
rand_phase = 0; %Change to 1 if random starting phases are desired for each component of the complex tone
phase = 0; %Default, randomized if desired below

%Allocate time vector
t = (0:(1/Fs):dur); 

%Other parameters and error messages for filters
rectangular_filt = 0; %1=rectangular filter (fir) made through firls with impulse response length 'el'; 0=gradual filter made through fir1 with order=5000.
order = 5000; %order of fir1 filter
el=1024; %length of impulse response for firls
switch length(Wn_freq)
    case 1
        switch filter_type
            case 0
            case 1
            case 2
            case 3
                error('In order to make the requested bandpass filter, Wn_freq must be of length 2')
            case 4
                error('In order to make the requested band-reject filter, Wn_freq must be of length 2')
        end
        fCo = Wn_freq; %Cutoff frequency for high and lowpass filters
        fCo_rel_nyq = fCo/(Fs/2); %Express fCo as a percentage of Nyquist rate (this is the type of input the filter syntax requires)
    case 2
        switch filter_type
            case 0
            case 1
                error('In order to make the requested low-pass filter, input only one frequency for Wn_freq')
            case 2
                error('In order to make the requested high-pass filter, input only one frequency for Wn_freq')
            case 3
            case 4
        end
        Lower = Wn_freq(1); %Cutoff frequencies for bandpass and bandreject filters
        Upper = Wn_freq(2);
        Wn1=Lower/(Fs/2);
        Wn2=Upper/(Fs/2);
end

%% Determine frequencies to generate
int_multiplier = [1:compnt_step:ncomponents]; 
freq_array = int_multiplier*f0;
pin = 0;

%% Omit fundamental if instructed
switch include_fundmntl
    case 0
        freq_array = freq_array(2:end);
    case 1
end

%% Step through each frequency component in the complex tone.
for f = [1:length(freq_array)];
    %Randomly vary the starting phase of each component - default turned off
    switch rand_phase
        case 0
        case 1
            phase = 2*pi * rand(1,1); %random starting phase
    end
    pin = pin + cos(2 * pi * freq_array(f) *t + phase);
end

%% Optional filters to apply to complex tone
switch rectangular_filt
    case 0 %gradually sloping (steepness is determined by filter order)
        switch filter_type
            case 0 %none
            case 1 %lowpass
                b=fir1(order,fCo_rel_nyq,'low');
                pin=filter(b,1,pin);
            case 2 %highpass
                b=fir1(order,fCo_rel_nyq,'high');
                pin=filter(b,1,pin);
            case 3 %bandpass
                b=fir1(order,[Wn1,Wn2],'bandpass');
                pin=filter(b,1,pin);
            case 4 %band reject
                b=fir1(order,[Wn1,Wn2],'stop');
                pin=filter(b,1,pin);
        end
        
    case 1 %rectangular (very steep)
        switch filter_type
            case 0 %none
            case 1 %lowpass
                f = [0 Wn1 Wn1 1];
                m = [1  1   0  0];
                b= firls(el,f,m);
                pin=filter(b,1,pin);
            case 2 %highpass
                f = [0 Wn1 Wn1 1];
                m = [0  0   1  1];
                b= firls(el,f,m);
                pin=filter(b,1,pin);
            case 3 %bandpass
                f = [0 Wn1 Wn1 Wn2 Wn2 1];
                m = [0  0   1   1   0  0];
                b= firls(el,f,m);
                pin=filter(b,1,pin);
            case 4 %band reject
                f = [0 Wn1 Wn1 Wn2 Wn2 1];
                m = [0  0   1   1   0  0];
                b= firls(el,f,m);
                pin=filter(b,1,pin);
        end
end

%% Gate and scale final complex tone
pin = tukeywin(length(pin), 2*rampdur/dur)' .* pin; % apply on/off ramps
pin = 20e-6 * 10.^(stimdB/20) * pin/rms(pin); % scale signal into Pascals
end