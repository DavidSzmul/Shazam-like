%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% The exercice consists in :
%%%
%%% - Realize a time-frequency representation of signal (spectrogram)
%%%
%%% - Extract the partition of a music:
%%%     - Determine frequency of notes
%%%     - Determine occurences
%%%
%%% - Try to regenerate the music using partition
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;
% close all;
% restoredefaultpath;
addpath([cd,'\lib']);

%% User Options
NameFile = 'lac_des_cygnes_extrait.mat';

%%%%%%% Parameters
%%% Lowpass is useful to remove harmonics around 4kHz that are more
%%% difficult to determine accurately
flag_add_lowpass = 0; %%% DEFAULT -> 0
%%% Compression of frequency
flag_compression_freq = 1;%%% DEFAULT -> 1

%%%%%%% Displays for verification
flag_display_temporal = 0;
flag_display_frequential = 0;
flag_hear_song = 1;
flag_display_spectrogram = 1;
flag_display_video = 0;

%% Load File
load(NameFile);
signal_save = signal;
t=0:1/Fe:(length(signal)-1)/Fe;

if flag_add_lowpass % Remove HF (to remove harmonics but not necessary at the end)
    Fc = 2500;
    N_order = 6;
    [b,a]=butter(N_order,Fc/(Fe/2));
    signal = filter(b,a,signal);
end

%% Test Display
%%% Annotations obtained from manual temporal/frequential analysis
%%%% - dt_notes = 0.3s, useful to estimate the time resolution of the spectogramm
t_new_notes = [0.1, 0.7, 1.3, 1.6, 1.9, 2.23, 2.57]; % Timing of major notes
f_notes = [1015,1057, 1347, 2015]; % Frequency of major notes

%%% Temporal visualization
if flag_display_temporal
    figure;
    plot(t,signal);
    grid on; hold on;
    plot([t_new_notes;t_new_notes],repmat([min(signal);max(signal)], 1,length(t_new_notes)),'r--');
    xlabel('s');
    legend('Signal', 'Manual Annotation');
    title('Temporal Domain');
end
%%% Frequential visualization
if flag_display_frequential
    % -> Frequency domain of sound between 20Hz and 20kHz (Fe enables to go
    % until 11kHz)
    df_resolution = 0.01; % (In Hz): Minimum resolution demanded for FFT or SFFT
    nb_mean = 5;
    type_axis = 'lin';
    [~, ~, nff] = optimize_parameters_SFFT(signal,Fe, df_resolution, 0);
    figure;
    plot_fft(signal, Fe, nff, nb_mean,type_axis);
    grid on; hold on;
    legend('Raw FFT', 'Smoothed FFT');
    title('Frequential Domain');
end

%% Spectrogramm (or SFFT)
df_accuracy = 50; % (In Hz): Minimum resolution demanded for the FFT of SFFT
dt_study = 0.05; % (In s): Approximation of delta time of spectogramme
[nsc, nov, nff] = optimize_parameters_SFFT(signal,Fe, df_accuracy, dt_study);

%%% Use Hamming window in order to improve performance by reducing border effects
win_hamm = hamming(nsc);
[~,F,T,P] = spectrogram(signal,win_hamm,nov,nff,Fe);
fprintf(['Frequency accuracy: ',num2str(F(2)-F(1)),' Hz\n',...
    'Time Window: ',num2str(nsc/Fe),' s\n',...
    'Delta Time: ',num2str(T(2)-T(1)),' s\n']);
real_freq_resolution = F(2)-F(1);
real_time_resolution = nsc/Fe;

%% Detection of notes
thr_min_DSP_dB = -50; % Minimum DSP power: below is considered as noise (obtained after analysing spectogram)
thr_min_DSP = 10^(thr_min_DSP_dB/10); %%% According to Matlab (why not 20?)

%%% Detect contours that contains notes (one contour is considered as one
%%% note)
P_analysis = P;
P_analysis(P<thr_min_DSP)=0;

%%%% ONLY if Image Processing Toolbox available
% [Boundary,Groups] = bwboundaries(P_analysis,'noholes');
%%%% Homemade equivalent function
[Boundary,Groups] = detect_groups(P_analysis);


nb_group = length(Boundary);
P_reshape = P_analysis(:);
frequency_group = zeros(1,nb_group); % Frequency of notes
time_group = zeros(1,nb_group); % Timing of notes
power_group = zeros(1,nb_group); % Power of notes

%%% Extract more pertinent frequency for each group
for i = 1:length(frequency_group)
    %%% Determine maximum signal inside group
    index_reshape = find(Groups==i);
    [~, index_in_frequency_reshape] = max(P_reshape(index_reshape));
    index_frequency_reshape = index_reshape(index_in_frequency_reshape);
    
    %%% Deduce frequency by 'deshaping'
    index_frequency = mod(index_frequency_reshape-1, size(Groups,1))+1;
    frequency_group(i)= F(index_frequency);
    power_group(i) = max(P_reshape(index_reshape));
    
    %%% Minimum time of group is the temporal start of the contour
    index_time = (index_frequency_reshape-index_frequency)/size(Groups,1)+1;
    time_group(i) = T(index_time)-real_time_resolution/2;
end

%%% IMPROVEMENTS TODO
% %%% Determine for each group the number of note
% %%% (it may happen that multiple notes of same frequency appears at
% %%% different times in a single group)
% for i=1:length(frequency_group)
% end


%% Compression frequency
%%% This part is not mandatory, it is only done because some frequencies
%%% are close each other, counted as a different note
%%% and the exercice is to determine each occurence of one same note
if flag_compression_freq
    tolerance_freq = 2*real_freq_resolution;
    disp(['Tolerance of estimated frequency around: ',num2str(tolerance_freq),' Hz']);
    frequency_group = floor(frequency_group/tolerance_freq)*tolerance_freq;
end

%% Generation of Partition and Encoded signal
%%% Outputs required for the exercice
Partition.frequency = sort(unique(frequency_group)); % Frequency of notes (In Hz)
Partition.occurence = cell(size(Partition.frequency));% Timing of notes (In s)
for i=1:length(Partition.frequency)
    Partition.occurence{i} = sort(time_group(frequency_group==Partition.frequency(i)));    
end
%%% Output as confirmation of good encoding
encoded_signal = generate_music_v2(time_group,frequency_group,power_group,Fe,real_freq_resolution, length(signal));

%% Save outputs
NameFile_new = 'lac_des_cygnes_encoded';
save([NameFile_new,'.mat'], 'Partition', 'encoded_signal', 'Fe');
audiowrite([NameFile_new,'.wav'],encoded_signal,Fe);
disp('Encoded audio saved');


%% Displays
%%% Hear song (enables to identify which frequency domain can be considered as noise)
if flag_hear_song
    %%% Original Sound
    sound(signal,Fe);
    pause(t(end)+0.5);
    
    %%% Auto generated song based on the partition deduced by algorithm
    sound(encoded_signal,Fe);
    pause(t(end)+0.5);
    
    %%% Compare Both signals
    figure;
    plot(t,signal);
    grid on; hold on;
    plot(t,encoded_signal);
    xlabel('s');
    legend('Original Signal', 'Encoded Signal');
    title('Encoding of signal');    
end

if flag_display_spectrogram
    figure;
    sp_spec = tight_subplot(2,1,[.01 .03],[.1 .1],[.1 .1]);
    axes(sp_spec(1)); % Without threshold
    spectrogram(signal,hamming(nsc),nov,nff,Fe,'yaxis');
    hold on; plot([t_new_notes;t_new_notes],repmat([0;Fe/2], 1,length(t_new_notes)),'r--');
    
    axes(sp_spec(2)); % With threshold
    spectrogram(signal,hamming(nsc),nov,nff,Fe,'MinThreshold',thr_min_DSP_dB,'yaxis');
    hold on; plot([t_new_notes;t_new_notes],repmat([0;Fe/2], 1,length(t_new_notes)),'r--');
    linkaxes(sp_spec,'xy');
    
    for k = 1:length(Boundary) %%% Plot detected notes
        T_boundary = T(Boundary{k}(:,2));
        F_boundary = F(Boundary{k}(:,1));
        plot(T_boundary, F_boundary/1000, 'g--', 'LineWidth', 1); % In kHz in spectogramm
        hold on;
    end
    
end
if flag_display_video
    dt_pause = 0.2;
    figure;
    sp(1)=subplot(2,1,1);sp(2)=subplot(2,1,2);
    axes(sp(1));
    plot(t,signal);
    grid on; hold on;
    p_line = plot([0,0], [min(signal);max(signal)],'g');
    plot([t_new_notes;t_new_notes],repmat([min(signal);max(signal)], 1,length(t_new_notes)),'r--');
    
    axes(sp(2));
    p_fft =plot(F,P(:,1));
    grid on; hold on;
    plot([F(1),F(end)],[thr_min_DSP, thr_min_DSP],'r--');
    %     ylim([0,max(max(P))/4]);
    ylim([0,thr_min_DSP*4]);
    legend('SFFT signal', 'Threshold Min');
    
    for i=1:length(T)
        set(p_line,'XData',[T(i),T(i)]);
        set(p_fft,'YData',P(:,i));
        pause(dt_pause);
    end
end