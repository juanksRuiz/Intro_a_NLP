clear
close all
clc

%% Setting the interfaz to acquire de audio data

fs = 8000 ; % Sampling Frequency
nBits = 16 ; % Number of bits per audio sample
nChannels = 1 ; % Number of Channels
ID = -1; % default audio input device 
recObj = audiorecorder(fs,nBits,nChannels,ID);
ns = 5; % Seconds of recording

%% Acquiring audio data

disp('Start speaking.') 
recordblocking(recObj, ns);
disp('End of Recording.');

%% Extracting audio

y = getaudiodata(recObj); % Extracting the audio data
y = y./std(y);
sound(y,fs) % Playing the recorded sound

%% Computing the spectogram

len_win = fix(0.02*fs); % Defining the samples for a 20ms window
window = hanning(len_win); % Computing the window length to compute the spectrogram.
len_over = fix(0.01*fs); % Defining number of samples for an overlap of 10ms between consecutive windows.
nfft = 512; % Number of points for the FFT

[s,f,t] = spectrogram(y,window,len_over,nfft,fs); % Computing the spectrogram 

Sw = abs(s); % Computing the absolute value of the spectrogram.
F = f*ones(1,length(t)); % frequency matrix for plotting
T = ones(length(f),1)*t; % time matrix for plotting

% plotting the spectrogram
figure, surf(T,F,Sw,'EdgeColor','none')
view([0,90])
xlabel('Time (s)')
ylabel('Frequency (Hz)')

% Computing the Mel-Spectrogram
n_coeff = 32;
[Sm,fm,tm] = melSpectrogram(y,fs,'Window',window,'OverlapLength',len_over, 'FFTLength',nfft,'NumBands',n_coeff,'FrequencyRange',[20,fs/2]);

Tm = tm*ones(1,length(fm)); % frequency matrix for plotting
Fm = ones(length(tm),1)*fm; % time matrix for plotting

% ploting the Mel-Spectrogram
figure, surf(Tm',Fm',Sm,'EdgeColor','none')
view([0,90])
xlabel('Time (s)')
ylabel('Frequency (Hz)')

%% Detecting changes in the Phonems

% using the spectrogram
D_Sw = sqrt(sum( (Sw(:,2:end)-Sw(:,1:end-1)).^2 ));
% using the Mel-spectrogram
D_Sm = sqrt(sum( (Sm(:,2:end)-Sm(:,1:end-1)).^2 ));

%% Extracting elements for the different fonems

% Using the Spectrogram

figure, plot(D_Sw)
a_s = 38:85;
e_s = 135:183;
i_s = 231:277;
o_s = 325:366;
u_s = 416:449;

% Using the Mel Spectrogram

figure, plot(D_Sm)
a_m = 19:42;
e_m = 67:89;
i_m = 115:137;
o_m = 163:181;
u_m = 208:221;

% Plots for the signature of the diferent vocal phonems using the
% spectrogram
figure, 
subplot(6,4,[1 2 5 6])
plot(f,mean(Sw(:,a_s+1),2),'k','LineWidth',2)
ylabel('Phonem a')
xlabel('Frequency [Hz]')
grid on
hold on
plot(f,mean(Sw(:,a_s+1),2)+std(Sw(:,a_s+1),[],2),'k--')
plot(f,mean(Sw(:,a_s+1),2)-std(Sw(:,a_s+1),[],2),'k--')

subplot(6,4,[3 4 7 7])
plot(f,mean(Sw(:,e_s+1),2),'k','LineWidth',2)
ylabel('Phonem e')
xlabel('Frequency [Hz]')
grid on
hold on
plot(f,mean(Sw(:,e_s+1),2)+std(Sw(:,e_s+1),[],2),'k--')
plot(f,mean(Sw(:,e_s+1),2)-std(Sw(:,e_s+1),[],2),'k--')

subplot(6,4,[9 10 13 14])
plot(f,mean(Sw(:,i_s+1),2),'k','LineWidth',2)
ylabel('Phonem i')
xlabel('Frequency [Hz]')
grid on
hold on
plot(f,mean(Sw(:,i_s+1),2)+std(Sw(:,i_s+1),[],2),'k--')
plot(f,mean(Sw(:,i_s+1),2)-std(Sw(:,i_s+1),[],2),'k--')

subplot(6,4,[11 12 15 16])
plot(f,mean(Sw(:,o_s+1),2),'k','LineWidth',2)
ylabel('Phonem o')
xlabel('Frequency [Hz]')
grid on
hold on
plot(f,mean(Sw(:,o_s+1),2)+std(Sw(:,o_s+1),[],2),'k--')
plot(f,mean(Sw(:,o_s+1),2)-std(Sw(:,o_s+1),[],2),'k--')

subplot(6,4,[18 19 22 23])
plot(f,mean(Sw(:,u_s+1),2),'k','LineWidth',2)
ylabel('Phonem u')
xlabel('Frequency [Hz]')
grid on
hold on
plot(f,mean(Sw(:,u_s+1),2)+std(Sw(:,u_s+1),[],2),'k--')
plot(f,mean(Sw(:,u_s+1),2)-std(Sw(:,u_s+1),[],2),'k--')

% Plots for the signature of the diferent vocal phonems using the
% Mel-spectrogram
figure, 
subplot(6,4,[1 2 5 6])
plot(fm,mean(Sm(:,a_m+1),2),'k','LineWidth',2)
ylabel('Phonem a')
xlabel('Frequency [Hz]')
grid on
hold on
plot(fm,mean(Sm(:,a_m+1),2)+std(Sm(:,a_m+1),[],2),'k--')
plot(fm,mean(Sm(:,a_m+1),2)-std(Sm(:,a_m+1),[],2),'k--')

subplot(6,4,[3 4 7 7])
plot(fm,mean(Sm(:,e_m+1),2),'k','LineWidth',2)
ylabel('Phonem e')
xlabel('Frequency [Hz]')
grid on
hold on
plot(fm,mean(Sm(:,e_m+1),2)+std(Sm(:,e_m+1),[],2),'k--')
plot(fm,mean(Sm(:,e_m+1),2)-std(Sm(:,e_m+1),[],2),'k--')

subplot(6,4,[9 10 13 14])
plot(fm,mean(Sm(:,i_m+1),2),'k','LineWidth',2)
ylabel('Phonem i')
xlabel('Frequency [Hz]')
grid on
hold on
plot(fm,mean(Sm(:,i_m+1),2)+std(Sm(:,i_m+1),[],2),'k--')
plot(fm,mean(Sm(:,i_m+1),2)-std(Sm(:,i_m+1),[],2),'k--')

subplot(6,4,[11 12 15 16])
plot(fm,mean(Sm(:,o_m+1),2),'k','LineWidth',2)
ylabel('Phonem o')
xlabel('Frequency [Hz]')
grid on
hold on
plot(fm,mean(Sm(:,o_m+1),2)+std(Sm(:,o_m+1),[],2),'k--')
plot(fm,mean(Sm(:,o_m+1),2)-std(Sm(:,o_m+1),[],2),'k--')

subplot(6,4,[18 19 22 23])
plot(fm,mean(Sm(:,u_m+1),2),'k','LineWidth',2)
ylabel('Phonem u')
xlabel('Frequency [Hz]')
grid on
hold on
plot(fm,mean(Sm(:,u_m+1),2)+std(Sm(:,u_m+1),[],2),'k--')
plot(fm,mean(Sm(:,u_m+1),2)-std(Sm(:,u_m+1),[],2),'k--')