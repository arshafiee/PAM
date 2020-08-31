%>>>>>>>>>> Binary Digital Communications <<<<<<<<<<%

clc;
clear;

%-----------

length_of_transmiting_signal = 10;
M = 100;
T = 0.001;
t = 0:T/M:(length_of_transmiting_signal*T)-T/M;     
R = 1/T;                                                          %bit_rate

select = input('Please choose a carrier signal :   1--->sinusoidal  0--->pulse\n');
sNr = input('Please input SNR :\n');

t1= 0:T/M:T-T/M;
S_T = 1*cos(2*pi/T*t1*select);                              %carrier signal

%>>>>>>>>>> Binary information at transmitter <<<<<<<<<<%

X=rand(1,length_of_transmiting_signal);
X=(X < 0.5); 
disp('Binary information at transmitter :');
disp(X);

%>>>>>>>>>> representation of transmitting binary information as digital
%signal <<<<<<<<<<%

bit=repmat(X,M,1);
bit=bit(1:end);                                            
figure();
subplot(3,1,1);
plot(t,bit,'LineWidth',2.5);
grid on;
axis([0 T.*length_of_transmiting_signal -0.5 1.5]);
ylabel('amplitude(volt)');
xlabel(' time(sec)');
title('transmitting information as a pulese');

%>>>>>>>>>> Binary modulation <<<<<<<<<<%

bit=2*repmat(X,M,1)-1; 
bit=S_T'.*bit;
Carrier=bit(1:end);
subplot(3,1,2);
plot(t,Carrier,'LineWidth',2.5);
grid on;
axis([0 T*length_of_transmiting_signal min(Carrier)-2 max(Carrier)+2]);
ylabel('amplitude(volt)');
xlabel(' time(sec)');
title('waveform for binary AM modulation coresponding binary information');

%>>>>>>>>>> AWGN (Additive White Gaussian Noise <<<<<<<<<<%

Noisy_Carrier = awgn(Carrier,sNr);
subplot(3,1,3);
plot(t,Noisy_Carrier,'LineWidth',2.5);
grid on;
axis([0 T*length_of_transmiting_signal min(Noisy_Carrier)-2 max(Noisy_Carrier)+2]);
ylabel('amplitude(volt)');
xlabel(' time(sec)');
title('AWGN noise is now added');

%>>>>>>>>>> Binary Demodulation <<<<<<<<<<%

u = @(t) t>=0;
Demod_func = S_T.*(u(t1) - u(t1-T));                  %Demodulator function
Demod_signal = conv(Demod_func,Noisy_Carrier);            %Demoduled signal
figure();
subplot(3,1,1);
plot(t1,Demod_func,'LineWidth',2.5);
grid on;
 axis([0 T*length_of_transmiting_signal min(S_T)-1 max(S_T)+1]);
ylabel('amplitude(volt)');
xlabel(' time(sec)');
title('Demodulator function');

subplot(3,1,2);
plot(0:T/M:(length(Demod_signal)*(T/M)-T/M),Demod_signal,'LineWidth',2.5);
grid on;
ylabel('amplitude(volt)');
xlabel(' time(sec)');
title('Demodulated signal(using matched-filter)');

th=.8*max(conv(Demod_func,Demod_func));
Y=zeros(1,length_of_transmiting_signal);
for i=1:length_of_transmiting_signal
    
     Y(i)=Demod_signal(i*(M))>th;

end

bit1=repmat(Y,M,1);
bit1=bit1(1:end);                                            
subplot(3,1,3);
plot(t,bit1,'LineWidth',2.5);
grid on;
axis([0 T.*length_of_transmiting_signal -0.5 1.5]);
ylabel('amplitude(volt)');
xlabel(' time(sec)');
title('recived information ');
disp('comparison beween transmited and recived data');
disp([Y' X']'); 
disp('BER(Bit Error Rate)');
disp(1-mean(Y==X));

%%
%>>>>>>>>>> Binary Digital Communications <<<<<<<<<<%

clc;
clear;

%--------------------

length_of_transmiting_signal=1000;
M=100;
T = 0.001;
t = 0:T/M:(length_of_transmiting_signal*T)-T/M;     
R = 1/T;

t1=0:T/M:(length_of_transmiting_signal*T/length_of_transmiting_signal)-T/M;
S_T = 10;                                                   %Carrier Signal

%-----------

for sNr=-80:20
    
%>>>>>>>>>> Binary information at transmitter <<<<<<<<<<%

X=rand(1,length_of_transmiting_signal);
X=(X<0.5); 

%>>>>>>>>>> representation of transmitting binary information as digital 
%signal <<<<<<<<<<%

bit=repmat(X,M,1);
bit=bit(1:end); 

%>>>>>>>>>> Binary modulation <<<<<<<<<<%

bit=2*repmat(X,M,1)-1; 
bit=S_T'.*bit;
Carrier=bit(1:end);

%>>>>>>>>>> AWGN (Additive White Gaussian Noise) <<<<<<<<<<%

Noisy_Carrier = awgn(Carrier,sNr);

%>>>>>>>>>> Binary Demodulation <<<<<<<<<<%

u = @(t) t>=0;
Demod_func = S_T.*(u(t1) - u(t1-T));                  %Demodulator function
Demod_signal = conv(Demod_func,Noisy_Carrier);            %Demoduled signal

th=.9*max(conv(Demod_func,Demod_func));
Y=zeros(1,length_of_transmiting_signal);
    for i=1:length_of_transmiting_signal
         Y(i)=Demod_signal(i*(M))>th;
    end
    
[Y' X']; %#ok
bre(sNr+81)=1-mean(Y==X); %#ok
end

figure();
stem(-80:20,bre);
xlabel('SNR(signal to noise ratio)');
ylabel('BER(Bit Error Rate)');
title('BER for different SNRs');

