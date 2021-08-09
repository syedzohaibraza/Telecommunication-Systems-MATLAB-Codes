% Shows time domain and spectral domain properties of time-discrete and time-continuous pseudorandom sequences
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear 
clf


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Simulation parameters
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

 L = 4;                                                      % Length of shift register
 start_val = '1000';                                         % Initial values of state memories
 fb = [1 , 4];                                               % Feedback paths


% L = 8;
% start_val = '10000000';                                   % Start value sequence 3
% fb = [3 , 5, 6, 8];

% L = 8;
% start_val = '11111110';                                   % Start value sequence 4
% fb = [1 , 7, 8];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% No changes after this line !!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

val = start_val;
M = 2^L-1;                                              % Number of possible states

for k=1:M
          
    for n = 1:length(fb)    
        xor_in(n) = str2num(val(fb(n)));				% Select input bits for xor operation                                              
    end;    
    
    result_xor = xor(xor_in(1),xor_in(2));              % 1st xor operation
    for n =1:length(fb)-2                               % Additional xor operations
        result_xor = xor(result_xor, xor_in(n+2));
    end;    

    y_mls(k) = result_xor;                              % Output value	


    for n = L:-1:2                                      % Shift operation
        val(n)=val(n-1);                               
    end;
    val(1)=num2str(result_xor);
    
    x = 1;
end;

y_mls = (y_mls-0.5)*2;
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Simulation of time-continous rectangular signal
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lr = 100;
f_T = 2000;                                               % Carrier frequency at receiver
Ts=.0001;                                                 % Symboldauer entspricht                                             

y_mls_rect = zeros(M*lr,1);                               % Rectangular signal  

for k = 1:M    
    y_mls_rect((k-1)*lr+1:k*lr) = y_mls(k)*ones(lr,1);        
end;    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xaxis =  [1:M*lr]/lr;

figure(1);stem((0:M-1),y_mls,'LineWidth',3); grid;
h = gca;
set(h,'FontSize',[14]) 
set(h,'YTick',[-1 -0.5 0 0.5  1])
xlabel('k \rightarrow');
title('Pseudorandom Sequence y(k)');
axis([-1 M, -1.1 1.1])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Spectrum of MLS
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Y = fft([y_mls y_mls y_mls y_mls y_mls y_mls y_mls y_mls]);

figure(2);plot(abs(Y),'Color',[1 0 0],'LineWidth',3);grid;

h = gca;
set(h,'FontSize',[14])
xlabel(' Frequency bin m \rightarrow');    
title(' Magnitude Spectrum |Y(m)|')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Spectrum of modulated rectangular signal
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xaxis =  [1:M*lr]/lr;

figure(3);plot(xaxis,y_mls_rect,'LineWidth',3); grid;
h = gca;
set(h,'FontSize',[14]) 
set(h,'YTick',[-1 -0.5 0 0.5  1])
xlabel('k \rightarrow');
title('Pseudorandom sequence y(t)');
axis([0 M, -1.1 1.1])

y_mls_rect = y_mls_rect.';
y_mls_rect_per = [y_mls_rect y_mls_rect y_mls_rect y_mls_rect y_mls_rect y_mls_rect y_mls_rect y_mls_rect y_mls_rect];

N=length(y_mls_rect_per);                                                % Number of symbols to be simulated

time=Ts*(N-1); t=0:Ts:time;                                              % Sampling interval and time vectors


z=real(y_mls_rect_per.*exp(j*(2*pi*f_T*t)));                            % Modulation by complex carrier

N=length(z);                                                            % length of the signal x
t=Ts*(1:N);                                                             % define a time vector
ssf=(ceil(-N/2):ceil(N/2)-1)/(Ts*N);                                    % frequency vector
fx=fft(z(1:N));                                                         % do DFT/FFT
fxs=fftshift(fx);                                                       % shift it for plotting

figure(4);plot(ssf,20*log10(abs(fxs)),'LineWidth',2,'Color',[1 0 0]);grid   % plot magnitude spectrum

h = gca;
set(h,'FontSize',[14])

xlabel('Frequency in Hz   \rightarrow'); title('Magnitude Spectrum |Y(f)|')       

if L == 4
    axis([1500 2500 0 80])
end;

if L == 8
    axis([1850 2150 0 80])
end;



