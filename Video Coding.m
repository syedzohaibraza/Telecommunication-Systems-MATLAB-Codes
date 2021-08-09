clear

% close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Simulation Parameters
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


N = 8;                                                              % N x N is the size of the unitary transform
                                                                    
Qscale = 1;                                                                                                                                      
                                                                                  
image = imread('Samples\lena512.tif');

QstepsY =  ...
[16 11 10 16 24 40 51 61;                                           % JPEG quantizer steps for 8x8 DCT
12 12 14 19 26 58 60 55;
14 13 16 24 40 57 69 56;
14 17 22 29 51 87 80 62;
18 22 37 56 68 109 103 77;
24 35 55 64 81 104 113 92;
49 64 78 87 103 121 120 101;
72 92 95 98 112 100 103 99];

% Alternative quantitzation schemes

% QstepsY =  ...
% [  5     8    13    20    29    40    53    68;                                                                                      
%    8    11    16    23    32    43    56    71;
%   13    16    21    28    37    48    61    76;
%   20    23    28    35    44    55    68    83;
%   29    32    37    44    53    64    77    92;
%   40    43    48    55    64    75    88   103;
%   53    56    61    68    77    88   101   116;
%   68    71    76    83    92   103   116   131;];
% 

% QstepsY =  ...
% [  5     8    13    20   29    255   255    255;
%    8    11    16    23   255   255   255    255;
%   13    16    21    255  255   255   255    255;
%   20    23    255   255  255   255   255    255;
%   29    255   255   255  255   255   255    255;
%   255   255   255   255  255   255   255    255;
%   255   255   255   255  255   255   255    255;
%   255   255   255   255  255   255   255    255;];



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% No changes beyond this line !!

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[Height,Width,Depth] = size(image);

if mod(Height,N) ~= 0 Height = floor(Height/N)*N; end;              % Adapt image size to DCT length
if mod(Width,N) ~= 0 Width = floor(Width/N)*N; end
x = image(1:Height,1:Width,:);

N_blocks = Height / N;                                              % number of DCT blocks in x direction (width)
M_blocks = Width / N;                                               % number of DCT blocks in y direction (height)

x = double(image(:,:,1));
x2 = zeros(size(x));                                                % reconstructed image

qx = zeros(Height,Width);
acBitsY = 0;
dcBitsY = 0;

% Compute the bits for the Y component
for m = 1:N:Height
    for n = 1:N:Width
        t = x(m:m+N-1,n:n+N-1) - 128;
        y_dct = dct_2(t,8);                                             % N x N 2D DCT of input image      
        
        temp = floor(y_dct./(Qscale*QstepsY) + 0.5);                    % quantize the DCT coefficients
        
        % Calculate bits for the DC difference
        if n == 1
            dc = temp(1,1);
            dcBitsY = dcBitsY + jpeg_bits(dc,'DC','Y');
        else
            dc = temp(1,1) - dc;
            dcBitsY = dcBitsY + jpeg_bits(dc,'DC','Y');
            dc = temp(1,1);
        end
                
        ACblkBits = jpeg_bits(temp,'AC','Y');                       % Calculate the bits for the AC coefficients
        acBitsY = acBitsY + ACblkBits;
        
        y_dct_quant(m:m+N-1,n:n+N-1) = temp .* (Qscale*QstepsY);
        
        qx(m:m+N-1,n:n+N-1)= idct_2((temp .* (Qscale*QstepsY)),8)+ 128;  % dequantize & IDCT the DCT coefficients
        
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computation of SNR and average bit rate

diff = x - qx;                                                      % Computation of SNR
mse = std_2(diff);
snr = 20*log10(std(x(:))/mse);
snr = sprintf('SNR = %4.2f dB',snr);

TotalBits = acBitsY + dcBitsY;                                      % Average bit rate
av_bitrate = TotalBits/(Height*Width);
av_bitrate = sprintf('bit rate: %4.2f bpp',av_bitrate)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computation of coding gain

y_cg = dct_2(x,8);

for i = 1:N
    for j = 1:N
        t = y_cg(i:N:Height,j:N:Width);
        st = std_2(t);
        CoefVar(i,j) = st * st;
    end
end

MeanVar = mean(CoefVar(:));
P1 = CoefVar .^ (1/(N^2));
GeoMeanVar = prod(P1(:));                                           % geometric mean of the coef. variances
G = 10*log10(MeanVar/GeoMeanVar)                                    % coding gain
sG = sprintf('Coding gain = %4.2f dB',G);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(1),imagesc(x); colormap(gray)
title(['Original bit rate: 8 bpp'])

figure(2),imagesc(qx); colormap(gray)
title(['JPEG compressed ', num2str(av_bitrate) ',  '  num2str(snr)] )

figure(3),imagesc(diff); colormap(gray)
title(['Difference: Original - coded image'])
