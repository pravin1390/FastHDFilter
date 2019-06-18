%% Example non-local means code

clc;
close all;
clear;

%% input
% [filename, user_canceled] = imgetfile;
filename='eyes_closeup.png';
I  =  double(imread(filename));
[m,n,d]=size(I);
sigma=0.08*255;
S=10;K=3;
fast_flag=1;
Iact=I;
I=I+sigma*randn(m,n,d);
I2=I./256;

error2 = reshape(Iact-I, [d*m*n,1]);
mse = (sum(error2.^2)/(d*m*n));
maxvalue=255;
max(I(:));
peaksnrvalue=20*log10(maxvalue)-10*log10(mse);
fprintf('\n The Peak-SNR value for noisy image is %0.4f (db) \n', peaksnrvalue);
%% Kmeans filtering

% Done in two steps : Clustering and Filtering
tic,
Cluster=31;
pcadim=25;

% Clustering
Apca=compute_pca(I2, K, pcadim);
pcadim=size(Apca,3);
Apcares=imresize(Apca,[256 256]);
Ares=reshape(Apcares,size(Apcares,1)*size(Apcares,2),pcadim);
Centre=kmeans_recursive(Ares,Cluster);

% Filtering
spatialtype='box';     
convmethod='O1'; % 'matlab' for matlab convolutions and 'O1' for O(1) convolutions
Ikmean=fastKmeansfiltapproxinter(I2,S,3.5*sigma/256,Centre,spatialtype,convmethod,fast_flag,Apca);      % nlm kmeans
Ikmean=Ikmean.*255;
Ikmean(Ikmean>=255)=255;
Ikmean(Ikmean<=0)=0;
toc
Tkmeans=toc;
fprintf('non local means by Kmeans complete with %d clusters \n',Cluster);
fprintf('time for non local means (ms)=%3.0f \n',Tkmeans*1000);

%% Displaying noisy and filtered image
figure;
imshow(uint8(Iact)),%title('Original image');
figure;
imshow(uint8(I)),%title('input image');
figure;
imshow(uint8(Ikmean));%title(['fast non local means with ',num2str(Cluster),' clusters']);

%% PSNR calculation
error2 = reshape(Iact-Ikmean, [d*m*n,1]);
mse = (sum(error2.^2)/(d*m*n));
maxvalue=255;
peaksnrvalue=20*log10(maxvalue)-10*log10(mse);
ssimval=ssimcalculate(Iact,Ikmean);
fprintf('\n mean sq error=%f \t The Peak-SNR value for filtered image is %0.4f (db) \t SSIM value is %0.4f \n',sqrt(mse), peaksnrvalue, ssimval);

