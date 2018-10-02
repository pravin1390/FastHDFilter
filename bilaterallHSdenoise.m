%% Example hyperspectral denoising code

clc;
clear;
close all;

%% Input Pavia
A=load('PaviaU.mat'); %% Download Pavia dataset from http://lesun.weebly.com/hyperspectral-data-set.html
Pavia=A.paviaU;
[m,n,d]=size(Pavia);
Y_scale=scaleForSVM(reshape(Pavia,m*n,d));

Y=reshape(Y_scale,m,n,d);
OriData3=Y;
I1=zeros(m,n,d);
noiselevel = 0.1;     
for i =1:d
     I1(:,:,i)=Y(:,:,i)  + noiselevel*randn(m,n);
end

sigs=3;
sigr=255*0.4;
%% Kmeans filtering

% Done in two steps : Clustering and Filtering

% Clustering
tic,
Cluster=32;
if m>256 && n>256
I2=imresize(I1,[256,256]);
elseif m>100 && n>100
I2=imresize(I1,[100,100]);   
else
I2=I1;
end
Ares=reshape(I2,size(I2,1)*size(I2,2),d);
Centre=kmeans_recursive(Ares,Cluster);
 
% Filtering
spatialtype='gaussian';     
convmethod='O1'; % 'matlab' for matlab convolutions and 'O1' for O(1) convolutions
Ikmean=fastKmeansfiltapproxinter(I1,sigs,sigr/255,Centre,spatialtype,convmethod,0);      % bilateral kmeans
Tkmeans=toc;
fprintf('Fast bilateral filter by Kmeans complete with %d clusters \n',Cluster);
fprintf('time for fast bilateral(ms)=%3.0f \n',Tkmeans*1000);

%% output
PSNRvector=zeros(1,d);
SSIMvector=zeros(1,d);
for k=1:1:d
    J=255*OriData3(:,:,k);
    I=255*Ikmean(:,:,k);
    error2 = reshape(J-I, [m*n,1]);
    mse = (sum(error2.^2)/(m*n));
    maxvalue=255;
    PSNRvector(1,k)=20*log10(maxvalue)-10*log10(mse);    
    SSIMvector(1,k)=ssimcalculate(J,I);  
end
MPSNR = mean(PSNRvector);
MSSIM = mean(SSIMvector);
OUTimag(:,:,1)=Ikmean(:,:,6);
OUTimag(:,:,2)=Ikmean(:,:,45);
OUTimag(:,:,3)=Ikmean(:,:,d);
Inimag(:,:,1)=OriData3(:,:,6);
Inimag(:,:,2)=OriData3(:,:,45);
Inimag(:,:,3)=OriData3(:,:,d);
Innoisyimag(:,:,1)=I1(:,:,6);
Innoisyimag(:,:,2)=I1(:,:,45);
Innoisyimag(:,:,3)=I1(:,:,d);
figure;
imshow(uint8(Inimag(:,:,1:3).*255))
figure;
imshow(uint8(Innoisyimag(:,:,1:3).*255))
figure;
imshow(uint8(OUTimag(:,:,1:3).*255))
fprintf('The Peak-SNR value for filtered image is %0.4f (db) \t SSIM value is %0.4f \n', MPSNR, MSSIM);


