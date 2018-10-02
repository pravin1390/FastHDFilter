function out = gauss_filter(in, sigma, fast_flag)
% Gaussian Smoothing at standard deviation euqal to sigma.
%
% fast_flag: 1 for standard separable Gaussian convolution
%            2 for image rescaling before Gaussian convolution for speedup.

if (fast_flag==0)
    out=in;
    for i=1:size(in,3)
         out(:,:,i)=young(in(:,:,i),sigma);
    end
else
    % Down sampling factor
    DSfaktor=round(min(max(1,2*floor((sigma+1)/4)),6));
    
    [hh,ww,~]=size(in);
    in=imresize(in,1/DSfaktor,'bilinear');
    out=in;
    for i=1:size(in,3)
        out(:,:,i)=young(in(:,:,i),sigma/DSfaktor);
    end
    out=imresize(out,[hh,ww]);
end
