function output_im=bilateral_filter(input_im,sigma_x,sigma_r)

%% THe source of the code is http://www.cvc.uab.es/LAMP/joost/BilateralFiltering/.
%   This is the code of following paper
%   Title: Global Color Sparseness and a Local Statistics Prior for Fast Bilateral Filtering,
%   Author: Mozerov, Mikhail G and Van De Weijer, Joost,
%   Journal: IEEE Transactions on Image Processing,
%

% standard bilateral filterinput_img
% input:
%       input_im:   input image can be grey or color
%       sigma_x :   spatial sigma
%       sigma_r :   range (tonal) sigma
% output:
%       output_im:  output image



filter_size = floor(3*sigma_x+0.5);
[hh,ww,zz]=size(input_im);
output_im=zeros(hh,ww,zz);

[yy,xx] = ndgrid(-filter_size:filter_size,-filter_size:filter_size);
gaussF = exp(-(xx.^2+yy.^2)/(2*sigma_x^2));   % normalization not necessary

CC=1/(2*sigma_r^2);
if(zz==1)
    for yy=1:hh
        for xx=1:ww
            
            x_s=max(xx-filter_size,1);
            x_e=min(xx+filter_size,ww);
            y_s=max(yy-filter_size,1);
            y_e=min(yy+filter_size,hh);
            
            patch=input_im(y_s:y_e,x_s:x_e);
            weights=exp(-(patch-input_im(yy,xx)).^2*CC).*gaussF((y_s:y_e)- yy + filter_size + 1,(x_s:x_e)- xx + filter_size + 1 );
            output_im(yy,xx)=sum(weights(:).*patch(:))./sum(weights(:));                        
        end
    end
else
    for yy=1:hh
        for xx=1:ww
            
            x_s=max(xx-filter_size,1);
            x_e=min(xx+filter_size,ww);
            y_s=max(yy-filter_size,1);
            y_e=min(yy+filter_size,hh);
            
            patch=input_im(y_s:y_e,x_s:x_e,:);
            dR=(patch(:,:,1)-input_im(yy,xx,1)).^2;
            dG=(patch(:,:,2)-input_im(yy,xx,2)).^2;
            dB=(patch(:,:,3)-input_im(yy,xx,3)).^2;
            
            weights=exp(-sum(dR+dG+dB,3)*CC).*gaussF((y_s:y_e)- yy + filter_size + 1,(x_s:x_e)- xx + filter_size + 1 );
            div=sum(weights(:));
            output_im(yy,xx,1)=sum(sum(weights.*patch(:,:,1)))/div;
            output_im(yy,xx,2)=sum(sum(weights.*patch(:,:,2)))/div;
            output_im(yy,xx,3)=sum(sum(weights.*patch(:,:,3)))/div;
        end
    end
end