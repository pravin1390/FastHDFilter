function out = box_filter(in, S, fast_flag)
% Routine for box convolutions
% fast_flag: 1 for standard separable box convolution
%            2 for image rescaling before box convolution for speedup.

if (fast_flag==0)
    out=in;
    for i=1:size(in,3)
        out(:,:,i)=box2D(in(:,:,i),S);
    end
else
    % Down sampling factor
    DSfaktor=4;   
    [hh,ww,~]=size(in);
    in=imresize(in,1/DSfaktor,'bilinear');
    out=in;
    for i=1:size(in,3)
        out(:,:,i)=box2D(in(:,:,i),ceil(S/DSfaktor));       
    end
    out=imresize(out,[hh,ww]);
end
