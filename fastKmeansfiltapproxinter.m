function A=fastKmeansfiltapproxinter(A,S,h,Centre,spatialkernel,convmethod,fast_flag,Aguide)
% Main high-dimensional filtering code

if ~exist('Aguide','var')
     % guideimage is same as inputimage
     Aguide=A;
end
if ~exist('fast_flag','var')
     % guideimage is same as inputimage
     fast_flag=1;
end
guided=size(Aguide,3);
[m,n,~]=size(A);
B=zeros(size(A));
Cluster=size(Centre,1);

%% Forming intermediate images and coefficients
C1=zeros(Cluster,Cluster);
for i=1:Cluster
    C1(i,i)=1;
    for j=i+1:Cluster
        C1(i,j)=exp(-sum((Centre(i,:)-Centre(j,:)).^2,2)/(2*h*h));
        C1(j,i)=C1(i,j);
    end
end
C1chan=pinv(C1);
W=zeros(m,n,Cluster);
for i=1:Cluster
    W(:,:,i)=sum((Aguide-reshape(Centre(i,:),1,1,guided)).^2,3);   
end
W=exp(-bsxfun(@rdivide,W,(2*(h^2))));
Wb=zeros(m,n);

%% Filtering using matlab command for convolutions

if strcmp(convmethod,'matlab')
    if strcmp(spatialkernel,'box')
        filt     = ones(2*S+1,2*S+1);       
    elseif strcmp(spatialkernel,'gaussian')       
        w  = round(6*S); if (mod(w,2) == 0); w  = w+1; end
        filt     = fspecial('gaussian', [w w], S);
    else
    end        
    for i=1:Cluster
        Wt=sum(W.*reshape(C1chan(i,:),1,1,Cluster),3);
        Wb=Wb+bsxfun(@times,Wt,imfilter(W(:,:,i),filt));
        B=B+bsxfun(@times,Wt,imfilter(bsxfun(@times,W(:,:,i),A),filt));              
    end
end

%% Filtering using O(1) convolutions
if strcmp(convmethod,'O1')
    for i=1:Cluster
        Wt=sum(W.*reshape(C1chan(i,:),1,1,Cluster),3);
        if strcmp(spatialkernel,'box')
            Wb=Wb+bsxfun(@times,Wt,box_filter(W(:,:,i),S,fast_flag));
            B=B+bsxfun(@times,Wt,box_filter(bsxfun(@times,W(:,:,i),A),S,fast_flag));          
        elseif strcmp(spatialkernel,'gaussian')
            Wb=Wb+bsxfun(@times,Wt,gauss_filter(W(:,:,i),S,fast_flag));
            B=B+bsxfun(@times,Wt,gauss_filter(bsxfun(@times,W(:,:,i),A),S,fast_flag));            
        else
        end 
    end
end 
A=bsxfun(@rdivide,B,Wb);

end

