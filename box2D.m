function out=box2D(in,S)
% Box filter of a 2D image.
    [m,n]=size(in);
    outtemp=in;
    out=in;
    for i=1:m
        outtemp(i,:)=box1D(in(i,:),S);
    end
    for j=1:n
        out(:,j)=box1D(outtemp(:,j),S);
    end
end

function out=box1D(in,S)
    m=length(in);
    out=zeros(size(in));
    out(1)=sum(in(1:1+S));
    for i=2:S+1
        out(i)=out(i-1)+in(i+S);
    end
    for i=S+2:m-S
        out(i)=out(i-1)+in(i+S)-in(i-(S+1));
    end    
    for i=m-S+1:m
        out(i)=out(i-1)-in(i-(S+1));
    end
end
