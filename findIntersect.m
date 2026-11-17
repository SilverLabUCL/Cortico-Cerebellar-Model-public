function tintersect = findIntersect(x,y1,y2,xlogflag,ylogflag)
% function is a cross between answers from:
% https://uk.mathworks.com/matlabcentral/answers/1987333-how-to-find-the-intersection-of-two-curves-with-the-input-data-being-two-vectors

L=1000;

if(all(isnan(y1)))
    tintersect=NaN(numel(x),1);
else
    % remove nans
    y1(isnan(y1))=max(y1);
    y2(isnan(y2))=max(y2);

    if(xlogflag==1)
        xint=logspace(log10(min(x)),log10(max(x)),L);
    else
        xint=linspace(min(x),max(x),L);
    end
    if(ylogflag==1)
        y1int = interp1(x,log(y1),xint);
        y2int = interp1(x,log(y2),xint);
    else
        y1int = interp1(x,y1,xint);
        y2int = interp1(x,y2,xint);
    end

    % find the x values where difference of functions changes sign
    % ixv=find(diff(sign(y1 - y2)));
    % L=numel(x);
    ixv=find(diff(sign(y1int - y2int)));
    
    tintersect=zeros(numel(ixv),1);
    nsteps=6;
    for k = 1:numel(ixv)
        idxrng = max(1, ixv(k)-nsteps) : min(L,ixv(k)+nsteps);
        xrange=[xint(idxrng(1)) xint(idxrng(end))];
        % find intersection between pointwise curves
        tintersect(k) = fzero(@zerofun,xrange);
    end
    if(isempty(tintersect))
        tintersect=NaN(numel(x),1);
    end
end

% define function which will equal zero when curves intersect
% interpolate to find values between defined points
    function val = zerofun(xq)
        % val = interp1(x,y1,xq) - interp1(x,y2,xq);
        val = interp1(xint,y1int,xq) - interp1(xint,y2int,xq);
    end
end