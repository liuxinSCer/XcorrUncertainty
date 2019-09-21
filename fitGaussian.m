function xparam=fitGaussian(svec, df, freqval)

sig0=0.025; % in Hz
[A0,indmax]=max(svec);

f0init=freqval(indmax);

    function fun=gauss(x0)
        %f0,sigma,A
        f0=x0(1);
        sigma = x0(2);
        A=x0(3);
        fun = A*exp(-(freqval-f0).^2/(2*sigma*sigma));
        fun=fun(:);
    end

misfit =  @(x0) norm( gauss(x0) - svec );
        
[xnew,fval] = fminunc(misfit,[f0init,sig0,A0]);
xparam=xnew;
end


