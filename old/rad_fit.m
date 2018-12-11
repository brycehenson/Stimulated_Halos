function fit_params=rad_fit(xdata,ydata,mu_guess,width_guess,amp_guess,isverbose)
%this function is designed to fit a gaussian to the correlation plots and
%return fit params as estimate('amp','mu','sig','off','grad'),SE



fo = statset('TolFun',10^-6,...
    'TolX',10^-6,...
    'MaxIter',10^6,...
    'UseParallel',0);
%amp*exp(damp*x1)*sin(2*pi*freq*x1 + phase)+offset'
%  'y~amp*exp((x1-mu)^2/(2*c^2)+1',...
fitobject=fitnlm(xdata,ydata,...
    'y~amp*exp(-1*((x1-mu)^2)/(2*sig^2))+off+x1*grad',...
   [amp_guess,mu_guess,width_guess,0,0],...
    'CoefficientNames',{'amp','mu','sig','off','grad'},'Options',fo);
xvalues=linspace(min(xdata),max(xdata),300);
if isverbose
    hold on
    plot(xvalues,feval(fitobject,xvalues),'r')
    %legend('data','fit')
    hold off
end

fit_params=[fitobject.Coefficients.Estimate,fitobject.Coefficients.SE];



end
