function fit_params=fit_sine(xdata,ydata,amp_guess,off_guess,phase_guess,freq_guess,isverbose)
%this function is designed to fit a gaussian to the correlation plots and
%return fit params as estimate{'amp','freq','phase','off'},SE

fo = statset('TolFun',10^-40,...
    'TolX',10^-40,...
    'MaxIter',10^6,...
    'UseParallel',1);
%amp*exp(damp*x1)*sin(2*pi*freq*x1 + phase)+offset'
%'y~amp*(freq*x1+phase)+off+x1*grad',...
fitobject=NonLinearModel.fit(xdata,ydata,...
    'y~amp*sin(2*3.14159265359*freq*x1+phase)+off',...
   [amp_guess,freq_guess,phase_guess,off_guess],...
    'CoefficientNames',{'amp','freq','phase','off'},'Options',fo);
%y val fit [2.5*10^-3,0, 40*2*pi, 4*10^-3, -1]
%x val [4.9*10^-3,-0.02, 134.1, 4*10^-3, -3.9, 0]

xvalues=linspace(min(xdata),max(xdata),300);
if isverbose
    hold on
    plot(xvalues,feval(fitobject,xvalues),'r')
    %legend('data','fit')
    hold off
end

fit_params=[fitobject.Coefficients.Estimate,fitobject.Coefficients.SE];

end
