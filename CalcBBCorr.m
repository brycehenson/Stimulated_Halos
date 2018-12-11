function [g2_normed,corr_fit_params] =CalcBBCorr(halo_centered_cells,corr,isverbose)
%this function takes 
    %a cell array, with each cell containing the hits for a particular halo.
    %structure corr with yy,dx,dt
    %is verbose
%and returns
    %g2 
    %gausian fit params (Y,X,Z)X('amp','mu','sig','off')X(estimate,SE)

%TO BE IMPROVED
%the normalization needs tweaking/improving

%output gives NAN as last entry ?


%halo_centered_cells=halo_centered_cells(1:100);

updates=100; %number of updates in parfor

if isverbose
        disp('Calculating Correlation For Each Halo')
        update_frac=length(halo_centered_cells)/updates;
        parfor_progress_imp(updates)
        tic
end
    
G2_corr=zeros(length(halo_centered_cells),3,size(corr.yy,2)); %initialize the G2_corr
for n=1:length(halo_centered_cells)
    dd= halo_centered_cells{n};
    counts=size(dd,1);
    %check if the 
    if ~isempty(halo_centered_cells{n})
        t=dd(:,1)';
        x=dd(:,2)';
        y=dd(:,3)';
        %here i assign to a temp variable so that parfor is happy
        G2_temp=[]; %this is to make things parfor friendly
        G2_temp(1,:)=find_pairs_back_to_back_c(corr.dx,corr.dt,corr.yy,y,x,t)./counts;
        G2_temp(2,:)=find_pairs_back_to_back_c(corr.dx,corr.dt,corr.yy,x,t,y)./counts;
        G2_temp(3,:)=find_pairs_back_to_back_c(corr.dx,corr.dt,corr.yy,t,y,x)./counts;
        G2_corr(n,:,:)=G2_temp;
    end%is cell empty
    if isverbose
            %dont want the parfor updating every cycle because it can actualy
            %take a decent bit of time
            if rand<(1/update_frac)
                parfor_progress_imp;
            end
    end
end% loop over cells
clear G2_temp

if isverbose
    parfor_progress(0);
    toc
end

sum_G2_corr=squeeze(sum(G2_corr,1));

counts_per_uncorr=500;
progress_scaling=1;
if corr.norm
     g2_normed=[];
    

    % need to find the little g to norm the correlations
    G2_uncorr=[];
    dd=vertcat(halo_centered_cells{:});
    %now im going to shuffle the order of the counts
    %then if i split i will have a lower likelyhood of having correlated
    %atoms
    dd=dd(randperm(size(dd,1)),:);
    %now to split into roughly equal bins,
    %as the size is unlikley divisible
    full_size=size(dd,1);%number of counts in all the halos
    
    %if the number of counts is equal or too small define a manual bin
    %note that we start at zero as otherwise the non overlab when using the
    %edges wont work
    if floor(full_size/counts_per_uncorr)==0 || full_size==counts_per_uncorr
        edges=[0,full_size];
    else
    %we normaly define our edges by rouunding a linspace    
    edges=round(linspace(0,full_size,ceil(full_size/counts_per_uncorr)));
    end
    
    
    if isverbose
        disp(['calculating normalized correlations'])
        disp(['processing ',num2str(full_size),' counts in ',num2str(counts_per_uncorr),' count chunks'])
        parfor_progress(numel(edges)-1);
        tic
    end

    %then we itterate over these bins
    dd_split=[];
    counts=[];
    for n=1:(numel(edges)-1)
        
        %temporary dd for the loop
        dd_sub=dd(edges(n)+1:edges(n+1),:);
        counts=size(dd_sub,1);
        G2_uncorr_temp=[]; %make things parfor friendly
        t=dd_sub(:,1)';
        x=dd_sub(:,2)';
        y=dd_sub(:,3)';
        G2_uncorr_temp(1,:)=find_pairs_back_to_back_c(corr.dx,corr.dt,corr.yy,y,x,t)/counts;
        G2_uncorr_temp(2,:)=find_pairs_back_to_back_c(corr.dx,corr.dt,corr.yy,x,t,y)/counts;
        G2_uncorr_temp(3,:)=find_pairs_back_to_back_c(corr.dx,corr.dt,corr.yy,t,y,x)/counts;
        G2_uncorr(n,:,:)=G2_uncorr_temp;
        if isverbose
                parfor_progress;
        end
    end
    if isverbose
        parfor_progress(0);
        toc
    end
    
    sum_G2_uncorr=squeeze(sum(G2_uncorr,1));
    %divide by the combined correlation
    g2_normed=sum_G2_corr./sum_G2_uncorr(:,:);

else
    g2_normed=sum_G2_corr;
end



%clean the nans
% g2_normed(isnan( g2_normed))=0;

 corr_fit_params=[];
 
figure(3)
subplot(2,2,1)
plot(corr.yy,  g2_normed(3,:))
if corr.fit
     corr_fit_params(1,:,:)=corr_fit(corr.yy, g2_normed(3,:),0,0.004); 
end
xlabel('Sep in Y')
ylabel('Counts')
title('Y')
subplot(2,2,2)             
plot(corr.yy,  g2_normed(2,:))    
if corr.fit
     corr_fit_params(2,:,:)=corr_fit(corr.yy, g2_normed(2,:),0,0.002);
end
xlabel('Sep in X')
ylabel('Counts')
title('X')
subplot(2,2,3)
plot(corr.yy,  g2_normed(1,:))
if corr.fit
     corr_fit_params(3,:,:)=corr_fit(corr.yy, g2_normed(1,:),0,0.003);
end
xlabel('Sep in Z')
ylabel('Counts')
title('T')
set(gcf,'Position',[400 100 600 600])
set(gcf,'Color',[1 1 1]);


%do this in the main file



% if isverbose && find_correlation
%     if ~isempty(pairs_corr) %protects from the case that no pairs are fround
%         disp([int2str(sum(pairs_corr(1,:))),' pairs in t'])
%         disp([int2str(sum(pairs_corr(2,:))),' pairs in x'])
%         disp([int2str(sum(pairs_corr(3,:))),' pairs in y'])
%     end
% end




end

function fit_params=corr_fit(xdata,ydata,mu_guess,width_guess)
%this function is designed to fit a gaussian to the correlation plots and
%return fit params as estimate(amp,mu,c),SE


hold on
fo = statset('TolFun',10^-10,...
    'TolX',10^-10,...
    'MaxIter',10^6,...
    'UseParallel',0);
%amp*exp(damp*x1)*sin(2*pi*freq*x1 + phase)+offset'
%  'y~amp*exp((x1-mu)^2/(2*c^2)+1',...
fitobject=fitnlm(xdata,ydata,...
    'y~amp*exp(-1*((x1-mu)^2)/(2*sig^2))+off',...
   [4,mu_guess,width_guess,1],...
    'CoefficientNames',{'amp','mu','sig','off'},'Options',fo);
xvalues=linspace(min(xdata),max(xdata),300);
plot(xvalues,feval(fitobject,xvalues),'r')
%legend('data','fit')
hold off

fit_params=[fitobject.Coefficients.Estimate,fitobject.Coefficients.SE];



end

