%data=[[1.1 2.2 3.6];[1.5 2.3 45];[7.8 9.2 5.4]];
data=vertcat(halo_centered_cells{1:40});
size(data,1)
rad_bins=linspace(50*10^-6,0.03,150);
tic
G2=CalcCorrRad(rad_bins,data);
rad_bins_cen=mean([rad_bins(1:end-1);rad_bins(2:end)]);
plot(rad_bins_cen,G2)
toc
