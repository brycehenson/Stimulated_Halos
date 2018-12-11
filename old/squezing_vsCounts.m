%function out=squezing(halo_centered)

%this takes the centered halo data as cells for each shot of z,x,y 
%converts it to radial cord and then finds the normalized number vairance
%the traditional way (http://arxiv.org/pdf/1008.0845.pdf)of 
%showing this squezing data is to plot the disjointed pairs and the jointed
%pairs, this sees crude and arbitrary( the x axis has no information)
%instead i plot using the angle between the bin centers when this is pi
%then you expect decreased num var above it not



%TO FIX
%the double loop just after disp('Fiding Norm Var (rearanging)') takes the
%majority of the time, by doing some smarter matrix stuff it may be sped up

%bins that wrap across 0,2pi lose counts



% ref system
% bec's are at the poles
% azimuthal angle is arround the equator (zero to 2 pi)
% inclination is angle from the poles (which is zero inclination)(zero to pi)
% this was chosen to reduce the wrapping problems arround zero/2 pi 
% 
% START USER INPUT
% define bins


%romans halos
% radial_min=0.0219;%used to exclude crap in middle of circle
% radial_max=0.0311;
%mag sens halos
radial_min=0.0147;%used to exclude crap in middle of circle
radial_max=0.03;

mirror_azm=1;

%range_azm=[0 2]*pi;
range_azm=[0.25 0.75]*pi;
range_elev=[0.25 0.75]*pi;
steps_azm=4;
steps_elev=1;
%defining the width seprately allows for under or over sampling bins
bin_width_azm=range(range_azm)/steps_azm;
bin_width_elev=range(range_elev)/steps_elev;


%END USER INPUT

%remove empty bins
halo_centered_cells=halo_centered_cells(~cellfun('isempty',halo_centered_cells));

% create the bins
%bit hacky to create bins that dont wrap across 0,2pi
bin_centers_azm=wrapTo2Pi(linspace(range_azm(1),range_azm(2),steps_azm+1)-range(range_azm)/(2*steps_azm));
bin_centers_azm=bin_centers_azm(2:end);
bin_pairs_azm=transpose([wrapTo2Pi(bin_centers_azm-bin_width_azm/2) ; wrapTo2Pi(bin_centers_azm+bin_width_azm/2)]);
%bin_centers_azm/pi
%bin_pairs_azm/pi

if mirror_azm
    bin_centers_azm=[bin_centers_azm, bin_centers_azm+pi];
    bin_pairs_azm=[bin_pairs_azm;bin_pairs_azm+pi];
end

bin_centers_elev=wrapTo2Pi(linspace(range_elev(1),range_elev(2),steps_elev+1)-range(range_elev)/(2*steps_elev));
bin_centers_elev=bin_centers_elev(2:end);
bin_pairs_elev=transpose([wrapTo2Pi(bin_centers_elev-bin_width_elev/2) ; wrapTo2Pi(bin_centers_elev+bin_width_elev/2)]);
%bin_centers_elev/pi
%bin_pairs_elev/pi


% mask=mask_rad & halo_azm>window_azm(1) & halo_azm<window_azm(2);
% num_origin=sum(mask);

disp('binning')
halo_counts=[];
file_bin_counts={}; %the cell array that will store the counts in each bin for each file
halo_centered_spherical={};
for n=1:size(halo_centered_cells,1)
    angle_counts_file={};
    single_halo=halo_centered_cells{n};
    halo_radial=sqrt(single_halo(:,1).^2+single_halo(:,2).^2+single_halo(:,3).^2);
    halo_azm=atan2(single_halo(:,2),single_halo(:,3))+pi;
    halo_elev=acos(single_halo(:,1)./halo_radial);
    %halo_centered_spherical{n} =[halo_radial,halo_azm,halo_elev];
    %[halo_radial,halo_azm/pi,halo_elev/pi]
    mask_rad=halo_radial>radial_min & halo_radial<radial_max;
    halo_counts(n)=sum(mask_rad);
    %launch into the nested loop of bining
    for m=[1:size(bin_pairs_azm,1)]
        for p=[1:size(bin_pairs_elev,1)]
        counts=sum(mask_rad & halo_azm>bin_pairs_azm(m,1) & halo_azm<bin_pairs_azm(m,2)  & halo_elev>bin_pairs_elev(p,1) & halo_elev<bin_pairs_elev(p,2));
            angle_counts_file{m,p}=[n,bin_centers_azm(m),bin_centers_elev(p),counts];        
        end
    end
    file_bin_counts{n}=vertcat(angle_counts_file{:});
end


figure(7)
disp(['mean in halo after radial mask', num2str(mean(halo_counts))])
hist(halo_counts,100)
title('Radialy Cleaned Halo')
xlabel('Number of Counts')
ylabel('Number of Files')

count_bins_max=2000;
count_bins_min=3;
count_bins=30;
count_bins_edge=linspace(count_bins_min,count_bins_max,count_bins+1);
count_bins_cen=linspace(count_bins_min,count_bins_max,count_bins+1);
count_bins_cen=count_bins_cen(1:end-1)+count_bins_max/(count_bins*2);

%now sort them into the bin ranges
num_bined_file_bined_counts={}; %horrible name, this will contan the halo_bins for each count range


for n=[1:size(count_bins_cen,2)]
    counts_min=count_bins_edge(n)
    counts_max=count_bins_edge(n+1)
    num_bined_file_bined_counts{n}={file_bin_counts{logical(halo_counts<counts_max & halo_counts>counts_min)}};
    %return
end

clear mask_tem
clear counts_min
clear counts_max

%all possible combinations of bins is given by
index_combs=nchoosek(1:size(file_bin_counts{1},1),2);
disp(['number of combinations of bins ',int2str(size(index_combs,1))]);
bin_pairs=zeros(size(halo_centered_cells,1),size(index_combs,1),2);


%want to calculate the angle between the bin_centers 
%https://stackoverflow.com/questions/5188561/signed-angle-between-two-3d-vectors-with-same-origin-within-the-same-plane
%https://math.stackexchange.com/questions/243142/what-is-the-general-formula-for-calculating-dot-and-cross-products-in-spherical
%r1r2(sin?1sin?2cos(?1??2)+cos?1cos?2)
%acos(sin(elev1)sin(elev2)cos(azm1?azm2)+cos(elev1)cos(elev2))
angle_pairs=zeros(1,size(index_combs,1));
for n=[1:size(index_combs,1)]
    azm1=file_bin_counts{1}(index_combs(n,1),2);
    elev1=file_bin_counts{1}(index_combs(n,1),3);
    azm2=file_bin_counts{1}(index_combs(n,2),2);
    elev2=file_bin_counts{1}(index_combs(n,2),3);
    angle_pairs(n)= acos(sin(elev1)*sin(elev2)*cos(azm1-azm2)+cos(elev1)*cos(elev2));
end


mean_squeze_halo_counts=[];
min_squeze_halo_counts=[];
disp('Fiding Norm Var for counts bins')
for i=[1:size(count_bins_cen,2)]

    

    angle_counts_mat=[];
    for n=1:size(num_bined_file_bined_counts{i},2)
        angle_counts_mat(:,:,n)=num_bined_file_bined_counts{i}{n};
    end
    %this selects the appropriate bins so we have a matix of filex*index in
    %comb list*(bin1,bin2)
    bin_pairs=[angle_counts_mat(index_combs(:,1),4,:), angle_counts_mat(index_combs(:,2),4,:)];
    %need to permute these get desired format
    bin_pairs=permute(bin_pairs,[3 1 2]);


    %the above is little more that rearangeing values, below is where the norm
    %var is actualy calculated
    diffs=(bin_pairs(:,:,1)-bin_pairs(:,:,2));
    diffs2=(diffs).^2;
    norm_var=(mean(diffs2,1)-mean(diffs,1).^2)./(mean(bin_pairs(:,:,1))+mean(bin_pairs(:,:,2)));


    %now i want to sort through the angle and norm var values, average the
    %norm_var values that have the same angle_pairs
    angle_tol=0.0001;
    uniq_angles=uniquetol(angle_pairs,0.0001);
    angle_var_sd=zeros(size(uniq_angles,2),3);
    %here i also remove zeros
    for n=1:size(uniq_angles,2)
        matching=norm_var((abs(angle_pairs-uniq_angles(n))<angle_tol) & norm_var~=0);
        angle_var_sd(n,:)=[uniq_angles(n)/pi, mean(matching),std(matching)/sqrt(size(matching,2))];
    end



%     figure(10)
%     errorbar(angle_var_sd(:,1),angle_var_sd(:,2),angle_var_sd(:,3),'x')
%     title('Squeezing')
%     set(gcf,'Color',[1 1 1]);
%     xlabel('Angle Between Bins /Pi')
%     ylabel('Normalised Variance')
%     line([-2 2], [1 1],'Color','red');
% 
%     set(gca,'Xlim',[-0.1 1.1])

    %also plot the more traditional type of plot where the x axis is the bin
    %comb number, two types of points one for the opst and the other for all
    %others

    norm_var_opst=norm_var((abs(angle_pairs-pi)<angle_tol) & norm_var~=0);
    norm_var_rest=norm_var(~(abs(angle_pairs-pi)<angle_tol) & norm_var~=0);

    norm_var_avg=angle_var_sd(:,1);
 %   disp(['mean opst bin' ,num2str(angle_var_sd(angle_var_sd(:,1)==1,2)),'±',num2str(angle_var_sd(angle_var_sd(:,1)==1,3))])
 %   disp(['min opst bin ',num2str(min(norm_var(angle_pairs==pi)))])

    mean_squeze_halo_counts(i,:)=[angle_var_sd(angle_var_sd(:,1)==1,2),angle_var_sd(angle_var_sd(:,1)==1,3)];
    min_squeze_halo_counts(i)=min(norm_var(angle_pairs==pi));
%     figure(11)
%     plot(1:size(norm_var_opst,2),norm_var_opst,'x',1:size(norm_var_rest,2),norm_var_rest,'+')
%     xlabel('Bin Pair Number')
%     ylabel('Normalised Variance')
%     line([0 size(norm_var_rest,2)], [1 1],'Color','red');

end

size(count_bins_cen)
size(mean_squeze_halo_counts)
figure(8)
plot(count_bins_cen,mean_squeze_halo_counts(:,1)') 
errorbar(count_bins_cen,mean_squeze_halo_counts(:,1),mean_squeze_halo_counts(:,2))
 hold on
plot(count_bins_cen,min_squeze_halo_counts,'r')
 hold off
 
 title('')
xlabel('Number of Counts')
ylabel('Norm. Var. Opposite bins')


%end