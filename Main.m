%======================================90char=======================================================
%+++++++++++++++++++++++++++++++++++++++ABOUT+++++++++++++++++++++++++++++++++++++++++++++
%Halo analysis for investigating the spontanious/ stimulated transistion

%This program will import the halo data using HaloImportData then
%plot
%for everything inc the bec's (turn on with save_all_points)
%   the 3d dist
%   tof dist
%the Time integrated histogram
%the radial distributionsin the halo
%the halo 3d sidtribution (plot3d_halo)

%then bin will bin  the halos by the number of counts for each bin we calc
%then correlation function (turn on with find_correlation)
%correlation
%squezing
%radial width

%then plot the output of these as a function of number bin

%++++++++++++++++++++++++++++++++++++To Be Improved+++++++++++++++++++++++++++++++++++++++

%look at pos correction as it seems to be dancing arround a fair bit

%take radial fft

%radial correlations

%calculate center of mass of each halo, look for trend with number in halo

%should find the halo radius from the distance between bec and condensate
%and compare with fit to get the velocity right on


%sqz
    %use proper expresssion for the uncert
    
%move to arrays for xy limits

%move to array for spherical cordinates instead of sep variables

%-----------------------------------START user var----------------------------------------


%filepath='C:\Tmp\dld_output\mag sens halos low atom num\d';
%filepath='C:\Tmp\dld_output\mag sens halos low atom num\d';


use_txy=1;                  %if false will remake the txy_forc files
use_saved_halos=1;          %if false will remake the halo files (usefull for any change to the halo cut)

%import settings (will only change things if use_txy =0)
files.count_min=1000;       %min mumber of counts to read in 
files.rot_angle=0.64;       %rotation anle from the txy data 
files.do_pos_correction=1;	%find the halo pos from the ceneter of bragg orders (position correction)
                            %if false will not radialy mask
clean_radialy=0;            %only keep the halo values within a certian radius
            
%what to do with imported halos                                                     
isverbose=1;                %print progress ect
    progress_scaling=50;    %what fraction of the time the progress bar should update
find_sqz=1;
find_bb_correlation=0;         %find the correlation in x,y,z
    corr.norm=1;            %normalize the correlations
    corr.fit=0;             %fit the correlaton with a gaussian
    %params for correlations
    corr.yy=linspace(-0.005,0.005,50);    %define the bin size
    corr.dx=0.0005;
    corr.dy=0.004;                       %define the condition in other axes
    corr.dt=corr.dx;
find_cl_correlation=0;
    est_corr_width=0.5e-3;
    corr.window=[est_corr_width,est_corr_width,est_corr_width]/2;
    corr.redges=linspace(0,est_corr_width*4,40);
    corr.xedges=linspace(-est_corr_width*5,est_corr_width*5,40);
    corr.chunk_factor=100;

files.save_all_points=0;    %create array that off all the points usefull for intial investigation and plot TOF and 3d
plot_sph_dist=1;            %plot the spherical distibution hitograms( radial, azm ,elev)
    fit_rad_dist=0;         %fit the radial distribution with a gaussian +offset & linear or just est from avg and sd
    azm_bins=300;            %number of bins for histogram and azm bin
    radial_width_azm_dep=0; %plot the radial width as a function of azm angle
        radial_width_plots=0;%plot each azm bin and delay by 1s to show
    fit_sine_azm=0;          %fit a sine wave to the counts and width
plot3d_hist=0;              %plot a 3d histogram of the halo colapsed in time
plot2d_hist=0;              %plot a 2d histogram(image) of the halo colapsed in time
    plot2d_hist_binw=0.0001; %bin width for hist in meters
    plot2d_hist_gauss=0.000;%apply gaussian blur in meters, for off set to zero otherwise specify blur radius
    plot2d_hist_each_shot=0;%plot for each shot indiv. to look for phase grains, will save
plot3d_halo=1;              %display all the combined halos as a 3d plot
plot_counts_dist=0;         %plot a histogram of the counts in the halo

movies3d=0;                 %make movies of the 3d plots


split_by_halocounts=0;
    halocounts_min=3;
    halocounts_max=1000;%375;
    halocounts_bins=5;
    plots_for_each_bin=1;
  
files.path='V:\TDC_user\ProgramFiles\my_read_tdc_gui_v1.0.1\dld_output\d'; 
files.numstart=1;           %start file num
files.numtoimport=3500;       %number of files to import

windows.bec.tmin=2.442;
windows.bec.tmax=2.445;

windows.bragg.tmin=2.455;
windows.bragg.tmax=2.458;

windows.all.tmin=2.44;
windows.all.tmax=2.46;

windows.halo.tmin=2.446;%.608
windows.halo.tmax=2.454;
windows.halo.rmin=0.0165;
windows.halo.rmax=0.02725;

windows.reflections=0; %allows reflections throguh the mask
windows.halo.reflecrmin=0.010;


windows.all.xmin=-0.03;
windows.all.xmax=0.035;
windows.all.ymin=-0.03;
windows.all.ymax=0.035;

files.velocity=9.8*0.430;    
    
% files.path='C:\Tmp\HOM v3\d'; 
% files.numstart=1;           %start file num
% files.numtoimport=8500;       %number of files to import 
% windows.bragg.tmin=4.812;
% windows.bragg.tmax=4.815;
% windows.all.tmin=4.822;
% windows.all.tmax=4.824; 
% windows.halo.tmin=4.8147;%.608
% windows.halo.tmax=4.8217;
% windows.halo.rmin=0.018;
% windows.halo.rmax=0.022;
% windows.reflections=0; %allows reflections throguh the mask
% windows.halo.reflecrmin=0.010;
% windows.all.xmin=-0.03;
% windows.all.xmax=0.035;
% windows.all.ymin=-0.03;
% windows.all.ymax=0.035;
% files.velocity=9.8*0.430;
  
%files.path='C:\Tmp\bragg 200us with mirror\d';
% files.path='C:\Tmp\halo with bs amp scan 3\d'; 
% files.numstart=1;           %start file num
% files.numtoimport=100;       %number of files to import
% 
% windows.bragg.tmin=4.812;
% windows.bragg.tmax=4.815;
% 
% windows.bec.tmin=4.822;
% windows.bec.tmax=4.824;
% 
% windows.halo.tmin=4.8145;%.608
% windows.halo.tmax=4.8217;
% windows.halo.rmin=0.0175;
% windows.halo.rmax=0.022;
% 
% windows.reflections=1; %allows reflections throguh the mask
% windows.halo.reflecrmin=0.007;
% 
% windows.all.xmin=-0.03;
% windows.all.xmax=0.035;
% windows.all.ymin=-0.03;
% windows.all.ymax=0.035;

% %for mag insensitve halos
%files.path='C:\Tmp\dld_output\mag insenstive pop scan\d';
%files.path='C:\Tmp\bragg 600us halo\d';
% files.path='C:\Tmp\bragg 200us with mirror\d'
% files.numstart=1;           %start file num
% files.numtoimport=200;       %number of files to import
% 
% windows.bragg.tmin=4.812;
% windows.bragg.tmax=4.815;
% 
% windows.bec.tmin=4.822;
% windows.bec.tmax=4.824;
% 
% windows.halo.tmin=4.815;%.608
% windows.halo.tmax=4.8204;
% windows.halo.rmin=0.0177;
% windows.halo.rmax=0.025;
% 
% windows.all.xmin=-0.03;
% windows.all.xmax=0.035;
% windows.all.ymin=-0.03;
% windows.all.ymax=0.035;
    
    
%for mag sensitive halo low num
%  files.path='C:\Tmp\dld_output\mag sens halos low atom num\d';
%  files.numstart=1;           %start file num
%  files.numtoimport=19000;       %number of files to import
% %for mag sensitive halo vary pop
% % files.path='C:\Tmp\dld_output\Mag Sens Halos Varying Num\run 3\d';
% % files.numstart=1;           %start file num
% % files.numtoimport=14000;       %number of files to import
% windows.bragg.tmin=0.60;
% windows.bragg.tmax=.604;
% 
% windows.bec.tmin=.612;
% windows.bec.tmax=.62;
% 
% windows.halo.tmin=.606;%.608
% windows.halo.tmax=.611;
% windows.halo.rmin=0.015;
% windows.halo.rmax=0.025;
% 
% windows.all.xmin=-0.025;
% windows.all.xmax=0.025;
% windows.all.ymin=-0.025;
% windows.all.ymax=0.025;

%for multi bragg
% windows.bragg.tmin=0.215;
% windows.bragg.tmax=.22;
% 
% windows.bec.tmin=.23;
% windows.bec.tmax=.235;
% 
% windows.halo.tmin=.221;%.608
% windows.halo.tmax=.229;
% 
% windows.all.xmin=-0.04;
% windows.all.xmax=0.03;
% windows.all.ymin=-0.035;
% windows.all.ymax=0.035;

%should be 2*9.8*0.6
%files.velocity=sqrt(2*9.8*0.7);%*50;


files.path='V:\TDC_user\ProgramFiles\my_read_tdc_gui_v1.0.1\dld_output\d';
files.numstart=1;           %start file num
files.numtoimport=560;       %number of files to import

windows.bragg.tmin=0.403;
windows.bragg.tmax=0.4071;

windows.bec.tmin=0.416;
windows.bec.tmax=0.420;

windows.halo.tmin=0.407;%.608
windows.halo.tmax=0.416;
windows.halo.rmin=0.023;
windows.halo.rmax=0.031;

windows.reflections=0; %allows reflections throguh the mask
windows.halo.reflecrmin=0.007;

windows.all.xmin=-0.035;
windows.all.xmax=0.035;
windows.all.ymin=-0.035;
windows.all.ymax=0.035;


%for the all points the time windows can be different, this is usefull for
%initialy finding where the bec is in time
tmin_allpoints=windows.bragg.tmin;
tmax_allpoints=windows.bec.tmax;



%-------------------------------------END user var----------------------------------------

tic

%initialize
close all;
halo_centered=cell(files.numtoimport,1);

if files.save_all_points
    all_points=cell(files.numtoimport,1);
    %as i do not save all the counts when saving the halo will have to
    %reprocess, this is fine as the use case is rare
    use_saved_halos=0;
end

%import the halo data
%[halo_centered_cells,halo_centered,bec_bragg,all_points]
[halo_centered_cells,halo_centered,bec_bragg,all_points]=HaloReflecImportData(files,windows,use_txy,use_saved_halos,1,progress_scaling);

%if all the points were saved then a TOF and 3d plot will be done
if files.save_all_points
    all_points=vertcat(all_points{:}); 
    
    %mask the dat 
    figure(4)
    mask=all_points(:,1)>tmin_allpoints & all_points(:,1)<tmax_allpoints;
    all_points=all_points(mask,:);
    %plot TOF
    hist(all_points(:,1),10000)
    %set(gcf,'Position',[400 100 600 600])
    set(gcf,'Color',[1 1 1]);
    
    %here i multiply by velocity for time plot comment this out
    all_points=all_points.*repmat([files.velocity 1 1],size(all_points,1),1);
    
    %plot the halo and bec
    figure(5)
    scatter3(all_points(:,1),all_points(:,2),all_points(:,3),1,'k.','SizeData',1)
    axis equal;
    axis vis3d;
    set(gcf,'Position',[200 50 800 600])
    set(gcf,'Color',[1 1 1]);
    xlabel('time*vel')
    ylabel('X(m)')
    zlabel('y(m)')
    if movies3d
        OptionZ.FrameRate=15;OptionZ.Duration=25;OptionZ.Periodic=true;
        CaptureFigVid([-20,10;-110,10;-190,80;-290,10;-380,10], 'halo_and_bec',OptionZ)
    end
end


%plots for all the count in the halo (all shots combined)
%plot the halo
if plot3d_halo
    figure(1)
    scatter3(halo_centered(:,1),halo_centered(:,2),halo_centered(:,3),1,'k.')
    %set(gcf,'Position',[200 50 800 600])
    set(gcf,'Color',[1 1 1]);
    axis vis3d;
    axis equal;
    xlabel('time*vel (m)')
    ylabel('X(m)')
    zlabel('y(m)')
    if movies3d
        OptionZ.FrameRate=15;OptionZ.Duration=25;OptionZ.Periodic=true;
        CaptureFigVid([-20,10;-110,10;-190,80;-290,10;-380,10], 'halo',OptionZ)
    end
end


if plot3d_hist
    figure(2)
    %hist3(halo_centered(:,2:3),[100 100])
    hist3(halo_centered(:,2:3),[100 100])
    set(get(gca,'child'),'FaceColor','interp','CDataMode','auto');
    %set(gcf,'Position',[400 100 600 600])
    set(gcf,'Color',[1 1 1]);
    xlabel('X(m)')
    ylabel('Y(m)')
    zlabel('Counts/Files')
end




if plot2d_hist==1 && plot2d_hist_each_shot==0
    figure(11)
    
    XEdges=windows.all.xmin:plot2d_hist_binw:windows.all.xmax;
    YEdges=windows.all.ymin:plot2d_hist_binw:windows.all.ymax;

    [counts,centers]=hist3(halo_centered(:,2:3),'edges',{XEdges,YEdges});
    if  ~plot2d_hist_gauss==0
        filter=fspecial('gaussian',round(10*plot2d_hist_gauss/plot2d_hist_binw),plot2d_hist_gauss/plot2d_hist_binw);
    	counts=imfilter(counts, filter, 'replicate');
    end
    %imagesc seems to plot the wrong way round so we transpose here
    imagesc(centers{1},centers{2},transpose(counts))
    set(gcf,'Color',[1 1 1]);
    title(['Halo Count Dist Combined Blur=',num2str(plot2d_hist_gauss)])
    
elseif plot2d_hist==1 && plot2d_hist_each_shot==1
    %calc the filter kernel out of the loop
    if  ~plot2d_hist_gauss==0
        filter=fspecial('gaussian',round(10*plot2d_hist_gauss/plot2d_hist_binw),plot2d_hist_gauss/plot2d_hist_binw);
    end
    XEdges=windows.all.xmin:plot2d_hist_binw:windows.all.xmax;
    YEdges=windows.all.ymin:plot2d_hist_binw:windows.all.ymax;
    figure(11)
    for n=1:length(halo_centered_cells)
        current_file_str = num2str(files.numstart+n-1);
        temp_halo=halo_centered_cells{n};

        if length(temp_halo)>300
            [counts,centers]=hist3(temp_halo(:,2:3),'edges',{XEdges,YEdges});
            if  ~plot2d_hist_gauss==0
                counts=imfilter(counts, filter, 'replicate');
            end
            imagesc(centers{1},centers{2},transpose(counts))
            set(gcf,'Color',[1 1 1]);
            title(['Halo Count Dist ',int2str(n),' Blur=',num2str(plot2d_hist_gauss),'m'])
            saveas(gcf,[files.path,current_file_str,'_2d_hist.jpg'])
        else
            disp(['file number ',current_file_str,' has too few counts to plot 2dhist'])
        end
    end
    clear temp_halo
    clear filter
end
clear XEdges
clear YEdges
clear counts
clear centers

            

%fing the number in each halo
halo_counts=cellfun(@(x) size(x,1),halo_centered_cells);

if plot_counts_dist
    figure(12)
    set(gcf,'Color',[1 1 1]);
    %find the size of the matrix in each cell
    hist(halo_counts,100)
    xlabel('Halo Counts')
    ylabel('Number of shots')
    title('Halo Count Dist')
end

%then for the corrrelation and the squezing we split into bins grouped by
%the number in the halo (radial masked)
if split_by_halocounts
	halo_centered_cells_count_bined={};
    disp('binning by num in halo')
  

    %define the bins
    count_bins_edge=linspace(halocounts_min,halocounts_max,halocounts_bins+1);
    count_bins_cen=linspace(halocounts_min,halocounts_max,halocounts_bins+1);
    count_bins_cen=count_bins_cen(1:end-1)+halocounts_max/(halocounts_bins*2);
    %now sort them into the bin ranges
    
    halo_centered_cells_count_bined={}; %this will contan the halo_bins for each count range
    halo_files_bined_counts=[];         
    for n=[1:size(count_bins_cen,2)]
        counts_min=count_bins_edge(n);
        counts_max=count_bins_edge(n+1);
        %select those files that are in the range
        mask=logical(halo_counts<counts_max & halo_counts>counts_min);
        halo_centered_cells_count_bined{n}=halo_centered_cells(mask);
        
    end
    clear mask
    clear counts_min
    clear counts_max
    verbose_corr=0;
else
    halo_centered_cells_count_bined={halo_centered_cells};
    count_bins_cen=[NaN];
    verbose_corr=1;
end



if find_bb_correlation   
    disp('Caclulating correlation for each bin')
    
    parfor_progress(size(halo_centered_cells_count_bined,2));
    corr_params=[];
    for n=1:size(halo_centered_cells_count_bined,2)
        [~,corr_params(n,:,:,:)]=CalcBBCorr(halo_centered_cells_count_bined{n}, corr, verbose_corr);
        saveas(gcf,[files.path,'_Correlations_Bin',num2str(count_bins_cen(n),'%5.0f'),'_counts.jpg'])
        saveas(gcf,[files.path,'_Correlations_Bin',num2str(count_bins_cen(n),'%5.0f'),'_counts.fig'])
        if verbose_corr==0
            parfor_progress;
        end
    end
    if verbose_corr==0
        parfor_progress(0);
    end
    if size(halo_centered_cells_count_bined,2)>2 && corr.fit
        figure(100)
        set(gcf,'Color',[1 1 1]);
        subplot(2,1,1)
        errorbar(count_bins_cen,corr_params(:,1,1,1),corr_params(:,1,1,2))
        hold on;
        errorbar(count_bins_cen,corr_params(:,2,1,1),corr_params(:,2,1,2),'r')
        errorbar(count_bins_cen,corr_params(:,3,1,1),corr_params(:,3,1,2),'g')
        hold off;
        xlabel('Counts')
        ylabel('Corr Amp')
        legend('Y','X','Z')
        
        subplot(2,1,2)
        errorbar(count_bins_cen,abs(corr_params(:,1,3,1)),corr_params(:,1,3,2))
        hold on;
        errorbar(count_bins_cen,abs(corr_params(:,2,3,1)),corr_params(:,2,3,2),'r')
        errorbar(count_bins_cen,abs(corr_params(:,3,3,1)),corr_params(:,3,3,2),'g')
        hold off;
        xlabel('Counts')
        ylabel('Corr Width (m)')
        legend('Y','X','Z')
        
    end
end

if find_cl_correlation
    %not set up to do with number bins
    CalcCorr(halo_centered_cells,corr)
    
end



if find_sqz
    disp('Caclulating sqz for bins')
    parfor_progress(size(halo_centered_cells_count_bined,2));
    sqz_params={};
    for n=1:size(halo_centered_cells_count_bined,2)
        sqz_params{n}=squezing(halo_centered_cells_count_bined{n},0);
        saveas(gcf,[files.path,'_Sqz_Bin_',num2str(count_bins_cen(n),'%5.0f'),'_counts.jpg'])
        saveas(gcf,[files.path,'_Sqz _Bin',num2str(count_bins_cen(n),'%5.0f'),'_counts.fig']) 
        parfor_progress;

    end
    parfor_progress(0);
    
    if size(halo_centered_cells_count_bined,2)>2
        figure(101)
        set(gcf,'Color',[1 1 1]);
        errorbar(count_bins_cen,cellfun(@(x) x{2}(1,1),sqz_params),cellfun(@(x) x{2}(1,2),sqz_params))
        xlabel('Counts')
        ylabel('Mean Norm Var (opst Bins)')
        title('Sqz')
        hold on
        plot(count_bins_cen,cellfun(@(x) x{2}(1,3),sqz_params),'r')
    end

end


if plot_sph_dist==1 
    disp('ploting radial distibution')
    rad_fit_params=[];
    for n=1:size(halo_centered_cells_count_bined,2)
        halo_centered_cells_count_bined_single=halo_centered_cells_count_bined{n};
        %here we plot the radial distribution of all halos combined
        halo_centered_cells_count_bined_single_comb=vertcat(halo_centered_cells_count_bined_single{:});
        bin_counts=length(halo_centered_cells_count_bined_single_comb);
        [halo_radial,halo_azm,halo_elev]=ConvToSph(halo_centered_cells_count_bined_single_comb);
        
        figure(6)
        set(gcf,'Color',[1 1 1]);
        %set(gcf,'Position',[400 100 600 600])
        
        subplot(2,2,1)
        hist(halo_elev/pi,100)
        title('Elevation')
        set(gcf,'Color',[1 1 1]);
        xlabel('Angle (Rad)/pi')
        ylabel('Counts')
        
        subplot(2,2,2)
        [azmcounts,azmedges]=histcounts(halo_azm/pi,azm_bins);
        %from https://au.mathworks.com/matlabcentral/answers/89845-how-do-i-create-a-vector-of-the-average-of-consecutive-elements-of-another-vector-without-using-a-l
        azmcenters=mean([azmedges(1:end-1);azmedges(2:end)]);
        plot(azmcenters,azmcounts)
        title('Azm Bin Counts')
        set(gcf,'Color',[1 1 1]);
        xlabel('Angle (Rad)/pi')
        ylabel('Counts')
        if fit_sine_azm
            disp('fit sine')
            %fit_params=fit_sine(xdata,ydata,amp_guess,phase_guess,freq_guess,isverbose)
            azm_counts_fit_params(n,:,:)=fit_azm_sine(azmcenters,azmcounts,range(azmcounts),mean(azmcounts),0,1,1);
        end
        
        subplot(2,2,3)
        %hist(halo_radial,100)
        [radcounts,radcenters]=hist(halo_radial,100);
        %from https://au.mathworks.com/matlabcentral/answers/89845-how-do-i-create-a-vector-of-the-average-of-consecutive-elements-of-another-vector-without-using-a-l
        plot(radcenters,radcounts)
        title('Com. Radial Width')
        set(gcf,'Color',[1 1 1]);
        xlabel('Distance From Cen. (m)')
        ylabel('Counts')
        if fit_rad_dist
            %coef*SE of {'amp','mu','sig','off','grad'}=rad_fit(xdata,ydata,mu_guess,width_guess,amp_guess)
            %estimate the fit params by assuming that it is gaussian like
            rad_fit_params(n,:,:)=rad_fit(radcenters,radcounts,mean(halo_radial),std(halo_radial),max(radcounts),1);
        end
        %want to calc the radial width as a function of the azm width
        %to start with will just calc the std radialy then maybe do a fit
        %if there are enough counts
        
        if radial_width_azm_dep
            halo_radial_std_azmbin=zeros(1,length(azmedges)-1);
            for m=1:(length(azmedges)-1)
                mask=halo_azm/pi>azmedges(m) & halo_azm/pi<azmedges(m+1);
                halo_radial_azmbin=halo_radial(mask);
                %can either find the sd of the radial for this azm bin or
                %better yet do a fit
                if fit_rad_dist  
                    [radcounts,radcenters]=hist(halo_radial_azmbin,100);
                    if radial_width_plots
                        figure(7)
                        plot(radcenters,radcounts,'+')   
                    end
                    %coef*SE of {'amp','mu','sig','off','grad'}=rad_fit(xdata,ydata,mu_guess,width_guess,amp_guess,isverbose)
                    %estimate the fit params by assuming that it is gaussian like
                    rad_fit_azmbin_params(n,m,:,:)=rad_fit(radcenters,radcounts,mean(halo_radial_azmbin),std(halo_radial_azmbin),max(radcounts),radial_width_plots);
                    if radial_width_plots
                        pause(0.1) 
                    end
                else
                   rad_fit_azmbin_params(n,m,:,1)=[max(radcounts),mean(halo_radial_azmbin),std(halo_radial_azmbin),0,0];
                    %should set params to best gueses from 
                    %rad_fit_azmbin_params(n,m,:,:)
                    %halo_radial_std_azmbin(m)=std(halo_radial_azmbin);
                end               

            end
            figure(6)
            subplot(2,2,4)
            errorbar(azmcenters,rad_fit_azmbin_params(n,:,3,1),rad_fit_azmbin_params(n,:,3,2))
            title('Azm. Bin Rad. Width')
            xlabel('Angle(rad)/pi')
            ylabel('Width')
            if fit_sine_azm
                disp('fit sine to rad width of azm bins')
                %fit_params=fit_sine(xdata,ydata,amp_guess,phase_guess,freq_guess,isverbose)
                azm_rad_fit_params(n,:,:)=fit_azm_sine(azmcenters,rad_fit_azmbin_params(n,:,3,1),range(rad_fit_azmbin_params(n,:,3,1)),mean(rad_fit_azmbin_params(n,:,3,1)),0,1,1);
            end
        end
        
        saveas(gcf,[files.path,'_RadDist_Bin_',num2str(count_bins_cen(n),'%5.0f'),'_counts.jpg'])
        saveas(gcf,[files.path,'_RadDist _Bin',num2str(count_bins_cen(n),'%5.0f'),'_counts.fig']) 
    end%loop over halo count bins
    
    clear halo_centered_cells_count_bined_single
    clear halo_centered_cells_count_bined_single_comb
    clear bin_counts
    clear mask
    if size(halo_centered_cells_count_bined,2)>2 && fit_rad_dist
        figure(8)
        set(gcf,'Color',[1 1 1]);
        subplot(2,1,1)
        
        errorbar(count_bins_cen,rad_fit_params(:,3,1),rad_fit_params(:,3,2),'x');
        xlabel('Counts')
        ylabel('Halo Radial Width (m)')
        title('Halo Width')
        subplot(2,1,2)
        fracwidth=rad_fit_params(:,3,1)./rad_fit_params(:,2,1);
        fracunc=fracwidth.*sqrt( (rad_fit_params(:,3,2)./rad_fit_params(:,3,1)).^2 +...
            (rad_fit_params(:,2,2)./rad_fit_params(:,2,1)).^2);
        errorbar(count_bins_cen,fracwidth,fracunc,'x');
        xlabel('Counts')
        ylabel('Halo Radial Width Fraction')
        
        clear fracunc
    end
    
end    %plot_sph_dist==1 
    
if isverbose
   disp('all done')
end
toc

%the halo width as frac of the halo radius
%rad_fit_params(:,3,1)./rad_fit_params(:,2,1)

%the bec width is given by
%bec_bragg(1,1,5) in x
%bec_bragg(1,1,6) in y
%mean(bec_bragg(:,1,5))
%mean(bec_bragg(:,1,6))

%as fraction of halo radius
%mean(bec_bragg(:,1,5))/rad_fit_params(:,2,1)
%mean(bec_bragg(:,1,6))/rad_fit_params(:,2,1)


