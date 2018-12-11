function [halo_centered_cells,halo_centered,bec_bragg,all_points]=HaloImportData(files,windows,use_txy,use_saved_halos,isverbose,progress_scaling)
%======================================90char=============================================
%+++++++++++++++++++++++++++++++++++++++ABOUT+++++++++++++++++++++++++++++++++++++++++++++
%This function reads in txy data for files from files.numstart to files.numstart+files.numtoimport
%(will convert DLD data if TXY does not exist) checks if there are at least files.count_min
%counts in the file (returns empty otherwise). The bec and bragg is
%windowed and the average is used to find the center of the halo. The halo
%counts have this center subtracted and multiplied by the velocity in thetime axis to 
%convert to position. The halo counts are then radialy selected to recuce
%the background. These counts and the window/file parameters are then saved
%to path+_all_halos_saved.mat which can just be imported the next time.

%INPUTS


%OUTPUTS
%bec_bragg files*2*7 matrix which contians for each file contains bragg avg (t,x,y)
%sd(t x y) numcounts as flat vecotr and then the same for bec


%TO IMPROVE: 
%if files =0 import all
%should pass a vector of paths for multi folder import
%should import params file for scans

%measure widths and number in BEC, bragg

%================================ START USER INPUT=====================================
min_counts_in_bragg_or_BEC=5;

%===================================END USER INPUT========================================

%check if saved halo data exists
if use_saved_halos && ~fileExists([files.path,'_all_halos_saved.mat'])
    use_saved_halos=0;
    disp('Halo data does not exist.')
end

if use_saved_halos
    %we will load the data and check that the same winsdows and files were
    %used to make it as what were requested
    
    %TO IMPROVE: this negects the case when the requested files are a
    %small subset of what has been converted. Could just always import
    %everything and then select after import...
    
    windows_req=windows;
    files_req=files;
    load([files.path,'_all_halos_saved.mat'],'halo_centered_cells','windows','files','bec_bragg')
    %compare
    if ~isequal(windows,windows_req) || ~isequal(files,files_req)
        if isverbose
            disp('Saved data windows/files does not agree with requested.')
            disp('Will reprocess.') 
        end
        use_saved_halos=0;
        files=files_req;
        windows=windows_req;
        
    end
        
    %clean up
    clear windows_req
    clear files_req
end
    
if ~use_saved_halos
%initalize variables
halo_centered_cells={};
%initialize the parfor compatable counters    
lowcountfiles=zeros(files.numtoimport,1);
missingfiles=zeros(files.numtoimport,1);
filestotxy=zeros(files.numtoimport,1);
bec_or_bragg_zero=zeros(files.numtoimport,1);
importokfiles=zeros(files.numtoimport,1);

if isverbose
    disp('Processing Data')
    parfor_progress(round(files.numtoimport/progress_scaling));
end

    bec_bragg=zeros(files.numtoimport,2,7);
    for n=1:files.numtoimport
        current_file_str = num2str(files.numstart+n-1);

        %check if the file exists
        if ~fileExists([files.path,current_file_str,'.txt'])
            missingfiles(n)=1;
        else
            %if the txy_forc does not exist or the use_txy flag is low (re) make it
            if ~fileExists([files.path,'_txy_forc',current_file_str,'.txt']) || ~use_txy
                filestotxy(n)=1;
                dld_raw_to_txy(files.path,files.numstart+n-1,files.numstart+n-1);
            end
            
            %three_channel_output=importdata([files.path,'_txy_forc',current_file_str,'.txt'],',');
            %three_channel_output=load('-ascii',[files.path,'_txy_forc',current_file_str,'.txt']);
            three_channel_output=txy_importer(files.path,current_file_str);

            if size(three_channel_output,1)<files.count_min
                %disp('selected shot has too few counts')
                lowcountfiles(n)=1;
                
            else
%                 if sum(sum(isnan(three_channel_output)))>0
%                     disp('NAN values found')
%                     return
%                 end

                three_channel_output_rot=zeros(size(three_channel_output));
                %rotate the spatial cord
                sin_theta = sin(files.rot_angle);
                cos_theta = cos(files.rot_angle);
                three_channel_output_rot(:,1) = three_channel_output(:,1);
                three_channel_output_rot(:,2) = three_channel_output(:,2)*cos_theta...
                    - three_channel_output(:,3)*sin_theta;
                three_channel_output_rot(:,3) = three_channel_output(:,2)*sin_theta...
                    + three_channel_output(:,3)*cos_theta;

                
                %which counts in the file are within the xy limits
                mask_XY=three_channel_output_rot(:,2)>windows.all.xmin &...
                three_channel_output_rot(:,2)<windows.all.xmax & ... 
                three_channel_output_rot(:,3)>windows.all.ymin & ...
                three_channel_output_rot(:,3)<windows.all.ymax;
                
                %save all points is mainly just usefull for a diagnostic
                if files.save_all_points
                    all_points{n}=three_channel_output_rot(mask_XY,:); 
                end

                %windows in time for bragg and BEC
                %TO IMPROVE: little point finding whats mask and bragg for
                %pos correction off,but this is a rare use case..
                mask_BEC=three_channel_output_rot(:,1)>windows.bec.tmin &...
                three_channel_output_rot(:,1)<windows.bec.tmax &...
                mask_XY;

                mask_bragg=three_channel_output_rot(:,1)>windows.bragg.tmin &...
                    three_channel_output_rot(:,1)<windows.bragg.tmax &... 
                    mask_XY;
                %if the bragg or bec is empty give up and set
                %halo_centered_cells_temp to empty
                if  (sum(mask_bragg)<min_counts_in_bragg_or_BEC || sum(mask_BEC)<min_counts_in_bragg_or_BEC) ...
                        && files.do_pos_correction
                    bec_or_bragg_zero(n)=1;
                    %write a halo file to prevent it getting here again
                    halo_centered_temp=[];
                else
                    three_channel_output_rot_bragg=three_channel_output_rot(mask_bragg,:);
                    three_channel_output_rot_BEC=three_channel_output_rot(mask_BEC,:);
                    %find the center of the bragg and bec and then subtract
                    %this position off the couns in the halo
                    if files.do_pos_correction
                        mean_bec=mean(three_channel_output_rot_BEC);
                        mean_bragg=mean(three_channel_output_rot_bragg);
                        bec_bragg(n,:,:)=[[mean_bec,std(three_channel_output_rot_BEC),sum(mask_BEC)];[mean_bragg,std(three_channel_output_rot_bragg),sum(mask_bragg)]];
                        halo_center=(mean_bragg+mean_bec)/2;
                    else
                        %same as not subtracting pos
                         halo_center=[0 0 0];
                    end
                    
                    %mask the halo
                    mask_halo=three_channel_output_rot(:,1)>windows.halo.tmin &...
                        three_channel_output_rot(:,1)<windows.halo.tmax &...
                        mask_XY;
                    three_channel_output_rot_halo=three_channel_output_rot(mask_halo,:);

                    %subtract the t,x,y center from every row
                    %repmat is apparently the fastest http://au.mathworks.com/matlabcentral/newsreader/view_thread/267082
                    %then convert the time axis into z pos by multiplying time by
                    %files.velocity
                    halo_centered_temp=(three_channel_output_rot_halo-repmat(halo_center,size(three_channel_output_rot_halo,1),1)).*...
                        repmat([files.velocity 1 1],size(three_channel_output_rot_halo,1),1);
                    
                    %here i mask radialy by computing the radius and seeing if it
                    %is within rmax,windows.halo.rmin
                    %will not do if files.do_pos_correction=0 as the halo is not zero
                    %centered
                    importokfiles(n)=1;
                    if files.do_pos_correction
                        radial_temp=sqrt(halo_centered_temp(:,1).^2+halo_centered_temp(:,2).^2+halo_centered_temp(:,3).^2);
                        radial_mask=radial_temp<windows.halo.rmax & radial_temp>windows.halo.rmin;
                        halo_centered_temp=halo_centered_temp(radial_mask,:);
                    end
                end
                
                %save the halo counts
                %use cells so that the output can be ragged
                halo_centered_cells{n}=halo_centered_temp;
                
                %clean up variables
                

            end %file size condition
        end %file exists condition
        
        if isverbose
            %dont want the parfor updating every cycle because it can actualy
            %take a decent bit of time
            if rand(1)<(1/progress_scaling)
            parfor_progress;
            end
        end
    end %file loop
    
%save the imported data
%TO IMPROVE: need to save more data such as count

save([files.path,'_all_halos_saved.mat'],'halo_centered_cells','windows','files','bec_bragg','-v7.3')

%clean up, this is only important if the loop is serial not parfor
clear('mask_BEC','mask_bragg','mask_XY','mask_halo','mean_bec','mean_bragg')
clear('three_channel_output','three_channel_output_rot_BEC','three_channel_output_rot_bragg')
clear('three_channel_output_rot_halo','halo_centered_temp','radial_temp','three_channel_output_rot')



%now i sum up all the parfor compatable counters
lowcountfiles=sum(lowcountfiles);
missingfiles=sum(missingfiles);
filestotxy=sum(filestotxy);
bec_or_bragg_zero=sum(bec_or_bragg_zero);
importokfiles=sum(importokfiles);

if isverbose
        parfor_progress(0);
        disp('Import Done')
end


end%use precomputed halo data


%halo_centered_table=cell2table(halo_centered_cells);
%writetable( halo_centered_table,'tabledata.dat')

halo_centered=vertcat(halo_centered_cells{:});


%jsut a dummy output to keep matlab happy
if ~files.save_all_points
    all_points=[];
end

if isverbose 
    if use_saved_halos
        disp('data loaded from all_halos_saved')
        missingfiles=0;
        bec_or_bragg_zero=0;
    else
        disp([int2str(lowcountfiles),' files with low counts'])
        disp([int2str(missingfiles),' files missing'])
        disp([int2str(filestotxy),' files converted to txy'])
        disp([int2str(bec_or_bragg_zero),' files with no counts in bec or bragg'])
        disp([int2str(importokfiles),' imported ok'])
        
    end
    disp([int2str(size(halo_centered,1)),' total counts in halo'])
    disp([num2str(size(halo_centered,1)/(files.numtoimport-(missingfiles+bec_or_bragg_zero))),' counts per file in halo'])
    disp(['mean x,y,t pos halo ',num2str(mean(halo_centered(:,2))),' , ', ...
         num2str(mean(halo_centered(:,3))),' , ', num2str(mean(halo_centered(:,1)))])
    
end


end