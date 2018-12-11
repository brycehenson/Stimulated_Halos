%%% Use this to convert raw files from DLD to txy format
%%% which is used by C++ program such as
%%% g2_calc_norm_across_files_spatial_dy.exe to compute g2
%%% there is a 1:1 correspondance between input files and 
%%% output files with '_txy_forc'


function[] = dld_raw_to_txy(filename_raw,startfile,end_file)
%tic

missingfiles=zeros(end_file-startfile,1);
alreadytxy=zeros(end_file-startfile,1);
convertedfiles=zeros(end_file-startfile,1);
parfor_progress(end_file-startfile)
parfor i=1:end_file-startfile
    file_no = num2str(i+startfile);
    filename_read = [filename_raw,file_no];
    if fileExists([filename_raw,'_txy_forc',file_no,'.txt'])
        alreadytxy(i)=1;   
    elseif ~fileExists([filename_read,'.txt'])
        missingfiles(i)=1;
    else
        convertedfiles(i)=1;
        [hits_sorted] = dld_read_5channels_reconst_multi_imp(filename_read,1,0,1,0);
        filename_write = [filename_raw,'_txy_forc',file_no,'.txt'];
        dlmwrite(filename_write,hits_sorted,'precision',8);
        
    end
    parfor_progress;
end
parfor_progress(0);
missingfiles=sum(missingfiles);
alreadytxy=sum(alreadytxy);
convertedfiles=sum(convertedfiles);
fprintf('files missing %i \n',missingfiles)
fprintf('files already txy %i \n',alreadytxy)
fprintf('files converted %i \n',convertedfiles)