function G2=CalcCorrRad(rad_bins,data)
%======================================90char=============================================
%+++++++++++++++++++++++++++++++++++++++ABOUT+++++++++++++++++++++++++++++++++++++++++++++
%there a number of ways to approach counting the number of points within a
%spcified radius

%simplest approach
%calc the distance from each point to every other
%prevent duplicates
%then for each point calc the number within that bin using the rad

%should look at how seaun does this

%data=vertcat(halo_centered_cells{1:50});
G2=zeros((length(rad_bins)-1),1);
%this loop should be written to write directly to a vector to remove the
%reshape step
rad=zeros(size(data,1));
for n=1:size(data,1)
    point1=data(n);
    for m=(n+1):(size(data,1))
        rad(n,m)=sqrt(sum((point1-data(m)).^2));
    end
end
rad=reshape(rad,[],1);
rad=rad(rad~=0);

for n=1:(length(rad_bins)-1)
    G2(n)=sum(rad<rad_bins(n+1) & rad>rad_bins(n));
end
% 
% for n=1:(length(rad_bins)-1)
%     for m=1:size(rad,1)
%         if rad(m)<rad_bins(n+1) && rad(m)>rad_bins(n)
%             G2(n)= G2(n)+1;
%         end
%     end
% end


end