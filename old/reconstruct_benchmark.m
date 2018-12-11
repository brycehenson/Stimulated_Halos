itt=10;
tic;
for n=1:itt
out1=dld_read_5channels_reconst_multi_imp('C:\Tmp\dld_data1\multihalos_1',1,0,1,0);
end; 
run1=toc;
fprintf('%f\n',run1/itt)
tic;
 for n=1:itt
out2=dld_read_5channels_reconst_multi('C:\Tmp\dld_data1\multihalos_1',1,1,0);
end; 
run2=toc;
fprintf('%f\n',run2/itt)
fprintf('%f\n',run1/run2)
fprintf('%i\n',isequal(out1,out2))