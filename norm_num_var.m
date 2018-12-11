%as given in http://arxiv.org/pdf/1008.0845.pdf
%Here i want to simulate what a given QE does to the norm varaiance for a
%perfectrly squezed state or an anti squezed state

num_bina=20;
num_binb=num_bina;

var_bina=0;%sqrt(num_bina+20^2);
var_binb=var_bina;



shots=1000;
effval=0.5;

bina=round(var_bina*randn(1,shots)+num_bina);
binb=round(var_binb*randn(1,shots)+num_binb);
%bina=sum(rand(bina,shots) > (1-effval));
%binb=sum(rand(binb,shots) > (1-effval));

%binb=binb-bina*anti_sqz;
diff2=(bina-binb).^2;
diff=(bina-binb);
Normvar=(mean(diff2)-mean(diff)^2)/(mean(bina)+mean(binb))

 
 
 
 
 
%  
%  num_bina=2;
% num_binb=20;
% shots=1000;
% anti_sqz=1.96;
% effvals=0:0.01:1;
% 
% for n=1:length(effvals)
% bina=sum(rand(num_bina,shots) > (1-effvals(n)));
% binb=sum(rand(num_binb,shots) > (1-effvals(n)));
% %binb=binb-bina*anti_sqz;
% diff2=(bina-binb).^2;
% diff=(bina-binb);
% 
% 
% Normvar(n)=(mean(diff2)-mean(diff)^2)/(mean(bina)+mean(binb));
%  end
%  plot(effvals,Normvar,effvals,-effvals+1)