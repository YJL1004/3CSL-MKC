clear
clc
addpath(genpath('.'))

data_set = 'bbcsport2view_Kmatrix';
repeate = 3;
alpha = 3;
beta = 0.125;
anchor = 2;
load(data_set);
cluster_num = length(unique(Y));
sample_num = length(Y);
ker_num = size(KH,3);

f_num = max(100, cluster_num*4);
KH = kcenter(KH);
KH = knorm(KH);

H = zeros(f_num, sample_num, ker_num);
opt.disp = 0;
for i=1:ker_num
    [Hi, ~] = eigs(KH(:,:,i), f_num, 'la', opt);
    H(:,:,i) = Hi';
end

res = zeros(repeate,5);
for retry = 1:repeate
    tic;
    [S, anchor_num, ~] = MKC_3C(H,Y,alpha,beta,anchor);
    t = toc;
    results = myNMIACC(S,Y,length(unique(Y)));
    res(retry,:) = [retry,results,t];
end
storge_file = fullfile('args/', [datestr(now,'dd-mmm-yyyy-HH-MM'),'-', data_set, '.mat']);
save(storge_file, 'res');

mean_res = mean(res(:,2:5),1);
std_res = std(res(:,2:5),0,1);
fprintf('@ [dataset:%s] @ \n', data_set);
fprintf('@ ACC:%3.4f(%3.4f) / NMI:%3.4f(%3.4f) / Precision:%3.4f(%3.4f) / Time:%6.2f(%6.2f) \n', mean_res(1),std_res(1),mean_res(2),std_res(2),mean_res(3),std_res(3),mean_res(4),std_res(4));