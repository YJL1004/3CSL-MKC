function [S,anchor_num,obj] = MKC_3C(H,y,alpha,beta,anchor)
[~,sample_num,ker_num] = size(H);
cluster_num = length(unique(y));

if ~exist('anchor','var')
    anchor = 2;
end
p_num1 = max(50, cluster_num*2);
p_num2 = max(100, cluster_num*4);
parameters = [cluster_num, p_num1, p_num2];

anchor_num = parameters(anchor);
f_num = cluster_num * 4;

Si = zeros(sample_num, anchor_num, ker_num);
E = zeros(anchor_num, cluster_num, ker_num);
F = zeros(sample_num, cluster_num, ker_num);

for i=1:ker_num
    tSi = H(1:anchor_num,:,i); 
    tSi(tSi < 0) = 0; tSi(tSi > 1) = 1; 
    Si(:,:,i) = tSi';
end
S = mean(Si,3);
H = H(1:cluster_num*4,:,:);
%[S, ~] = eigs(mean(K, 3), anchor_num, 'la', opt); S(S < 0) = 0; S(S > 1) = 1;
[~, tB] = kmeans(mean(H, 3)', cluster_num); B = tB';
P = randn(sample_num, anchor_num); %P = orth(P);

mit = 20;
obj = [];
for iter = 1:mit
    %For P OK OK
    tP = zeros(sample_num, anchor_num);
    for i = 1 : ker_num
        tP = tP + H(:,:,i)'*(H(:,:,i)*Si(:,:,i)+beta*B*E(:,:,i)');
    end
    [U, ~, V] = svd(tP, 'econ');
    P1 = U * V';
    p_diff = norm(P1-P, 'fro');
    P = P1;
    clear tP;clear U;clear V;clear P1;
    
    %For B OK OK
    tB = zeros(f_num, cluster_num);
    for i = 1 : ker_num
        tB = tB + H(:,:,i)*(F(:,:,i)+P*E(:,:,i));
    end
    [U, ~, V] = svd(tB, 'econ');
    B = U * V';
    clear tB;clear U;clear V;
    
    %For Si OK OK
    for i = 1:ker_num
        A =H(:,:,i)'*(H(:,:,i)*P);
        tSi = (2 * alpha * S + A) / (2*alpha);
        tSi(tSi > 1) = 1; tSi(tSi < 0) = 0;
        %tSi = NormalizeFea(tSi);
        Si(:,:,i) = tSi;
        clear tSi;clear A;
    end
    
    
    %For Fi  OK OK
    for i = 1 : ker_num
        tFI = B'*H(:,:,i);
        [U, ~, V] = svd(tFI, 'econ');
        F(:,:,i) = V*U';
        F(:,:,i) = 0.5 * (F(:,:,i) + abs(F(:,:,i)));
        clear tFI;clear U;clear V;
    end
    
    
    %For Ei OK OK
    for i = 1 : ker_num
        tEI = B'*H(:,:,i)*P;
        [U, ~, V] = svd(tEI, 'econ');
        E(:,:,i) = V*U';
        E(:,:,i) = 0.5 * (E(:,:,i) + abs(E(:,:,i)));
        clear tEI;clear U;clear V;
    end
    
    %For S OK    
    S = mean(Si, 3);
    S(S < 0) = 0; S(S > 1) = 1;

    obj(iter) = p_diff;
    if iter > 2 && p_diff < 10^(-3)
        break;
    end  
%     imagesc(E(:,:,1));
%     pause(0.5);
end


% B = rand(f_num, cluster_num);
% B = orth(B);
%P =ones(sample_num, anchor_num)/sample_num;
%S = zeros(sample_num, anchor_num);

% Si = zeros(sample_num, anchor_num, ker_num);
% for i = 1:ker_num
%     A =H(:,:,i)'*(H(:,:,i)*P);
%     tSi = (2 * alpha * S + A) / (2*alpha);
%     tSi(tSi > 1) = 1;
%     tSi(tSi < 0) = 0;
%     Si(:,:,i) = tSi;
%     clear tSi; clear A;
% end
%S = mean(Si, 3);







    %         for i = 1 : sample_num
    %             S(i,:) = projection(S(i,:)',1)';
    %         end
