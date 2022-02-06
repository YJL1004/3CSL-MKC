function D = SampleDensity_V3(KH , k)
[SNum , ~ , KerNum] = size(KH);
D = zeros(KerNum , SNum);

for ii = 1 : KerNum
    TKH = KH(:,:,ii);
    TR = sort(TKH(:) , 'descend');
    ThresholdSimilarityIndex = round(SNum * k);
    ThresholdSimilarity = TR(ThresholdSimilarityIndex^2);
    Index = TKH > ThresholdSimilarity;
    D(ii , :) = ones(1 , SNum) * Index;
end
% Normalization
% NormDColumn = repmat(diag(D*D') , 1 , SNum);
NormDColumn = repmat(sum(D , 2) , 1 , SNum);
D = D ./ NormDColumn;

% % ����ÿһ��kernel���ܶ����Ĳ���������ͳ������
% TCount = zeros(KerNum , round(SNum * k/2) );
% for ii = 1 : KerNum
%     [~ , hehe] = sort(D(ii , :) , 'descend');
%     TCount(ii , :) = hehe(1 : round(SNum*k/2));
% end

end