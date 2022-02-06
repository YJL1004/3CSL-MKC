function data_save_print1(dataset, res,paras)


tm=datestr(now,'yyyy-mm-dd_HH_MM_SS');

ACC = res(1)*100;
MIhat = res(2)*100;
purity = res(3)*100;
TIME = res(4);
fprintf('@ %s \n', dataset);
fprintf('@ alpha:%5.4f / beta:%5.4f/ anchor:%d/ \n', paras.alpha, paras.beta, paras.anchor_num);
fprintf('@ ACC:%5.2f / NMI:%5.2f / Purity:%5.2f / \n', ACC,MIhat,purity);

fileID = fopen(['log/',dataset,'_RES.txt'],'a');
fprintf(fileID,'%s\n',tm);
fprintf(fileID,'@ alpha:%5.4f / beta:%5.4f/ anchor:%d/ \n', paras.alpha, paras.beta, paras.anchor_num);
fprintf(fileID,'ACC:%5.2f NMI:%5.2f Purity:%5.2f Time:%5.2f\n',ACC,MIhat,purity,TIME);
fprintf(fileID,'\n');
fclose(fileID);