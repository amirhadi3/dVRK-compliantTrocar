function ati = readForceData()
ati = table2array(readtable('atiData.txt'));
ati = ati-mean(ati(1:500,:));
% ati(:,1) = -ati(:,1);
% ati(:,3) = -ati(:,3);
% ati(:,6) = -ati(:,6);
% ati(:,4) = -ati(:,4);
end