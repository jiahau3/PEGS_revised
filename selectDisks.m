clear all;

frame='216';
id1=23;
cb1=0885;

% id3=24;
% cb3=0004;
% id4=19;
% cb4=0226;

disks = ['../s2/maskR0.50_mat/',frame,'.mat'];
load(disks);

diskAllcomb = ['../s2/maskR0.50_matAll/',frame,'.mat'];
load(diskAllcomb);
N = length(pres);
disk1 = presAll(id1,cb1);
% disk3 = presAll(id3,cb3);
% disk4 = presAll(id4,cb4);
% clearvars presAll;

pid2=24;
pcb2=01;
% pid4=12;
% pcb4=34;
% 
eachcomb = ['../s2_eachDisk/maskR0.50_matAll/',frame,'p',sprintf('%02d',pid2),'.mat'];
load(eachcomb);
disk2 = presAll(pid2, pcb2);
% 
% eachcomb2 = ['../s2_eachDisk/maskR0.50_matAll/',frame,'p',sprintf('%02d',pid4),'.mat'];
% load(eachcomb2);
% disk4 = presAll(pid4, pcb4);

pres(id1) = disk1;
pres(pid2) = disk2;
% pres(id3) = disk3;
% pres(id4) = disk4;

dirname = '../s2_picked/maskR0.50_mat/';
if ~exist(dirname, 'dir') mkdir(dirname); end
save(['../s2_picked/maskR0.50_mat/',frame,'.mat'],'pres');
