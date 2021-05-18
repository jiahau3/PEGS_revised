clear all;

frame='137';

id1=09;
cb1=0190;
% 
% id3=17;
% cb3=0509;
% 
% id5=21;
% cb5=0021;

% id7=03;
% cb7=0073;


disks = ['../s2/maskR0.50_mat/',frame,'.mat'];
load(disks);

diskAllcomb = ['../s2/maskR0.50_matAll/',frame,'.mat'];
load(diskAllcomb);
N = length(pres);
disk1 = presAll(id1,cb1);
% disk3 = presAll(id3,cb3);
% disk5 = presAll(id5,cb5);
% disk7 = presAll(id7,cb7);
% clearvars presAll;

pid2=19;
pcb2=01;
% 
% pid4=18;
% pcb4=01;

% pid6=18;
% pcb6=03;

eachcomb = ['../s2_eachDisk/maskR0.50_matAll/',frame,'p',sprintf('%02d',pid2),'.mat'];
load(eachcomb);
disk2 = presAll(pid2, pcb2);
% 
% eachcomb2 = ['../s2_eachDisk/maskR0.50_matAll/',frame,'p',sprintf('%02d',pid4),'.mat'];
% load(eachcomb2);
% disk4 = presAll(pid4, pcb4);

% eachcomb3 = ['../s2_eachDisk/maskR0.50_matAll/',frame,'p',sprintf('%02d',pid6),'.mat'];
% load(eachcomb3);
% disk6 = presAll(pid6, pcb6);

pres(id1) = disk1;
pres(pid2) = disk2;
% pres(id3) = disk3;
% pres(pid4) = disk4;
% pres(id5) = disk5;
% pres(pid6) = disk6;
% pres(id7) = disk7;

dirname = '../s2_picked/maskR0.50_mat/';
if ~exist(dirname, 'dir') mkdir(dirname); end
save(['../s2_picked/maskR0.50_mat/',frame,'.mat'],'pres');
