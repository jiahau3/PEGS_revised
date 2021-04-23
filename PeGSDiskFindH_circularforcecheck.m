function [ particle ] = PeGSDiskFindH_circularforcecheck( Rimg, Rlarge, SL, pxPerMeter, fsigma, Rcompensate, Rbound)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

%Disc detection for large particles happens here
% Rimg = Rimg > 50;
[centers, radii, metric] = imfindcircles(Rimg,Rlarge,'ObjectPolarity','dark','Method','TwoStage','Sensitivity',SL);
XoutL = find(centers(:,2) - Rbound < 0);
centers(XoutL,:) = [];
radii(XoutL,:) = [];
XoutR = find(centers(:,2) + Rbound > size(Rimg,1));
centers(XoutR,:) = [];
radii(XoutR,:) = [];
YoutT = find(centers(:,1) - Rbound < 0);
centers(YoutT,:) = [];
radii(YoutT,:) = [];
YoutB = find(centers(:,1) + Rbound > size(Rimg, 2));
centers(YoutB,:) = [];
radii(YoutB,:) = [];
Nlarge = length(centers) / 2;
particle(1:Nlarge) = struct('id',0,'x',0,'y',0,'r',0,'rm',0,'color','','fsigma',0,'z',0,'f',0,'g2',0,'forces',[],'betas',[],'alphas',[],'neighbours',[],'contactG2s',[],'forceImage',[]);
% display(Nlarge);
% display(centers);
%Bookkeeping
for n=1:Nlarge
    particle(n).id= n;
    particle(n).x = centers(n,1);
    particle(n).y = centers(n,2);
    particle(n).r = radii(n) + Rcompensate;
end

for n=1:Nlarge
    particle(n).rm = particle(n).r * pxPerMeter;
    particle(n).fsigma = fsigma;
    particle(n).color = 'r';
    end
end


