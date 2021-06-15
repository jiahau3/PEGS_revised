function [SL] = PeGSHoughDebug_cali(Rimg, Rlarge, SL, DS, NlargeExp, Rbound)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

Nlarge = 0;
countlimit = 0;
while Nlarge ~= NlargeExp
    [centers, radii, ~] = imfindcircles(Rimg,Rlarge,'ObjectPolarity','dark','Method','TwoStage','Sensitivity',SL);
    XoutL = find(centers(:,1) - Rbound < 0);
    centers(XoutL,:) = [];
    radii(XoutL,:) = [];
    XoutR = find(centers(:,1) + Rbound > size(Rimg,1));
    centers(XoutR,:) = [];
    radii(XoutR,:) = [];
    YoutT = find(centers(:,2) - Rbound < 0);
    centers(YoutT,:) = [];
    radii(YoutT,:) = [];
    YoutB = find(centers(:,2) + Rbound > size(Rimg, 2));
    centers(YoutB,:) = [];
    radii(YoutB,:) = [];
    Nlarge=round(length(centers) / 2);
    particle(1:Nlarge) = struct('id',0,'x',0,'y',0,'r',0,'rm',0,'metric',0,'color','','fsigma',0,'z',0,'f',0,'g2',0,'forces',[],'betas',[],'alphas',[],'neighbours',[],'contactG2s',[],'forceImage',[]);
    if Nlarge > NlargeExp
        if (countlimit > 20)
            break
        end
        SL = SL - DS;
        % fprintf('Too many large particles, lowering SL\n');
        countlimit = countlimit + 1;
        % display(countlimit);
    elseif Nlarge < NlargeExp
        SL = SL + DS;
        % fprintf('Not enough large particles, raising SL\n');
        countlimit = countlimit + 1;
    end
end

% Rose = 0;
% Nsmall = 0;
% while Nsmall ~= NsmallExp
%     [centers, ~, ~] = imfindcircles(Rimg,Rsmall,'ObjectPolarity','bright','Method','TwoStage','Sensitivity',SS);
%     Nsmall = length(centers);
%     N = Nsmall + Nlarge;
%     particle(Nlarge+1:N) = struct('id',0,'x',0,'y',0,'r',0,'rm',0,'metric',0,'color','','fsigma',0,'z',0,'f',0,'g2',0,'forces',[],'betas',[],'alphas',[],'neighbours',[],'contactG2s',[],'forceImage',[]);
%     if Nsmall > NsmallExp
%         if Rose == 1 % Rose variable used to prevent infinite loop
%             break
%         end
%         SS = SS - DS;
%         fprintf('Too many small particles, lowering SS\n');
%     elseif Nsmall < NsmallExp
%         SS = SS + DS;
%         fprintf('Not enough small particles, raising SS\n');
%         Rose = 1;
%     end
% end

end

