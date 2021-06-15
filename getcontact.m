% Particle Detector, Neighbour finder and Contact validator that creates input files for my Photoelastic Disk Solver
% Particle Detection and Neigbour Finding Adapted from my Earlier Script (joCentersNonMonodisperse.m as of 2016/05/03)
% Photoelastic Disk Solver inspired from peDiskSolve by James Puckett (Phd-Thesis 2012) http://nile.physics.ncsu.edu

% If you use this please cite the follwoing paper
% K.E. Daniels, J. E. Kollmer & J. G. Puckett, "Photoelastic force measurements in granular materials", Rev. Sci. Inst. (201X)
% DOI: XXXXXX

% last edit on 2018/08/09 by Joshua Miller (jsmille9@ncsu.edu)

function getcontact(imagepath, centerdir, outputdir, Dm, Dpx, g2guess, FS, DT, conR, cG2Thrsd, ctrstL, ctrstH, shift4calibration)

    % close all % Housekeeping
    % clear all % Housekeeping

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                           User defined values                           %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    verbose = true; %Generates lots of plots showing results

% imagepath = '../raw/001.png';
% centerdir = '../center/txt/';
% outputdir = '../check_g2/'; 
% MeterPerpx = 0.008 / 100;
% g2cal = 100;
% fsigma = 6.25e4*2e-3; 
% dtol = 5;     
% CR = 10; 
% contactG2Threshold = 0.0001; 
% contrast_percent = [0.005, 0.995];
% shift4calibration = 0;
    
    MeterPerpx = Dm / Dpx;
    g2cal = g2guess; %Calibration Value for the g^2 method
    fsigma = FS; %photoelastic stress coefficient
    dtol = DT; % How far away can the outlines of 2 particles be to still be considered Neighbours
    CR = conR; %radius around a contactact point that is checked for contact validation
    contactG2Threshold = cG2Thrsd*CR^2; %sum of g2 in a contact area larger than this determines a valid contact
    contrast_percent = [ctrstL, ctrstH];
    
%     if nargin > 3
%         MeterPerpx = Dm / Dpx;
%         g2cal = g2guess; %Calibration Value for the g^2 method
%         fsigma = FS; %photoelastic stress coefficient
%         dtol = DT; % How far away can the outlines of 2 particles be to still be considered Neighbours
%         CR = conR; %radius around a contactact point that is checked for contact validation
%         contactG2Threshold = cG2Thrsd*CR^2; %sum of g2 in a contact area larger than this determines a valid contact
%     else
%         MeterPerpx = 0.008 / 455;
%         g2cal = 4000;
%         fsigma = 228; 
%         dtol = 7;     
%         CR = 40; 
%         contactG2Threshold = 0.008*CR^2; 
%     end


    imagedir = imagepath(1:max(strfind(imagepath, '/')));
    imagefile = imagepath(max(strfind(imagepath, '/'))+1:end);
    directory4raw = imagedir;
    
    if ~exist(outputdir, 'dir')
        mkdir(outputdir);
    end

    % directory = '/mnt/nas2/tora/jhc/force_measurement/tmp/';

    % files = dir([directory4raw, 'center_1_',sprintf('%3d', weight),'g_e.txt']); %Which files are we processing?
    files = dir([directory4raw, imagefile]);
    %files = dir('Centers*0001.txt'); %Alternatively, centers files can be loaded. This requires that both particle detections be flagged false however.
    [~, index] = natsortfiles({files.name}); % Sorting files as increasing number
    files = files(index);


    nFrames = length(files); %How many files are we processing ?

    % Hough Transform Values
    
    doParticleDetectionH = false; %Detect particles using Hough Transform?
    HoughDebug = false; %Debugs Hough Sensitivities so particles are found "better"

    DS = 0.0005; % How much should we adjust sensitivity if wrong number of particles are found
    RlargeH = [30 42]; %What radius (in pixels) range do we expect for the central dark area of large discs?
    Rcompensate = 182; % From the edge of core to edge of ring, it affects the contact point 
    Rbound = 165; % Outer ring radius in pixel.
    % RsmallH = [45 55]; %What radius (in pixels) range do we expect for the small discs?
    SL = 0.9300; %Sensitivity of the Hough Transform disc detetcor, exact value is Voodo magic...
    % SS = 0.89; %Sensitivity of the Hough Transform disc detetcor, exact value is Voodo magic...

    % NsmallH = 15; %Number of small discs. Only used in Hough Debug.
    NlargeH = 3; %Number of large discs. Only used in Hough Debug.

    % Convolution Method Values

    doParticleDetectionC = false; %Detect particles using convolution method?
    ConvDebug = false;

    % RlargeC = 66; %What radius (in pixels) do we expect for the large discs?
    % RsmallC = 47; %What radius (in pixels) do we expect for the small discs? 
    %Note: The above can be input in a range ONLY if the ConvDebug is set to
    %Note: true. Otherwise, the program needs the radius that works.
    % NsmallC = 15; %Number of small discs. Needed for Convolution.
    % NlargeC = 14; %Number of large discs. Needed for Convolution.

    % Neighbour Finding Values

    findNeighbours = true;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                 User Input Not Required Below This Line                 %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



    for frame = 1:nFrames %Loops for total number of images


        if HoughDebug

            imageFile = [directory4raw,files(frame).name]; %input filename
            img = imread(imageFile);
            Rimg = img(:,:,1);
            [SL] = PeGSHoughDebug_cali(Rimg, RlargeH, SL, DS, NlargeH, Rbound);

        end


        if doParticleDetectionH || doParticleDetectionC

            imageFile = [directory4raw,files(frame).name]; %input filename
            img = imread(imageFile); %read a color image that has particles in red and forces in green channel
            Rimg = img(:,:,1); %particle image
            Gimg = img(:,:,2); %force image

        else 

            centersfile = [centerdir, files(frame).name(1:end-3), 'txt']; %input filename
            ImgFile = [directory4raw, files(frame).name];  %adjusted force image filename
    %         ImgFile = [directory, centersfile(42:end-3),'png'];
    %         rImgFile = [directory, centersfile(8:end-3),'png'];  %adjusted force image filename
            img = imread(ImgFile); 
            Rimg =  img(:,:,1); %particle image
            Gimg = img(:,:,2); %force image

        end

        Gimg = im2double(Gimg);
        Rimg = im2double(Rimg);
        Gimg = Gimg-0.30*Rimg;
        Gimg = Gimg.*(Gimg > 0);
        % Gimg = imadjust(Gimg,stretchlim(Gimg,[contrast_percent(1) contrast_percent(2)]));
    %     Rimg = imadjust(Rimg,[0.06 0.6],[]);
    %     Rimg = imgaussfilt(Rimg,0.5);
    %     Rimg = im2bw(Rimg, 0.6);

    %     if (verbose)
    %         
    %         figure(1); %Draw the particle Image
    %         imshow(Rimg);
    %         
    %         figure(2); %Draw the Force Image
    %         imshow(Gimg);
    %         
    %     end

        if doParticleDetectionH

            particle = PeGSDiskFindH_circularforcecheck(Rimg, RlargeH, SL, MeterPerpx, fsigma, Rcompensate, Rbound);

        elseif doParticleDetectionC

            particle = PeGSDiskFindC(img, RsmallC, NsmallC, RlargeC, NlargeC);

        else

            pData = dlmread(centersfile, '\t'); %Read Position data from centers file
            XoutL = find(pData(:,2) - pData(:,3) < 0);
            pData(XoutL,:) = [];            
            XoutR = find(pData(:,2) + pData(:,3) > size(Rimg,1));
            pData(XoutR,:) = [];            
            YoutT = find(pData(:,1) - pData(:,3) < 0);
            pData(YoutT,:) = [];            
            YoutB = find(pData(:,1) + pData(:,3) > size(Rimg,2));
            pData(YoutB,:) = [];
            N = size(pData,1);
            clear particle;
            particle(1:N) = struct('id',0,'x',0,'y',0,'r',0,'rm',0,'color','','fsigma',0,'z',0,'f',0,'g2',0,'sImgg2',0,'forces',[],'betas',[],'alphas',[],'neighbours',[],'contactG2s',[],'forceImage',[]);
            for n=1:N %Bookkeeping
                particle(n).id= n;
                particle(n).x = pData(n,1); %-xoffset;
                particle(n).y = pData(n,2); %-yoffset;
                particle(n).r = pData(n,3);
                particle(n).rm = particle(n).r * MeterPerpx;
                particle(n).fsigma = fsigma;
                particle(n).color = 'r';
            end

        end

        N = length(particle);

        if(verbose)
            %add some information about the particles to the plots
            f1 = figure('Visible', 'off'); %Draw the particle Image
%             d = linspace(min(Rimg(:)),max(Rimg(:)),256);
%             Rimg1 = uint8(arrayfun(@(x) find(abs(d(:)-x)==min(abs(d(:)-x))),Rimg));
            Rimg1 = uint8(255 * mat2gray(Rimg));
%             Rimg = uint8(Rimg);
            image(Rimg1); colormap(gray); axis off; hold on
            for n=1:N
                viscircles([particle(n).x; particle(n).y]', particle(n).r,'EdgeColor',particle(n).color); %draw particle outline
                hold on
                plot(particle(n).x,particle(n).y,'rx'); %Mark particle centers
                text(particle(n).x,particle(n).y,num2str(particle(n).id),'Color','w');
            end
            if ~exist([outputdir, 'Rchannel/'], 'dir')
                mkdir([outputdir, 'Rchannel/']);
            end
            reddir = [outputdir, 'Rchannel/'];
            saveas(f1, [reddir, files(frame).name(1:end-4), '.png']);
            f2 = figure('Visible', 'off'); 
%             d = linspace(min(Gimg(:)),max(Gimg(:)),256);
%             Gimg1 = uint8(arrayfun(@(x) find(abs(d(:)-x)==min(abs(d(:)-x))),Gimg));
            Gimg1 = uint8(255 * mat2gray(Gimg));
%             figure(2);
            image(Gimg1); colormap(gray); axis off; hold on
            for n=1:N
                viscircles([particle(n).x; particle(n).y]', particle(n).r,'LineWidth',1,'EdgeColor',particle(n).color); %draw particle outline
                hold on
                plot(particle(n).x,particle(n).y,'rx'); %Mark particle centers
                text(particle(n).x,particle(n).y,num2str(particle(n).id),'FontSize',10,'Color','g');
            end
            drawnow;
        end
        for n=1:N
            %create a circular mask
            % => Find a better way yo do this masking!
            r = particle(n).r;
            mask = abs(-r:r);
            mask = mask.^2 + mask.^2';
            mask1 = double(sqrt(mask) <= r);
            mask2 = double(sqrt(mask) <= r*0.9-1);

            %This crops out a particle
            cropXstart = round(particle(n).x-r);
            cropXstop = round(particle(n).x-r)+ size(mask1,1)-1;
            cropYstart = round(particle(n).y-r);
            cropYstop = round(particle(n).y-r)+ size(mask1,2)-1;
            if (cropXstart <= 0 || cropXstop <= 0 || cropYstart  <= 0 || cropYstop <= 0)
                continue
            end
            cimg = im2double(Gimg(cropYstart:cropYstop, cropXstart:cropXstop));
            particleImg = cimg.*mask1;
            particle(n).forceImage=particleImg;

            se = strel('disk', 2);
            sigma = 1;
            template = imadjust(particleImg, stretchlim(particleImg, [0.05,0.95]));
            template = imgaussfilt(template, sigma);
            [gx,gy] = gradient(mask2.*template);
            g2 = (gx.^2 + gy.^2);
            particle(n).g2 = sum(sum(g2));


            %create a circular mask with a radius that is one pixel smaller
            %for cropping out the relevant gradient

%             mask2 = double(sqrt(mask) <= (r-0.2*r)-1);

%             %Compute G^2 for each particle
%             se = strel('disk', 2);            
%             sigma = 2;
% %             particleImg = particleImg - 0.4;
% %             particleImg = imerode(particleImg, se);
%             particleImg = imdilate(particleImg.*mask2, se);
%             particleImg = imgaussfilt(particleImg, sigma);
%             particleImg = imerode(particleImg, se);
%             particleImg = imadjust(particleImg);
%             particleImg = imgaussfilt(particleImg, 2);
%             particleImg = imadjust(particleImg, stretchlim(particleImg, [0.1,0.9]));
%             [gx,gy] = gradient(particleImg);
%             g2 = (gx.^2 + gy.^2);
%             [gx,gy] = imgradient(particleImg);
%             g2 = (gx.^2 + gy.^2);
%             particle(n).g2 = sum(sum(g2));
%             particle(n).f = particle(n).g2/g2cal;
        end

        if findNeighbours

            particle = PeGSNeighbourFind(Gimg, contactG2Threshold, dtol, CR, verbose, particle, shift4calibration);

        end
        
%         beta = particle(n).betas;

        
%         for n=1:N
%             px = size(particle(n).forceImage,1);
%             cx=px/2;cy=px/2;ix=px;iy=px;
%             r = particle(n).r;
%             maskR = 0.5*r;
%             xyr=[(r-maskR)*cos(beta)', (r-maskR)*sin(beta)', zeros(z,1)+(maskR-maskR*0.2)];
%             xyr=permute(xyr, [3, 2, 1]);
%             [x,y]=meshgrid(-(cx-1):(ix-cx),-(cy-1):(iy-cy));
%             mask2=any(hypot(x-xyr(1,1,:), y-xyr(1,2,:)) <= xyr(1,3,:),3);
%             
%             se = strel('disk', 2);            
%             sigma = 2;
%             particleImg = imdilate(particleImg.*mask2, se);
%             particleImg = imgaussfilt(particleImg, sigma);
% %             particleImg = imerode(particleImg, se);
% %             particleImg = imadjust(particleImg);
% %             particleImg = imgaussfilt(particleImg, 2);
% %             particleImg = imadjust(particleImg, stretchlim(particleImg, [0.1,0.9]));
%             [gx,gy] = gradient(particleImg);
%             g2 = (gx.^2 + gy.^2);
% %             [gx,gy] = imgradient(particleImg);
% %             g2 = (gx.^2 + gy.^2);
%             particle(n).g2 = sum(sum(g2));
%             particle(n).f = particle(n).g2/g2cal;            
%         end

        if ~exist([outputdir, 'Gcontact/'], 'dir')
            mkdir([outputdir, 'Gcontact/']);
        end
        contactdir = [outputdir, 'Gcontact/'];
        saveas(f2, [contactdir, files(frame).name(1:end-4), '.png']);

        
        
%         figure(3);
        f3 = figure('Visible', 'off');
%         d = linspace(min(Gimg(:)),max(Gimg(:)),256);
%         Gimg1 = uint8(arrayfun(@(x) find(abs(d(:)-x)==min(abs(d(:)-x))),Gimg));
        Gimg1 = uint8(255 * mat2gray(Gimg));
        image(Gimg1); colormap(gray); axis off; hold on
        for n = 1:N
            z = particle(n).z; %get particle coordination number
            if (z>0) %if the particle does have contacts
                for m = 1:z %for each contact
                    %draw contact lines
                    lineX(1)=particle(n).x;
                    lineY(1)=particle(n).y;
                    lineX(2) = lineX(1) + particle(n).r * cos(particle(n).betas(m));
                    lineY(2) = lineY(1) + particle(n).r * sin(particle(n).betas(m));
                    cX = lineX(1) + (particle(n).r-CR) * cos(particle(n).betas(m));
                    cY = lineY(1) + (particle(n).r-CR) * sin(particle(n).betas(m));
                    hold on; % Don't blow away the image.
                    plot(lineX, lineY,'-y','LineWidth',1);hold on;
                end
            end
        end
        if ~exist([outputdir, 'center2c/'], 'dir')
             mkdir([outputdir, 'center2c/']);
        end
        center2contactdir = [outputdir, 'center2c/'];   
        saveas(f3, [center2contactdir, files(frame).name(1:end-4), '.png']);
    %     saveas(figure(2), [contactdir, files(frame).name(1:end-4), '.png']);
        %Save what we got so far

        if ~exist([outputdir, 'mat/'], 'dir')
            mkdir([outputdir, 'mat/']);
        end
        if ~exist([outputdir, 'csv/'], 'dir')
            mkdir([outputdir, 'csv/']);
        end
        save([[outputdir, 'mat/'], files(frame).name(1:end-4),'.mat'],'particle');
        writetable(struct2table(rmfield(particle, 'forceImage'),'AsArray',true), [[outputdir, 'csv/'], files(frame).name(1:end-4),'.csv']);
    end
end