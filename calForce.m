function calForce(matpath, outdir, rMask, tF1, tF2, tF3, tF4, g2_cal_a, g2_cal_b, exitX, exitY, D_arch2exit, given_force, calibrate, optimization) 
% matpath='../s1/mat/001.mat';
% outdir='../test_opt_1set3_weightMax/';
% rMask=0.45;
% tF1=3.0;tF2=0.0;tF3=0.0;tF4=0.0;
% g2_cal_a = 80/0.88; % [g^2]/[N]
% g2_cal_b = 0;
% exitX=427;
% exitY=475;
% D_arch2exit = 260;
% given_force = 0;
% calibrate = 0;
% optimization  = 1;

    

    tF=[tF1,tF2,tF3,tF4];   
    if (nnz(tF) == 1)
        tF = tF1;
    elseif (nnz(tF) == 2)
        tF = [tF1, tF2];
    elseif (nnz(tF) == 3)
        tF = [tF1,tF2,tF3];
    else
        tF=[tF1,tF2,tF3,tF4];
    end
    ntF=length(tF); tS=zeros(ntF^5,5);
    directory = matpath(1:max(strfind(matpath, '/')));
    matname = matpath(max(strfind(matpath, '/'))+1:end);
    if ~exist(outdir, 'dir')
        mkdir(outdir);
    end
    files = dir([directory, matname]); 

    
    [~, index] = natsortfiles({files.name}); % Sorting files as increasing number
    files = files(index);

    scaling = 1; %scale the image by this factor before doing the fit
    verbose = true; 
    fitoptions = optimoptions('lsqnonlin','Algorithm','levenberg-marquardt','MaxIter',100,'MaxFunEvals',400,'TolFun',0.001,'Display','final-detailed');

    nFrames = length(files); %how many files are we processing ?
    for frame = 1:nFrames %loop over these frames 
        fileName = [directory,files(frame).name]; %which file/frame are we processing now ?
        load(fileName); %load the particle data file to process
        pres=particle; %initialize a copy of the particle vector to mess around with
        N = length(particle); %number of particles in this frame
        clearvars presAll;
        presAll(N,ntF) = struct();
        dirname0 = [outdir, files(frame).name(1:end-4),'/'];
        if ~exist(dirname0, 'dir') mkdir(dirname0); end
        dirname1 = [dirname0, sprintf('maskR_%1.2f/', rMask)];
        if ~exist(dirname1, 'dir') mkdir(dirname1); end
        
%         disp(['processing file ',fileName, ' containing ' ,num2str(N), ' particles']); %status indicator
        for n=1:N
            if sqrt((particle(n).x - exitX)^2 + (particle(n).y - exitY)^2) < D_arch2exit
                clearvars pCompare;
                clearvars pres1;
    %             disp(['fitting force(s) to particle ',num2str(n)]); %status indicator
                if (particle(n).z > 0 )
                    fsigma = particle(n).fsigma;
                    rm = particle(n).rm;
                    template = particle(n).forceImage;
                    template = imresize(template,scaling);  
                    se = strel('disk', 2);
                    sigma = 1;
%                     template = template.*(template > 0.1);                    
                    template = imadjust(template, stretchlim(template, [0.05,0.95]));
%                     template = imerode(template, se);
%                     template = imdilate(template, se);
                    template = imgaussfilt(template, sigma);
                    px = size(template,1); %size of the force image
                    if verbose
                        % f1 = figure('Visible', 'off'); 
                        % template1 = uint8(255 * mat2gray(template));
                        % ax1 = subplot(1,2,1);
                        % image(template1); colormap(gray);
                        % pbaspect(ax1, [1 0.9 1]);
                        % axis off;
                    end
                    z = particle(n).z; NumOfCombination=ntF^z;  % disp(NumOfCombination); return;
                    forces = zeros(z,1);
                    cg2s = sum(particle(n).contactG2s);
                    dirname2 = [dirname1, sprintf('p%02d/', n)];
                    if ~exist(dirname2, 'dir')
                        mkdir(dirname2)
                    end          
                    k=1;
                    for t5=1:ntF for t4=1:ntF for t3=1:ntF for t2=1:ntF for t1=1:ntF
                        tS(k,1)=tF(t1);tS(k,2)=tF(t2);tS(k,3)=tF(t3);tS(k,4)=tF(t4);tS(k,5)=tF(t5);
                        k=k+1; 
                    end; end; end; end; end

                    tSTb = zeros(size(tS,1), size(tS,2)+1);
                    for i=1:size(tS,1)
                        tSTb(i,:)= [i, tS(i,:)]; 
                    end

                    dlmwrite([dirname2, 'cmbntable.txt'], tSTb(1:NumOfCombination,1:z+1), 'delimiter', '\t', 'precision', '%1.1f');

                    pCompare(1:NumOfCombination) = struct('NumComb',0,'fitError',0);
                    pres1(1) = struct();
                    
                    cx=px/2;cy=px/2;ix=px;iy=px;r=rMask*px; %create a circular mask
                    [x,y]=meshgrid(-(cx-1):(ix-cx),-(cy-1):(iy-cy));
                    c_mask=((x.^2+y.^2)<=r^2-1);
                    % [gx1,gy1] = gradient(c_mask.*template);
                    % g21 = (gx1.^2 + gy1.^2);
                    % g2_center = sum(sum(g21));
                    
                    for k=1:NumOfCombination
                        dirname3 = [dirname2, sprintf('%04d', k), '/'];
                        if ~exist(dirname3, 'dir') mkdir(dirname3); end
                        if ~exist([dirname3, 'img/'], 'dir') mkdir([dirname3, 'img/']); end
                        if ~exist([dirname3, 'csv/'], 'dir') mkdir([dirname3, 'csv/']); end 
                        if ~exist([dirname3, 'mat/'], 'dir') mkdir([dirname3, 'mat/']); end

                       
                        if (given_force == true)
                            for i=1:z
                                forces(i)=str2num(files(frame).name(1:end-4))*9.8e-3;
                            end
                        else
                            for i=1:z
                                forces(i) = tS(k,i)*(particle(n).g2-g2_cal_b)/g2_cal_a*particle(n).contactG2s(i)/cg2s;
                            end
                        end

                        alphas = zeros(z,1); %Initial guesses for the angles of attack of the forces
                        beta = particle(n).betas;
    %                     [alphas,forces] = forceBalance(forces,alphas,beta); %apply force balance to the initial guesses

                        % cx=px/2;cy=px/2;ix=px;iy=px;
                        % r = particle(n).r;
                        % maskR = maskradius2pxsize*r;[0.446940883540157,0.238812162754233,0.561448596751543,0.253990740185760]
                        % xyr=[(r-maskR)*cos(beta)', (r-maskR)*sin(beta)', zeros(z,1)+(maskR-5)];
                        % xyr=permute(xyr, [3, 2, 1]);
                        % [x,y]=meshgrid(-(cx-1):(ix-cx),-(cy-1):(iy-cy));
                        % c_mask=any(hypot(x-xyr(1,1,:), y-xyr(1,2,:)) <= xyr(1,3,:),3);

                        func = @(par) joForceImg(z, par(1:z),par(z+1:z+z), beta(1:z), fsigma, rm, px, verbose); %+par(2*z+1); %this is the function I want to fit (i.e. synthetic stres image), the fitting paramters are in vector par
                        err = @(par) real(sum(sum( ( c_mask.*(template-func(par)).^2) ))); %BUG: for some reason I sometimes get imaginary results, this should not happen
    %                     err = @(par) 1-ssim(c_mask.*(template), c_mask.*(func(par)));
                        p0(1:z) = forces; %Set up initial guesses
                        p0(z+1:2*z) = alphas;
                        if (optimization == true)
                            p=lsqnonlin(err,p0,[],[],fitoptions);
                            
                            forces = p(1:z); %get back the result from fitting
                            alphas = p(z+1:z+z);
                            fitError = err(p);
                        else
                            fitError = err(p0);
                        end

    %                     errstore = [errstore, fitError];
                        imgFit = joForceImg(z, forces, alphas, beta, fsigma, rm, px*(1/scaling), verbose); %generate an image with the fitted parameters
                        [gx,gy] = gradient(c_mask.*imgFit);
                        g2 = (gx.^2 + gy.^2);
                        [alphas,forces] = forceBalance(forces,alphas,beta);
                        forces = abs(forces);
                        pres(n).NumComb = k;                    
                        pres(n).fitError = fitError; %store the new information into a new particle vector 
                        pres(n).sImgg2 = sum(sum(g2));
                        pres(n).betas = beta;
                        pres(n).forces = forces;
                        pres(n).alphas= alphas;
                        pres(n).synthImg = imgFit;
                        % pres(n).g2 = g2_center;
                        pres1 = pres(n);

                        presAll(n,k).id = n;
                        presAll(n,k).NumComb = k; 
                        presAll(n,k).fitError = fitError; %store the new information into a new particle vector 
                        presAll(n,k).x = pres1.x;
                        presAll(n,k).y = pres1.y;
                        presAll(n,k).r = pres1.r;
                        presAll(n,k).rm = pres1.rm;
                        presAll(n,k).color = pres1.color;
                        presAll(n,k).fsigma = pres1.fsigma;
                        presAll(n,k).z = pres1.z;
                        presAll(n,k).f = pres1.f;
                        presAll(n,k).g2 = pres1.g2;
                        presAll(n,k).sImgg2 = pres1.sImgg2;
                        presAll(n,k).betas = pres1.betas;
                        presAll(n,k).forces = forces;
                        presAll(n,k).alphas= alphas;
                        presAll(n,k).neighbours = pres1.neighbours;
                        presAll(n,k).contactG2s = pres1.contactG2s;
                        presAll(n,k).forceImage = pres1.forceImage;
                        presAll(n,k).contactIs = pres1.contactIs;
                        presAll(n,k).synthImg = imgFit;

                        pCompare(k).NumComb = k;
                        pCompare(k).fitError = fitError;
                        pCompare(k).forces = forces;
                        pCompare(k).betas = beta;
                        pCompare(k).alphas = alphas;
                        pCompare(k).maskRsize = rMask;
                        pCompare(k).synthImg = imgFit;
                        pCompare(k).sImgg2 = pres1.sImgg2;

                        imwrite(particle(n).forceImage, [dirname3, 'img/raw.png']);
                        imwrite(c_mask.*template, [dirname3, 'img/processed.png']);
                        imwrite(imgFit, [dirname3, 'img/syn.png']);                       
                        writetable(struct2table(rmfield(pres1,{'forceImage', 'synthImg'}),'AsArray',true), [[dirname3, 'csv/'], sprintf('%04d.csv', k)]); 
                        save([[dirname3, 'mat/'], sprintf('%04d.mat', k)],'pres1');


                    end
                    writetable(struct2table(rmfield(pCompare,'synthImg'),'AsArray',true), [dirname2, 'pixelError.csv']);
                    [~, minI] = min([pCompare.fitError]);
                    pres(n).sImgg2 = pCompare(minI).sImgg2;
                    pres(n).forces = pCompare(minI).forces;
                    pres(n).alphas = pCompare(minI).alphas;
                    pres(n).NumComb = pCompare(minI).NumComb;                
                    pres(n).fitError = pCompare(minI).fitError;
                    pres(n).synthImg = pCompare(minI).synthImg;
                end
            end
            if (calibrate == true)
                if ~exist([outdir, sprintf('maskR%1.2f_img/', rMask)], 'dir')
                    mkdir([outdir, sprintf('maskR%1.2f_img/', rMask)])
                end
                imwrite(pres(n).synthImg, [outdir, sprintf('maskR%1.2f_img/', rMask), files(frame).name(1:end-4),'.png']);
            end
        end
        
        if ~exist([outdir, sprintf('maskR%1.2f_mat/', rMask)], 'dir')
            mkdir([outdir, sprintf('maskR%1.2f_mat/', rMask)])
        end
        save([[outdir, sprintf('maskR%1.2f_mat/', rMask)], files(frame).name(1:end-4), '.mat'],'pres');
        
        if ~exist([outdir, sprintf('maskR%1.2f_csv/', rMask)], 'dir')
            mkdir([outdir, sprintf('maskR%1.2f_csv/', rMask)])
        end       
        if (isfield(pres,'synthImg') == 1)
            writetable(struct2table(rmfield(pres,{'synthImg', 'forceImage'}),'AsArray',true), [outdir, sprintf('maskR%1.2f_csv/', rMask),files(frame).name(1:end-4), '.csv']);
        end
        
        if ~exist([outdir, sprintf('maskR%1.2f_matAll/', rMask)], 'dir')
            mkdir([outdir, sprintf('maskR%1.2f_matAll/', rMask)])
        end
        save([[outdir, sprintf('maskR%1.2f_matAll/', rMask)], files(frame).name(1:end-4), '.mat'],'presAll');
        
%         if ~exist([outdir, sprintf('maskR%1.2f_csvAll/', maskradius2pxsize)], 'dir')
%             mkdir([outdir, sprintf('maskR%1.2f_csvAll/', maskradius2pxsize)])
%         end   
%         if (isfield(presAll,'synthImg')== 1)
%             writetable(struct2table(rmfield(presAll,{'synthImg', 'forceImage'})), [outdir, sprintf('maskR%1.2f_csv/', maskradius2pxsize),files(frame).name(1:end-4), '.csv']);        
%         end
    end
end