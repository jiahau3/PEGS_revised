
%joForceAdjMat

%Compute adjacency matrices and some statistics on joDiscSolve output
%derived from the joFroceAdjMat script as of 2016/08/12
%wich was derived from the joForceStat.m script as of 2016/07/13

%works for the joDiscSolve output file format as of September 14 2016

%Written by Jonathan Kollmer
%Last update 2016/09/29

% function ForceAdjMat(matpath, imgdir, outdir)

matpath = '../s2_Optimized/maskR0.50_mat/*.mat';
imgdir = '../raw/';
outdir = '../s3_Optimized/';
    E = 1.24172e6; % Young's modulus
    nu = 0.5; % Poisson's ratio
    thickness = 2e-3;
    inputdir = matpath(1:max(strfind(matpath, '/')));
    matname = matpath(max(strfind(matpath, '/'))+1:end);

    % close all %Housekeeping
    % clear all %Housekeeping


    files = dir([inputdir,matname]); %which files are we processing ?
    [~, index] = natsortfiles({files.name}); % Sorting files as increasing number
    files = files(index);

    outputdir = outdir;
    if ~exist(outputdir, 'dir')
        mkdir(outputdir);
    end

    nFrames = length(files); %how many files are we processing ?

    %PARAMETERS NEEDED TO RUN THIS SCRIPT ARE SET HERE

    fmin = 0.001; %minimum force (in Newton) to consider a contact a valid contact
    fmax = 20; %maximum force (in Newton) to consider a contact a valid contact
    emax = 100000; %maximum fit error/residual to consider a contact a valid contact
    fs=16; %plot font size
    verbose = false; %make lots of plots as we go
    % jut in case your data coordinate system is offset from the image
    % coordinate system (i.e. you only processed a small part of a larger
    % image)
    % xoffset = 1000;
    % yoffset = 700;
    % xsize = 2470-xoffset;
    % ysize = 2260-yoffset;

    % workspacefilename = [];

    %Global Metrics will be stored in these structures
    allContacts = struct('frame',0,'id',0,'contact_neighbor',0,'x',0,'y',0,'fAbs',0,'fNorm',0,'fTan',0,'alpha',0,'beta',0,'contactX',0,'contactY',0,'error',0); %data structure to store information about contacts
    aID = 1; %Global contact counter over all contacts in all cycles.

    for cycle = 1:nFrames %loop over these cycles 
        clearvars particle;
        clearvars contact;

        %input filnames
        peInfilename = [inputdir,files(cycle).name]; %input filename 
        camImageFileName = [imgdir, files(cycle).name(1:end-4),'.png'];  %adjusted force image filename 
        %camImageFileName = ['frame1.png'];  %adjusted force image filename 

        %output filenames 
        if ~exist([outputdir, 'mat/'], 'dir') mkdir([outputdir, 'mat/']); end
        if ~exist([outputdir, 'png/'], 'dir') mkdir([outputdir, 'png/']); end
        if ~exist([outputdir, 'energy_img/'], 'dir') mkdir([outputdir, 'energy_img/']); end
        if ~exist([outputdir, 'energy_array/'], 'dir') mkdir([outputdir, 'energy_array/']); end        
        if ~exist([outputdir, 'adjmatri/'], 'dir') mkdir([outputdir, 'adjmatri/']); end
        if ~exist([outputdir, 'csv/'], 'dir') mkdir([outputdir, 'csv/']); end
        workspacefilename =  [[outputdir, 'mat/'], files(cycle).name(1:end-4),'.mat']; 
        synthImgFilename = [[outputdir, 'png/'], files(cycle).name(1:end-4),'.png'];  %output filename 
        energyImgFilename = [[outputdir, 'energy_img/'], files(cycle).name(1:end-4),'.png'];  %output filename 
        energyArrayFilename = [[outputdir, 'energy_array/'], files(cycle).name(1:end-4),'.txt'];  %output filename         
        adjMatAbsFilename = [[outputdir, 'adjmatri/'], files(cycle).name(1:end-4),'.txt'];  %output filename 


        % NO PARAMETERS SHOULD BE SET BY HAND BELOW THIS LINE

        %check if the data we want to read exists
        %if it does, load it, else abort
        if ~(exist(peInfilename, 'file') == 2) %if the file we try to open does not exist
            display(['File not Found:', peInfilename]); %complain about it
            return %and end the execution of this script
        else
            load(peInfilename); %read peDiscSolve ouput
            particle = pres;
            NN = length(particle);
        end

        
        % Gernate Combined Syntehtic Force Image
        img = imread(camImageFileName); %camera force image
        %img = imcrop(img,[xoffset, yoffset, xsize, ysize]); %particle image
        bigSynthImg = zeros(size(img,1),size(img,2)); %make an empty image with the same size as the camera image
        bigEnergyImg = zeros(size(img,1),size(img,2));        
        for n=1:NN %for all particles
            %display(['fitting force(s) to particle ',num2str(n)]); %status indicator
            if (particle(n).z > 0 )
                %Add the syntetic peImage for the particle to the
                %synthetic image of our whole packing 
                x = floor(pres(n).x); %interger rounded x coordinate of the current particle
                y = floor(pres(n).y); %interger rounded y coordinate of the current particle
                r = pres(n).r; %radius (in pixels) of the current particle
                sImg = pres(round(n)).synthImg; %synthetic force image for the current particle
                sx = size(sImg,1)/2; %width of the synthetic force image of the current particle
                sy = size(sImg,2)/2; %heights of the synthetic force image of the current partice
                bigSynthImg(round(y-sy+1):round(y+sy),round(x-sx+1):round(x+sx)) = bigSynthImg(round(y-sy+1):round(y+sy),round(x-sx+1):round(x+sx))+sImg; %Add the syntetic Force Image of the current particle to the appropriate location
            end
        end

        %DATA EVALUATION AND ANALYSIS STARTS HERE

        %particle(1:size(data,1)) = struct('id',0,'x',0,'y',0,'z',0,'fx',0,'fy',0); %data structure to store particle information
        contact = struct('id1',0,'id2',0,'x',0,'y',0,'fAbs',0,'fNorm',0,'fTan',0,'alpha',0,'beta',0,'contactX',0,'contactY',0,'error',0); %data structure to store information about contacts
        cID = 1; %contact counter
        %A = zeros(NN); %empty binary adjacency matrix
        W = zeros(NN); %empty force weighted adjacency matrix
        %N = zeros(NN); %empty normal force weighted adjacency matrix
        %T = zeros(NN); %empty tangential force weighted adjacency matrix

        for n = 1:NN %for each particle
            if (isfield(particle(n), 'fitError') == 0)
                continue
            end
            err = particle(n).fitError; %get fit error 
%             err = 1000;
            z = particle(n).z; %get coordination number
            r = particle(n).r; %get particle radius in pixel
            rm = particle(n).rm; %get particle radius in meters
            px = size(particle(n).forceImage,1);

            if(length(particle(n).neighbours) > 0 && length(particle(n).alphas) > 0) % particle is in contact
                contacts = particle(n).neighbours; %get IDs of all contacting particles
                betas = particle(n).betas; %get the beta angle (position of contact point) associated with each contact
                forces = particle(n).forces; %get the force associated with each contact
                alphas = particle(n).alphas; %get the alpha angle (direction of force) associated with each contact
                %Calculate stress ditribution and energy
                
                elastic_energy = stress_distribution(z,forces, alphas, betas, particle(n).fsigma, rm, px, verbose, E, nu,thickness);
                totEnergy = sum(sum(elastic_energy))*thickness;
                
                x = floor(particle(n).x); %interger rounded x coordinate of the current particle
                y = floor(particle(n).y); %interger rounded y coordinate of the current particle
                sx = size(particle(n).synthImg,1)/2; %width of the synthetic force image of the current particle
                sy = size(particle(n).synthImg,2)/2; %heights of the synthetic force image of the current partice
%                 rr=1*sx; %create a circular mask
%                 [a,b]=meshgrid(-(sx-1):(sx),-(sy-1):(sy));
%                 c_mask=((a.^2+b.^2)<=(rr)^2);      
                EImg = log(elastic_energy+(elastic_energy==0));
                EImg = fliplr(EImg);
                EImg = rot90(EImg,1);
                % bigEnergyImg(round(y-sy+1):round(y+sy),round(x-sx+1):round(x+sx)) = bigEnergyImg(round(y-sy+1):round(y+sy),round(x-sx+1):round(x+sx))+EImg;
                bigEnergyImg = EImg;
                
                for m=1:length(forces) %for each contact
                    
                    if(forces(m) > fmin && err < emax && forces(m) < fmax) %is this a valid contact ?
                                 
                        %put information about the first particle involved in this
                        %contact in the corresponding particle structure vector
                        
                        %ideally the accumulated fx and fy should be zero, that is
                        %the particle is in force balance
                        %particle(n).fx = particle(n).fx + forces(m) * cos(betas(m)-pi); %x component of total force vector %CHECK AGAIN IF THIS IS GEOMETRICALLY CORRECT
                        %particle(n).fy = particle(n).fy + forces(m) * sin(betas(m)-pi); %y component of total force vector %CHECK AGAIN IF THIS IS GEOMETRICALLY CORRECT
                        %particle(n).z = particle(n).z+1; %increment the real contact number for the current particle
                        
                        %put all the information about this contact
                        %into the contact struct vector
                        contact(cID).id1 = n; %first particle involved in this contact
                        contact(cID).id2 = contacts(m); %second particle involved in this contact
                        contact(cID).x = particle(n).x;
                        contact(cID).y = particle(n).y;
                        contact(cID).fAbs = forces(m); %absolute force
                        contact(cID).fNorm = forces(m)*cos(alphas(m)); %normal force (see Eq. 4.16)
                        contact(cID).fTan = forces(m)*sin(alphas(m)); %tangential force
                        contact(cID).alpha = alphas(m); %the alpha angle (direction of force) associated with this contact
                        contact(cID).beta = betas(m); %the beta angle (position of contact point) associated with this contact
                        contact(cID).contactX =  r * cos(betas(m)-pi); %x component of vector to contact point
                        contact(cID).contactY =  r * sin(betas(m)-pi); %y component of vector to contact point
                        contact(cID).error = err; %fit error for this particle (the first particle in the contact)
                        
                        cID = cID + 1; %increment contact counter
                        
                        allContacts(aID).frame = str2num(files(cycle).name(1:end-4));
                        allContacts(aID).id = n;
                        allContacts(aID).contact_neighbor = contacts(m);
                        allContacts(aID).x = particle(n).x;
                        allContacts(aID).y = particle(n).y;
                        allContacts(aID).fAbs = forces(m); %absolute force
                        allContacts(aID).fNorm = forces(m)*cos(alphas(m)); %normal force (see Eq. 4.16)
                        allContacts(aID).fTan = forces(m)*sin(alphas(m)); %tangential force
                        allContacts(aID).alpha = alphas(m);
                        allContacts(aID).beta = betas(m);
                        allContacts(aID).contactX = round(r * cos(betas(m)-pi) + particle(n).x);
                        allContacts(aID).contactY = round(r * sin(betas(m)-pi) + particle(n).y);
                        allContacts(aID).error = err;
                        allContacts(aID).energydisp = elastic_energy;
                        allContacts(aID).totenergy = totEnergy;
                        aID = aID+1;
                        
                        %build some adjacency matrices
                        if (contacts(m)>0) %correct for negative contact IDs in peDiscsolve, i.e.non-wall contacts only
                            %A(n,contacts(m)) = 1; %mark contact in the binary adjacency matrix
                            W(n,contacts(m)) = forces(m); %write the corrsponding force as a weight into an adjacency matrix
                            %N(n,contacts(m)) = forces(m)*cos(alphas(m)); %write the corrsponding normal force as a weight into an adjacency matrix
                            %T(n,contacts(m)) = forces(m)*sin(alphas(m)); %write the corrsponding tangential force as a weight into an adjacency matrix
                        end
                    end
                end
            end
        end


        %make sure our adjacency matrix is nice
        %discard singular contacts
            B = double((W.*W') > 0);
            W = (W.*B);
        %average reciprocal contacts so the matrix is fully symmetric
            W = (W+W')/2;

        %DATA OUTPUT STARTS HERE
        %write out the weighte adjacency matrix for the current frame
        dlmwrite(adjMatAbsFilename, W, '\t'); 
        %write out the synthetic force image for the current frame
        imwrite(bigSynthImg,synthImgFilename);

        f1 = figure('Visible','off');
        imagesc(bigEnergyImg); axis off; colormap(jet); colorbar;
        saveas(f1, energyImgFilename);
        
        dlmwrite(energyArrayFilename,bigEnergyImg, '\t');        

        %Plot what we have learned so far
        if verbose
            %set up a new figure environment and scale it appropriately
            close all;
            pFig = figure ;
            %set(pFig,'Position',[0 0 2000 1200]);
            %set(pFig,'paperorientation','landscape');

            %this suplot shows the original image, overlayed with the
            %detected contacts
            figure(1)
                %read and display the original image used as input to peDisc
                %img = imcrop(imread(camImageFileName),[xoffset, yoffset, xsize, ysize]); %force image
                img = imread(camImageFileName); %force image
                imshow(img); hold on;
                colormap(gray);
                %plot the centers of all particles associated with a contact
                plot([contact.x],[contact.y],'or')
                %plot arrows from the centers of all particles associated with a contact to
                %the contact point
                quiver([contact.x],[contact.y],[contact.contactX],[contact.contactY],0,'LineWidth',1.5)
                %set font sizes and labels
                set(gca,'FontSize',fs);
                title('camera image','FontSize',fs);
                hold off;
            figure(2)
                %read and display the original image used as input to peDisc
                %img = imcrop(imread(camImageFileName),[xoffset, yoffset, xsize, ysize]); %force image
                img = bigSynthImg; %force image
                imshow(img); hold on;
                colormap(gray);
                %plot the centers of all particles associated with a contact
                plot([contact.x],[contact.y],'or')
                %plot arrows from the centers of all particles associated with a contact to
                %the contact point
                quiver([contact.x],[contact.y],[contact.contactX],[contact.contactY],0,'LineWidth',1.5)
                %set font sizes and labels
                set(gca,'FontSize',fs);
                title('camera image','FontSize',fs);
                hold off;
            figure(3)
                imagesc(W); 
                colormap(jet);
        end
    end

    %save everything
    save(workspacefilename, 'allContacts');
    writetable(struct2table(rmfield(allContacts,'energydisp')), [[outputdir, 'csv/'], files(end).name(1:end-4), '.csv']);

% end
