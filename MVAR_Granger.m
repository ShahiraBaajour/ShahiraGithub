%Code by SHAHIRA BAAJOUR- Brain Imaging Research Division @ Wayne State Univerity
%--------------------------------------------------------------------------
%This code uses input fMRI timeseries data, # of subjects,
%models, conditions, onset times to generate multivariate autoregressive analysis
%consistent with principles of Granger Causality using the fitglm function.
%The output is a single matrix of MVAR coefficients representing mean 
%strength of directional interactions between regions.

%Date Last modified: March 3,2020
%Contact: baajour35@gmail.com
%--------------------------------------------------------------------------
%open data files
folderlist=dir();
for i=4:length(folderlist)
        cd (folderlist.name(i),C)%load timeseries data
        data_ts(i,:)=C;%assemble timeseries data [nsubs x nregs]
end

nsubs=size(data_ts,1);%# of Subjects n
nreg=size(data_ts,2);% # of regions involved
ncond=2;% 2 Conditions, e.g. Cooling, Warming
nsets=4;
nmod = zeros(1,nsets);

% 150s 
nmod(1)=57;
nmod(2)=57;
nmod(3)=57;
nmod(4)=57;
indy = zeros(max(nmod),nsets);
indX = zeros(max(nmod),nsets);
indX(1:nmod(1),1) = 150:206;
% ** model 1, y sequence, relative onset times
indy(1:nmod(1),1) = 151:207;% lag=1 introduced

indX(1:nmod(2),2) = 427:483;
indy(1:nmod(2),2) = 428:484;

indX(1:nmod(3),3) = 46:102;
indy(1:nmod(3),3) = 47:103;

indX(1:nmod(4),4) = 300:356;
indy(1:nmod(4),4) = 301:357;

npair = nreg*(nreg-1);
coeff = cell(nsubs,1);
% initialize 3D matrices for resultant MVAR Coeffs
for n=1:nsubs
     coeff{n,1} = zeros(nreg,nreg,nsets);
end 
for imageset=1:nsets
    % disp(imageset);
    nummodels = nmod(imageset);
    imsy = indy(1:nummodels,imageset);
    imsX = indX(1:nummodels,imageset);
    % imsy is the index of y sequence images for the current image set
    % imsX is the index of X sequence images for the current image set
    
    % loop over subjects
        % get data file name for each subject
        for sub=1:nsubs
            Z = data_ts(sub,:);
        
            for reg=1:nreg
            % left
                y = Z(imsy,reg,1);
                X = Z(imsX,:,1);
            
            % use fitglm function to create generalized linear model of X
            % and y 
                fit = fitglm(X,y);
                coeff{sub,1}(reg,:,imageset) = ...
                (table2array(fit.Coefficients(2:(nreg+1),1)))';%Store the Coeffs 
            end
        end 
end
s=cell(nsubs,ncond);%define matrix to combine all results
for i=1:nsubs
    for j=1:ncond
        s{i,j}=zeros(nreg,nreg,1);
    end 
end
for sub=1:nsubs
        %organize the combination of results in s
        s{sub,1}=mean(coeff{sub,1}(:,:,1:2),3);%sets 1:2= Cooling Condition
        s{sub,2}=mean(coeff{sub,1}(:,:,3:4),3);%sets 3:4 = Warming Condition
        
end 
for sub=1:nsubs
    for i=1:nreg
        for j=1:nreg
             if i==j
                 % ignore region interaction with themselves
                s{sub,1}(i,j,:) = NaN;
                s{sub,2}(i,j,:)=NaN;
             end
        end
    end
end
average=cell(1,2);%average out across all subjects, across sets of conditions
for i=1:2
    average{1,i}=zeros(nreg,nreg);
end 
for sub=1:nsubs
    for j=1:2
        average{1,j}=average{1,j}+s{sub,j};
    end 
end 
for n=1:2
    average{1,n}=average{1,n}./nsubs;%Final MVSR coeff matrix for both conditions
end 


