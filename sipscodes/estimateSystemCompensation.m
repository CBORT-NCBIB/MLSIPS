
function [out,histBin] = estimateSystemCompensation(S1,procStruct)

    modTrue = 0;
    isFilt = 0;
    fwx = 6;
    fwxProvided = 0;
    N = size(S1,3);
    % For Histograms
    Nel = 1001;
    Naz = 2001;

    fnames = fieldnames(procStruct);
    for ind = 1:numel(fnames)
        if strcmp(fnames{ind},'isMod')
            modTrue = procStruct.isMod;
        elseif strcmp(fnames{ind},'isFilt')
            isFilt = procStruct.isFilt;
        elseif strcmp(fnames{ind},'fwx')
            fwx = procStruct.fwx;
            fwxProvided = 1;
        end
    end
   
    % If modulating, take every other A-line
    if modTrue ~=0
        S1 = S1(1:2:end,:,:,:,:);
    end
    % If unfiltered stokes vectors added, filter before further processing
    if isFilt==0
        disp("Filtering Stokes Vectors...")
        if fwxProvided ==0
            disp("Lateral filtering width not provided...")
            disp("...:::... Proceeding with default 6 ...:::...")
            S1f = filterStokes(S1,fwx);
        else
            disp("Lateral filtering width provided...")
            disp("...:::... Proceeding with fwx = "+fwx+" ...:::...")
            S1f = filterStokes(S1,fwx);
        end
    else
        S1f = S1;
    end
   
   [cOut,histBin] = EstimateAsymmetricCompensation(S1f,Naz,Nel);

    % Get Symmetric Matrix
    CT = cOut.symMatrix;
    C = pagetranspose(CT);
    CvecSIPS = decomposeRot(C);

    disp("Aligning Q Matrices")
    for i = 1:size(S1,5)
        disp("Working on Slice #"+num2str(i)+"/"+num2str(size(S1,5))+"...")
        % Get Q matrix (alignment)
        qOut = EstimateSymmetricCompensation(C,N,S1f(:,:,:,:,i));
        tempSC = qOut.sysCompensation;

        % Stand in so averaging works properly
        tempSC.Hmat = zeros(4,4,N);
        tempSC.NH = 1;
        SC(i)=tempSC;
    end

    % Average system compensation
    SCTot = mean(SC(i),2);
    % Feed in C vector
    SCTot.symRotVec = CvecSIPS;
    out = SCTot;
end
function [S1fTot] = filterStokes(S1tot,fwx)   
fwz = 1;
for i = 1:size(S1tot,5)

    S1 = squeeze(S1tot(:,:,:,:,i));
    %% Filter & Get DOP
    dim = size(S1);

    if dim(4) == 3 % 3-component Stokes vector was provided
        S1 = cat(4,sqrt(dot(S1,S1,4)),S1);
    end 
    
    nx = (round(fwx*1.5)-1)/2;
    nx = linspace(-nx,nx,round(fwx*1.5))*2*sqrt(log(2))/fwx;
    h = exp(-nx.^2);
    if fwz>1
        nz = (round(fwaxial*1.5)-1)/2;
        nz = linspace(-nz,nz,round(fwaxial*1.5))*2*sqrt(log(2))/fwaxial;
        h = exp(-nz(:).^2)*h;
    end
    h = h/sum(h(:));
    S1f = imfilter(S1,h,'circular');
    S1fTot(:,:,:,:,i) = S1f;
end
end

