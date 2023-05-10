%systemCompensation is a matlab class to manage the various system
%compensation settings, provide read, save, and visualization routines, and
%manage the retrieval of the compensation parameters over multiple frames.

classdef systemCompensation
    properties
        symRotVec % rotation vector to make SO3 matrices D-transpose symmetric with a left-side multiplication; 3 x Nbins
        Hmat % H matrix computed during the symmetrization process
        NH %number of data points used for Hmat
        symErrInit
        symErr
        
        alignRotVec % rotation vector to align to central bin, D*R.'*D*J*R, where R = makeRot(alignRotVec)
        alignSumProd % data matrix for computing alignRotVec
        NSumProd %number of data points used for alignSumProd
        alignErrInit
        alignErr
    end
    properties (Transient = true)
        initialized %bool to check if obj is read and initialized
    end

    methods
        function obj = systemCompensation(varargin)
            % systemCompensation constructs an empty systemCompensation
            % object; if filename is provided, it will read the file.

            obj.initialized = false;
            if nargin>0
                out = obj.read(varargin{:});
                if ~isempty(out)
                    obj = out;
                else
                    obj.initialized = false;
                end
            end
        end
        
        function out = read(~,fileName,Nbins)
        % read(fileName) reads the systemCompensation parameters from a .mat file
            if isfile(fileName)
                load(fileName,'loc');
                out = loc;
                out.initialized = true;
            elseif isfolder(fileName) % assume fileName is a directory, and we are looking for a .mat file
                list = ls(fullfile(fileName,'*SysCom.mat'));
                inds = regexp(list,'.mat');
                for ind = 1:numel(inds)
                    load(fullfile(fileName,list(1:inds(1)+3)),'loc');
                    out(ind) = loc;
                    Ndim(ind) = size(out(ind).symRotVec,2);% retrieve the number of bins
                end
                if numel(out)>1 && nargin>2
                    ind = find(Ndim==Nbins,1,'first');
                    if ~isempty(ind)
                        out = out(ind);
                        out.initialized = true;
                    else
                        out = [];
                    end
                else
                    out = out(1);
                    out.initialized = true;
                end
            else
                out = [];
            end                
        end
        function write(obj,fileName)
        % write(fileName) writes the current setting to a .mat file
            loc = obj;
            save(fileName,'loc');
        end
        
        function obj = plus(obj1, obj2)
            obj = systemCompensation();
            obj.Hmat = (obj1.Hmat*obj1.NH + obj2.Hmat*obj2.NH)/(obj1.NH + obj2.NH);
            obj.NH = obj1.NH + obj2.NH;
         
            obj.alignSumProd = (obj1.alignSumProd*obj1.NSumProd + obj2.alignSumProd*obj2.NSumProd)/(obj1.NSumProd + obj2.NSumProd);
            obj.NSumProd = obj1.NSumProd + obj2.NSumProd;
        end
        
        function out = mean(objArray,dim)
            if nargin<2 
                if size(dim,1)==1
                    dim = 2;
                else
                    dim = 1;
                end
            end
            for ind = 1:size(objArray,setdiff(1:2,dim))
                Hmat = zeros(size(objArray(1).Hmat));
                NH = 0;
                alignSumProd = zeros(size(objArray(1).alignSumProd));
                NSumProd = 0;
                for jnd = 1:size(objArray,dim)
                    Hmat = Hmat + objArray(ind,jnd).Hmat*objArray(ind,jnd).NH;
                    NH = NH + objArray(ind,jnd).NH;
                    alignSumProd = alignSumProd + objArray(ind,jnd).alignSumProd*objArray(ind,jnd).NSumProd;
                    NSumProd = NSumProd + objArray(ind,jnd).NSumProd;
                end
                loc = systemCompensation;
                loc.Hmat = Hmat/NH;
                loc.NH = NH;
                loc.alignSumProd = alignSumProd/NSumProd;
                loc.NSumProd = NSumProd;
                
                [~,loc] = compensateSystem([],[],loc);
                out(ind) = loc;
            end
        end
        
        function flipBool = checkForFlip(obj1,obj2)
            % check if obj1 and obj2 are in the same S1/S2 order, or not
            loc = obj1;
            locf = obj1.flipS1S2;
            % we want to apply the correction of template to the data of loc
            loc.symRotVec = obj2.symRotVec;
            loc.alignRotVec = obj2.alignRotVec;
            locf.symRotVec = obj2.symRotVec;
            locf.alignRotVec = obj2.alignRotVec;

            [~,loc] = compensateSystem([],[],loc,true);
            [~,locf] = compensateSystem([],[],locf,true);
            if mean(loc.symErr)<mean(locf.symErr)
                flipBool = false;
            else
                flipBool = true;
            end
        end
        
        function obj = flipS1S2(obj)
            % in some instances, we have the C and L matrices, but the next
            % measurement has inverted S1 and S2 columns. Here we convert 
            % the existing correction to this setting, by modifying the
            % Hmat and alignSumProd matrices and recomputing the
            % corresponding symRotVec and alignRotVec.
        
            P = [0,1,0;1,0,0;0,0,-1];% This describes the switch of S1 and S2, and resulting flip of S2, Mp = M*P
            A = [1 0 0 1; 1 0 0 -1; 0 1 1 0; 0 1i -1i 0]/sqrt(2);% used to compute H-matrix
            FtoH = @(x)reshape(x([1,9,3,11,5,13,7,15,2,10,4,12,6,14,8,16]),[4,4]);% used to compute H-matrix
            
            % this is the former symmetrization matrix
            Cold = makeRot3x3(obj.symRotVec);
        
            % compute new Hmat, corresponding to right multiplying the
            % original M with P
            for indw = 1:size(obj.Hmat,3)
                Hloc = obj.Hmat(:,:,indw);
                Hloc = FtoH(FtoH(Hloc)*A'*[1,0,0,0;zeros(3,1),P]*A);
                obj.Hmat(:,:,indw) = Hloc;
            end
            % re-compute the new M
            obj.symRotVec = [];
            [~,obj] = makeSymmetric([],[],obj);
            %[~,loc] = compensateSystem([],[],loc);

            Cnew = makeRot3x3(obj.symRotVec);
            % modify alignSumProd, which involves first correcting for the
            % effect of the former Cold, applying P, and the new Cnew
            for indw = 1:size(obj.Hmat,3)
                obj.alignSumProd(:,:,indw) = kron(P,Cnew(:,:,indw)*Cold(:,:,indw).')*obj.alignSumProd(:,:,indw)*kron(P,Cnew(:,:,indw)*Cold(:,:,indw).').';
            end
            obj.alignRotVec = [];
            [~,obj] = alignToCentralBin([],[],obj);
%            [~,obj] = compensateSystem([],[],obj);                  
        end
        
        function alignLocalRotVec = convertToLocalAlignment(obj)
            % use the sym and alignRotVec to generate the alignment
            % correction used when averaging spectral bins after taking the
            % depth-derivative in the local metric.
            C = makeRot3x3(obj.symRotVec);
            RT = pagetranspose(makeRot3x3(obj.alignRotVec)).*[1,1,-1;1,1,-1;-1,-1,1];
            RC = pagemtimes(RT,C);
            indwc = ceil(size(RC,3)/2);
            alignLocalRotVec = decomposeRot(pagemtimes(RC(:,:,indwc).',RC));% remove effect of correction on central bin
            alignLocalRotVec(:,indwc) = zeros(3,1);
        end
        
        function visualizeCompensation(obj,fh,lineStyle)
            %visualizeCompensation(figureHandle,lineStyle) opens a figure
            %(new, if no handle provided) and displays the system
            %compensation elements.
            if nargin<2
                figure;
                lineStyle = '-';
                holdonBool = false;
            elseif nargin<3
                figure(fh);
                clf
                holdonBool = false;
                lineStyle = '-';
            else
                figure(fh);
                holdonBool = true;
            end
           
            cm = lines;
            
            
            subplot(2,2,1)
            if holdonBool
                hold on
            end
            plot(obj.symRotVec(1,:),lineStyle,'color',cm(1,:))
            hold on
            plot(obj.symRotVec(2,:),lineStyle,'color',cm(2,:))
            plot(obj.symRotVec(3,:),lineStyle,'color',cm(3,:))
            xlabel('Spectral bins')
            ylabel('[rad]')
            title('Symmetrization rotation vector')
            
            subplot(2,2,2)
            if holdonBool
                hold on
            end
            plot(obj.alignRotVec(1,:),lineStyle,'color',cm(1,:))
            hold on
            plot(obj.alignRotVec(2,:),lineStyle,'color',cm(2,:))
            plot(obj.alignRotVec(3,:),lineStyle,'color',cm(3,:))
            xlabel('Spectral bins')
            ylabel('[rad]')
            title('Alignment rotation vector')

            subplot(2,2,3)
            if holdonBool
                hold on
            end
            plot(obj.symErrInit,lineStyle,'color',cm(1,:))
            hold on
            plot(obj.symErr,lineStyle,'color',cm(2,:))
            xlabel('Spectral bins')
            ylabel('Error per pixel')
            title('Symmetrization error')

            subplot(2,2,4)
            if holdonBool
                hold on
            end
            plot(obj.alignErrInit,lineStyle,'color',cm(1,:))
            hold on
            plot(obj.alignErr,lineStyle,'color',cm(2,:))
            xlabel('Spectral bins')
            ylabel('Error per pixel')
            title('Alignment error')
        end
        
    end
end

