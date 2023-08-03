%systemCompensation is a matlab class to manage the various system
%compensation settings, provide read, save, and visualization routines, and
%manage the retrieval of the compensation parameters over multiple frames.

classdef systemCompensation
    properties
        symRotVec % rotation vector to make SO3 matrices D-transpose symmetric with a left-side multiplication; 3 x Nbins
        
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
            obj.alignSumProd = (obj1.alignSumProd*obj1.NSumProd + obj2.alignSumProd*obj2.NSumProd)/(obj1.NSumProd + obj2.NSumProd);
            obj.NSumProd = obj1.NSumProd + obj2.NSumProd;
            obj.symRotVec = decomposeRot(euclideanRotation(makeRot3x3(obj1.symRotVec)+makeRot3x3(obj2.symRotVec)));
            [~,obj] = compensateSystem([],[],obj);
        end
        
        function out = mean(objArray,dim)
            if nargin<2 
                if size(objArray,1)==1
                    dim = 2;
                else
                    dim = 1;
                end
            end
            for ind = 1:size(objArray,setdiff(1:2,dim))
                alignSumProd = zeros(size(objArray(1).alignSumProd));
                NSumProd = 0;
                for jnd = 1:size(objArray,dim)
                    symRotVec(:,:,jnd,:) = makeRot3x3(objArray(ind,jnd).symRotVec);
                    alignSumProd = alignSumProd + objArray(ind,jnd).alignSumProd*objArray(ind,jnd).NSumProd;
                    NSumProd = NSumProd + objArray(ind,jnd).NSumProd;
                end
                loc = systemCompensation;
                loc.alignSumProd = alignSumProd/NSumProd;
                loc.NSumProd = NSumProd;
                loc.symRotVec = decomposeRot(euclideanRotation(mean(symRotVec,3)));
                [~,loc] = compensateSystem([],[],loc);
                out(ind) = loc;
            end
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
            
            
            subplot(1,3,1)
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
            
            subplot(1,3,2)
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

            subplot(1,3,3)
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

