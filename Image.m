classdef Image
    %IMAGE Representation of the image:
    %   voxel matrices
    %   image dimensions
    %   voxel size
    %   time vector
    %   image origin (point 0,0,0)
    
    properties (Constant)
        DIR_X = 2;
        DIR_Y = 1;
        DIR_Z = 3;
        DIR_T = 4;
    end
    
    properties (Access = public)
        time; % 1-D matrix containing time values in seconds
        dim; % 1-D 4 element matrix containing dimensions: [nx ny nz nt]
        voxelSize; % 1-D 3 element matrix containg voxel sizes in direction: [dx dy dz]
        imageOrigin; %  1-D 3 element matrix containg location of the image origin: [ox oy oz] skalowanie
        voxels; % 4-D matrix containg all image voxels: [time z x y];
    end
    
    methods (Access = public)
        function obj = Image(nx,ny,nz,nt)
            %IMAGE Constructs an empty image of the given dimensions.
            obj.voxels = zeros(ny,nx,nz,nt);
            obj.dim = [ny nx nz nt];
            obj.voxelSize = [1 1 1];
            obj.imageOrigin = [0 0 0];
            obj.time = 1:nt;
        end
        
        function obj = read2DImageFromScript(obj,filePath,fileName,isFilter,filterSize)
            %IMAGE Read image from file selected by user ////
            % Wykonaj skrypt z pliku
            run(fullfile(filePath,fileName))
            macierz = BB_data_decim;
            %dim
            obj = Image(1,1,1,1);
            [ny, nx, nt] = size(BB_data_decim);       %zm
            obj.dim(obj.DIR_X) = nx;
            obj.dim(obj.DIR_Y) = ny;
            obj.dim(obj.DIR_Z) = 1;
            obj.dim(obj.DIR_T) = nt;
            %voxelSize
            obj.voxelSize(obj.DIR_X) = voxS(1); %x
            obj.voxelSize(obj.DIR_Y) = voxS(2); %y
            if obj.dim(obj.DIR_Z)==1
                obj.voxelSize(obj.DIR_Z)=1;
            else
                obj.voxelSize(obj.DIR_Z) = voxS(3);
            end
            %time
            obj.time = linspace(1,60,obj.dim(obj.DIR_T));
            %imageOrgin
            obj.imageOrigin(obj.DIR_X) = round(obj.dim(obj.DIR_X)/2);
            obj.imageOrigin(obj.DIR_Y) = round(obj.dim(obj.DIR_Y)/2);
            obj.imageOrigin(obj.DIR_Z) = 1;
            %voxels
            tmpVoxels = zeros(obj.dim(obj.DIR_T),obj.dim(obj.DIR_Z),obj.dim(obj.DIR_Y),obj.dim(obj.DIR_X));
            for t = 1:obj.dim(obj.DIR_T)
                for z = 1:obj.dim(obj.DIR_Z)
                    if isFilter==TRUE
                        tmpVoxels(t,z,:,:) = imgaussfilt(macierz(:,:,t),filterSize);
                    else
                        tmpVoxels(t,z,:,:) = macierz(:,:,t);
                end
            end
            obj.voxels = permute(tmpVoxels,[3 4 2 1]);
        end
        
        function [ox, oy, oz] = getRealFromMatrix(obj,ix, iy, iz)
            ox = obj.voxelSize(obj.DIR_X)*ix-obj.voxelSize(obj.DIR_X)*obj.imageOrigin(obj.DIR_X);
            oy = obj.voxelSize(obj.DIR_Y)*iy-obj.voxelSize(obj.DIR_Y)*obj.imageOrigin(obj.DIR_Y);
            oz = obj.voxelSize(obj.DIR_Z)*iz-obj.voxelSize(obj.DIR_Z)*obj.imageOrigin(obj.DIR_Z);
        end
        
        function [ix, iy, iz] = getMatrixFromReal(obj,ox, oy, oz)
            ix = (ox/obj.voxelSize(obj.DIR_X))+obj.imageOrigin(obj.DIR_X);
            iy = (oy/obj.voxelSize(obj.DIR_Y))+obj.imageOrigin(obj.DIR_Y);
            iz = (oz/obj.voxelSize(obj.DIR_Z))+obj.imageOrigin(obj.DIR_Z);
        end
        
        function obj = ImageReflectionDown(obj)
            ReflectedImage = zeros(obj.dim(obj.DIR_T),obj.dim(obj.DIR_Z),obj.dim(obj.DIR_Y)*2,obj.dim(obj.DIR_X));
            mat = permute(obj.voxels,[4 3 1 2]);
            for t = 1:obj.dim(obj.DIR_T)
                for z = 1:obj.dim(obj.DIR_Z)
                    matrix = squeeze(mat(t,z,:,:));
                    ReflectedImage(t,z,1:obj.dim(obj.DIR_Y),:) = matrix;
                    fmatrix = flipud(matrix);
                    ReflectedImage(t,z,obj.dim(obj.DIR_Y)+1:2*obj.dim(obj.DIR_Y),:) = fmatrix;
                end
            end
            ReflectedImage = permute( ReflectedImage,[3 4 2 1]);
            obj.dim(obj.DIR_Y)=2*obj.dim(obj.DIR_Y);
            obj.imageOrigin(obj.DIR_Y) = round(obj.dim(obj.DIR_Y)/2);
            obj.voxels = ReflectedImage;
        end
        
        function obj = ImageReflectionUp(obj)
            ReflectedImage = zeros(obj.dim(obj.DIR_T),obj.dim(obj.DIR_Z),obj.dim(obj.DIR_Y)*2,obj.dim(obj.DIR_X));
            mat = permute(obj.voxels,[4 3 1 2]);
            for t = 1:obj.dim(obj.DIR_T)
                for z = 1:obj.dim(obj.DIR_Z)
                    matrix = squeeze(mat(t,z,:,:));
                    fmatrix = flipud(matrix);
                    ReflectedImage(t,z,1:obj.dim(obj.DIR_Y),:) = fmatrix;
                    ReflectedImage(t,z,obj.dim(obj.DIR_Y)+1:2*obj.dim(obj.DIR_Y),:) = matrix;
                end
            end
            ReflectedImage = permute( ReflectedImage,[3 4 2 1]);
            obj.dim(obj.DIR_Y)=2*obj.dim(obj.DIR_Y);
            obj.imageOrigin(obj.DIR_Y) = round(obj.dim(obj.DIR_Y)/2);
            obj.voxels = ReflectedImage;
        end
        
        function obj = ImageReflectionRight(obj)
            ReflectedImage = zeros(obj.dim(obj.DIR_T),obj.dim(obj.DIR_Z),obj.dim(obj.DIR_Y),obj.dim(obj.DIR_X)*2);
            mat = permute(obj.voxels,[4 3 1 2]);
            for t = 1:obj.dim(obj.DIR_T)
                for z = 1:obj.dim(obj.DIR_Z)
                    matrix = squeeze(mat(t,z,:,:));
                    ReflectedImage(t,z,:,1:obj.dim(obj.DIR_X)) = matrix;
                    fmatrix = fliplr(matrix);
                    ReflectedImage(t,z,:,obj.dim(obj.DIR_X)+1:2*obj.dim(obj.DIR_X)) = fmatrix;
                end
            end
            ReflectedImage = permute( ReflectedImage,[3 4 2 1]);
            obj.dim(obj.DIR_X)=2*obj.dim(obj.DIR_X);
            obj.imageOrigin(obj.DIR_X) = round(obj.dim(obj.DIR_X)/2);
            obj.voxels = ReflectedImage;
        end
        
        function obj = ImageReflectionLeft(obj)
            ReflectedImage = zeros(obj.dim(obj.DIR_T),obj.dim(obj.DIR_Z),obj.dim(obj.DIR_Y),obj.dim(obj.DIR_X)*2);
            mat = permute(obj.voxels,[4 3 1 2]);
            for t = 1:obj.dim(obj.DIR_T)
                for z = 1:obj.dim(obj.DIR_Z)
                    matrix = squeeze(mat(t,z,:,:));
                    fmatrix = fliplr(matrix);
                    ReflectedImage(t,z,:,1:obj.dim(obj.DIR_X)) = fmatrix;
                    ReflectedImage(t,z,:,obj.dim(obj.DIR_X)+1:2*obj.dim(obj.DIR_X)) = matrix;
                end
            end
            ReflectedImage = permute( ReflectedImage,[3 4 2 1]);
            obj.dim(obj.DIR_X)=2*obj.dim(obj.DIR_X);
            obj.imageOrigin(obj.DIR_X) = round(obj.dim(obj.DIR_X)/2);
            obj.voxels = ReflectedImage;
        end
        
    end
end
