classdef Image
    %IMAGE Representation of the image:
    %   voxel matrices
    %   image dimensions
    %   voxel size
    %   time vector
    %   image origin (point 0,0,0)
    
    properties
        time; % 1-D matrix containing time values in seconds
        dim; % 1-D 4 element matrix containing dimensions: [nx ny nz nt]
        voxelSize; % 1-D 3 element matrix containg voxel sizes in direction: [dx dy dz]
        imageOrigin; %  1-D 3 element matrix containg location of the image origin: [ox oy oz]
        voxels; % 4-D matrix containg all image voxels: [time z x y];
    end
    
    methods
        function obj = Image(nx,ny,nz,nt)
            %IMAGE Constructs an empty image of the given dimensions.
            obj.voxels = zeros(nt,nz,nx,ny);
            obj.dim = [nx ny nz nt];
            obj.voxelSize = [1 1 1];
            obj.imageOrigin = [0 0 0];
            obj.time = 1:nt;
        end
        
        function obj = readImage(imageFile)
            %IMAGE Read image from file selected by user
            [fileName,filePath] = uigetfile('*.*','Select and image');
            if  isequal(fileName,0)
                disp('User selected Cancel');
            else
                disp(['User selected ', fullfile(filePath,fileName)]);
                imageFile = imread(fileName);
            end 
        end
        
        function [ox oy oz] = getRealFromMatrix(ix, iy, iz)
            % TODO for Kasia: convert matrix coordinates into milimeters
        end;
        
        function [ox oy oz] = getMatrixFromReal(ix, iy, iz)
            % TODO for Kasia: convert matrix coordinates into milimeters
        end;
    end
end

