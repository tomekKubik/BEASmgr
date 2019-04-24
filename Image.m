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
        imageOrigin; %  1-D 3 element matrix containg location of the image origin: [ox oy oz] skalowanie
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
        
        function [obj, macierz] = read2DImageFromScript(obj)
            %IMAGE Read image from file selected by user
            [fileName,filePath] = uigetfile('*.m','Select an image file');
            if  isequal([filePath fileName],0)
                disp('User selected Cancel');
            else
                disp(['User selected ', fullfile(filePath,fileName)]);
                % Wykonaj skrypt z pliku
                run(fullfile(filePath,fileName));
                macierz = BB_data_decim;
%dim
                [nx, ny, nt] = size(BB_data_decim);       
                obj.dim(1) = nx;
                obj.dim(2) = ny;
                obj.dim(3) = 1;
                obj.dim(4) = nt;
%voxelSize
                obj.voxelSize(1) = 100/obj.dim(1);
                obj.voxelSize(2) = 100/obj.dim(2);
                if obj.dim(3)==1
                    obj.voxelSize(3)=1;
                else
                    obj.voxelSize(3) = 100/obj.dim(3);  
                end
%time
                obj.time = linspace(1,60,obj.dim(4));
%imageOrgin
                obj.imageOrigin(1) = round(obj.dim(1)/2);
                obj.imageOrigin(2) = round(obj.dim(2)/2);
                obj.imageOrigin(3) = 1;
%voxels 
                obj.voxels = zeros(obj.dim(4),obj.dim(3),obj.dim(1),obj.dim(2));
                for t = 1:obj.dim(4)
                    for z = 1:obj.dim(3)
                        obj.voxels(t,z,:,:) = macierz(:,:,t);
                    end
                end 
            end
        end
        
        function [ox, oy, oz] = getRealFromMatrix(obj,ix, iy, iz)
                ox = obj.voxelSize(1)*ix-obj.voxelSize(1)*obj.imageOrgin(1);
                oy = obj.voxelSize(2)*iy-obj.voxelSize(2)*obj.imageOrgin(2);
                oz = obj.voxelSize(3)*iz-obj.voxelSize(3)*obj.imageOrgin(3);
        end
        
        function [ix, iy, iz] = getMatrixFromReal(obj,ox, oy, oz)
                ix = int(ox/obj.voxelSize(1))+obj.imageOrigin(1);
                iy = int(oy/obj.voxelSize(2))+obj.imageOrigin(2);
                iz = int(oz/obj.voxelSize(3))+obj.imageOrigin(3);
        end
        
        function showFilm(obj)
%                 fig = uifigure('Dane', Data);
%                 pnl = uipanel(fig);
%                 sld = uislider(pnl,'Position',[50 50 50 3]); % [left bottom width height]
%                 sld.Limits  = [0 obj.voxels(1)];
              for t = 1:obj.dim(4)
                    imagesc(abs(squeeze(obj.voxels(t,1,:,:))));
                    hold on 
                    plot(obj.imageOrigin(2),obj.imageOrigin(1),'r+', 'MarkerSize', 10);
                    pause(0.03)  
              end
        end
        
        
        function displayImageApp(obj)
            %TODO: Application
        end
    end
end
%uislider
%showImage - funkcja (scrolowanie po czasie, krzy�yk pokazuj�cy origina,
%wczytanie)
%gui + voxels - todo
%nie wczytywane s� do struktury obrazu warto�ci z obrazu