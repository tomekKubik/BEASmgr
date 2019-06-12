classdef Model
    %MODEL Segmentation model
    %   Class represents the segmentation model streched over radially
    %   equal nodes and approximated by Bsplines.
    
    properties (Access = public)
        middle; % Middle point - should be 1D array.
        % Model coordinates:
        thetas;
        phis;
        rs;
        image Image
        app app1
    end
    
    methods (Access = public)
        function obj = Model(x,y,z,ntheta,nphi)
            %MODEL Construct an instance of this class
            %   Model is initialized with its middle coordinates and with the specific radial
            %   spacing (ntheta = 0 means model 2D).
            obj.middle = [x y z];
            obj.thetas = linspace(0,180,ntheta+1);
            obj.phis = linspace(0,360,nphi+1);
            obj.rs = ones(ntheta+1,nphi+1);
        end
        
        function [xunit, yunit, obj] = create2DModelBasedOnElipse(obj,image,x,y,z,ra,rb,ntheta,nphi)
            obj = Model(x,y,z,ntheta,nphi);
            obj.image = image;
      %      obj = Model (1,1,1,1,1);   
            [ix, iy, iz] = getMatrixFromReal(obj.image, x, y, z);
            obj.middle(1) = ix;
            obj.middle(2) = iy;
            
            if ntheta==1 
                obj.thetas = pi/2;
            else
                obj.thetas = linspace(0,pi,ntheta);
            end
            
            obj.phis = linspace(0,2*pi,nphi+1);
            
            obj.rs = ones(ntheta,nphi+1);
            b = (ra^2*(sin(obj.phis).^2)+rb^2*(cos(obj.phis).^2));
            for i = 1:nphi+1
                obj.rs(i) = sqrt((ra^2*rb^2)/b(i));
            end
            
            xunit = (ra * cos(obj.phis) / image.voxelSize(2)) + ix;
            yunit = (rb * sin(obj.phis) / image.voxelSize(1)) + iy;
             % TODO: Geometrycne s³upki
            
            
        end
        
        function [mid] = getMiddle(obj)
            mid = zeros(1,3);
            mid(1) = obj.middle(1);
            mid(2) = obj.middle(2);
            mid(3) = obj.middle(3);
       
        end
        
        function obj = setR(obj,nth,nph,r)
            rs(nth,nph) = r;
        end
        
        function r = getR(obj,nth,nph)
            rs(nth,nph);
        end
        
        function testModel(obj)
            obj
            obj.rs
        end
    end
end

