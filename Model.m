classdef Model
    %MODEL Segmentation model
    %   Class represents the segmentation model streched over radially
    %   equal nodes and approximated by Bsplines.
    
    properties
        middle; % Middle point - should be 1D array.
        % Model coordinates:
        thetas;
        phis;
        rs;
    end
    
    methods
        function obj = Model(x,y,z,ntheta,nphi)
            %MODEL Construct an instance of this class
            %   Model is initialized with its middle coordinates and with the specific radial
            %   spacing (ntheta = 0 means model 2D).
            obj.middle = [x y z];
            obj.thetas = linspace(0,180,ntheta+1);
            obj.phis = linspace(0,360,nphi+1);
            obj.rs = ones(ntheta+1,nphi+1);
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
    end
end

