classdef Model
    %MODEL Segmentation model
    %   Class represents the segmentation model streched over radially
    %   equal nodes and approximated by Bsplines.
    
    properties (GetAccess = private)
        middle; % Middle point - should be 1D array.
        r;
    end
    
    methods
        function obj = Model(x,y,z)
            %MODEL Construct an instance of this class
            %   Model is initialized with its middle coordinates.
            obj.middle = [x y z];
        end
        
        function [mid] = getMiddle(obj)
            mid = zeros(1,3);
            mid(1) = obj.middle(1);
            mid(2) = obj.middle(2);
            mid(3) = obj.middle(3);
        end
    end
end

