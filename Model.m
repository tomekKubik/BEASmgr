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
        bspline;
        image Image
        app app1
        neighbourhood
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
            obj.neighbourhood = [5,5];
        end
        
        function [xunit, yunit, obj] = create2DModelBasedOnElipse(obj,image,x,y,z,ra,rb,ntheta,nphi,valueInMM)
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
            
            knots = [obj.phis];
            knots = [-knots(2) knots knots(length(knots))+knots(2)];
            weights = obj.rs;
            obj.bspline = fastBSpline(knots,weights);
            obj.neighbourhood = createNeighbourhood(obj, valueInMM);
            test = 1;
        end
       
        function [Xm, Ym] = calculateCoordinatesOfTheNode2D(obj,n)
                r = evalAt(obj.bspline,obj.phis(n));
                Xm = r*cos(obj.phis(n));
                Ym = r*sin(obj.phis(n));
        end
        
        function [R, Phi] = calculateCartesianToPolar(obj, Xm,Ym)
            R = sqrt((Xm-obj.middle(1)).^2 + (Ym-obj.middle(2)).^2);
            Phi = atan((Ym-obj.middle(2))/(Xm-obj.middle(1))); 
        end
             
        function location = isPixelInOrOut(obj, x, y)
            epsi = 10^(-5);
            [R, Phi] = calculateCartesianToPolar(obj, x,y);
            if R == 0
                Phi = 0;   
            end
            modelR = evalAt(obj.bspline,Phi);
            location = 1/2*(1+2/pi*atan((R-modelR)/epsi));          
        end
        
        function [u,v] = calculateAvrVOxelsIntensiti(obj, n, neighbourhood) %s¹siedztwo + pêtla po ca³oœci 
            [xn, yn] = calculateCoordinatesOfTheNode2D(obj,n);
            sumInside = 0;
            sumOutside = 0;
            weightsInside = 0;
            weightsOutside = 0;
            xn=round(xn);
            yn=round(yn);
            for x=xn-neighbourhood(2):1:xn+neighbourhood(2)
                for y=yn-neighbourhood(1):1:yn+neighbourhood(1)
                    if x<1 || y<1 || x>obj.image.dim(1) || y>obj.image.dim(2)
                        continue
                    end
                    hv = isPixelInOrOut(obj,x,y);
                    sumInside = sumInside + obj.image.voxels(1,1,y,x)*hv;
                    sumOutside = sumOutside + obj.image.voxels(1,1,y,x)*(1-hv);
                    weightsInside = weightsInside + hv;
                    weightsOutside = weightsOutside + (1-hv);
                end
            end
            u = sumInside/weightsInside; 
            v = sumOutside/weightsOutside;
        end
         
        function fInside = fIn(obj,n,u, neighbourhood)
            [xn, yn] = calculateCoordinatesOfTheNode2D(obj,n);
            fInside = 0;
            xn = round(xn);
            yn = round(yn);
            for x=xn-neighbourhood(2):1:xn+neighbourhood(2)
                for y=yn-neighbourhood(1):1:yn+neighbourhood(1)
                     if x<1 || y<1 || x>obj.image.dim(1) || y>obj.image.dim(2)
                        continue
                    end
                    fInside = fInside + (obj.image.voxels(1,1,y,x) - u)^2;
                end
            end
        end
            
           
        function fOutside = fOut(obj,n, v, neighbourhood)
            [xn, yn] = calculateCoordinatesOfTheNode2D(obj,n);
            xn = round(xn);
            yn = round(yn);
            fOutside = 0;
            for x=xn-neighbourhood(2):1:xn+neighbourhood(2)
                for y=yn-neighbourhood(1):1:yn+neighbourhood(1)
                     if x<1 || y<1 || x>obj.image.dim(1) || y>obj.image.dim(2)
                        continue
                    end
                    fOutside = fOutside + (obj.image.voxels(1,1,y,x) - v)^2;
                end
            end
        end
        
        function energy = energyOfModel(obj)   %w pêtli po wêz³ach
            energy = 0;
            for n = 1:length(obj.phis)
                [u, v] = calculateAvrVOxelsIntensiti(obj, n, obj.neighbourhood);
                fInside = fIn(obj,n,u, obj.neighbourhood);
                fOutside = fIn(obj,n,v, obj.neighbourhood);
                energy = energy + fInside + fOutside;
            end
            disp(['Energy of the model is: ', num2str(energy)])
        end
        
        function neighbourhood = createNeighbourhood(obj, valueInMM)
            neighbourhood = [round(valueInMM/obj.image.voxelSize(2)), round(valueInMM/obj.image.voxelSize(1))];  
            
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
            obj.phis
        end
        
        function segResult = getModelImage(obj)
            
        end
    end
end

