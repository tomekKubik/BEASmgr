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
        neighbourhood
    end
    
    methods (Access = public)
        function obj = Model(x,y,z,ntheta,nphi)
            %MODEL Construct an instance of this class
            %   Model is initialized with its middle coordinates and with the specific radial
            %   spacing (ntheta = 0 means model 2D).
            obj.middle = [y x z];
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
            obj.middle(Image.DIR_X) = ix;
            obj.middle(Image.DIR_Y) = iy;
            
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
            
            xunit = (ra * cos(obj.phis) / image.voxelSize(Image.DIR_X)) + ix;
            yunit = (rb * sin(obj.phis) / image.voxelSize(Image.DIR_Y)) + iy;
            
            knots = [obj.phis];
            knots = [-knots(2) knots knots(length(knots))+knots(2)];
            weights = obj.rs;
            obj.bspline = fastBSpline(knots,weights);
            obj.neighbourhood = createNeighbourhood(obj, valueInMM);
        end
        
        function [Xm, Ym] = calculateCoordinatesOfTheNode2D(obj,n)
            r = evalAt(obj.bspline,obj.phis(n));
            Xm = obj.middle(Image.DIR_X)+r*cos(obj.phis(n))/obj.image.voxelSize(Image.DIR_X);
            Ym = obj.middle(Image.DIR_Y)+r*sin(obj.phis(n))/obj.image.voxelSize(Image.DIR_Y);
        end
        
        function [R, Phi] = calculateCartesianToPolar(obj, Xm,Ym)
            diffX = (Xm-obj.middle(Image.DIR_X))*obj.image.voxelSize(Image.DIR_X);
            diffY = (Ym-obj.middle(Image.DIR_Y))*obj.image.voxelSize(Image.DIR_Y);
            if diffY>=0 && diffX>=0
                shift=0;
            elseif diffX<0
                shift=pi;
            else
                shift=2*pi;
            end
            R = sqrt(diffX^2 + diffY^2);
            Phi = atan(diffY/diffX)+shift;
        end
        
        function location = isPixelInOrOut(obj, x, y)
            epsi = 10^(-5);
            [R, Phi] = calculateCartesianToPolar(obj, x,y);
            if R == 0
                Phi = 0;
            end
            modelR = evalAt(obj.bspline,Phi);
            location = 1/2*(1+2/pi*atan((modelR-R)/epsi));
        end
        
        function [u,v] = calculateAvrVOxelsIntensiti(obj, n) %s�siedztwo + p�tla po ca�o�ci
            [xn, yn] = calculateCoordinatesOfTheNode2D(obj,n);
            sumInside = 0;
            sumOutside = 0;
            weightsInside = 0;
            weightsOutside = 0;
            xn=round(xn);
            yn=round(yn);
            for x=xn-obj.neighbourhood(Image.DIR_X):1:xn+obj.neighbourhood(Image.DIR_X)
                for y=yn-obj.neighbourhood(Image.DIR_Y):1:yn+obj.neighbourhood(Image.DIR_Y)
                    if x<1 || y<1 || x>obj.image.dim(Image.DIR_X) || y>obj.image.dim(Image.DIR_Y)
                        continue
                    end
                    hv = isPixelInOrOut(obj,x,y);
                    sumInside = sumInside + obj.image.voxels(y,x,1,1)*hv;
                    sumOutside = sumOutside + obj.image.voxels(y,x,1,1)*(1-hv);
                    weightsInside = weightsInside + hv;
                    weightsOutside = weightsOutside + (1-hv);
                end
            end
            u = sumInside/weightsInside;
            v = sumOutside/weightsOutside;
        end
        
        function fInside = fIn(obj,n,u)
            [xn, yn] = calculateCoordinatesOfTheNode2D(obj,n);
            fInside = 0;
            xn = round(xn);
            yn = round(yn);
            for x=xn-obj.neighbourhood(Image.DIR_X):1:xn+obj.neighbourhood(Image.DIR_X)
                for y=yn-obj.neighbourhood(Image.DIR_Y):1:yn+obj.neighbourhood(Image.DIR_Y)
                    if x<1 || y<1 || x>obj.image.dim(Image.DIR_X) || y>obj.image.dim(Image.DIR_Y)
                        continue
                    end
                    fInside = fInside + (obj.image.voxels(y,x,1,1) - u)^2;
                end
            end
        end
        
        
        function fOutside = fOut(obj,n, v)
            [xn, yn] = calculateCoordinatesOfTheNode2D(obj,n);
            xn = round(xn);
            yn = round(yn);
            fOutside = 0;
            for x=xn-obj.neighbourhood(Image.DIR_X):1:xn+obj.neighbourhood(Image.DIR_X)
                for y=yn-obj.neighbourhood(Image.DIR_Y):1:yn+obj.neighbourhood(Image.DIR_Y)
                    if x<1 || y<1 || x>obj.image.dim(Image.DIR_X) || y>obj.image.dim(Image.DIR_Y)
                        continue
                    end
                    fOutside = fOutside + (obj.image.voxels(y,x,1,1) - v)^2;
                end
            end
        end
        
        function energy = energyOfModel(obj)   %w p�tli po w�z�ach
            energy = 0;
            for n = 1:length(obj.phis)
                [u, v] = calculateAvrVOxelsIntensiti(obj, n);
                fInside = fIn(obj,n,u);
                fOutside = fIn(obj,n,v);
                energy = energy + fInside + fOutside;
            end
            disp(['Energy of the model is: ', num2str(energy)])
        end
        
        function neighbourhood = createNeighbourhood(obj, valueInMM)
            neighbourhood = [round(valueInMM/obj.image.voxelSize(Image.DIR_Y)), round(valueInMM/obj.image.voxelSize(Image.DIR_X))];
            
        end
        
        function gradient = gradientOfModel(obj,n)
            [u, v] = calculateAvrVOxelsIntensiti(obj, n);
            fInside = fIn(obj,n,u);
            fOutside = fIn(obj,n,v);
            gradient = fInside - fOutside;
        end
        
        function testModel(obj)
            obj.rs
            obj.phis
        end
        
        function segResult = getModelImage(obj)
             segResult = zeros(obj.image.dim(Image.DIR_T),obj.image.dim(Image.DIR_Z),obj.image.dim(Image.DIR_Y),obj.image.dim(Image.DIR_X));
             for x = 1:obj.image.dim(Image.DIR_X)
                 for y = 1:obj.image.dim(Image.DIR_Y)
                     location = isPixelInOrOut(obj, x, y);
                     if location <= 0.5
                         segResult(:,:,y,x) = 1;
                     else
                         segResult(:,:,y,x) = 0;
                     end
                 end
             end
        end
        
        function nIter = runSegmentation(obj,iter,wIter,lambda,lambdaP,lambdaN)
            startEnergy = energyOfModel(obj);
            nWrongIter = 0;
            nIter = 0;
            for i=1:iter
                disp(['Iteration ' num2str(i)]);
                dsplmc = zeros(1,length(obj.rs));
                for n=1:length(dsplmc)
                    gradient = gradientOfModel(obj,n);
                    dsplmc(n) = gradient*lambda;
                    obj.rs(n) = obj.rs(n) + dsplmc(n);
                end
                    
                newEnergy = energyOfModel(obj);
                if newEnergy < startEnergy
                    startEnergy = newEnergy;
                    lambda = lambda*lambdaP;
                    nWrongIter = 0;
                else
                    obj.rs(n) = obj.rs(n) - dsplmc(n);
                    lambda = lambda/lambdaN;
                    nWrongIter = nWrongIter + 1;
                    if nWrongIter == wIter
                        break
                    end
                end
            nIter = nIter + 1;
            end
        end
    end
end



