classdef Model
    %MODEL Segmentation model
    %   Class represents the segmentation model streched over radially
    %   equal nodes and approximated by Bsplines.
    
    properties (Access = public)
        middle; % Middle point - should be 1D array.
        % Model coordinates:
        thetas;
        phis;
        spacingTheta;
        spacingPhi;
        rs;
        bspline;
        image Image
        neighbourhood
        alfa = 1;
        beta = 1;
        gamma = 1;
    end
    
    methods (Access = public)
        function obj = Model(x,y,z,ntheta,nphi)
            %MODEL Construct an instance of this class
            %   Model is initialized with its middle coordinates and with the specific radial
            %   spacing (ntheta = 0 means model 2D).
            obj.middle = [y x z];
            obj.thetas = linspace(0,180,ntheta+1);
            obj.phis = linspace(0,360,nphi+1);
            obj.spacingTheta = pi/ntheta;
            obj.spacingPhi = 2*pi/nphi;
            obj.rs = ones(ntheta+1,nphi+1);
            obj.neighbourhood = [5,5];
        end
        
        function [xunit, yunit, obj, firstPhis] = create2DModelBasedOnElipse(obj,image,x,y,z,ra,rb,ntheta,nphi,valueInMM)
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
            firstPhis = linspace(0,2*pi,nphi+1);
            
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
        
        function obj = updateBsplineNodes(obj)
            obj.bspline = fastBSpline(obj.bspline.knots,obj.rs);
        end
        
        function obj = setWeights(obj, alfa, beta, gamma)
            obj.alfa = alfa;
            obj.beta = beta;
            obj.gamma = gamma;
        end
                
        function [Xm, Ym] = calculateCoordinatesOfTheNode2D(obj,n)
            r = evalAt(obj.bspline,obj.phis(n));
            Xm = obj.middle(Image.DIR_X)+r*cos(obj.phis(n))/obj.image.voxelSize(Image.DIR_X);
            Ym = obj.middle(Image.DIR_Y)+r*sin(obj.phis(n))/obj.image.voxelSize(Image.DIR_Y);
        end
        
        function [Xm, Ym] = calculateCoordinatesOfTheR(obj,n,r)
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
        
        function [u,v] = calculateAvrVOxelsIntensiti(obj, n, middle) %s�siedztwo + p�tla po ca�o�ci
            [xn, yn] = calculateCoordinatesOfTheR(obj,n,middle);
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
            if sumInside==0 && weightsInside==0
                u = 0;
            else
                u = sumInside/weightsInside;
            end
            if sumOutside==0 && weightsOutside==0
                v = 0;
            else
                v = sumOutside/weightsOutside;
            end
        end
        
        function fInside = fIn(obj,n,middle,u)
            [xn, yn] = calculateCoordinatesOfTheR(obj,n,middle);
            fInside = 0;
            xn = round(xn);
            yn = round(yn);
            for x=xn-obj.neighbourhood(Image.DIR_X):1:xn+obj.neighbourhood(Image.DIR_X)
                for y=yn-obj.neighbourhood(Image.DIR_Y):1:yn+obj.neighbourhood(Image.DIR_Y)
                    if x<1 || y<1 || x>obj.image.dim(Image.DIR_X) || y>obj.image.dim(Image.DIR_Y)
                        continue
                    end
                    location = isPixelInOrOut(obj, x, y);
                    fInside = fInside + location*((obj.image.voxels(y,x,1,1) - u)^2);
                end
            end
        end
        
        
        function fOutside = fOut(obj,n,middle,v)
            [xn, yn] = calculateCoordinatesOfTheR(obj,n,middle);
            xn = round(xn);
            yn = round(yn);
            fOutside = 0;
            for x=xn-obj.neighbourhood(Image.DIR_X):1:xn+obj.neighbourhood(Image.DIR_X)
                for y=yn-obj.neighbourhood(Image.DIR_Y):1:yn+obj.neighbourhood(Image.DIR_Y)
                    if x<1 || y<1 || x>obj.image.dim(Image.DIR_X) || y>obj.image.dim(Image.DIR_Y)
                        continue
                    end
                    location = isPixelInOrOut(obj, x, y);
                    fOutside = fOutside + (1-location)*(obj.image.voxels(y,x,1,1) - v)^2;
                end
            end
        end
        
        function elasticEnergy = elasticEnergyOfModel(obj)
            ee = zeros(1,length(obj.phis));
            for n = 1:length(obj.phis)
                nl = n-1;
                if nl<1
                    nl = length(obj.rs)-1;
                end
                nr = n+1;
                if nr>length(obj.rs)
                    nr = 2;
                end
                ee = (((obj.rs(nr)-2*obj.rs(n) +(obj.rs(nl)))/obj.spacingPhi^2)^2)/2;
            end
%            ee = normalizeEnergy(obj,ee);
            elasticEnergy = sum(ee);
        end
        
        function stiffnessEnergy = stiffnessEnergyOfModel(obj)
            se = zeros(1,length(obj.phis));
            for n = 1:length(obj.phis)
                nn = n-1;
                if nn<1
                    nn = length(obj.rs)-1;
                end
                se(n) = ((obj.rs(n)-obj.rs(nn))/obj.spacingPhi)^2;
            end
%            se = normalizeEnergy(obj,se);
            stiffnessEnergy = sum(se);
        end
        
        function imageEnergy = imageEnergyOfModel(obj,middleRs)   %w p�tli po w�z�ach
            ie = zeros(1,length(obj.phis));
            for n = 1:length(obj.phis)
                [u, v] = calculateAvrVOxelsIntensiti(obj,n,middleRs(n));
                fInside = fIn(obj,n,middleRs(n),u);
                fOutside = fOut(obj,n,middleRs(n),v);
                ie(n) = fInside + fOutside;
            end
%            ie = normalizeEnergy(obj,ie);
            imageEnergy = sum(ie);
        end
        
        function wholeEnergy = energyOfModel(obj, middleRs)
            energy = [0 0 0];
            energy(1) = imageEnergyOfModel(obj,middleRs);
            energy(2) = stiffnessEnergyOfModel(obj);
            energy(3) = elasticEnergyOfModel(obj);
            wholeEnergy = obj.gamma*energy(1)+obj.alfa*energy(2)+obj.beta*energy(3);
            disp(['Energy of the model is: ', num2str(wholeEnergy)])
        end
        
        function ne = normalizeEnergy(obj,energy)
            maxEnergy = max(energy);
            if energy~=0
                ne = energy/maxEnergy;
            else
                ne = maxEnergy;
            end
        end
        
        function neighbourhood = createNeighbourhood(obj, valueInMM)
            neighbourhood = [round(valueInMM/obj.image.voxelSize(Image.DIR_Y)), round(valueInMM/obj.image.voxelSize(Image.DIR_X))];
        end
        
        function gradients = modelGradient(obj,middle)
            gradients = zeros(1,length(obj.rs));
            dsplmc = zeros(3,length(obj.rs));
            for n=1:length(obj.rs)
                dsplmc(1,n) = imageGradient(obj,n, middle);
                dsplmc(2,n) = stiffnessGradient(obj,n);
                dsplmc(3,n) = elasticGradient(obj,n);
            end
            for e=1:3
                minD = min(dsplmc(e,:));
                maxD = max(dsplmc(e,:));
                maxVal = max([abs(minD),abs(maxD)]);
                if maxD~=0
                    dsplmc(e,:) = dsplmc(e,:)./maxVal;
                end
            end
            for n=1:length(obj.rs)
                gradients(n) = obj.gamma*dsplmc(1,n) + obj.alfa*dsplmc(2,n) + obj.beta*dsplmc(3,n);
            end
        end
        
        function gradient = imageGradient(obj,n,middle)
            [u, v] = calculateAvrVOxelsIntensiti(obj,n,middle(n));
            fInside = fIn(obj,n,middle,u);
            fOutside = fOut(obj,n,middle,v);
            gradient = fOutside - fInside;
        end
        
        function stiffGradient = stiffnessGradient(obj,n)
            nl = n-1;
            if nl<1
                nl = length(obj.rs)-1;
            end
            nr = n+1;
            if nr>length(obj.rs)
                nr = 2;
            end
            stiffGradient = ((obj.rs(nr) - obj.rs(n))/obj.spacingPhi) - ((obj.rs(n)-obj.rs(nl))/obj.spacingPhi);
        end
        
        function elastGradient = elasticGradient(obj,n)
            nl = n-1;
            if nl<1
                nl = length(obj.rs)-1;
            end
            nr = n+1;
            if nr>length(obj.rs)
                nr = 2;
            end
            elastGradient = (obj.rs(nr)-4*obj.rs(nr)+6*obj.rs(n)-4*obj.rs(nl)+obj.rs(nl))/(obj.spacingPhi)^4;
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
        
        function obj = runSegmentation(obj,iter,wIter,lambda,lambdaP,lambdaN)
            disp('Segmentation started.');
            startLambda = lambda;
            startEnergy = energyOfModel(obj,obj.rs);
            nWrongIter = 0;
            nIter = 0;
            for i=1:iter
                nIter = nIter + 1;
                disp(['Iteration = ' num2str(i) ', wrong iterations = ' num2str(nWrongIter) ', lambda = ' num2str(lambda)]);
                oldRs = obj.rs;
                dsplmc = modelGradient(obj,oldRs);
                for n=1:length(dsplmc)
                    move = dsplmc(n)*lambda;
                    obj.rs(n) = obj.rs(n) + move;
                    if obj.rs(n)<0
                        obj.rs(n)=0;
                    end
                end
                obj.rs = medfilt1(obj.rs);
                obj = updateBsplineNodes(obj);
                newEnergy = energyOfModel(obj,oldRs);
                if newEnergy < startEnergy
                    startEnergy = newEnergy;
                    if( lambda<startLambda )
                        lambda = startLambda;
                    else
                        lambda = lambda*lambdaP;
                    end
                    nWrongIter = 0;
                else
                    obj.rs = oldRs;
                    obj = updateBsplineNodes(obj);
                    if( lambda>startLambda )
                        lambda = startLambda;
                    else
                        lambda = lambda/lambdaN;
                    end
                    nWrongIter = nWrongIter + 1;
                    if nWrongIter == wIter
                        break
                    end
                end
            end
            disp(['Segmentation finished after ' num2str(nIter) ' iterations.']);
        end
    end
end



