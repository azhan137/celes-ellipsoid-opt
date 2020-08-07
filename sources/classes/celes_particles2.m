classdef celes_particles2
    properties
       %type of particle (sphere or ellipsoid supported)
       type = 'sphere'
       
       %positions of all the particles in cartesian coordinates (x,y,z)
       %format [x(:),y(:),z(:)] nx3 array
       positionArray = zeros(1,3)
       
       %parameter array for particles
       %spheres nx2 array
       %    [:,1] = radii
       %    [:,2] = refractive indices
       
       %cylinder nx3 array
       %    [:,1] = radii
       %    [:,2] = height
       %    [:,3] = refractive index
       
       %ellipsoid nx5 array
       %    [:,1] = axis a
       %    [:,2] = axis b
       %    [:,3] = axis c
       %    [:,4] = rotation array
       %    [:,5] = refractive index
       
       parameterArray = ones(1,2)
       
       
       
       %maximum particle spacing
       maxParticleDistance
    end
    
    properties (Dependent)
        %total number of particles
        number
        
        %refractiveIndices of the particles
        refractiveIndexArray        
        
        %set of all unique particles defined in parameterArray
        %ordered by radius (sphere) or axis a (ellipsoid)
        uniqueParticles
        
        %number of unique particles defined in parameterArray
        numUniqueParticles
        
        %particle array indexed by the uniqueParticles
        particleArrayIndex
        
        %unique values of the particle array index
        uniqueParticleArrayIndex
        
        singleParticleArrayIndex
        
        singleUniqueParticleArrayIndex
    end
    
    methods
        % ======================================================================
        %> @brief Set method for type
        % ======================================================================
        function obj = set.type(obj,value)
            switch value
                case 'sphere'
                    obj.type = value;
                case 'ellipsoid'
                    obj.type = value;
                case 'cylinder'
                    obj.type = value;
                otherwise
                    error('this particle type is at the moment not implemented')
            end
        end
        
        % ======================================================================
        %> @brief Set method for positionArray
        % ======================================================================
        function obj = set.positionArray(obj,value)
            if length(value(1,:))==3
                obj.positionArray = single(value);
            else
                error('illegal position array')
            end
        end
        
        % ======================================================================
        %> @brief Set method for parameterArray
        
        %parameter array for particles
        %spheres nx2 array
        %    [:,1] = radii
        %    [:,2] = refractive indices
       
        %ellipsoid nx5 array
        %    [:,1] = axis a (x)
        %    [:,2] = axis b (y)
        %    [:,3] = axis c (z)
        %    [:,4] = rotation array
        %    [:,5] = refractive index
        
        % ======================================================================
        function obj = set.parameterArray(obj,value)
            obj.parameterArray = single(value);
        end
        
        % ======================================================================
        %> @brief Get method for particle number
        % ======================================================================
        function value = get.number(obj)
            value=length(obj.positionArray(:,1));
        end
        
        % ======================================================================
        %> @brief Get method for particle number
        % ======================================================================
        function value = get.refractiveIndexArray(obj)
            value=obj.parameterArray(:,end);
        end
        
        % ======================================================================
        %> @brief Get method for unique particle values, returns ordered
        %vector of unique particle parameters
        % ======================================================================
        function value = get.uniqueParticles(obj)
            value=unique(obj.parameterArray,'rows');
        end
        
        % ======================================================================
        %> @brief Get method for the number of unique particles
        % ======================================================================
        function value = get.numUniqueParticles(obj)
            value=length(obj.uniqueParticles(:,1));
        end  
        
        % ======================================================================
        %> @brief Get method for particle array in terms of indices given
        %by a four pair function. pairs a,b,c, and rotation angles
        % ======================================================================
        function value = get.particleArrayIndex(obj)
            switch obj.type
                case 'sphere'
                    value = dsearchn(obj.uniqueParticles(:,1),obj.parameterArray(:,1));
                case 'cylinder'
                    rs = dsearchn(obj.uniqueParticles(:,1),obj.parameterArray(:,1));
                    hs = dsearchn(obj.uniqueParticles(:,2),obj.parameterArray(:,2));
                    value = obj.pairIndices(rs,hs);
                case 'ellipsoid'
                    as = dsearchn(obj.uniqueParticles(:,1),obj.parameterArray(:,1));
                    bs = dsearchn(obj.uniqueParticles(:,2),obj.parameterArray(:,2));
                    cs = dsearchn(obj.uniqueParticles(:,3),obj.parameterArray(:,3));
                    rs = dsearchn(obj.uniqueParticles(:,4),obj.parameterArray(:,4));
                    ab = obj.pairIndices(as,bs);
                    cr = obj.pairIndices(cs,rs);
                    value = obj.pairIndices(ab,cr);
                otherwise
                    error('particle type not currently supported, sphere or ellipsoid only');
            end
        end 
        
        % ======================================================================
        %> @brief Get method for particle array in terms of indices given
        % to find unique particle indices
        % ======================================================================
        function value = get.uniqueParticleArrayIndex(obj)
            value = unique(obj.particleArrayIndex,'rows');
        end
        
        % ======================================================================
        %> @brief Get method for all particles by location in array
        % ======================================================================
        function value = get.singleParticleArrayIndex(obj)
            value = dsearchn(obj.uniqueParticleArrayIndex,obj.particleArrayIndex);
        end
        
        % ======================================================================
        %> @brief Get method for unique particles by location in array
        % used for calculation of mie coefficients/T-matrices
        % ======================================================================        
        function value = get.singleUniqueParticleArrayIndex(obj)
            value = unique(obj.singleParticleArrayIndex);
        end

        % ======================================================================
        %> @brief Set the maximalParticleDistance attribute to the correct value
        % ======================================================================
        function obj = compute_maximal_particle_distance(obj)
            %value=max(pdist(obj.positionArray));  pdist part of statistics and machine learning toolbox and might not be available
            obj.maxParticleDistance=0;
            for jp1=1:obj.number
                diffs=bsxfun(@plus,obj.positionArray((jp1+1):end,:),-obj.positionArray(jp1,:));
                dists2 = diffs(:,1).^2+diffs(:,2).^2+diffs(:,3).^2;
                if max(dists2)>obj.maxParticleDistance^2
                    obj.maxParticleDistance=sqrt(max(dists2));
                end
            end
        end        
        
        %function for pairing indices (cantor method)
        function cantorPairedIndex = pairIndices(~,index1,index2)
            cantorPairedIndex = 1/2*(index1+index2).*(index1+index2+1)+index2;
        end
    end
end
       