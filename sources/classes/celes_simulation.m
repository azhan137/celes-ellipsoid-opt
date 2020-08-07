%  Copyright (c) 2017, Amos Egel (KIT), Lorenzo Pattelli (LENS)
%                      Giacomo Mazzamuto (LENS)
%  All rights reserved.
%
%  Redistribution and use in source and binary forms, with or without
%  modification, are permitted provided that the following conditions are met:
%
%  * Redistributions of source code must retain the above copyright notice, this
%    list of conditions and the following disclaimer.
%
%  * Redistributions in binary form must reproduce the above copyright notice,
%    this list of conditions and the following disclaimer in the documentation
%    and/or other materials provided with the distribution.
%
%  * Neither the name of the copyright holder nor the names of its
%    contributors may be used to endorse or promote products derived from
%    this software without specific prior written permission.
%
%  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
%  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
%  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
%  ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
%  LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
%  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
%  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
%  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
%  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
%  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
%  POSSIBILITY OF SUCH DAMAGE.

%> @file celes_simulation.m
% ======================================================================
%> @brief Central data structure of the celes software
%
%> The simulation class contains all input, intermediate results and 
%> output for one calculation.
% ======================================================================
classdef celes_simulation
    
    properties
        %> celes_input object which contains the parameters that specify
        %> the simulation geometry and initial field
        input
        
        %> celes_numerics object which contains the numerical settings
        numerics
        
        %> celes_tables object which contains lookup tables and other
        %> intermediate results
        tables
        
        %> celes_output object which contains the results of the
        %> simulation
        output
        
    end
    
    properties (Dependent)
        %> single array which contains a grid of distances used for the
        %> lookup of the spherical hankel function in the particle coupling
        lookupParticleDistances
    end
    
    methods
        % ======================================================================
        %> @brief get method for dependent property lookupParticleDistances
        % ======================================================================
        function value = get.lookupParticleDistances(obj)
            value = [0,0:obj.numerics.particleDistanceResolution:(obj.input.particles.maxParticleDistance + 3*obj.numerics.particleDistanceResolution)];  % add two zeros at beginning to allow for interpolation also in first segment
        end
        
        % ======================================================================
        %> @brief Evaluate the Mie coefficients
        %>
        %> @return celes_simulation object with updated mieCoefficients
        % ======================================================================
        function obj = computeMieCoefficients(obj)
            %fprintf(1,'compute Mie coefficients ...');
            switch obj.input.particles.type
                case 'sphere'
                    obj.tables.mieCoefficients = zeros(obj.input.particles.numUniqueParticles,obj.numerics.nmax,'single');
                    for u_i=1:obj.input.particles.numUniqueParticles
                        for tau=1:2
                            for l=1:obj.numerics.lmax
                                for m=-l:l
                                    jmult = multi2single_index(1,tau,l,m,obj.numerics.lmax);
                                    obj.tables.mieCoefficients(u_i,jmult) = T_entry(tau,l,obj.input.k_medium,obj.input.k_particle(1),obj.input.particles.uniqueParticles(u_i,1));
                                end
                            end
                        end
                    end
                    obj.tables.singleParticleArrayIndex = obj.input.particles.singleParticleArrayIndex;
                    obj.tables.particleType = 'sphere';
                case 'ellipsoid'
                    obj.tables.mieCoefficients = zeros(obj.input.particles.numUniqueParticles,obj.numerics.nmax,obj.numerics.nmax,'single');
                    for u_i=1:obj.input.particles.numUniqueParticles
                        [T,dT] = compute_T(obj.numerics.lmax,40,40,obj.input.particles.uniqueParticles(u_i,1:4),obj.input.mediumRefractiveIndex,obj.input.particles.refractiveIndexArray(1),obj.input.wavelength);
                        obj.tables.mieCoefficients(u_i,:,:) = T;
                        obj.tables.gradMieCoefficients(u_i,:,:,:) = dT;
                    end
                    obj.tables.singleParticleArrayIndex = obj.input.particles.singleParticleArrayIndex;
                    obj.tables.particleType = 'ellipsoid';
                otherwise
                    error('particle type not implemented')
            end
            fprintf(1,' done\n');
        end

        
        % ======================================================================
        %> @brief Evaluate the Mie coefficients (parallel)
        %>
        %> @return celes_simulation object with updated mieCoefficients
        % for the case of ellipsoids the T-matrix derivatives are also
        % returned
        % ======================================================================
        function obj = computeParallelMieCoefficients(obj)
            %fprintf(1,'compute Mie (Parallel) coefficients ...');
            switch obj.input.particles.type
                case 'sphere'
                    lmax = obj.numerics.lmax;
                    mieCoefficients = zeros(obj.input.particles.numUniqueParticles,obj.numerics.nmax,'single');
                    k_medium = obj.input.k_medium;
                    k_particle = obj.input.k_particle(1);
                    uniqueParticles = obj.input.particles.uniqueParticles;
                    nmax = obj.numerics.nmax;
                    parfor u_i=1:obj.input.particles.numUniqueParticles
                        mieSlice = ones(1,nmax);
                        for tau=1:2
                            for l=1:lmax
                                for m=-l:l
                                    jmult = multi2single_index(1,tau,l,m,lmax);
                                    mieSlice(1,jmult) = T_entry(tau,l,k_medium,k_particle,uniqueParticles(u_i,1));
                                end
                            end
                        end
                        mieCoefficients(u_i,:) = mieSlice;
                    end
                    obj.tables.mieCoefficients = mieCoefficients;
                    obj.tables.singleParticleArrayIndex = obj.input.particles.singleParticleArrayIndex;
                    obj.tables.particleType = 'sphere';
                case 'cylinder'
                    lmax = obj.numerics.lmax;
                    mediumIndex = obj.input.mediumRefractiveIndex;
                    particleIndex = obj.input.particles.refractiveIndexArray(1);
                    wavelength = obj.input.wavelength;
                    uniqueParticles = obj.input.particles.uniqueParticles(:,1:3);
                    mieCoefficients = zeros(obj.input.particles.numUniqueParticles,obj.numerics.nmax,obj.numerics.nmax);
                    gradMieCoefficients = zeros(obj.input.particles.numUniqueParticles,obj.numerics.nmax,obj.numerics.nmax,2);
                    
                    cylParams = uniqueParticles(:,1:2);
                    
                    parfor u_i = 1:obj.input.particles.numUniqueParticles
                        [T,dT] = compute_T_cyl(lmax,100,cylParams(u_i,:),mediumIndex,particleIndex,wavelength);
                    
                        mieCoefficients(u_i,:,:) = -T;
                        gradMieCoefficients(u_i,:,:,:) = -dT;
                    end
                    
                    obj.tables.mieCoefficients = mieCoefficients;
                    obj.tables.gradMieCoefficients = gradMieCoefficients;
                    obj.tables.singleParticleArrayIndex = obj.input.particles.singleParticleArrayIndex;
                    obj.tables.particleType = 'cylinder';
                 
                case 'ellipsoid'
                    lmax = obj.numerics.lmax;
                    mediumIndex = obj.input.mediumRefractiveIndex;
                    particleIndex = obj.input.particles.refractiveIndexArray(1);
                    wavelength = obj.input.wavelength;
                    uniqueParticles = obj.input.particles.uniqueParticles(:,1:4);
                    mieCoefficients = zeros(obj.input.particles.numUniqueParticles,obj.numerics.nmax,obj.numerics.nmax);
                    gradMieCoefficients = zeros(obj.input.particles.numUniqueParticles,obj.numerics.nmax,obj.numerics.nmax,4);
                    
                    parfor u_i = 1:obj.input.particles.numUniqueParticles                        
                        [T,dT] = compute_T(lmax,60,60,uniqueParticles(u_i,1:3),mediumIndex,particleIndex,wavelength);
                        [T,dT] = axial_rotation(lmax,T,dT(:,:,1:3),uniqueParticles(u_i,4));
                        mieCoefficients(u_i,:,:) = -T;
                        gradMieCoefficients(u_i,:,:,:) = -dT;
                    end
                    
                    obj.tables.mieCoefficients = mieCoefficients;
                    obj.tables.gradMieCoefficients = gradMieCoefficients;
                    obj.tables.singleParticleArrayIndex = obj.input.particles.singleParticleArrayIndex;
                    obj.tables.particleType = 'ellipsoid';
                otherwise
                    error('particle type not implemented')
            end
        end
        
        % ======================================================================
        %> @brief Prepare a lookup for the a and b coefficients for
        %> particle coupling
        %>
        %> @return celes_simulation object with updated translationTable
        % ======================================================================
        function obj = computeTranslationTable(obj)
            %fprintf(1,'compute translation table ...');
            obj.tables.translationTable = translation_table_ab(obj.numerics.lmax);
            %fprintf(1,' done\n');
        end
        
        % ======================================================================
        %> @brief Evaluate the initial field coefficients \f$a^S_{0,n}\f$
        %> of the initial field expansion around each particle:
        %> \f$\mathbf{E}_0=\sum_n a^S_{0,n}\mathbf{\Psi}^{(1)}_n\f$
        %>
        %> @return celes_simulation object with updated initialFieldCoefficients
        % ======================================================================
        function obj = computeInitialFieldCoefficients(obj)
            %fprintf(1,'compute initial field coefficients ...');
            if isfinite(obj.input.initialField.beamWidth) && obj.input.initialField.beamWidth
                fprintf(1,' Gaussian beam ...');
                if obj.input.initialField.normalIncidence
                    obj.tables.initialFieldCoefficients = initial_field_coefficients_wavebundle_normal_incidence(obj);
                else
                    error('this case is not implemented')
                end
            else % infinite or 0 beam width
                %fprintf(1,' plane wave ...');
                obj.tables.initialFieldCoefficients = initial_field_coefficients_planewave(obj);
            end
            %fprintf(1,' done\n');
        end
        
        % ======================================================================
        %> @brief Evaluate the power flux of the initial field
        %>
        %> @return celes_simulation object with updated initialFieldPower
        % ======================================================================
        function obj = computeInitialFieldPower(obj)
            %fprintf(1,'compute initial field power ...');
            if obj.input.initialField.normalIncidence
                obj.output.initialFieldPower = initial_power_wavebundle_normal_incidence(obj);
            else
                error('this case is not implemented')
            end
            %fprintf(1,' done\n');
        end
        
        % ======================================================================
        %> @brief Compute the scattered field coefficients b by iteratively 
        %> solving the linear system M*b=T*aI
        %>
        %> @param Optional: b0, initial guess for scattered field
        %> coefficients
        %> @return celes_simulation object with updated initialFieldPower
        % ======================================================================
        function obj = computeScatteredFieldCoefficients(obj,varargin)
            %fprintf(1,'compute scattered field coefficients ...');
            mmm = @(x) obj.masterMatrixMultiply(x);
            [b,convHist] = obj.numerics.solver.run(mmm,obj.tables.rightHandSide(:),varargin{:});
            obj.tables.scatteredFieldCoefficients = reshape(gather(b),size(obj.tables.rightHandSide));
            obj.output.convergenceHistory = convHist;
        end
        
        % ======================================================================
        %> @brief Compute the plane wave pattern of the scattered field
        %> (i.e., the expansion coefficients of the scattered field in plane
        %> vector wave functions)
        %>
        %> @return celes_simulation object with updated
        %> output.scatteredFieldPlaneWavePattern
        % ======================================================================
        function obj = computeScatteredFieldPWP(obj)
            %fprintf(1,'compute scattered field plane wave pattern: ');
            obj.output.scatteredFieldPlaneWavePattern = scattered_field_plane_wave_pattern(obj);
            %fprintf(1,' ... done\n');
        end
        
        % ======================================================================
        %> @brief Compute the plane wave pattern of the total field
        %> (i.e., the expansion coefficients of the scattered field in plane
        %> vector wave functions)
        %>
        %> @return celes_simulation object with updated
        %> output.totalFieldPlaneWavePattern
        % ======================================================================
        function obj = computeTotalFieldPWP(obj)
            %fprintf(1,'compute total field coefficients table ...');
            obj.output.totalFieldPlaneWavePattern=cell(2,1);
            pwpScat = obj.output.scatteredFieldPlaneWavePattern;
            pwpIn = initial_field_plane_wave_pattern(obj);
            for pol=1:2
                obj.output.totalFieldPlaneWavePattern{pol} = pwpScat{pol};
                obj.output.totalFieldPlaneWavePattern{pol} = obj.output.totalFieldPlaneWavePattern{pol}.addTo(pwpIn{pol});
            end
            %fprintf(1,' done\n');
        end
        
        % ======================================================================
        %> @brief Evaluate the power flux of the total field, both in
        %> forward and in backward direction
        %>
        %> @return celes_simulation object with updated
        %> output.totalFieldForwardPower and output.totalFieldBackwardPower
        % ======================================================================
        function obj = computeTotalFieldPower(obj)
            %fprintf(1,'compute total field power ...');
            obj.output.totalFieldForwardPower = gather(pwp_power_flux(obj.output.totalFieldPlaneWavePattern{1},obj,'forward') + pwp_power_flux(obj.output.totalFieldPlaneWavePattern{2},obj,'forward'));
            obj.output.totalFieldBackwardPower = gather(pwp_power_flux(obj.output.totalFieldPlaneWavePattern{1},obj,'backward') + pwp_power_flux(obj.output.totalFieldPlaneWavePattern{2},obj,'backward'));
            %fprintf(1,' done\n');
        end
        
        % ======================================================================
        %> @brief First prepare the scattered and total field's plane wave
        %> pattern, then evaluate the power flux
        %>
        %> @return celes_simulation object with updated
        %> scatteredFieldPlaneWavePattern, totalFieldPlaneWavePattern,
        %> output.totalFieldForwardPower and output.totalFieldBackwardPower
        % ======================================================================
        function obj = evaluatePower(obj)
            tpow = tic;
            obj = obj.computeScatteredFieldPWP;
            obj = obj.computeTotalFieldPWP;
            obj = obj.computeTotalFieldPower;
            obj.output.powerEvaluationTime = toc(tpow);
            %fprintf(1,'power flux evaluated in %.1f seconds.\n',obj.output.powerEvaluationTime);
        end
        
        % ======================================================================
        %> @brief Evaluate the initial (near)field at the positions 
        %> specified in the input. The field can then be plotted.
        %>
        %> @return celes_simulation object with updated 
        %> output.InitialField
        % ======================================================================
        function obj = evaluateInitialField(obj)
            %fprintf(1,'evaluate initial field ...');
            obj.output.initialField = compute_initial_field(obj);
            %fprintf(1,' done\n');
        end
        
        % ======================================================================
        %> @brief Evaluate the scattered (near)field at the positions 
        %> specified in the input. The field can then be plotted.
        %>
        %> @return celes_simulation object with updated 
        %> output.scatteredField
        % ======================================================================
        function obj = evaluateScatteredField(obj)
            %fprintf(1,'evaluate scattered field ...');
            obj.output.scatteredField = compute_scattered_field(obj);
            %fprintf(1,' done\n');
        end
        
        % ======================================================================
        %> @brief Evaluate the internal (near)field at the positions 
        %> specified in the input. The field can then be plotted.
        %>
        %> @return celes_simulation object with updated 
        %> output.internalField
        % ======================================================================
        function obj = evaluateInternalField(obj)
            %fprintf(1,'evaluate internal field ...');
            [obj.output.internalField,obj.output.internalIndices] = compute_internal_field(obj);
            %fprintf(1,' done\n');
        end
        
        % ======================================================================
        %> @brief Evaluate both the initial and the scattered (near)field 
        %> at the positions specified in the input. The field can then be
        %> plotted.
        %>
        %> @return celes_simulation object with updated 
        %> output.initialField and output.scatteredField
        % ======================================================================
        function obj = evaluateFields(obj)
            tfld = tic;
            obj = obj.evaluateInitialField;
            obj = obj.evaluateScatteredField;
%            obj = obj.evaluateInternalField;
            obj.output.fieldEvaluationTime = toc(tfld);
            fprintf(1,'fields evaluated in %.1f seconds.\n',obj.output.fieldEvaluationTime);
        end
        
        % ======================================================================
        %> @brief Multiply the master matrix M=1-T*W to some vector x
        %>
        %> @param Vector x of incoming field SVWF coefficients
        %> @param verbose (logical, optional): If true (default), display detailed timing information
        %> @return Vector M*x
        % ======================================================================
        function Mx = masterMatrixMultiply(obj,value,varargin)
            value=value(:);
            
            if isempty(varargin)
                verbose=true;
            else
                verbose=varargin{1};
            end
            
            Wx=coupling_matrix_multiply(obj,value,varargin{:});
            
            if verbose
                %fprintf('apply T-matrix ... ')
            end
            %t_matrix_timer = tic;
            
            Wx=reshape(Wx,obj.input.particles.number,obj.numerics.nmax);
            switch obj.input.particles.type
                case 'sphere'
                    TWx = obj.tables.mieCoefficients(obj.input.particles.singleParticleArrayIndex,:).*Wx;
                case 'ellipsoid'
                    TWx = zeros(size(Wx));
                    mieCoefficients = obj.tables.mieCoefficients(obj.input.particles.singleParticleArrayIndex,:,:);
                    for u_i = 1:obj.input.particles.number
                        TWx(u_i,:) = squeeze(mieCoefficients(u_i,:,:))*gather(Wx(u_i,:)).';
                    end
                case 'cylinder'
                    TWx = zeros(size(Wx));
                    mieCoefficients = obj.tables.mieCoefficients(obj.input.particles.singleParticleArrayIndex,:,:);
                    for u_i = 1:obj.input.particles.number
                        TWx(u_i,:) = squeeze(mieCoefficients(u_i,:,:))*gather(Wx(u_i,:)).';
                    end
                otherwise
                    disp('particle type not supported, sphere or ellipsoid');
            end
            Mx = gather(value - TWx(:));
            
            %t_matrix_time = toc(t_matrix_timer);
            if verbose
                %fprintf(' done in %f seconds.\n', t_matrix_time)
            end
        end

        % ======================================================================
        %> @brief Multiply the master matrix M^T=1-W*T^T to some vector x
        %>
        %> @param Vector x of incoming field SVWF coefficients
        %> @param verbose (logical, optional): If true (default), display detailed timing information
        %> @return Vector M^T*x
        % ======================================================================
        function Mtx = masterMatrixMultiplyTranspose(obj,value)
            value = value(:);
            
            switch obj.input.particles.type
                case 'sphere'
                    mieCoeff = obj.tables.mieCoefficients(obj.input.particles.singleParticleArrayIndex,:);
                    WtTtx = coupling_matrix_multiply_T(obj,mieCoeff(:).*value(:));
                case 'cylinder'
                    mieCoefficients = obj.tables.mieCoefficients(obj.input.particles.singleParticleArrayIndex,:,:);
                    value = reshape(value,obj.input.particles.number,obj.numerics.nmax);
                    Ttx = zeros(size(value));
                    for u_i = 1:obj.input.particles.number
                        Ttx(u_i,:) = squeeze(mieCoefficients(u_i,:,:)).'*value(u_i,:).';
                    end
                    WtTtx = coupling_matrix_multiply_T(obj,Ttx(:));
                case 'ellipsoid'
                    mieCoefficients = obj.tables.mieCoefficients(obj.input.particles.singleParticleArrayIndex,:,:);
                    value = reshape(value,obj.input.particles.number,obj.numerics.nmax);
                    Ttx = zeros(size(value));
                    for u_i = 1:obj.input.particles.number
                        Ttx(u_i,:) = squeeze(mieCoefficients(u_i,:,:)).'*value(u_i,:).';
                    end
                    WtTtx = coupling_matrix_multiply_T(obj,Ttx(:));
                otherwise
                    disp('particle type not supported, spheres and ellipsoids only');
            end
           
            Mtx = gather(value(:)-WtTtx);
        end

        % ======================================================================
        %> @brief Evaluate the gradient of Mie coefficients with respect to
        %> radius of spheres
        %>
        %> @return celes_simulation object with updated gradMieCoefficients
        % ======================================================================
        function obj = computeGradMieCoefficients(obj)
            fprintf(1,'compute Mie coefficients derivatives ...');
            switch obj.input.particles.type
                case 'sphere'
                    obj.tables.gradMieCoefficients = zeros(obj.input.particles.numUniqueParticles(:,1),obj.numerics.nmax,'single');
                    for u_i=1:obj.input.particles.numUniqueParticles
                        for tau=1:2
                            for l=1:obj.numerics.lmax
                                for m=-l:l
                                    jmult = multi2single_index(1,tau,l,m,obj.numerics.lmax);
                                    obj.tables.gradMieCoefficients(u_i,jmult) = gradT_entry(tau,l,obj.input.k_medium,obj.input.k_particle(1),obj.input.particles.uniqueParticles(u_i,1));
                                end
                            end
                        end
                    end
                case 'ellipsoid'
                    disp('already computed');
                otherwise
                    error('particle type not implemented')
            end
            fprintf(1,' done\n');
        end
        
        function obj = computeParallelGradMieCoefficients(obj)
            fprintf(1,'compute Mie coefficients derivatives parallel ...');
            switch obj.input.particles.type
                case 'sphere'
                    gradMieCoefficients = zeros(obj.input.particles.numUniqueParticles(:,1),obj.numerics.nmax,'single');
                    lmax = obj.numerics.lmax;
                    k_medium = obj.input.k_medium;
                    k_particle = obj.input.k_particle(1);
                    uniqueParticles = obj.input.particles.uniqueParticles(:,1);  
                    nmax = obj.numerics.nmax;
                    parfor u_i=1:obj.input.particles.numUniqueParticles
                        gradMieSlice = ones(1,nmax);
                        for tau=1:2
                            for l=1:lmax
                                for m=-l:l
                                    jmult = multi2single_index(1,tau,l,m,lmax);
                                    gradMieSlice(1,jmult) =  gradT_entry(tau,l,k_medium,k_particle,uniqueParticles(u_i));
                                end
                            end
                        end
                        gradMieCoefficients(u_i,:) = gradMieSlice;
                    end
                    obj.tables.gradMieCoefficients = gradMieCoefficients;
                case 'ellipsoid'
                    disp('already computed');
                otherwise
                    error('particle type not implemented')
            end
            fprintf(1,'done\n');
        end
                    

        % ======================================================================
        %> @brief Compute the adjoint field coefficients lambda by iteratively 
        %> solving the linear system M^T*lambda=dI/db_i,n
        %>
        %> @param:
        %> rhs: computed rhs of the adjoint equation
        %> Optional: varargin, initial guess for scattered field
        %> coefficients
        %> 
        %> @return celes_simulation object with updated adjointFields
        % ======================================================================
        function obj = computeAdjointCoefficients(obj,rhs,varargin)
            fprintf(1,'compute adjoint coefficients ...');
            mmm = @(x) obj.masterMatrixMultiplyTranspose(x);
            [b,convHist] = obj.numerics.inverseSolver.run(mmm,rhs(:),varargin{:});
            obj.tables.adjointFields = reshape(gather(b),size((obj.tables.rightHandSide).'));
            obj.output.convergenceHistory = convHist;
        end      
        
        % ======================================================================
        %> @brief Compute the adjoint field coefficients lambda by iteratively 
        %> solving the linear system M^T*lambda=dI/db_i,n
        %>
        %> @param Optional: b0, initial guess for scattered field
        %> coefficients
        %> @return celes_simulation object with updated adjointFields
        %> 
        %> new method that calls 'compute_adjoint_rhs' -> more suited for
        %> every situation
        % ======================================================================
        %%designed for more than one point max sum_i(|E(r_i)|^2
        function obj = computeAdjointCoefficients2(obj,E1,points,varargin)
            fprintf(1,'compute adjoint coefficients ...');
            mmm = @(x) obj.masterMatrixMultiplyTranspose(x);
            rhs = compute_adjoint_rhs(obj,E1,points).';
            [b,convHist] = obj.numerics.inverseSolver.run(mmm,rhs(:),varargin{:});
            obj.tables.adjointFields = reshape(gather(b),size((obj.tables.rightHandSide).'));
            obj.output.convergenceHistory = convHist;
        end 
        
        % ======================================================================
        %> @brief Compute the adjoint field coefficients lambda by iteratively 
        %> solving the linear system M^T*lambda=dI/db_i,n
        %>
        %> @param Optional: b0, initial guess for scattered field
        %> coefficients
        %> @return celes_simulation object with updated adjointFields
        %> 
        %> old method that calls 'compute_bessel_value' -> slow and
        %> inefficient memory usage.
        %>
        %> two wavelengths (E1,E2, runs two forward simulations, and two inverse
        %> simulations 
        % ======================================================================
        function obj = computeChromaticAdjointCoefficients(obj,E1,E2,points,varargin)
            fprintf(1,'compute adjoint coefficients ...');
            mmm = @(x) obj.masterMatrixMultiplyTranspose(x);
            rhs = compute_adjoint_chromatic_derivative(obj,E1,E2,points).';
            [b,convHist] = obj.numerics.inverseSolver.run(mmm,rhs(:),varargin{:});
            obj.tables.adjointFields = reshape(gather(b),size((obj.tables.rightHandSide).'));
            obj.output.convergenceHistory = convHist;
        end
        
        % ======================================================================
        %> @brief Compute the adjoint field coefficients lambda by iteratively 
        %> solving the linear system M^T*lambda=dI/db_i,n
        %>
        %> @param Optional: b0, initial guess for scattered field
        %> coefficients
        %> @return celes_simulation object with updated adjointFields
        %> 
        %> new method that calls 'compute_adjoint_rhs' -> more suited for
        %> every situation
        % ======================================================================
        %minimizing squared error sum_i(I_0(r_i) - I_k(r_i))^2, where I_0
        %is some desired image (intensity), and I_k is the iterate
        function obj = computeImageAdjointCoefficients(obj,image,E,points,varargin)
            fprintf(1,'compute adjoint coefficients ...');
            mmm = @(x) obj.masterMatrixMultiplyTranspose(x);
            rhs = compute_image_rhs(obj,image,E,points).';
            [b,convHist] = obj.numerics.inverseSolver.run(mmm,rhs(:),varargin{:});
            obj.tables.adjointFields = reshape(gather(b),size((obj.tables.rightHandSide).'));
            obj.output.convergenceHistory = convHist;
        end  
        
	% ======================================================================
        %> @brief Compute the adjoint field coefficients lambda by iteratively 
        %> solving the linear system M^T*lambda=dI/db_i,n for polarization multiplexing
        %>
        %> @param Optional: b0, initial guess for scattered field
        %> coefficients
        %> @return celes_simulation object with updated adjointFields
        %> 
        %> new method that calls 'compute_adjoint_rhs' -> more suited for
        %> every situation
        % ======================================================================
        %minimizing squared error sum_i(I_0(r_i) - I_k(r_i))^2, where I_0
        %is some desired image (intensity), and I_k is the iterate
        function obj = computePolMultiplexImageAdjointCoefficients(obj,image_TE,E_TE,points_TE,image_TM,E_TM,points_TM,varargin)
            fprintf(1,'compute adjoint coefficients ...');
            mmm = @(x) obj.masterMatrixMultiplyTranspose(x);
            rhs = compute_image_rhs(obj,image_TE,E_TE,points_TE).'+...
                compute_image_rhs(obj,image_TM,E_TM,points_TM).';
            [b,convHist] = obj.numerics.inverseSolver.run(mmm,rhs(:),varargin{:});
            obj.tables.adjointFields = reshape(gather(b),size((obj.tables.rightHandSide).'));
            obj.output.convergenceHistory = convHist;
        end  
        
        %%same as above, but uses SSIM image similarity index (probably
        %%unstable)
        function obj = computeSSIMAdjointCoefficients(obj,img_ref,mean_ref,var_ref,E_iter,mean_iter,var_iter,covar,img_pts,varargin)
            fprintf(1,'compute adjoint coefficients ...');
            mmm = @(x) obj.masterMatrixMultiplyTranspose(x);
            rhs = compute_SSIM_rhs(obj,img_ref,mean_ref,var_ref,E_iter,mean_iter,var_iter,covar,img_pts).';
            [b,convHist] = obj.numerics.inverseSolver.run(mmm,rhs(:),varargin{:});
            obj.tables.adjointFields = reshape(gather(b),size((obj.tables.rightHandSide).'));
            obj.output.convergenceHistory = convHist;
        end          
        
        % ======================================================================
        %> @brief Compute the adjoint field coefficients lambda by iteratively 
        %> solving the linear system M^T*lambda=dI/db_i,n
        %>
        %> @param Optional: b0, initial guess for scattered field
        %> coefficients
        %> @return celes_simulation object with updated adjointFields
        %> 
        %> old method that calls 'compute_bessel_value' -> slow and
        %> inefficient memory usage.
        %>
        %> two wavelengths (E1,E2, runs two forward simulations, and two inverse
        %> simulations 
        %>
        %> attempts to optimize two different images at two different
        %> 
        % ======================================================================
        function obj = computeChromaticImageAdjointCoefficients(obj,simul2,E1,E2,points1,points2,varargin)
            fprintf(1,'compute adjoint coefficients ...');
            mmm = @(x) obj.masterMatrixMultiplyTranspose(x);
            rhs = compute_adjoint_chromatic_image_derivative(obj,simul2,E1,E2,points1,points2).';
            [b,convHist] = obj.numerics.inverseSolver.run(mmm,rhs(:),varargin{:});
            obj.tables.adjointFields = reshape(gather(b),size((obj.tables.rightHandSide).'));
            obj.output.convergenceHistory = convHist;
        end
      
        %>
        function obj = computeAbsVectorAdjointCoefficients(obj,E_0,E,points,fieldType,varargin)
            fprintf(1,'compute adjoint coefficients ...');
            mmm = @(x) obj.masterMatrixMultiplyTranspose(x);
            rhs = compute_abs_vector_rhs(obj,E_0,E,points,fieldType).';
            [b,convHist] = obj.numerics.inverseSolver.run(mmm,rhs(:),varargin{:});
            obj.tables.adjointFields = reshape(gather(b),size((obj.tables.rightHandSide).'));
            obj.output.convergenceHistory = convHist;
        end   
        
        
        % ======================================================================
        %> @brief Compute the adjoint field coefficients lambda by iteratively 
        %> solving the linear system M^T*lambda=dI/db_i,n
        %>
        %> @param Optional: b0, initial guess for scattered field
        %> coefficients
        %> @return celes_simulation object with updated adjointFields
        %> 
        %> new method that calls 'compute_adjoint_rhs' -> more suited for
        %> every situation
        % ======================================================================
        %minimizing squared error sum_i(I_0(r_i) - I_k(r_i))^2, where I_0
        %is some desired image (intensity), and I_k is the iterate
        function obj = computeEfieldImageAdjointCoefficients(obj,image,E,points,varargin)
            fprintf(1,'compute adjoint coefficients ...');
            mmm = @(x) obj.masterMatrixMultiplyTranspose(x);
            %rhs = compute_image_e_rhs(obj,image,E,points).';
            rhs = compute_image_field_rhs(obj,image,E,points).';
            [b,convHist] = obj.numerics.inverseSolver.run(mmm,rhs(:),varargin{:});
            obj.tables.adjointFields = reshape(gather(b),size((obj.tables.rightHandSide).'));
            obj.output.convergenceHistory = convHist;
        end  
        
%         % ======================================================================
%         %> @brief calculate right hand side of optimization equation using
%         %> stored scattered field coefficients
%         %> 
%         %> @return Vector grad_T*W*b+grad_T*aI
%         % ======================================================================     
%         function optRHS = rightHandSideOpt(obj,particleNumber)
%             scattCoeff = obj.tables.scatteredFieldCoefficients;
%             Wb=coupling_matrix_multiply(obj,scattCoeff(:));
%             Wb=reshape(Wb,obj.input.particles.number,obj.numerics.nmax);
%             tempGradMie = zeros(obj.input.particles.number,obj.numerics.nmax);
%             tempGradMie(particleNumber,:) = obj.tables.gradMieCoefficients(particleNumber,:);
%             TWb = tempGradMie.*Wb;
%             optRHS = TWb+obj.tables.rightHandSideInitial(particleNumber);
%         end        

%       % ======================================================================
%       % > @brief calculate right hand side of optimization equation using
%       % > stored scattered field coefficients
%       % > 
%       % > @return Vector grad_T*W*b+grad_T*aI
%       % ======================================================================     
        function optRHS = rightHandSideOpt(obj,particleNumber,particleIndex,varargin)
            Wb=gather(obj.tables.coupledScatteringCoefficients);
            switch obj.input.particles.type
                case 'sphere'
                    tempGradMie = zeros(obj.input.particles.number,obj.numerics.nmax);
                    tempGradMie(particleNumber,:) = obj.tables.gradMieCoefficients(particleIndex,:);
                    TWb = tempGradMie.*Wb;
                    optRHS = TWb+gather(obj.tables.rightHandSideInitial(particleNumber,particleIndex));
                case 'ellipsoid'
                    tempGradMie = squeeze(obj.tables.gradMieCoefficients(particleIndex,:,:,varargin{1}));
                    TWb = zeros(obj.input.particles.number,obj.numerics.nmax);
                    TWb(particleNumber,:) = tempGradMie*Wb(particleNumber,:).';
                    optRHS = TWb+gather(obj.tables.rightHandSideInitial(particleNumber,particleIndex,varargin{1}));
                case 'cylinder'
                    tempGradMie = squeeze(obj.tables.gradMieCoefficients(particleIndex,:,:,varargin{1}));
                    TWb = zeros(obj.input.particles.number,obj.numerics.nmax);
                    TWb(particleNumber,:) = tempGradMie*Wb(particleNumber,:).';
                    optRHS = TWb+gather(obj.tables.rightHandSideInitial(particleNumber,particleIndex,varargin{1}));
                otherwise
                    disp('Particle type not supported, spheres and ellipsoids only');
            end
            
        end        
        
        % ======================================================================
        %> @brief compute coupled scattering coefficients and store them
        %>
        %> @return celes_simulation object with updated coupledScatteringCoefficients
        %> represented as W*b in equation
        % ======================================================================
        function obj = computeCoupledScatteringCoefficients(obj)
            scattCoeff = obj.tables.scatteredFieldCoefficients;
            obj.tables.coupledScatteringCoefficients = reshape(coupling_matrix_multiply(obj,scattCoeff(:)),obj.input.particles.number,obj.numerics.nmax);
        end
        
        
        % ======================================================================
        %> @brief Compute the scattered field coefficients b by iteratively 
        %> solving the linear system M*b=T*aI
        %>
        %> @param Optional: b0, initial guess for scattered field
        %> coefficients
        %> @return celes_simulation object with updated initialFieldPower
        
        % for NONADJOINT calculations!!!
        % ======================================================================
        function obj = computeGradScatteredFieldCoefficients(obj,varargin)
            fprintf(1,'compute grad scattered field coefficients ...');
            mmm = @(x) obj.masterMatrixMultiply(x);
            gradScatFieldCoef = zeros(obj.input.particles.number,obj.numerics.nmax,obj.input.particles.number);
            particleIndices = obj.input.particles.radiusArrayIndex;
            particleIndices = particleIndices(:);
            for p_i = 1:obj.input.particles.number
                particleIndex = particleIndices(p_i);
                rhs = obj.rightHandSideOpt(p_i,particleIndex);
                [grad_b,convHist] = obj.numerics.solver.run(mmm,rhs(:),varargin{:});
                grad_b = reshape(gather(grad_b),size(obj.tables.rightHandSide));
                obj.output.convergenceHistory = convHist;
                gradScatFieldCoef(:,:,p_i) = grad_b;
            end
            
            obj.tables.gradScatteredFieldCoefficients = gradScatFieldCoef;
        end        
        % ======================================================================
        %> @brief Run the simulation.
        %> 
        %> A simulation run includes:
        %> - computation of initial field power
        %> - computation of Mie coefficients
        %> - computation of the translation table
        %> - computation of the maximal distance between pairs of particles
        %> - preparation of the particle partitioning (if blockdiagonal
        %>   preconditioner is active)
        %> - preparation of the blockdiagonal preconditioner (if active)
        %> - computation of initial field coefficients
        %> - solution of linear system
        %> 
        %> @param Optional: Initial guess b0 for the scattered field
        %> coefficients vector b
        %> @return celes_simulation object with various fields updated
        % ======================================================================
        function obj = run(obj,varargin)
            print_logo
            print_parameters(obj)
            tcomp=tic;
            % cuda_compile(obj.numerics.lmax);
            fprintf(1,'starting simulation.\n');
            obj = obj.computeInitialFieldPower;
            obj = obj.computeMieCoefficients;
            obj = obj.computeTranslationTable;
            fprintf(1,'compute maximal particle distance ...');
            obj.input.particles = obj.input.particles.compute_maximal_particle_distance;
            fprintf(1,' done\n');
            tprec=tic;
            if strcmp(obj.numerics.solver.preconditioner.type,'blockdiagonal')
                fprintf(1,'make particle partition ...');
                partitioning = make_particle_partion(obj.input.particles.positionArray,obj.numerics.solver.preconditioner.partitionEdgeSizes);
                obj.numerics.solver.preconditioner.partitioning = partitioning;
                fprintf(1,' done\n');
                obj = obj.numerics.solver.preconditioner.prepare(obj);
            end
            obj.output.preconiditionerPreparationTime = toc(tprec);
            fprintf(1,'preconditioner prepared in %.1f seconds.\n',obj.output.preconiditionerPreparationTime);
            tsolv=tic;
            obj = obj.computeInitialFieldCoefficients;
            obj = obj.computeScatteredFieldCoefficients(varargin{:});
            obj.output.solverTime = toc(tsolv);
            fprintf(1,'solver terminated in %.1f seconds.\n',obj.output.solverTime);
%             obj.numerics.solver.preconditioner.factorizedMasterMatrices = []; % clear memory intensive fields
%             obj.numerics.solver.preconditioner.masterMatrices = [];
            obj.output.runningTime = toc(tcomp);
            fprintf(1,'simulation ran in %.1f seconds.\n',obj.output.runningTime);
        end
    end
end
