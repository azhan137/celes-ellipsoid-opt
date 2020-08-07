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

%> @file celes_numerics.m
% ======================================================================
%> @brief Parameters describing the numerical settings used in the simulation
% ======================================================================

classdef celes_numerics_v2
    properties
        %> maximal polar multipole order
        lmax
        
        %> array of polar angles (between 0 and pi) that specifies the
        %> discretization of angles in the plane wave patterns
        polarAnglesArray
        
        %> array of azimuthal angles (between 0 and 2*pi) that specifies 
        %> the discretization of angles in the plane wave patterns
        azimuthalAnglesArray
        
        %> shall the GPU be used in computations of, e.g., the
        %> preconditioner?
        %> NOTE: In the current version, the GPU is alway used for the
        %> matrix vector products. Future versions might offer the
        %> opportunity to run the whole calculation on the CPU
        gpuFlag = true
        
        %> celes_solver object that contains all the information about
        %> how the linear system is to be solved
        solver = celes_solver
        
        %> celes_solver object that contains information about how the
        %> optimization linear system is solved
        inverseSolver = celes_solver
        
        %> resolution of the radial particle distance in the lookup tables
        %> of the translation operator
        particleDistanceResolution = single(10);
        
        %> for type blockdiagonal: edge sizes of the cuboids that define
        %> the partition of the particles (i.e. the blocks)
        %> format: [dx,dy,dz]
        partitionEdgeSizes
        
        %> for type blockdiagonal: cell array of particle indices that are
        %> grouped into the same partition
        partitioning
        
        %> for type blockdiagonal: cell array of SVWF coefficient indices
        %> that belong to the corresponding partitioning
        partitioningIdcs
        
        %> for type blockdiagonal: cell array of blocks of coupling matrix
        %> W that corresponds to coupling matrix of one partitioning
        couplingMatrices
    end
    
    properties (Dependent)
        %> number of unknowns per particle
        nmax
    end
    
    methods
        % ======================================================================
        %> @brief prepare the coupling matrix
        %>
        %> @param simulation object
        %> @return W coupling matrix
        % ======================================================================
        function simul = prepareW(obj,simul)
            switch obj.solver.preconditioner.type
                case 'blockdiagonal'
                    fprintf(1,'prepare blockdiagonal preconditioner ...\n');
                    msg = '';
                    l_max = simul.numerics.lmax;
                    k = simul.input.k_medium;
                    
                    for jp=1:length(obj.partitioning)
                        spherArr=obj.partitioning{jp};
                        NSi = length(spherArr);
                        
                        Idcs= [];
                        for n=1:simul.numerics.nmax
                            Idcs = [Idcs;spherArr+simul.input.particles.number*(n-1)];
                        end
                        simul.numerics.partitioningIdcs{jp} = Idcs;
                        
                        fprintf(1,repmat('\b',[1,length(msg)]));
                        msg = sprintf('partition %i of %i: %i particles -- compute master matrix ...',jp,length(obj.partitioning),length(spherArr));
                        fprintf(1,msg);
                        
                        W = zeros(NSi*simul.numerics.nmax,'single');
                        
                        [x2,x1]=meshgrid(simul.input.particles.positionArray(spherArr,1));
                        [y2,y1]=meshgrid(simul.input.particles.positionArray(spherArr,2));
                        [z2,z1]=meshgrid(simul.input.particles.positionArray(spherArr,3));
                        x1mnx2 = x1-x2;
                        y1mny2 = y1-y2;
                        z1mnz2 = z1-z2;
                        dTab = sqrt(x1mnx2.^2+y1mny2.^2+z1mnz2.^2);
                        ctTab = z1mnz2./dTab;
                        stTab = sqrt(1-ctTab.^2);
                        phiTab = atan2(y1mny2,x1mnx2);
                        Plm = legendre_normalized_trigon(ctTab,stTab,2*l_max);
                        sphHank = zeros(2*l_max,NSi,NSi);
                        for p=0:2*l_max
                            sphHank(p+1,:,:) = sph_bessel(3,p,k*dTab);
                        end
                        
                        for tau1=1:2
                            for l1=1:l_max
                                for m1=-l1:l1
                                    n1=multi2single_index(1,tau1,l1,m1,l_max);
                                    n1S1Arr=(1:NSi)+(n1-1)*NSi;
                                    for tau2=1:2
                                        for l2=1:l_max
                                            for m2=-l2:l2
                                                n2=multi2single_index(1,tau2,l2,m2,l_max);
                                                n2S2Arr=(1:NSi)+(n2-1)*NSi;
                                                for p=max(abs(m1-m2),abs(l1-l2)+abs(tau1-tau2)):l1+l2
                                                    Wpn1n2 = squeeze(sphHank(p+1,:,:)).*Plm{p+1,abs(m1-m2)+1}.*simul.tables.translationTable.ab5(n2,n1,p+1).*exp(1i*(m2-m1)*phiTab);
                                                    s1eqs2=logical(eye(NSi));
                                                    Wpn1n2(s1eqs2(:))=0;
                                                    W(n1S1Arr,n2S2Arr) = W(n1S1Arr,n2S2Arr)+Wpn1n2;
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                        end
                        
                        fprintf(1,repmat('\b',[1,length(msg)]));
                        msg = sprintf('done');
                        fprintf(1,msg);
                        
                        simul.numerics.couplingMatrices{jp}=W;
                    end
                    fprintf(' done\n')
                case 'none'
            end
        end
        
        % ======================================================================
        %> @brief Get method for nmax
        % ======================================================================
        function value=get.nmax(obj)
            value=jmult_max(1,obj.lmax);
        end
        
        % ======================================================================
        %> @brief Returns a gpuArray if gpuFlag, and a usual array
        %> otherwise
        %>
        %> @param Array x
        %> @return gpuArray(x) if gpuFlag, otherwise x
        % ======================================================================
        function arr_out = deviceArray(obj,arr_in)
            if obj.gpuFlag
                arr_out=gpuArray(arr_in);
            else
                arr_out=arr_in;
            end
        end
        
    end
end

