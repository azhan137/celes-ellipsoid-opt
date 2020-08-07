function Mt = prepareMt(simul)

lmax=simul.numerics.lmax;
nmax=simul.numerics.nmax;
k = simul.input.k_medium;

spherArr=1:simul.input.particles.number;
spherArr=spherArr(:);
NSi = length(spherArr);

Mt = eye(NSi*simul.numerics.nmax,'single');

[x2,x1]=meshgrid(simul.input.particles.positionArray(spherArr,1));
[y2,y1]=meshgrid(simul.input.particles.positionArray(spherArr,2));
[z2,z1]=meshgrid(simul.input.particles.positionArray(spherArr,3));
x1mnx2 = -x1+x2;
y1mny2 = -y1+y2;
z1mnz2 = -z1+z2;
dTab = sqrt(x1mnx2.^2+y1mny2.^2+z1mnz2.^2);
ctTab = z1mnz2./dTab;
stTab = sqrt(1-ctTab.^2);
phiTab = atan2(y1mny2,x1mnx2);
Plm = legendre_normalized_trigon(ctTab,stTab,2*lmax);
particleArrayInd = simul.input.particles.singleUniqueArrayIndex;
for p=0:2*lmax
    sphHank = sph_bessel(3,p,k*dTab);
    for tau1=1:2
        for l1=1:lmax
            for m1=-l1:l1
                n1=multi2single_index(1,tau1,l1,m1,lmax);
                n1S1Arr=(1:NSi)+(n1-1)*NSi;
                for tau2=1:2
                    for l2=1:lmax
                        for m2=-l2:l2
                            if abs(m1-m2)<=p
                                n2=multi2single_index(1,tau2,l2,m2,lmax);
                                n2S2Arr=(1:NSi)+(n2-1)*NSi;
                                TWn1n2 = simul.tables.mieCoefficients(particleArrayInd(spherArr(:)),n2)*simul.tables.translationTable.ab5(n1,n2,p+1) .* Plm{p+1,abs(m1-m2)+1} .* sphHank .* exp(1i*(m1-m2)*phiTab) ;
                                s1eqs2=logical(eye(NSi));
                                TWn1n2(s1eqs2(:))=0; % jS1=jS2
                                Mt(n1S1Arr,n2S2Arr) = Mt(n1S1Arr,n2S2Arr) - TWn1n2;
                            end
                        end
                    end
                end
            end
        end
    end
end
end