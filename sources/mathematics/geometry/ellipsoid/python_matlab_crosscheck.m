python_data_path = 'C:/users/azhan/documents/github/cylinder_t_matrix/cylinder_python_test.mat';

load(python_data_path);

radiusp = double(radiusp);
heightp = double(heightp);
Nthetap = double(Nthetap);
lmaxp = double(lmaxp);

[J11,J12,J21,J22,dJ11,dJ12,dJ21,dJ22] = compute_J_cyl_celes_posm(lmaxp,Nthetap,radiusp,heightp,nip,nsp,wavelengthp,3);
[rJ11,rJ12,rJ21,rJ22,drJ11,drJ12,drJ21,drJ22] = compute_J_cyl_celes_posm(lmaxp,Nthetap,radiusp,heightp,nip,nsp,wavelengthp,1);

[Q, dQ] = compute_Q_cyl(lmaxp,Nthetap,radiusp,heightp,nip,nsp,wavelengthp,3);
[rQ, drQ] = compute_Q_cyl(lmaxp,Nthetap,radiusp,heightp,nip,nsp,wavelengthp,1);

[T0, dT0] = compute_T_cyl(lmaxp,Nthetap,[radiusp,heightp],nip,nsp,wavelengthp);

DJ11 = abs(J11p-J11);
DJ12 = abs(J12p-J12);
DJ21 = abs(J21p-J21);
DJ22 = abs(J22p-J22);
DdJ11 = abs(dJ11p-dJ11);
DdJ12 = abs(dJ12p-dJ12);
DdJ21 = abs(dJ21p-dJ21);
DdJ22 = abs(dJ22p-dJ22);

DQ = abs(Q-Qp);
DrQ = abs(rQ-rQp);
DdQ = abs(dQ-dQp);

DT = abs(T0-T0p);
DdT = abs(dT0-dT0p);

[Tmp, dTmp] = compute_T_p(Qp, rQp, dQp, drQp);

DTmp = abs(T0-Tmp);
DdTmp = abs(dT0-dTmp);

invQ = inv(Q);
invQmp = inv(Qp);

DinvQ = abs(invQ-invQp);

T1 = rQp*invQp;
T2 = rQp*invQmp;