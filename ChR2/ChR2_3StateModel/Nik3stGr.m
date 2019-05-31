function Out3stdot = Nik3stGr(O, D, P, Gd, Gr)

%Nikolic 3-state model
Odot = P*(1-O-D) - Gd*O;
Ddot = Gd*O - Gr*D;


Out3stdot = [Odot Ddot];