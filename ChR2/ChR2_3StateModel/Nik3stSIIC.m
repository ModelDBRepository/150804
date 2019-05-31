function Out3stdot = Nik3stSIIC(t,Out3st,P)

% this function is integrating the 3-state model 

global Gd Gr

O = Out3st(1); D = Out3st(2);

% the 3-state model
Odot = P*(1-O-D) - Gd*O;
Ddot = Gd*O - Gr*D;


Out3stdot = [Odot Ddot];