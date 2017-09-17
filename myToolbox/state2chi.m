function chi = state2chi(Rot,v,x,varargin)
if(nargin == 4)
    PosAmers = varargin{1};
    NbAmers = size(PosAmers,2);
    chi = eye(5+NbAmers);
    chi(1:3,1:5) = [Rot v x];
    if NbAmers > 0
        chi(1:3,6:end) = PosAmers;
    end
else
    chi = [Rot v x; zeros(2,3) eye(2)];
end
end
