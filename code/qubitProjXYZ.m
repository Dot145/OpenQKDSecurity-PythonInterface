function [op] = qubitProjXYZ(basis, index, sender, bra)
% Gives the projection operator for the X (1), Y (2), or Z (3) basis.
% Note that for the Y basis, the sender and receiver have opposite states.
% For example, Alice having |R><R| means Bob gets |L><L|.

% Select the appropriate transformation matrix
U = eye(2);
switch basis
    case 1
        % X basis, apply Hadamard
        U = [1,1;1,-1]/sqrt(2);
    case 2
        % Y basis, check if this is sender or receiver
        if sender
            U = [1,1;-1i,1i]/sqrt(2); %U|1> = |L>, U|2> = |R>
        else
            U = [1,1;1i,-1i]/sqrt(2); %U|1> = |R>, U|2> = |L>
        end
    case 3
        % Z basis, no transformation needed
    otherwise
        error('Basis must be 1, 2, or 3 (corresponding to X, Y, or Z). Entered: %d',basis);
end
op = U*zket(2, index);
op = op';
if(~exist('bra', 'var'))
    op = op'*op;
end
end