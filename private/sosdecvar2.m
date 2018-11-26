function s = sosdecvar2(cstr,z)
% Returns SOS decision variable with minimal affine subspace.

z2 = kron(z,z);
if length(monomials(z2))/length(z2) > 0.125
    % use SOS decision variable as long as the number of additional 
    % components is reasonable small
    s = sosdecvar(cstr,z);
else
    % switch to minimal affine subspace decision variable
    % for large problems
    s = polydecvar(cstr,monomials(z2));
end

end