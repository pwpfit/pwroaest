function s = sosmdecvar(cstr, z, k)
% Returns k-by-1 vector of SOS decision variables.

s = polynomial(zeros(k,1));
for i=1:k
    s(i) = sosdecvar(sprintf('%s__%d',cstr,i), z);
end

end