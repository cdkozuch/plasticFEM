function ind = globElemInd(elements,elemID,ndofn)
%%GLOBELEMIND obtains global indices of given element's dofs
%Input:
%   elements (nel x nnel): element connectivity
%   elemID (scalar): element number of element in question
%   ndofn (scalar): dofs per node
%Output:
%   ind (1 x nnel*ndofn): global indices

nnel = size(elements,2); %nodes per element
ind = zeros(1,nnel*ndofn); %preallocate

for j=1:nnel
    r = (elements(elemID,j)-1)*ndofn + 1; %global index of first dof
    s = (j-1)*ndofn + 1;
    ind(s:s+ndofn-1) = r:r+ndofn-1;
end

end