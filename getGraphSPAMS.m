function graph = getGraphSPAMS(sizeImg,sizeBatch)

m = sizeImg(1);
n = sizeImg(2);
a = min(sizeBatch(1),m);
b = min(sizeBatch(2),n);

numX = (m-a+1);
numY = (n-b+1);
numGroup = numX*numY;

graph.eta_g=ones(1,numGroup);
graph.groups=sparse(numGroup,numGroup);
graph.groups_var=sparse(m*n,numGroup);

for i = 1:numGroup 
    indiMatrix = zeros(m,n);
    indX = mod(i-1,numX)+1;
    indY = ceil(i/numX);
    indiMatrix(indX:indX+a-1,indY:indY+b-1)=1;
    graph.groups_var(:,i)=indiMatrix(:);
end