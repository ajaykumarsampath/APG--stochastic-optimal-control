function [ sys ] = Normalise_constraints(sys)
% This function  normalize the constraints of the matrix 

%sys=sys_no_precond;
K=size(sys.Ft,1);

%Normalising the constraints in sys.F

nc=size(sys.F,1);
for i=1:nc
    if(abs(sys.g(i))>0)
        sys.F(i,:)=sys.F(i,:)/sys.g(i,:);
        sys.G(i,:)=sys.G(i,:)/sys.g(i,:);
        sys.g(i,1)=1;
    end
end

for i=1:K
    nc=size(sys.Ft{i},1);
    for j=1:nc
        if(abs(sys.gt{i}(j))>0)
            sys.Ft{i}(j,:)=sys.Ft{i}(j,:)/sys.gt{i}(j);
            sys.gt{i}(j,1)=1;
        end
    end
end
end

