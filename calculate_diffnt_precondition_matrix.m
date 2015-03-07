function [ sys_pre,Hessian ] = calculate_diffnt_precondition_matrix( sys,V,Tree,ops)
% This function calculate the precondion matirix.
% The hessian of the system is H=kron(prob(1:non_leaves,1),[F G]diag(Q,R)[F G]'
% ;Ft_{i}*prob(leaves)*Ft{i})
% The precondtion matrix is H^(-1/2).

sys_pre=sys;
Nd=length(Tree.stage);
Ns=length(Tree.leaves);
H1=sys.F*(V.Q\sys.F')+sys.G*(V.R\sys.G');
if(ops.use_cell==1)
    sys_pre.F=cell(Nd-Ns,1);
    sys_pre.G=cell(Nd-Ns,1);
    sys_pre.g=cell(Nd-Ns,1);
    if(ops.use_hessian==0)
        for i=1:Nd-Ns
            Hessian.H=sqrt(diag(diag(H1)));
            %sys_pre.F{i}=sqrt(Tree.prob(i))*(Hessian.H\sys.F);
            %sys_pre.G{i}=sqrt(Tree.prob(i))*(Hessian.H\sys.G);
            %sys_pre.g{i}=sqrt(Tree.prob(i))*(Hessian.H\sys.g);
            sys_pre.F{i}=sqrt(Tree.prob(i))*(sys.F);
            sys_pre.G{i}=sqrt(Tree.prob(i))*(sys.G);
            sys_pre.g{i}=sqrt(Tree.prob(i))*(sys.g);
        end
    else
        for i=1:Nd-Ns
            Hessian.H=sqrt(diag(diag(H1)));
            sys_pre.F{i}=Hessian.H\sys.F;
            sys_pre.G{i}=Hessian.H\sys.G;
            sys_pre.g{i}=Hessian.H\sys.g;
        end
    end
    
else
    Hessian.H=sqrt(diag(diag(H1)));
    sys_pre.F=Hessian.H\sys.F;
    sys_pre.G=Hessian.H\sys.G;
    sys_pre.g=Hessian.H\sys.g;
end
%[VB,VD]=eig(H1);
%Hessian.H=((VB'\sqrt(VD))*VB')';
%norm(H1-Hessian.H*Hessian.H)
for i=1:Ns
    H1=sys.Ft{i}*(V.Vf{i}\sys.Ft{i}');
    Hessian.Ht{i}=sqrt(diag(diag(H1)));
    %[VB,VD]=eig(H1);
    %Hessian.Ht{i}=((VB'\sqrt(VD))*VB')';
    %norm(H1-Hessian.Ht{i}*Hessian.Ht{i})
end

if(ops.use_hessian==0)
    for i=1:Ns
        %sys_pre.Ft{i}=sqrt(Tree.prob(Tree.leaves(i)))*(Hessian.Ht{i}\sys_pre.Ft{i});
        %sys_pre.gt{i}=sqrt(Tree.prob(Tree.leaves(i)))*(Hessian.Ht{i}\sys_pre.gt{i});
        sys_pre.Ft{i}=sqrt(Tree.prob(Tree.leaves(i)))*(sys_pre.Ft{i});
        sys_pre.gt{i}=sqrt(Tree.prob(Tree.leaves(i)))*(sys_pre.gt{i});
    end
else
    for i=1:Ns
        sys_pre.Ft{i}=Hessian.Ht{i}\sys_pre.Ft{i};
        sys_pre.gt{i}=Hessian.Ht{i}\sys_pre.gt{i};
    end
end
end


