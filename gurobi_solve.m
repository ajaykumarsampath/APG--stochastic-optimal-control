function [ model] = gurobi_solve( sys,V,Tree)

% This functions solves the scenario tree based optimization problem
% Inputs:       sys         Input system dynamics
%              tree         contains the tree structure
%               ops         options of the gurobi solver.

%sys=sys_no_precond;
non_leaf_node=size(Tree.children,1);
K=size(Tree.leaves,2);
nvar=sys.nx+sys.nu;
nc=size(sys.F,1);
nc_t=size(sys.Ft{1},1);
model.Q=sparse(nvar*non_leaf_node+K*sys.nx,nvar*non_leaf_node+K*sys.nx);
model.A=sparse((nc+sys.nx)*non_leaf_node+(nc_t+sys.nx)*K,nvar*non_leaf_node+K*sys.nx);
model.rhs=zeros((nc+sys.nx)*non_leaf_node+(nc_t+sys.nx)*K,1);

%model.varnames = cell(1);

Qsps=sparse(V.Q);
Rsps=sparse(V.R);
PolyConstsps=sparse([sys.F sys.G]);

DynamicSps=sparse((non_leaf_node+K)*sys.nx,nvar*non_leaf_node+K*sys.nx);
krow=nc*non_leaf_node+nc_t*K;
%DynamicSps(1:sys.nx,1:sys.nx)=speye(sys.nx);

node_pos=0;
%JJ=zeros((child+K)*sys.nx,nvar*child+K*sys.nx);
for i=1:non_leaf_node
    %{
    for j=1:sys.nx
        model.varnames{(i-1)*nvar+j} = ['x_' num2str(i) ',' num2str(j)];
    end
    for j=1:sys.nu
        model.varnames{(i-1)*nvar+sys.nx+j}=['u_' num2str(i) ',' num2str(j)];
    end
    %}
    model.Q((i-1)*nvar+1:(i-1)*nvar+sys.nx,(i-1)*nvar+1:(i-1)*nvar+sys.nx)=Tree.prob(i)*Qsps;
    model.Q((i-1)*nvar+sys.nx+1:i*nvar,(i-1)*nvar+sys.nx+1:i*nvar)=Tree.prob(i)*Rsps;
    nchild=length(Tree.children{i});
    model.A((i-1)*nc+1:i*nc,(i-1)*nvar+1:i*nvar)=PolyConstsps;
    %if(nchild>1)
    for j=1:nchild
        DynamicSps((node_pos+j-1)*sys.nx+1:(node_pos+j)*sys.nx,(i-1)*nvar+1:i*nvar)...
            =sparse([sys.A{node_pos+j} sys.B{node_pos+j}]);
        model.rhs(krow+(node_pos+j-1)*sys.nx+1:krow+(node_pos+j)*sys.nx,1)=-Tree.value(node_pos+j,:)';
        if(i<=non_leaf_node-K)
            DynamicSps((node_pos+j-1)*sys.nx+1:(node_pos+j)*sys.nx,(node_pos+j)*nvar+1:...
                (node_pos+j)*nvar+sys.nx)=-speye(sys.nx);
        else
            DynamicSps((node_pos+j-1)*sys.nx+1:(node_pos+j)*sys.nx,non_leaf_node*nvar+(node_pos-non_leaf_node+j)*sys.nx+1:...
                non_leaf_node*nvar+(node_pos-non_leaf_node+j+1)*sys.nx)=-speye(sys.nx);
        end

    end
    
    node_pos=node_pos+nchild;
end
model.rhs(1:nc*non_leaf_node)=kron(ones(non_leaf_node,1),sys.g);
DynamicSps((non_leaf_node+K-1)*sys.nx+1:(non_leaf_node+K)*sys.nx,1:sys.nx)=speye(sys.nx);
%JJ=full(model.A);
for i=1:K
    %{
    for j=1:sys.nx
        model.varnames{(non_leaf_node+i-1)*nvar+j} = ['x_' num2str(non_leaf_node+i) ',' num2str(j)];
    end
    %}
    model.Q(nvar*non_leaf_node+(i-1)*sys.nx+1:nvar*non_leaf_node+i*sys.nx,...
        nvar*non_leaf_node+(i-1)*sys.nx+1:nvar*non_leaf_node+i*sys.nx)=Tree.prob(Tree.leaves(i))*sparse(V.Vf{i});
    model.A(non_leaf_node*nc+(i-1)*nc_t+1:non_leaf_node*nc+i*nc_t,nvar*non_leaf_node+(i-1)*sys.nx+1:nvar*non_leaf_node+i*sys.nx)=...
        sparse(sys.Ft{i});
    model.rhs(nc*non_leaf_node+(i-1)*nc_t+1:nc*non_leaf_node+i*nc_t)=sys.gt{i};
end
model.A(nc*non_leaf_node+nc_t*K+1:(nc+sys.nx)*non_leaf_node+(nc_t+sys.nx)*K,1:nvar*non_leaf_node+sys.nx*K)=DynamicSps;
model.sense=[repmat('<',nc*non_leaf_node+nc_t*K,1);repmat('=',(non_leaf_node+K)*sys.nx,1)];
model.obj=zeros(non_leaf_node*nvar+K*sys.nx,1);
model.lb=-inf*ones(non_leaf_node*nvar+K*sys.nx,1);
end

