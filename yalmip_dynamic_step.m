function [ yalmip_dp,yalmip_const_time] = yalmip_dynamic_step(sys,V,Tree,Y,ops)
%This a that creates an yalmip optimizer object to solve the optimization
%on a scenario tree. 
%
% Syntax:  [ yalmip_tree,yalmip_const_time] = yalmip_dynamic_step(sys,V,Tree,ops)
% 
%  Input:             sys= Details of the system including the A,B,
%                          constraints and Terminal constraints.
%                       V=cost function.
%                    Tree=Descripition of the tree structure which includes
%                         Children, ancestor and value at each node.
%                     ops=Containt the options for the yalmip solver.
%                         Default is to use Gurobi.
%
% Output:     yalmip_tree= Output yalmip optimizer.
%       yalmip_const_time= Time take to construct the yamip objective.


default_options=sdpsettings('solver','gurobi','verbose',0,'cachesolvers',1);
sdp_options=default_options;
if nargin==5, % User-provided options
    sdp_options = ops;
    flds = fieldnames(default_options);
    for i=1:numel(flds),
        if ~isfield(sdp_options, flds(i))
            sdp_options.(flds{i})=default_options.(flds{i});
        end
    end
end
    
Ns=length(Tree.leaves);

Xyal=sdpvar(sys.nx,Tree.leaves(end));
Uyal=sdpvar(sys.nu,Tree.ancestor(Tree.leaves(end)));
Xinit=sdpvar(sys.nx,1);

J_yal=0;
const_yal=(Xyal(:,1)==Xinit);
tic;
for i=1:Tree.ancestor(Tree.leaves(end))
    pl=Tree.prob(i);
    J_yal=J_yal+pl*(Xyal(:,i)'*V.Q*Xyal(:,i)+Uyal(:,i)'*V.R*Uyal(:,i)); 
    J_yal=J_yal+Y.y(i,:)*(sys.F{i}*Xyal(:,i)+sys.G{i}*Uyal(:,i)-sys.g{i});
    for j=1:length(Tree.children{i})
        const_yal=const_yal+(Xyal(:,Tree.children{i}(j))==sys.A{Tree.children{i}(j)}*Xyal(:,i)+...
            sys.B{Tree.children{i}(j)}*Uyal(:,i)+Tree.value(Tree.children{i}(j),:)');
    end
end
for i=1:Ns
    pl=Tree.prob(Tree.leaves(i));
    J_yal=J_yal+pl*(Xyal(:,Tree.leaves(i))'*V.Vf{i}*Xyal(:,Tree.leaves(i)));
    J_yal=J_yal+Y.yt{i}*(sys.Ft{i}*Xyal(:,Tree.leaves(i))-sys.gt{i});
end
details.yalmip_const=toc;
% building a optimizer 
%sdp_options = sdpsettings('solver','gurobi','verbose',0,'cachesolvers',1);
tic
yalmip_dp=optimizer(const_yal,J_yal,sdp_options,{Xinit},{Xyal,Uyal,J_yal});
yalmip_const_time=toc;



end




