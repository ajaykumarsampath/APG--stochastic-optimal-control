function [Ptree] = GPAD_dynamic_formulation_precondition( sys,V,Tree)
% This function calculate the off-line elements for the dynamic programming
% step of the GPAD algorithm


Ptree=struct('P',cell(1,1),'c',cell(1,1),'d',cell(1,1),'f',cell(1,1),'h',cell(1,1),...
    'Phi',cell(1,1),'Theta',cell(1,1),'sigma',cell(1,1));
for i=1:length(Tree.leaves)
    %pchld=Tree.prob(Tree.ancestor(Tree.leaves(i)))/Tree.prob(Tree.leaves(i));
    Ptree.P{Tree.leaves(i),1}=Tree.prob(Tree.leaves(i))*V.Vf{i};
    %Ptree.P{Tree.leaves(i),1}=V.Vf{i};
end

for i=sys.Np+1:-1:1
    if(i==sys.Np+1)
    else
        nodes_stage=find(Tree.stage==i-1);
        total_nodes=length(nodes_stage);
        for j=1:total_nodes
            Pchld=Tree.prob(Tree.children{nodes_stage(j)})/Tree.prob(nodes_stage(j));
            no_child=length(Pchld);
            Pbar=zeros(sys.nx);
            neta=zeros(sys.nx,1);
            for k=1:no_child
                Pbar=Pbar+Ptree.P{Tree.children{nodes_stage(j)}(k)};% \bar{P}
                neta=neta+Ptree.P{Tree.children{nodes_stage(j)}(k)}...
                    *Tree.value(Tree.children{nodes_stage(j)}(k),:)';%phi_{k-1}^{(i)}
            end
            
            %terms in the control u_{k-1}^{\star (i)}
            Rbar=2*(Tree.prob(nodes_stage(j))*V.R+sys.B'*Pbar*sys.B);
            %Rbar=2*(V.R+sys.B'*Pbar*sys.B);
            Rbar_inv=Rbar\eye(sys.nu);
            Ptree.K{nodes_stage(j)}=-2*Rbar_inv*(sys.B'*Pbar*sys.A);%K_{k-1}^{(i)}
            Ptree.Phi{nodes_stage(j)}=-Rbar_inv*sys.G{nodes_stage(j)}';%\Phi_{k-1}^{(i)}
            if(i==sys.Np)
                Ptree.Theta{nodes_stage(j)}=-Rbar_inv*sys.B'*sys.Ft{j,1}';%\Theta_{k-1}^{(i)}
            else
                Ptree.Theta{nodes_stage(j)}=-Rbar_inv*sys.B';%\Theta_{k-1}^{(i)}
            end
            Ptree.sigma{nodes_stage(j)}=-2*Rbar_inv*(sys.B'*neta);%sigma_{k-1}^{(i)}
            
            %terms in the linear cost
            Ptree.c{nodes_stage(j)}=2*(sys.A+sys.B*Ptree.K{nodes_stage(j)})'*neta;%c_{k-1}^{(i)}
            Ptree.d{nodes_stage(j)}=sys.F{nodes_stage(j)}+sys.G{nodes_stage(j)}*Ptree.K{nodes_stage(j)};%d_{k-1}^{(i)}
            if(i==sys.Np)
                Ptree.f{nodes_stage(j)}=sys.Ft{j,1}*(sys.A+sys.B*Ptree.K{nodes_stage(j)});%f_{k-1}^{(i)}
            else
                Ptree.f{nodes_stage(j)}=(sys.A+sys.B*Ptree.K{nodes_stage(j)});%f_{k-1}^{(i)}
            end
            
            Ptree.h{nodes_stage(j)}=(2*sys.A'*Pbar*sys.B+Ptree.K{nodes_stage(j)}'*Rbar')';%h_{k-1}^{(i)}
            
            %Quadratic cost
            if(i==sys.Np)
                Ptree.P{nodes_stage(j)}=Tree.prob(nodes_stage(j))*(V.Q+Ptree.K{nodes_stage(j)}'*V.R*Ptree.K{nodes_stage(j)})...
                    +(sys.A+sys.B*Ptree.K{nodes_stage(j)})'*Pbar*(sys.A+sys.B*Ptree.K{nodes_stage(j)});
            else
                Ptree.P{nodes_stage(j)}=Tree.prob(nodes_stage(j))*(V.Q+Ptree.K{nodes_stage(j)}'*V.R*Ptree.K{nodes_stage(j)})...
                    +Ptree.f{nodes_stage(j)}'*Pbar*Ptree.f{nodes_stage(j)};
            end
            
            %Ptree.P{nodes_stage(j)}=(V.Q+Ptree.K{nodes_stage(j)}'*V.R*Ptree.K{nodes_stage(j)})...
            %+Ptree.f{nodes_stage(j)}'*Pbar*Ptree.f{nodes_stage(j)};
        end
    end
    
end

end


