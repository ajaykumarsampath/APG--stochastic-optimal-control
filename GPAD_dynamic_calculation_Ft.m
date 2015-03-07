function [ Z,Q] = GPAD_dynamic_calculation_Ft( sys,Ptree,Tree,Y,xinit)
%This function calculate the solution of the dynamic programming step on 
%the tree. All the off-line terms are calculated and passed as input to the
%function. The inital x, and dual varaibles (y) and terminal F_N are
%passed as inputs to this function
% Z is the output and containing 

Z.X=zeros(sys.nx,Tree.leaves(end));
Z.U=zeros(sys.nu,Tree.ancestor(Tree.leaves(end)));
S=zeros(sys.nu,Tree.ancestor(Tree.leaves(end)));

q=zeros(sys.nx,Tree.ancestor(Tree.leaves(end)));
qt=cell(1,length(Tree.leaves));
for i=1:length(Tree.leaves)
    %q(:,Tree.leaves(i))=sys.Ft{i,1}'*Y.yt{i,:}';
    qt{1,i}=Y.yt{i,:}';
end
for i=sys.Np:-1:1
    nodes_stage=find(Tree.stage==i-1);
    total_nodes=length(nodes_stage);
    for j=1:total_nodes
        no_child=length(Tree.children{nodes_stage(j)});
        if(no_child>1)
            sum_q=zeros(sys.nx,1);
            for k=1:no_child
                sum_q=sum_q+q(:,Tree.children{nodes_stage(j)}(k)-1);
            end
            S(:,nodes_stage(j))=Ptree.Phi{nodes_stage(j)}*Y.y(nodes_stage(j),:)'...
                +Ptree.Theta{nodes_stage(j)}*sum_q+Ptree.sigma{nodes_stage(j)};
        else
            if(i==sys.Np)
                %sum_q=sys.Ft{j,1}'*qt{1,j};
                sum_q=qt{1,j};
                S(:,nodes_stage(j))=Ptree.Phi{nodes_stage(j)}*Y.y(nodes_stage(j),:)'...
                   +Ptree.Theta{nodes_stage(j)}*qt{1,j}+Ptree.sigma{nodes_stage(j)};  
            else
                sum_q=q(:,Tree.children{nodes_stage(j)});
                S(:,nodes_stage(j))=Ptree.Phi{nodes_stage(j)}*Y.y(nodes_stage(j),:)'...
                    +Ptree.Theta{nodes_stage(j)}*q(:,Tree.children{nodes_stage(j)})+Ptree.sigma{nodes_stage(j)};
            end
        end 
        q(:,nodes_stage(j))=Ptree.c{nodes_stage(j)}+Ptree.d{nodes_stage(j)}'*Y.y(nodes_stage(j),:)'+...
            Ptree.f{nodes_stage(j)}'*sum_q+Ptree.h{nodes_stage(j)}'*S(:,nodes_stage(j));
    end
end
Z.X(:,1)=xinit;
for i=1:Tree.ancestor(Tree.leaves(end))
    Z.U(:,i)=Ptree.K{i}*Z.X(:,i)+S(:,i);
    for j=1:length(Tree.children{i})
        Z.X(:,Tree.children{i}(j))=sys.A*Z.X(:,i)+sys.B*Z.U(:,i)+Tree.value(Tree.children{i}(j),:)';
    end
end
Q.q=q;
Q.qt=qt;
end

