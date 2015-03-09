% Where to store the header file:
%header_filename = '/media/ajay/New\ Volume/work/gpadsys_data.h';
header_filename = 'gpadsys_simpleheader.h';
int_type = 'uint_t';

%create header
tic;
value=Tree.value';
size_of_FN = zeros(length(sys.Ft),1);
size_of_FN_cum = zeros(size(size_of_FN));
for kk=1:length(sys.Ft),
    size_of_FN(kk) = length(sys.gt{kk});
    if (kk>1), size_of_FN_cum(kk) = size_of_FN_cum(kk-1)+size_of_FN(kk-1); end
end
FN_numel = sum(size_of_FN)*sys.nx;
GN_numel = sum(size_of_FN);

% Total number of children
N_children = 0;
Tree_Num_Children = zeros(length(Tree.stage),1);
Tree_Num_Children_cum =  zeros(length(Tree.stage),1);
for i=1:length(Tree.children),
    N_children = N_children + length(Tree.children{i});
    Tree_Num_Children(i) = length(Tree.children{i});
    Tree_Num_Children_cum(i) = N_children;
end
Tree_Num_Children_cum(i+1:end)=N_children;

Tree_nodes_per_stage = zeros(sys.Np,1);
Tree_nodes_per_stage_cum = zeros(sys.Np,1);
for i=0:sys.Np,
    Tree_nodes_per_stage(i+1) = sum(Tree.stage==i);
    if (i<sys.Np+1), Tree_nodes_per_stage_cum(i+2) = Tree_nodes_per_stage_cum(i+1) + Tree_nodes_per_stage(i+1); end
end

n_child = zeros(length(Tree.children),1);
for i=1:length(Tree.children),
    n_child(i) = length(Tree.children{i});
end

nc=size(sys.F{1},1);
if(test_cuda)
    wt=zeros(1,nc*Nd+GN_numel);
    for i=1:Nd
        wt(1,(i-1)*nc+1:i*nc)=W.y(i,:);
    end
    if(Ns==1)
        wt(nc*Nd+1:nc*Nd+size_of_FN(Ns))=W.yt{Ns};
    else
        for i=1:Ns-1
            wt(1,nc*Nd+size_of_FN_cum(i)+1:nc*Nd+size_of_FN_cum(i+1))=W.yt{i};
        end
        wt(1,nc*Nd+size_of_FN_cum(Ns)+1:nc*Nd+size_of_FN_cum(Ns)+size_of_FN(Ns))=W.yt{Ns};
    end
else
    wt=[];
    Zdgpad.X=[];
    Zdgpad.U=[];
end
if(details_solvers.method==0)
    Z_yalmip=cell(1,3);
end
    
f = fopen(header_filename,'w+');
fprintf(f, '/* Auto-generated file */\n');
fprintf(f, '/* Header file : sys_simpleheader.h */\n\n\n');
fprintf(f, '#ifndef __SIMPLE_GPAD_HEADER_\n');
fprintf(f, '#define __SIMPLE_GPAD_HEADER_\n\n\n');
fprintf(f,'typedef int uint_t; \n');
fprintf(f,'typedef float real_t;\n');


Dimensions = {
    {'NX              ', sys.nx, 'State dimension'}
    {'NU              ', sys.nu, 'Input dimension'}
    {'N               ',  sys.Np, 'Prediction horizon'}
    {'NC              ', size(sys.F{1},1), 'Number of mixed state-input constraints'}
    {'K               ',  length(Tree.leaves), 'Number of scenarios (leaf nodes)'}
    {'FN_NUMEL        ', FN_numel, 'Total number of elements to be stored in FN'}
    {'GN_NUMEL        ', GN_numel, 'Total number of elements to be stored in gN'}
    {'N_NODES         ',length(Tree.prob),'Number of nodes of the tree'}
    {'N_CHILDREN_TOT  ', N_children, 'Total number of children'}
    {'N_NONLEAF_NODES ', length(Tree.children), 'Number of non-leaf nodes in the tree'}
    {'DIM_GPAD_K_GAIN ', sys.nu*sys.nx*length(Tree.children),'Dimension of K (total number of elements)'}
    {'test_size',        Test_points,'Number of test cases'}
    {'test_cuda',        test_cuda,'Check is it in testing mode'}
    {'multi_uncertanity', ops_system.sys_uncert, 'The system contain multiplicative uncertanity'}
    };




Details = {
    {'TREE DATA'}
    {int_type,  'TREE_STAGES',          Tree.stage,             'N_NODES',          'The stage of each node in the tree (std. node enumeration)'}
    {int_type,  'TREE_NODES_PER_STAGE', Tree_nodes_per_stage,   'N+1',              'Number of nodes at each stage 0,...,N'}
    {int_type,  'TREE_NODES_PER_STAGE_CUMUL', Tree_nodes_per_stage_cum,   'N+2',      'Cumulative counterpart of TREE_NODES_PRE_STAGE'}
    {int_type,  'TREE_LEAVES',          Tree.leaves,            'K',                'Indices of the leaf nodes of the tree'}
    {int_type,  'TREE_CHILDREN',        Tree.children,          'N_CHILDREN_TOT',   'Children indices (look-up array)'}
    {int_type,  'TREE_NUM_CHILDREN',    n_child,                'N_NONLEAF_NODES',  'Number of children of each node'}
    {int_type,  'TREE_ANCESTOR',        Tree.ancestor,          'N_NODES',          'Ancestors of all nodes'}
    {int_type,  'TREE_N_CHILDREN_CUMUL',Tree_Num_Children_cum,  'N_NODES',          'Cumulative number of children'}
    {'real_t',  'TREE_PROB',            Tree.prob,              'N_NODES',          'Probability of every node'}
    {int_type,   'iterate',             3000,           '1'}
    {'SYSTEM DATA'}
    {int_type,  'FN_ROWS',              size_of_FN,             'K',                'Sizes of the terminal sets (# inequalities)'}
    {int_type,  'FN_ROWS_CUMUL',        size_of_FN_cum,         'K',                'FN_ROWS cumulative'}
    {'real_t',  'Q',                    V.Q,                       'NX*NX',            'Cost function Q'}
    {'real_t',  'R',                    V.R,                       'NX*NU',            'Cost function R'}
    {'SOLVER DATA'}
    {'real_t',  'hessian',              ops_GPAD.alpha,         '1'}
    {'Test data for the GPAD'}
    {'real_t',   'dgpad_x0_test',       x_rand,             'NX*test_size'}
    {'real_t',   'dgpad_ures_test'      result.u{1},           'NU*test_size'}
    {'Test DP step of GPAD'}
    {'real_t',  'dp_x_test',           Zdgpad.X,           'NX*N_NODES'}
    {'real_t',  'dp_u_test',           Zdgpad.U,           'NU*N_NONLEAF_NODES'}
    {'real_t',  'dp_y_test',                 wt,           'NC*N_NONLEAF_NODES+GN_NUMEL'}
    {'Test GPAD algorithm'}
    %{'real_t'    'gpad_u_test',      Z_gpad_pre.U,       'NU*N_NONLEAF_NODES'}
    {'real_t',    'gpad_u_test'      Z_yalmip{1,2},      'NU*N_NONLEAF_NODES'}
    };

if(test_cuda)
    len_details=length(Details);
else
    len_details=length(Details);
end

fprintf(f, '/* Dimensions */\n');
for kk=1:length(Dimensions),
    fprintf(f,'#define\t\t%s\t\t%d', Dimensions{kk}{1}, Dimensions{kk}{2});
    if length(Dimensions{kk})==3,
        fprintf(f, '\t\t/**< %s */', Dimensions{kk}{3});
    end
    fprintf(f, '\n');
end
fprintf(f,'\n\n');



fprintf(f, '/* Basic problem data */\n');


for kk=1:len_details,
    m_matrix = Details{kk};
    if (length(Details{kk})==1),
        fprintf(f, '\n\n/**** %s ****/\n\n', m_matrix{1});
    else
        if (length(Details{kk})==5),
            fprintf(f,'/** %s */\n', m_matrix{5});
        else
            fprintf(f,'/** Matrix %s */\n', m_matrix{2});
        end
        fprintf(f,'%s %s[%s] = {',m_matrix{1}, m_matrix{2}, m_matrix{4});
        
        m = m_matrix{3};
        if isnumeric(m), % m is a matrix
            for i=1:size(m,2),
                for j=1:size(m,1),
                    if strcmp(m_matrix{1},'int')==1,
                        fprintf(f,'%d', m(j,i));
                    else
                        fprintf(f,'%g', m(j,i));
                    end
                    if (~(i==size(m,1) && j==size(m,2))),
                        fprintf(f,', ');
                    end
                end
            end
            fprintf(f, '};\n\n');
        elseif iscell(m), % m is a cell
            for s=1:length(m),
                mm=m{s};
                for i=1:size(mm,2),
                    for j=1:size(mm,1),
                        if strcmp(m_matrix{1},'int')==1,
                            fprintf(f,'%d', mm(j,i));
                        else
                            fprintf(f,'%g', mm(j,i));
                        end
                        if (~(i==size(mm,1) && j==size(mm,2) && s==length(m))),
                            fprintf(f,', ');
                        end
                    end
                end
            end
            fprintf(f, '};\n\n');
        end
    end
end
fprintf(f,'\n\n');
fprintf(f, '#endif /* __SIMPLE_GPAD_HEADER_ */\n');
fprintf(f, '/* File generated automatically in %5.3fs */', toc);
fclose(f);


