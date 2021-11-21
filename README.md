# SLNMF
Overview:

This is code to do structural learning non-negative matrix factorization clustering of single-cell RNA-seq data given in the "experiment" section of the paper: 

Wenming Wu, Xiaoke Ma*. "A Network-based Structural Learning Nonnegative Matrix Factorization Algorithm for Clustering of scRNA-seq Data"  

The coding here is a generalization of the algorithm given in the paper. SLNMF is written in the MATLAB programming language. **To use, please download the SLNMF folder and follow the tutorial provided in the “README.doc”. Running the "main_slnmf.m" file to realize SLNMF experiment.** 

Files description:

slnmf.m - The main function.

main_slnmf.m - A script with a real scRNA-seq data to show how to run the code.

Biase_data.mat - A real scRNA-seq data used in the cell type clustering example.  We retain the genes that are expressed in at least 10% cells for the dataset. The Biase dataset contains 49 mouse fetal brain cells sequenced using SMAR T-Seq platform, which consists of three cell types, zygote cells, Two-cell cells and Four-cell cells.  

constructW.m - Compute adjacent matrix W.

PMI.m - Cell-cell network construction.

bestMap.m - permute labels of L2 to match L1 as good as possible.

ARI.m - Program for calculating the Adjusted Rand Index ( Hubert & Arabie) between two clusterings.

hungarian.m - Solve the Assignment problem using the Hungarian method.

data$Biase.expr.csv - Biase dataset original expression data. 

data$Biase.celltype.csv - The cell type of Biase dataset. 

Example:

**Follow the steps below to run SLNMF（also contained in the " main_slnmf.m" file）. Here use a real scRNA-seq data (Biase_data) set as an example.
**
clear all;

clc;

load('Biase_data.mat')%Loading the normalized scRNA-seq data

X=Data;%%%%%%%Rows are genes, columns are cell sample

    %==============Constructing a weight matrix==============
    
    %Preset value before constructing adjacency matrix
    
options = [];

option.Metric = 'Cosine';

options.NeighborMode = 'KNN';%KNN

options.k =5;%5 nearest neighbors

options.WeightMode = 'Cosine';%Weights are 0 or 1, it can eplace with 'HeatKernel', 'Euclidean' 

W = constructW(X',options);


M=PMI(W,1);%%% High-order cell-cell PMI network construction

k=3;lambda=10;iter=150;

[B,F,nn,error]=slnmf(k,lambda,iter,M);  %% Call the main function to solve the variables

%%%%%%%%%%% Clustering cell type label

for e=1:size(B,1)

    v=B(e,:);
    
    ma=max(v);
    
    [s,t]=find(v==ma);
    
    l(e)=t(1);
    
end

[newl] = bestMap(real_label,l);%%%  the label originally identified by the authors

ari = ARI(real_label,max(real_label),newl,max(newl))%% Calculating the Adjusted Rand Index (ARI)

pre_label =newl;

if ~isempty(real_label) 

exact = find(pre_label == real_label);

accuracy = length(exact)/length(newl) %% Calculating the accuracy

else

accuracy = []

end


Contact:

Please send any questions or found bugs to Xiaoke Ma xkma@xidian.edu.cn 
