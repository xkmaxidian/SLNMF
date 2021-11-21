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

