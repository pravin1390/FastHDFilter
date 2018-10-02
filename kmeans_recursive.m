function [Centre, minCentre]=kmeans_recursive(X,K0)

%   kmeans_recursive    Bisecting k-means clustring
%   [Centre, minCentre] = kmeans_recursive(X, Cluster) partititions the N x
%   P data matrix X into K0 clusters through a fully vectorized algorithm, where N is the number
%   of data points and P is the number of dimensions (variables). First whole data is bisected into clusters.
%   Chooses the cluster with point-centre distance is maximum and bisect it again to two clusters
%   minCentre is the returned N x 1 vector contains the cluster indices of each point.
%   Centre is the K0 cluster centroid locations in the K0 x P matrix.
%
K=1;                                                                % Initializing cluster count
[minCentre,Centretemp,newdistvect]=kmeans_cluster(X);       % Performing K means clustering

d=size(Centretemp,2);
Centre=zeros(K0,d);
var=zeros(K0,1);

K=K+1;

Centre(1,:)=Centretemp(1,:); 
Centre(K,:)=Centretemp(2,:);   % Incrementing cluster count  

var(1)=newdistvect(1); 
var(K)=newdistvect(2);
  
[~,maxindex]=max(var);                                    % Calculating maximum point-centre distance and cluster to be bisected
while (K<K0)   
    [label,Centretemp,newdistvect,newclustvect]=kmeans_cluster(X(minCentre==maxindex,:));      % Performing K means clustering   
   if (newclustvect(1) == 0)
   var(maxindex)=0; 
   elseif (newclustvect(2) == 0)
   var(maxindex)=0; 
   else
    K=K+1;                                                       % Incrementing cluster count
        
    % Populating new cluster centres
    Centre(maxindex,:)=Centretemp(1,:); 
    Centre(K,:)=Centretemp(2,:);   
  
    % Populating new cluster indices vector
    label=(label==1).*maxindex+(label==2).*K;
    minCentre(minCentre==maxindex)=label;   
        
    var(maxindex)=newdistvect(1); 
    var(K)=newdistvect(2);
    
   end
   if(any(var)==0) 
       break; 
   end
    [~,maxindex]=max(var);                                % Calculating maximum point-centre distance and cluster to be bisected    end
end
end
