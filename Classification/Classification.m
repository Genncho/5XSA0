clear all;
close all;
clc;

data = load("2D_points_for_Kmeans.txt");
data_new = permute(reshape(data.', 2,1,[]), [1 2 3]);
ClusterCenters = [ 
    [1 2]; 
    [3 4]; 
    [5 6] 
 ];
[centroids_new, iteration] = K_Cluster_Algorithm(data, 3);

function[centroids_new, numIterations] = K_Cluster_Algorithm(dataPoints, ClustersAmount, ClusterCenters)
    % properly index data
    data = permute(reshape(dataPoints.', 2,1,[]), [1 2 3]);

    centroids = zeros(2,1,ClustersAmount);
    % create random cluster points
    if nargin == 2 % rng them, otherwise direct assign
        ClusterCenters = rand(ClustersAmount, 2) .* 10;
    end
    centroids = permute(reshape(ClusterCenters.', 2, 1, []), [1 2 3]);

    % recursion part
    [centroids_new, numIterations] = RecursionPart(data, centroids, 1);
end

function[centroids_new, numIterations] = RecursionPart(dataPoints, centroids, iteration)
    % save point indx in cluster arrays
    % actually I can just compare the old and new clusters to check if
    % there has been a change of a point parenthood, aka if it has converge
    % in case of no change - points are indexed anyways so yes
    % ORRRR if there is almost no change in the centroids then it has also
    % converged?

    for i = 1:size(centroids, 3)
        clusters{i} = [];
    end

    % find nearest centroid for each point and assing to cluster i+1
    for point_idx = 1:size(dataPoints, 3)
        minDist = realmax;  % max big number
        minClusterIdx = -1;
        for clusterIdx = 1:size(centroids, 3)
            %x = dataPoints(1,:,point_idx) - centroids(1,:,clusterIdx);
            %y = dataPoints(2,:,point_idx) - centroids(2,:,clusterIdx);
            dist = pointDist2D(dataPoints(:,:,point_idx), centroids(:,:,clusterIdx));  %sqrt( x^2 + y^2 );
            if dist < minDist
                minDist = dist;
                minClusterIdx = clusterIdx;
            end
        end
        % now I have the closest centroid - assign
        clusters{minClusterIdx}(end+1) = point_idx;
    end
    plotClusters(dataPoints, centroids, clusters, iteration);
    pause(1);

    % calculate new centroids - by means of clustered points
    for clusterIdx = 1:size(centroids, 3)
        pointIndices = clusters{clusterIdx};
        if ~isempty(pointIndices)
            clusterPoints = dataPoints(:,:,pointIndices);  
            centroids_new(:,:,clusterIdx) = mean(clusterPoints, 3);
        else
            warning("Empty Cluster Detected. Initializing new clusters randomly!");
            randIdx = randi(size(dataPoints, 3));
            centroids_new(:,:,clusterIdx) = dataPoints(:,:,randIdx);
        end
    end

    threshhold = 0.01;
    % compare with old centroids and check for convergence

    converged = true;
    for i = 1:size(centroids, 3)
        dist = pointDist2D(centroids(:,:,i), centroids_new(:,:,i));
        if dist > threshhold
            converged = false;
            break;
        end
    end
    
   if converged
        numIterations = iteration;
        return;
    else
        [centroids_new, numIterations] = RecursionPart(dataPoints, centroids_new, iteration + 1);
    end
end

function d = pointDist2D(p1, p2)
    % p1 and p2 are 2x1 vectors (x; y)
    d = sqrt(sum((p1 - p2).^2));
end

function plotClusters(dataPoints, centroids, clusters, iteration)
    colors = lines(numel(clusters));  % unique colors per cluster

    figure(1); clf; hold on;
    title(['K-Means Clustering - Iteration ', num2str(iteration)]);
    
    % Plot each cluster
    for i = 1:numel(clusters)
        if ~isempty(clusters{i})
            clusterData = dataPoints(:,:,clusters{i});
            x = squeeze(clusterData(1,1,:));
            y = squeeze(clusterData(2,1,:));
            scatter(x, y, 36, colors(i,:), 'filled');
        end
    end

    % Plot centroids
    for i = 1:size(centroids, 3)
        cx = centroids(1,1,i);
        cy = centroids(2,1,i);
        plot(cx, cy, 'kx', 'MarkerSize', 15, 'LineWidth', 2);  % black X
    end

    axis equal; grid on;
    drawnow;
end

%%
% Algorithm outline:
% Place random centorids
% Repeat following:
%   for each point -> find nearest centroid then assing to cluster j
%   for each cluster -> find new centroid by the mean of all points
%       in the cluster
%   repeat until no points change cluster -> convergence


