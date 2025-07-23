function [DTnew, theta] = filterAcuteTriangles(DT, obtuseAngleDegrees)
% FILTERACUTETRIANGLES Filters triangles in a Delaunay triangulation that
% contain edges forming acute angles with neighboring edges.
%
% Input:
%   DT - delaunayTriangulation object
%
% Output:
%   DTnew - New delaunayTriangulation object with filtered connectivity
%
% Written with the assistance of ChatGPT on 12/17/2024
%   Joe Prisaznuk added theta, making the thresholding logic simpler.
%   We just have to check if any angle is overly obtuse, if so, it is no
%   longer a valid triangle. Having theta as an output also makes it easier
%   to diagnose how many triangles would be eliminated.

    % Extract points and current triangulation connectivity
    roiPts = DT.Points;
    connectivityList = DT.ConnectivityList;

    % Initialize array to keep valid triangles
    validTriangles = true(size(connectivityList, 1), 1);
    theta = zeros(size(connectivityList, 1),3);

    % Loop through all triangles
    for ii = 1:size(connectivityList, 1)
        % Get the points for the current triangle
        triPts = roiPts(connectivityList(ii, :), :);

        % Compute edge vectors for the triangle
        v1 = triPts(2, :) - triPts(1, :); % Edge 1-2
        v2 = triPts(3, :) - triPts(2, :); % Edge 2-3
        v3 = triPts(1, :) - triPts(3, :); % Edge 3-1

        % Normalize the edge vectors
        v1n = v1 / norm(v1);
        v2n = v2 / norm(v2);
        v3n = v3 / norm(v3);

        % Compute cosines of the angles using the dot product
        cosTheta(1) = dot(-v3n, v1n); % Angle at vertex 1
        cosTheta(2) = dot(-v1n, v2n); % Angle at vertex 2
        cosTheta(3) = dot(-v2n, v3n); % Angle at vertex 3
        theta(ii,:) = acos(cosTheta);

        thresholdAngle = obtuseAngleDegrees * pi/180; % angle in degrees

        if any(theta(ii,:) > thresholdAngle)
            validTriangles(ii) = false;
        end
    end

    % Rebuild the new connectivity list with only valid triangles
    newConnectivityList = connectivityList(validTriangles, :);

    % Create a new delaunayTriangulation object with the updated triangles
    DTnew = triangulation(newConnectivityList, roiPts);
end
