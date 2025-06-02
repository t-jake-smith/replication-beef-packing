function d = distances(row_lat, row_lng, col_lat, col_lng)
% Function - compute distance between given two matrixes of [lat,lng]

% Radius of the earth in miles; needed for distance calc
R = 3959;

% Define constants for the haversine calculation
phi1 = col_lat*pi()/180; % Latitude of the column vector in radians
phi2 = row_lat*pi()/180; % Latitude of the row vector in radians
deltaPhi = phi1-phi2; % Difference between column and row latitudes in radians
deltaLambda = (col_lng-row_lng)*pi()/180; % Difference between column and row longitudes in radians

a = sin(deltaPhi/2).*sin(deltaPhi/2) + cos(phi1).*cos(phi2).*sin(deltaLambda/2).*sin(deltaLambda/2);
c = 2*atan2(sqrt(a),sqrt(1-a));

d = R*c; % row x col matrix of distances, where d(r,c) = distance between r and c
end