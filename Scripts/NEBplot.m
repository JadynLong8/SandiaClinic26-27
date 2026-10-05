% Matrix B contains specific points from A to be used as markers
% Format: [Distance (A), Energy (H)]
B1 = [
0.0000   0.00000000      0.00000000
0.2246   9.21009588      0.00255950
0.3128  12.82724833      0.01067922
0.3605  14.78131585      0.02217416
0.4054  16.62395085      0.02356506
0.4500  18.45299729      0.02147992
0.4991  20.46632497      0.02083793
0.5496  22.53973064     -0.00062689
1.0000  41.00774345     -0.05891633
];
B2 = [
0.0000   0.00000000      0.00000000
0.2251   8.61211513      0.00860505
0.2737  10.47373280      0.06745207
0.3223  12.33462953      0.06514408
0.3725  14.25263785      0.00754025
0.6133  23.46864129      0.00188289
1.0000  38.26719971     -0.01274898
];
B3 =[
0.0000   0.00000000      0.00000000
0.2650  10.93845440      0.00192600
0.4082  16.84473941      0.00498986
0.4852  20.02621120      0.01077655
0.5421  22.37391881      0.01518321
0.5989  24.71717631      0.00071386
1.0000  41.27033300     -0.01767690
];
B4 = [
0.0000   0.00000000      0.00000000
0.2283   9.79124623      0.04089276
0.3042  13.04912441      0.10083027
0.3756  16.11219606      0.10718705
0.4472  19.18355123      0.09468473
0.5293  22.70319300      0.05660755
0.6323  27.12068019      0.07810998
0.7355  31.54738456      0.03733856
1.0000  42.89457827      0.02092683
];

% Shift Matrix B (Column 2: Distance, Column 3: Energy)
B2(:,2) = B2(:,2) + B1(end,2);
B2(:,3) = B2(:,3) + B1(end,3);

B3(:,2) = B3(:,2) + B2(end,2);
B3(:,3) = B3(:,3) + B2(end,3);

B4(:,2) = B4(:,2) + B3(end,2);
B4(:,3) = B4(:,3) + B3(end,3);

% Combine 
B = [B1; B2; B3; B4];

% Convert Hartree to eV
B(:,3) = B(:,3) .* 27.211386245981; 

% --- RESCALE DISTANCE TO [0, 1] ---
% Find the maximum distance in the combined matrix
max_distance = max(B(:,2));
% Divide all distance values by the maximum distance
B(:,2) = B(:,2) / max_distance;

% --- INTERPOLATION STEP ---
% 1. Remove duplicate stitched points for interpolation
[x_unique, idx] = unique(B(:,2));
y_unique = B(idx, 3);

% 2. Create a dense set of X points for a smooth curve
x_interp = linspace(min(x_unique), max(x_unique), 1000);

% 3. Calculate the interpolated Y values globally using a cubic spline
y_interp = interp1(x_unique, y_unique, x_interp, 'pchip');

% --- PLOTTING WITH DIFFERENT COLORS ---
figure;
hold on;

% Define 3 base colors from MATLAB's colormap
base_colors = lines(3); 

% Assign the colors to the 4 segments (B2 and B3 share the 2nd color)
segment_colors = [base_colors(1,:); ... % B1
                  base_colors(2,:); ... % B2
                  base_colors(2,:); ... % B3 (Same as B2)
                  base_colors(3,:)];    % B4

% Calculate the exact row indices where each segment ends in the combined matrix 'B'
b_idx = [1, ...
         size(B1,1), ...
         size(B1,1) + size(B2,1), ...
         size(B1,1) + size(B2,1) + size(B3,1), ...
         size(B,1)];

% Array to store handles for the legend
h_lines = zeros(1, 4);

for i = 1:4
    % Define start and end indices for the markers of the current segment
    start_idx = b_idx(i);
    if i > 1
        % Start one index later to avoid plotting the overlapping stitch marker twice
        start_idx = b_idx(i) + 1; 
    end
    end_idx = b_idx(i+1);
    
    % Define the X-coordinate boundaries for masking the interpolated line
    x_start = B(b_idx(i), 2);
    x_end = B(b_idx(i+1), 2);
    
    % Create a logical mask to isolate the interpolated points for this segment
    mask = (x_interp >= x_start) & (x_interp <= x_end);
    
    % Plot the interpolated segment
    h_lines(i) = plot(x_interp(mask), y_interp(mask), 'Color', segment_colors(i,:), 'LineWidth', 1.5);
    
    % Plot the markers for this segment on top
    plot(B(start_idx:end_idx, 2), B(start_idx:end_idx, 3), 'o', ...
        'MarkerSize', 4, 'MarkerEdgeColor', segment_colors(i,:), 'MarkerFaceColor', segment_colors(i,:));
end

% --- ADD BLACK VERTICAL LINES AT TRANSITIONS ---
% The boundaries between segments occur at indices 2, 3, and 4 in the b_idx array
transition_x_coords = [B(b_idx(2), 2), B(b_idx(4), 2)];

for t_x = transition_x_coords
    % Draw a solid black vertical line ('k-'). Use HandleVisibility 'off' so it avoids the legend.
    % If you prefer dashed lines, simply change 'k-' to 'k--'.
    xline(t_x, 'k-', 'LineWidth', 1, 'HandleVisibility', 'off');
end

% Formatting
xlabel('Reaction Coordinate');  
ylabel('Energy (eV)'); 
title('tBuPA-Pt Catalyst for Chalk-Harrod Mechanism Step 1-3');
grid on;

% Only pass handles 1, 2, and 4 to the legend to group B2 & B3
legend([h_lines(1), h_lines(2), h_lines(4)], {'Step 1', 'Step 2', 'Step 3'}, 'Location', 'best');

hold off;

% Save the plot
exportgraphics(gcf, 'NEB_Energy_Profile.png', 'Resolution', 300);