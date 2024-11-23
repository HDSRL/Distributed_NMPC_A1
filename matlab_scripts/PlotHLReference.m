clc; clear all; close all;
% Define file paths
addpath(genpath(pwd))
pathFile = 'HLPath.txt';
velocityFile = 'HLVelocity.txt';
waypointsFile = 'Waypoints.txt';

% Load data
Pr_refined = load(pathFile);
Prd_refined = load(velocityFile);
Pstart = load('Pstart.txt');
Pobs = load("Pobs.txt");

% Load waypoints from CSV
waypointsData = readmatrix(waypointsFile);

% Assuming there are 4 agents and each has 2 rows (one for X and one for Y)
numAgents = size(Pr_refined, 1) / 2; % number of agents

% Parse waypoints into a cell array
waypoints = cell(numAgents, 1);
for i = 1:size(waypointsData, 1)
    agentID = waypointsData(i, 1) + 1; % MATLAB indices start from 1
    waypoint = waypointsData(i, 2:3);
    waypoints{agentID} = [waypoints{agentID}; waypoint];
end

% Plot positions for each agent
% figure(1);
% hold on;
% for i = 1:numAgents
%     % Assuming the first row is X and the second row is Y for each agent
%     plot(Pr_refined(2*i-1, 1:end-10), Pr_refined(2*i, 1:end-10), 'LineWidth', 2, 'DisplayName', ['Agent ' num2str(i) ' Position']);
%     % Plot waypoints
%     plot(waypoints{i}(:,1), waypoints{i}(:,2), 'o', 'LineWidth', 2, 'DisplayName', ['Agent ' num2str(i) ' Waypoints']);
% end
% % Plot square obstacles
% obstacleSize = 0.2;
% halfSize = obstacleSize / 2;  % Half size to calculate lower left corner from center
% for j = 1:size(Pobs, 2)
%     rectangle('Position', [Pobs(1, j) - halfSize, Pobs(2, j) - halfSize, obstacleSize, obstacleSize], ...
%               'FaceColor', 'y', 'EdgeColor', 'k');
% end
% 
% title('Agent Positions');
% xlabel('X Position');
% ylabel('Y Position');
% xlim([-1 10])
% ylim([-5.5 5.5])
% axis square
% legend('show');
% grid on
% 
% hold off;

% Plot velocities for each agent
figure(2);
hold on;
for i = 1:numAgents
    % Assuming the first row is X velocity and the second r`ow is Y velocity for each agent
    % Calculate velocity magnitude as sqrt(vx^2 + vy^2)
    vx = Prd_refined(2*i-1, 1:end-10);
    vy = Prd_refined(2*i, 1:end-10);
    velocityMagnitude = sqrt(vx.^2 + vy.^2);
    plot(1:length(velocityMagnitude), velocityMagnitude, 'LineWidth', 2, 'DisplayName', ['Agent ' num2str(i) ' Velocity']);
end
title('Agent Velocities');
xlabel('Time Step');
ylabel('Velocity Magnitude');
legend('show');
hold off;

path = computePath(Pr_refined, Prd_refined);

autoArrangeFigures();

% Assuming computePath has been called and 'path' matrix is available
numAgents = size(path, 1) / 12; % Assuming 12 states per agent
numSteps = size(path, 2);

agentWidth = 0.6;
agentHeight = 0.4;

figure(3);
hold on;

for i = 1:numAgents
    % Assuming the first row is X and the second row is Y for each agent
    plot(Pr_refined(2*i-1, 1:end-10), Pr_refined(2*i, 1:end-10), 'LineWidth', 2, 'DisplayName', ['Agent ' num2str(i) ' Position']);
    % Plot waypoints
    plot(waypoints{i}(:,1), waypoints{i}(:,2), 'o', 'LineWidth', 2, 'DisplayName', ['Agent ' num2str(i) ' Waypoints']);
end
% Plot square obstacles
obstacleSize = 0.2;
halfSize = obstacleSize / 2;  % Half size to calculate lower left corner from center
for j = 1:size(Pobs, 2)
    rectangle('Position', [Pobs(1, j) - halfSize, Pobs(2, j) - halfSize, obstacleSize, obstacleSize], ...
              'FaceColor', 'y', 'EdgeColor', 'k');
end

axis equal;
grid on;
xlabel('X Position');
ylabel('Y Position');
title('Agent Path Animation');
autoArrangeFigures();

% Calculate XY limits based on the entire path

xlim([-1 10])
ylim([-5.5 5.5])
axis square

% Initialize path plots and agent patches
pathHandles = gobjects(1, numAgents);
agentHandles = gobjects(1, numAgents);
colors = ['r', 'g', 'b', 'k', 'm', 'y', 'c']; % Extend as needed

for agent = 1:numAgents
    colorIndex = mod(agent, length(colors)) + 1;
    pathHandles(agent) = plot(nan, nan, 'Color', colors(colorIndex), 'LineStyle', '--');
    agentHandles(agent) = drawRotatedRect([0, 0], agentWidth, agentHeight, 0, colors(colorIndex));
end

% Animate agents
for i = 1:numSteps
    for agent = 1:numAgents
        % Update path
        set(pathHandles(agent), 'XData', path(1 + (agent-1)*12, 1:i), 'YData', path(2 + (agent-1)*12, 1:i));
        
        % Update agent position and heading
        agentCenter = [path(1 + (agent-1)*12, i), path(2 + (agent-1)*12, i)];
        agentAngle = path(7 + (agent-1)*12, i);
        updateRotatedRect(agentHandles(agent), agentCenter, agentWidth, agentHeight, agentAngle);
    end
    
    drawnow;
    pause(0.01); % Adjust for desired animation speed
end

%% Animation

% Define the source directory with the wildcard to specify .txt files

% DO NOT UPDATE FILES
if (0)
    sourceDirHL = '/home/basit/workspace/SRB_NMPC/build/tmp/*.txt';
    sourceDirLL = '/home/basit/workspace/SRB_NMPC/matlab_scripts/LL/*.txt';
    
    % Define the destination directory
    destinationDir = '/home/basit/workspace/SRB_NMPC/matlab_scripts/cpp_output/';
    
    % Copy all .txt files from source to destination
    status = copyfile(sourceDirHL, destinationDir);
    status = copyfile(sourceDirLL, destinationDir);
    
    
    % Check if the operation was successful
    if status
        disp('Files were copied successfully.');
    else
        disp('Error in copying files.');
    end
end

plotfromcpp = 1;

if (plotfromcpp)
    clearvars; close all;
    cd cpp_output/
    agent_id = 1;
    feet = readmatrix("feet_" + agent_id + ".txt");
    uout = readmatrix("uout_" + agent_id + ".txt");
    tout = readmatrix("tout_" + agent_id + ".txt");
    xout = readmatrix("xout_" + agent_id + ".txt");
    traj = readmatrix("traj_" + agent_id + ".txt");
    cbf = readmatrix("cbf_" + agent_id + ".txt");
    NMPC_solve_time = readmatrix("NMPC_solve_time_" + agent_id + ".txt");
    LL = parseData("agent_" + agent_id + ".txt");
    LL0 = parseData("agent_0.txt");
    LL1 = parseData("agent_1.txt");
    cd ..
end



TF = nnz(NMPC_solve_time);
% TF = 100;
cbf_ = cbf(1:TF-1);
feet = feet(1:TF-1, :);
uout = uout(1:TF-1, :);
tout = tout(1:TF-1, :);
xout = xout(1:TF-1, :);
traj = traj(1:TF-1, :);
solve_time_avg = mean(NMPC_solve_time(1:TF-1))


PlotGRFs(tout,uout)
% PlotCOMPos(tout,traj,xout, LL)
PlotCOMVel(tout,traj,xout, LL, cbf_)
% PlotEulerAngls(tout,traj,xout)
% PlotDerEulerAngls(tout,traj,xout)
% PlotCBF(tout, cbf_)

% autoArrangeFigures();

% Assuming parseData is a function that reads the data and returns a struct
% LL0 = parseData("agent_0.txt");
% LL1 = parseData("agent_1.txt");

% Extract positions and time
x0 = LL0.pos(:, 1); % x position of 0th agent
y0 = LL0.pos(:, 2); % y position of 0th agent

x1 = LL1.pos(:, 1); % x position of 1st agent
y1 = LL1.pos(:, 2); % y position of 1st agent

time = LL0.time(:); % time

idxx = 2000:10000;

x0 = x0(idxx); y0 = y0(idxx);
x1 = x1(idxx); y1 = y1(idxx);
% Calculate the distance between agent 0 and agent 1
dist = sqrt((x1 - x0).^2 + (y1 - y0).^2);

% Calculate the LJ potential
epsilon = 0.5;
sigma = 0.8;
LJ_potential = 1e9 * 4 * epsilon * ( (sigma ./ dist).^12 - (sigma ./ dist).^6 );

% Plot the LJ potential over time
figure(8);
subplot(2, 1, 1);
plot(time(idxx), LJ_potential);
hold on;
xlabel('Time');
ylabel('LJ Potential');
title('Consensus Protocol Cost between Agent 0 and Agent 1');
subplot(2, 1, 2);
hold on;
plot(time(idxx), dist, 'LineWidth', 4);
ylim([0.5 1.5])
xlabel('Time');
ylabel('Distance (m)');
title('Distance between Agent 0 and Agent 1');
grid on;

%%
% ShowAnimation(tout,xout,uout,feet)
%% Supporting Functions
function h = drawRotatedRect(center, width, height, angle, color)
    % Draw a rotated rectangle and return the patch handle
    theta = linspace(0, 2*pi, 100);
    x = width/2 * cos(theta);
    y = height/2 * sin(theta);
    R = [cos(angle) -sin(angle); sin(angle) cos(angle)];
    coords = R * [x; y];
    xRot = coords(1, :) + center(1);
    yRot = coords(2, :) + center(2);
    h = fill(xRot, yRot, color);
end

function updateRotatedRect(h, center, width, height, angle)
    % Update the position and rotation of the rectangle
    theta = linspace(0, 2*pi, 100);
    x = width/2 * cos(theta);
    y = height/2 * sin(theta);
    R = [cos(angle) -sin(angle); sin(angle) cos(angle)];
    coords = R * [x; y];
    xRot = coords(1, :) + center(1);
    yRot = coords(2, :) + center(2);
    set(h, 'XData', xRot, 'YData', yRot);
end

function path = computePath(Pr_refined_, Prd_refined_)
    numAgents = size(Pr_refined_, 1) / 2; % Each agent has 2 rows (x and y)
    numSteps = size(Pr_refined_, 2); % Columns represent time steps
    path = zeros(12 * numAgents, numSteps); % Adjusted for multiple agents
    velocityChangeThreshold = 0.1; % Threshold for considering velocity change significant

    for agent = 1:numAgents
        % Calculate the row index offset for the current agent in the path matrix
        offset = (agent - 1) * 12;

        for i = 1:numSteps
            % Set position and velocity states
            path(1 + offset:6 + offset, i) = [Pr_refined_((agent - 1) * 2 + (1:2), i); 0; Prd_refined_((agent - 1) * 2 + (1:2), i); 0];

            if i > 1
                dx = Prd_refined_(1 + (agent - 1) * 2, i);
                dy = Prd_refined_(2 + (agent - 1) * 2, i);
                velocityChange = sqrt(dx^2 + dy^2);

                if velocityChange > velocityChangeThreshold
                    theta = atan2(dy, dx);
                else
                    % If the change is below the threshold, use the previous theta
                    % This avoids unnecessary theta changes for minor adjustments
                    theta = path(7 + offset, i-1);
                end
            else
                theta = 0; % Initial theta
            end
            path(7 + offset, i) = theta; % Update theta

            % Placeholder for gamma, phi, and their derivatives (set to zero)
            path(8 + offset:11 + offset, i) = 0; % gamma, phi, dot{gamma}, dot{phi}
        end
    end
end
