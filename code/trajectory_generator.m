function[pos,ori,vel_x,vel_y,yaw_rate,delta] = trajectory_generator(Fs,toa,wp)
    

    %% Generate the figure of 8 trajectory:
    % Specify waypoints, times of arrival, and sampling rate

    % Create trajectory. 
    l_r = 0.18;
    l_f = l_r;
    traj = waypointTrajectory(wp, toa, SampleRate=Fs);
%     [position,orientation,velocity,acceleration,angularVelocity] = traj();
    % Get position.
    k = 0:1/Fs:toa(end);
    [position,orientation,velocity,acceleration,angularVelocity] = lookupPose(traj, k);
    
    pos = position(:,1:2);
    ori = rotvecd(orientation)*pi/180;
    ori = unwrap(ori(:,3));
%     ori(ori > 0.5 *pi) =  ori(ori > 0.5 *pi) - 2*pi;
    
    vel_x = velocity(:,1).*cos(ori) + velocity(:,2).*sin(ori);
    vel_y = velocity(:,2).*cos(ori) - velocity(:,1).*sin(ori);
    yaw_rate = angularVelocity(:,3);
    delta =  unwrap(atan((l_f+l_r)*yaw_rate./vel_x));
%     delta = delta + (delta < 0)*2*pi - pi;

%     x_des = [pos';ori';vel_x';vel_y';yaw_rate'; delta'];
    
    % Plot.
    % figure
    % plot(pos(:,1), pos(:,2))

end
