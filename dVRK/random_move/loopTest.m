clc;clear all;close all;
%%
addpath '/home/david/catkin_ws/src/dvrk-ros/dvrk_matlab'
import zaber.motion.Library;
import zaber.motion.binary.Connection;
import zaber.motion.Units;

Library.enableDeviceDbStore();
%%
% The rest of your program goes here
connection = Connection.openSerialPort('/dev/ttyUSB2');
h = uicontrol('Style', 'pushbutton', 'String', 'Stop', ...
    'Callback', 'delete(gcf)');
drawnow()
try
    deviceList = connection.detectDevices();
    fprintf('Found %d devices.\n', deviceList.length);
    device = deviceList(1);
    device.home();
catch exception
    connection.close();
    rethrow(exception);
end

% --- Start Matlab Global Node
try
    rosnode list
catch exp   % Error from rosnode list
    rosinit  % only if error: rosinit
end
% ....

%%
name = 'PSM2';
% jacobian body
topic = strcat(name, '/insertion');
insertion_subscriber = ...
    rossubscriber(topic, rostype.std_msgs_Float64);

insertion_subscriber.NewMessageFcn = @(a, b, c)[];
%%
r = psm(name);
[position, velocity, effort, timestamp] = r.get_state_joint_current();
psm_transPos_init = position(3)
zaber_init = device.getPosition(Units.LENGTH_METRES)
while ishandle(h)
    [position, velocity, effort, timestamp] = r.get_state_joint_current();
    try
        position = insertion_subscriber.LatestMessage.Data;
    catch
        position = psm_transPos_init;
    end
    posCmd = position-psm_transPos_init + zaber_init;
    posCmd = max(0,posCmd);
    if (posCmd <= 0.15) && (posCmd >= 0)
        device.moveAbsolute(posCmd, Units.LENGTH_METRES);
    end
    pause(0.01);
end

connection.close();

