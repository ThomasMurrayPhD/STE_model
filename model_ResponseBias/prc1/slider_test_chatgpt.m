

close all; clear; clc;



% open screen
Screen('Preference', 'SkipSyncTests', 1);
PsychDefaultSetup(2);
screenNumber = max(Screen('Screens'));
[win, windowRect] = Screen('OpenWindow', screenNumber, [210,210,210]);
[screenWidth, screenHeight] = Screen('WindowSize', win);
Screen('BlendFunction', win, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA'); % Alpha blending
[xCentre, yCentre] = RectCenter(windowRect);
Screen('TextSize', win, 30);
fix_size=30;
fix_coords = [...
    -fix_size fix_size 0 0;...
    0 0 -fix_size fix_size];

scale_length=600; % how many pixels long
slider_height=30; % how many pixels high
scale_y = yCentre+20; % y position of slider
scale_rect = [xCentre-(scale_length/2), scale_y-(slider_height/2), xCentre+(scale_length/2), scale_y+(slider_height/2)];
scale_line = [scale_rect(1), scale_y, scale_rect(3), scale_y];
next_button_rect = [xCentre-50, scale_rect(4)+30, xCentre+50, scale_rect(4)+80];
value_rect = [xCentre-(scale_length/2), scale_y-100, xCentre+(scale_length/2), scale_y-50];

% Draw slider
starting_P = 30; % value to start slider at
starting_x = (scale_line(3)*(starting_P/100)) + (scale_line(1)*((100-starting_P)/100)); % x coordinate of starting position
Screen('DrawLine', win, [0,0,0], scale_line(1), scale_line(2), scale_line(3), scale_line(4), 5); % draw scale line
Screen('DrawLine', win, [255,0,0], starting_x, scale_rect(2), starting_x, scale_rect(4), 10); % draw slider 

x_slider = starting_x;

% Draw button
% Screen('FillRect', win, [128, 128, 128], next_button_rect);
% DrawFormattedText(win, 'Next', 'center', 'center', [255,255,255], [],[],[],[],[],next_button_rect);

% make sure mouse is still shown
ShowCursor('Arrow', win);


% Slider state
mouseWasDown = false;
P = nan;
next_trial = false;

max_secs = 5;
trial_start = Screen('Flip', win);

while ~next_trial

    % Timeout check
    if GetSecs - trial_start > max_secs
        timed_out = true;
        break;
    end

    % Get mouse once per frame
    [x, y, buttons] = GetMouse(win);
    mouseIsDown = buttons(1);

    % Check if mouse is inside slider
    within_slider = ...
        (x >= scale_rect(1)) && (x <= scale_rect(3)) && ...
        (y >= scale_rect(2)) && (y <= scale_rect(4));

    % Check if mouse is inside button
    within_button = ...
        (x >= next_button_rect(1)) && (x <= next_button_rect(3)) && ...
        (y >= next_button_rect(2)) && (y <= next_button_rect(4));

    % Slider dragging
    if mouseIsDown && within_slider
        % Restrict x to bounds of slider
        x = min(max(x, scale_rect(1)), scale_rect(3));

        % Convert to percentage
        P = 100 * (x - scale_rect(1)) / scale_length;
    end

    % Scale line
    Screen('DrawLine', win, [0 0 0], scale_line(1), scale_line(2), scale_line(3), scale_line(4), 5);

    % Slider handle
    Screen('DrawLine', win, [255,0,0], x_slider, scale_rect(2), x_slider, scale_rect(4), 10); 

    % Update slider handle x position (only after first movement)
    if ~isnan(P)
        x_slider = scale_rect(1) + (P/100)*scale_length;
        DrawFormattedText(win, sprintf('%d', round(P)), 'center', 'center', [255 255 255], [],[],[],[],[],value_rect);
        Screen('FillRect', win, [128 128 128], next_button_rect); % draw button
        DrawFormattedText(win, 'Next', 'center', 'center', [255 255 255], [],[],[],[],[], next_button_rect);
        Screen('Flip', win);
    end

    % Button click triggers next trial
    if mouseIsDown && ~mouseWasDown && within_button && ~isnan(P)
        next_trial = true;
        Screen('Flip', win);
    end

    % Store previous mouse state
    mouseWasDown = mouseIsDown;
end


WaitSecs(0.5);
sca;



