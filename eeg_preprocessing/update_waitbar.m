%15-may-2017    17:40
%Samran         Auckland

%updates waitbar w, with resolution, taking current position and the total
%number of steps, and returning the updated position

function [waitbar_step] = update_waitbar (w, waitbar_resolution, waitbar_step, ...
    waitbar_steps, message)

waitbar_step = waitbar_step + waitbar_resolution;

if (isempty(message) || strcmp(message, ''))
    waitbar (waitbar_step/waitbar_steps, w);
else
    waitbar (waitbar_step/waitbar_steps, w, message);
end

end