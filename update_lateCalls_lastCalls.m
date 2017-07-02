function [ lateCalls ] = update_lateCalls_lastCalls( lateCalls, R, updateCoord )
% PROBLEM: Both updates do not take into account the change in the "last"
% agent due to lack of information. So the new lastCalls should not be
% trusted.

% coordinates of added agent
i_p = updateCoord(1); % group
s_p = updateCoord(2); % shift
% coordinates of removed agent
i_m = updateCoord(3); % group
s_m = updateCoord(4); % shift

% % For each last call of the removed agent, move it to the lateCalls,
% % because it can no longer be picked up and needed to be queued.
% nLast = size(lastCalls,2);
% k = 1;
% while k <= nLast
%     if lastCalls(6,k) == i_m;
%        lateCalls = [lateCalls, lastCalls(1:5,k)];
%        lastCalls(:,k) = []; 
%        nLast = nLast - 1;
%     else
%         k = k + 1;
%     end
% end

% Remove the late call when additional agent of certain group is added.
nLate = size(lateCalls,2);
nextAvail = 0;
j = 1;
% for every late call that can be picked up by group i
if ~isempty(lateCalls)
    while j <= nLate && R(i_p, lateCalls(2,j)) < 1000
        % check calls in lateCalls that are is not expired or picked up
        % late yet (still in queue)
        if lateCalls(3,j) > nextAvail && lateCalls(4,j) > nextAvail
            % available after mean service time
            if ~isnan(lateCalls(5,j))
                nextAvail = lateCalls(1,j) + lateCalls(5,j);
            else
                nextAvail = lateCalls(1,j) + exprnd(meanST(lateCalls(2,j),i_p));
            end
            lateCalls(:,j) = [];
            nLate = nLate - 1;
        else
            j = j + 1;
        end
    end
end
end

