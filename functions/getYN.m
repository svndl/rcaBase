function resp = getYN(message)
% yes = getYN(message)
% Echo message to the command window and wait for the user to type a yY/nN.
% Return True on y or Y and False on n or N

% 3/3/06: changed to work with java vm
% 8/20/15: changed to use built-in matlab function, input

while 1
    fprintf(message);
    resp = input('\nEnter Y or N: ', 's');     % This two step gets
    fprintf('\n');         % fprintf to handle \n embedded in the message
    %fprintf('%s\n', resp);  
    if (resp == 'y' || resp == 'Y')
        resp = true;
        break
    elseif (resp == 'n' || resp == 'N')
        resp = false;
        break
    else
        fprintf('I do not understand that response. Please type y, Y, n, or N.');
    end
end
            