function strN = dbkNum2StrFixDecPlaces(N, numDigits)
    
    N = roundoff(N, numDigits);
    strN = num2str(N);
    
    if strfind(strN, '.')
        % There is currently a decimal
        % Count how many digits after the decimal place
        numCurrDigits = length(strN(strfind(strN, '.')+1:end));
        
        strN = [strN, repmat('0', 1, numDigits - numCurrDigits)];
    else
        % There is no decimal place in the number as it currently is
        strN = [strN, '.' repmat('0', 1, numDigits)];
    end

end