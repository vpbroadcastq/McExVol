function Vab = exvolExact(a, b)
     % Excluded volume between molecules a and b

    if numel(a) == 1 && numel(b) == 1
        % sphere-sphere
        Va = (4/3)*pi*(a(1))^3; Vb = (4/3)*pi*(b(1))^3; 
        Ha = a(1); Hb = b(1); 
        Sa = 4*pi*(a(1)^2); Sb = 4*pi*(b(1)^2); 
    elseif numel(a) == 1 && numel(b) == 2
        % sphere-spherocylinder
        Va = (4/3)*pi*(a(1))^3; Vb = (b(1)^3)*((4*pi)/3)*(1+3*b(2)/2);
        Ha = a(1); Hb = b(1)*(1 + b(2)/2);
        Sa = 4*pi*(a(1)^2); Sb = (4*pi*b(1)^2)*(1 + b(2));
    elseif numel(a) == 1 && numel(b) == 2
        % spherocylinder-sphere
        Vb = (4/3)*pi*(b(1))^3; Va = (a(1)^3)*((4*pi)/3)*(1+3*a(2)/2);
        Hb = b(1); Ha = a(1)*(1 + a(2)/2);
        Sb = 4*pi*(b(1)^2); Sa = (4*pi*a(1)^2)*(1 + a(2));
    elseif numel(a) == 2 && numel(b) == 2
        % spherocylinder-spherocylinder
        Va = (a(1)^3)*((4*pi)/3)*(1+3*a(2)/2); Vb = (b(1)^3)*((4*pi)/3)*(1+3*b(2)/2);
        Ha = a(1)*(1 + a(2)/2); Hb = b(1)*(1 + b(2)/2);
        Sa = (4*pi*a(1)^2)*(1 + a(2)); Sb = (4*pi*b(1)^2)*(1 + b(2));
    end
    
    Vab = Va + Vb + Ha*Sb + Hb*Sa;
end
