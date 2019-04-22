function cRGB = Pej_Color(Color)
if ~isnumeric(Color)
Color = strrep(lower(Color), ' ', '');
switch Color

    case 'blue'
        cRGB = [0 0.4470 0.7410];
    case 'gray'
        cRGB = [0.5 .5 .5];
    case 'red'
        cRGB = [0.8500 0.3250 0.0980];
    case 'yellow'
        cRGB = [0.9290    0.6940    0.1250];
    case 'purple'
        cRGB = [0.4940    0.1840    0.5560];
    case 'green'
        cRGB = [0.4660    0.6740    0.1880];
    case 'lightblue'
        cRGB = [0.3010    0.7450    0.9330];
    case 'darkred'
        cRGB = [0.6350    0.0780    0.1840];
end

else
    Cbox = [
        [0 0.4470 0.7410];
        [0.8500 0.3250 0.0980];
        [0.9290    0.6940    0.1250];
        [0.4940    0.1840    0.5560];
        [0.4660    0.6740    0.1880];
        [0.3010    0.7450    0.9330];
        [0.6350    0.0780    0.1840]];
   
    I = mod(Color, size(Cbox,1));
    I(I==0)= size(Cbox,1);
    cRGB = Cbox(I,:);
end
end