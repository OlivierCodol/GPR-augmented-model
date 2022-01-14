function str = date_n_time
% Get date and time and make it usable for a file name
str1 = char(datetime);


ix = (strfind(str1,'-'));
y = str1(  (ix(2)+1)  :  (strfind(str1,' ')-1)  );
m = month2num(str1((ix(1)+1):(ix(2)-1)));
d = str1(1:(ix(1)-1));

str2 = regexprep(str1,':','-');
str3 = str2((strfind(str2,' ')+1):end);

str = [y '-' m '-' d '_' str3];

end

function N = month2num(str)

switch str
    case 'Jan'
        n = 1;
    case 'Feb'
        n = 2;
    case 'Mar'
        n = 3;
    case 'Apr'
        n = 4;
    case 'May'
        n = 5;
    case 'Jun'
        n = 6;
    case 'Jul'
        n = 7;
    case 'Aug'
        n = 8;
    case 'Sep'
        n = 9;
    case 'Oct'
        n = 10;
    case 'Nov'
        n = 11;
    case 'Dec'
        n = 12;
end
N = num2str(n);
end