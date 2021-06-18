function jday = julian_day( date )
%@function. Calculate the specify data is the jday of that year.
%           For example, January 3th is the third day of the year.
%@param @date. The date in form of yyyymmdd.

    year = floor(date / 10000);
    month = mod(floor(date / 100), 100);
    day = mod(date, 100);

    Feb = 28;
    
    if (mod(year, 100) == 0)
        if (mod(year, 400) == 0)
            Feb = 29;
        end
    elseif (mod(year, 4) == 0)
        Feb = 29;
    end
    
    jday = 0;
    
    switch month
        case 1
            jday = day;
        case 2
            jday = 31 + day;
        case 3
            jday = 31 + Feb + day;
        case 4
            jday = 62 + Feb + day;
        case 5
            jday = 92 + Feb + day;
        case 6
            jday = 123 + Feb + day;
        case 7
            jday = 153 + Feb + day;
        case 8
            jday = 184 + Feb + day;
        case 9
            jday = 215 + Feb + day;
        case 10
            jday = 245 + Feb + day;
        case 11
            jday = 276 + Feb + day;
        case 12
            jday = 306 + Feb + day;
    end
   
    return;
    
    

