function outArr = removeElementByIndex(arr, i)
%Function to remove element of 1-D matrix by it's index
%   %arr - 1-D vector, i - index of removing element

%initialising with input value
outArr=arr;

if(i>length(arr))||(i<1)%if index exseeds limits
    return
end

if(i==length(arr))%deleting last element
    outArr=arr(1:(end-1));
    return
end

if(i==1)%deleting first element
    outArr=arr((i+1):end);
    return
end

%deleting element in the middle
outArr=[arr(1:(i-1)) arr((i+1):end)];
end

