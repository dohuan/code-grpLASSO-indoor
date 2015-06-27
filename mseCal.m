function mse = mseCal(predict,target)
    [m,~] = size(predict);
    err_tmp = predict - target;
    err = sum(err_tmp.^2,2);
    mse = sum(err)/m;
end