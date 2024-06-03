function [gcp]=gini(data)
    data_sort=sort(data,2);
    n=size(data,2);
    
    for i=1:size(data,1)
        x=data_sort(i,:);
        G_num=sum(((2*[1:n])-n-1).*x); % (2*[1:n])-n-1) is a vector stands for weight. eg.n=5,2*[1:5] - 5 - 1 生成 [-4, -2, 0, 2, 4], means the middle one has the lowest weight
        G_den=sum(x)*n;
        gc(i)=(G_num/G_den)*100;
    end
    for i=1:numel(gc)
        gcp(i)=prctile(data(i,:),gc(i));
    end
    gcp=gcp';
end 