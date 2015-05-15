function [ output_avgs ] = compute_array_averages_gen2( channel_data, offset, repeats, avg_type )


%%
k=1;
p=1;
datac1 = channel_data;
for j = 1:(size(datac1,1)/(repeats));
    for i = 1:(repeats-offset)
        datac1_1(p,:) = [j datac1(k+offset, :)];
        k=k+1;
        p=p+1;
    end
    k=k+offset;
end

%%

n= size(datac1_1);

n = n(:,1);
runtot = [];
avgs = [];
avg=[];

for j= 1:8
    %j=1;
    
    for i = 1:n
        if(datac1_1(i,1) == j)
            runtot = [runtot; datac1_1(i,:)];
        end
    end
    
    if(~isempty(runtot))
        
%something else
        runtot=[];
    end
end

output_avgs = avgs;

end

