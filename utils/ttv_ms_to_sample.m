function t_s = ttv_ms_to_sample(ms, times)
    t_s = [];
	for i = 1:length(ms)
        t_s = [t_s find(times == ms(i))];
	end
end