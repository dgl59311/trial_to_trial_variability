function split_data = find_in_cell(data, names)
    split_data = cell(length(names), 1);
    for i = 1:length(names)
        split_data{i} = data(find(contains({data.name}, names{i})));
    end
end


