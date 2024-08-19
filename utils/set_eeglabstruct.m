function [struct] = set_eeglabstruct(input, data)
    input.data = data;
    input.trials = 1;
    input.epoch = [];
    input.event= [];
    input.urevent = [];
    rmbase = MyRmBaseline(input, [-0.5 0]);
    struct = pop_reref(rmbase, []); 
end