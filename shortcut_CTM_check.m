% Shortcut summary goes here
sprintf('%s: %d %d %d %d','C',arrayfun(@(i) check_CTM(tensor_C{i}),1:4))
sprintf('%s: %d %d %d %d %d %d %d %d','T',...
    arrayfun(@(i) check_CTM(tensor_T{i}),1:8))