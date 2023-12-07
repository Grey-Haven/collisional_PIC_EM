function logical_coordinate = physical_to_logical_space(loc,padding,stagger,offset,diff)
    index_offset = 1;
    logical_coordinate = padding+index_offset+(loc-stagger+offset)./diff;
end