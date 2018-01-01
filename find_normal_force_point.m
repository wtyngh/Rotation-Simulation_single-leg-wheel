function normal_force_point = find_normal_force_point(leg_contour,landscape_table)

    contact_area_index_1 = leg_contour.leg_1.contour.y < lookup_table(landscape_table(1,:),landscape_table(2,:),leg_contour.leg_1.contour.x);
    contact_area_index_2 = leg_contour.leg_2.contour.y < lookup_table(landscape_table(1,:),landscape_table(2,:),leg_contour.leg_2.contour.x);

    leg_1_touch_gound = ~isempty( find(contact_area_index_1, 1) ) ;
    leg_2_touch_gound = ~isempty( find(contact_area_index_2, 1) ) ;

    if leg_1_touch_gound

        contact_area_x_1 = leg_contour.leg_1.contour.x( 1, contact_area_index_1 );

        [min_dif_1,min_dif_index_1] = ...
         min( leg_contour.leg_1.contour.y(1, contact_area_index_1) - lookup_table(landscape_table(1,:),landscape_table(2,:),contact_area_x_1));

        contact_point_x = contact_area_x_1(min_dif_index_1);

        normal_force_point.point_1 = ...
         [contact_point_x, lookup_table(landscape_table(1,:),landscape_table(2,:),contact_point_x),min_dif_1];
    else
        normal_force_point.point_1 = [];
    end


    if leg_2_touch_gound

        contact_area_x_2 = leg_contour.leg_2.contour.x( 1, contact_area_index_2 );

        [min_dif_2,min_dif_index_2] = ...
         min( leg_contour.leg_2.contour.y(1, contact_area_index_2) - lookup_table(landscape_table(1,:),landscape_table(2,:),contact_area_x_2));

        contact_point_x = contact_area_x_2(min_dif_index_2);

        normal_force_point.point_2 = [contact_point_x, lookup_table(landscape_table(1,:),landscape_table(2,:),contact_point_x) ,min_dif_2];
    else
        normal_force_point.point_2 = [];
    end


end



